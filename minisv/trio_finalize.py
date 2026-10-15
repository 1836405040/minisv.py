"""Unify catalog allele discovery and audited C3 evidence without truth input.

Local events use parent-relative sequence semantics, not reference VCF semantics.
Cached full-read distances are sufficient statistics for the error-rate sweep.
"""
import concurrent.futures
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def classify(local, child_key, delta, rate, min_reads=3, min_depth=10, margin=10, minimum=50):
    counts, novel, depth, names = {}, {}, {}, {}
    for role in ('C', 'F', 'M'):
        reads = local['full_reads'][role]
        if len({r['name'] for r in reads}) != len(reads):
            raise ValueError('Duplicate molecule in cached evidence')
        depth[role] = len(reads)
        target = [r for r in reads if r['distances'][child_key] == min(r['distances'].values())
                  and r['distances'][child_key] / max(1, r['length']) <= rate]
        counts[role] = len(target)
        specific = [r for r in target if min(v for k, v in r['distances'].items() if k.startswith(('F', 'M')))
                    - r['distances'][child_key] >= margin]
        novel[role] = len(specific)
        names[role] = [r['name'] for r in specific]
    reasons = []
    if abs(delta) < minimum:
        reasons.append(f'BELOW_{minimum}_BP')
    if novel['C'] < min_reads:
        reasons.append('CHILD_SPECIFIC_READS_LT_MIN')
    if min(depth['F'], depth['M']) < min_depth:
        reasons.append('PARENT_SPANNING_DEPTH_LT_MIN')
    if counts['F'] or counts['M']:
        reasons.append('PARENT_TARGET_SUPPORT')
    return dict(status='RULE_PASS' if not reasons else 'REVIEW', reasons=reasons,
                depth=depth, target_support=counts, novel_support=novel, novel_read_names=names)


def process_locus(job):
    sample, path, rates = job[:3]
    minimum = job[3] if len(job) > 3 else 50
    raw = Path(path).read_bytes()
    local = json.loads(raw)
    if local['child'] != sample:
        raise ValueError(f'Sample mismatch: {path}')
    if not local.get('full_reads'):
        raise ValueError(f'Missing full BAM evidence: {path}')
    m = re.fullmatch(r'(.+)_(\d+)_(\d+)_.*', local['trid'])
    if not m:
        raise ValueError(f'Invalid catalog coordinates: {path}')
    transmissions = local['transmissions']
    if not transmissions:
        return dict(events=[], input_hash=hashlib.sha256(raw).hexdigest(), path=path,
                    warning='NO_SUPPORTED_DIPLOID_TRANSMISSION')
    best = transmissions[0]['total_edit_distance']
    ties = [t for t in transmissions if t['total_edit_distance'] == best]
    groups = {}
    for t in ties:
        for ck, pk in [(t['paternal_child'], t['father']), (t['maternal_child'], t['mother'])]:
            cs, ps = local['models'][ck], local['models'][pk]
            delta = len(cs) - len(ps)
            if abs(delta) < minimum:
                continue
            key = (cs, delta)
            if key not in groups:
                groups[key] = dict(child_key=ck, child_sequence=cs, delta=delta, parents={})
            groups[key]['parents'][pk] = ps
    events = []
    for (seq, delta), g in groups.items():
        key = g['child_key']
        event_id = f'{sample}.{local["trid"]}.{digest(seq)[:16]}.{delta:+d}'
        event = dict(event_id=event_id, sample=sample, trid=local['trid'], chrom=m[1],
                     start=int(m[2]), end=int(m[3]), type='INS' if delta > 0 else 'DEL',
                     delta=delta, basis='PARENT_ALLELE', child_sequence=seq,
                     child_model=key, parent_sequences=g['parents'],
                     best_transmission_ties=len(ties), origin_ambiguous=len({k[0] for k in g['parents']}) > 1,
                     branches=['LOCAL_SEQUENCE'], c3_event_ids=[], c3_record_ids=[],
                     evidence={str(rate): classify(local, key, delta, rate, minimum=minimum) for rate in rates})
        events.append(event)
    return dict(events=events, input_hash=hashlib.sha256(raw).hexdigest(), path=path)


def related(c3, local):
    if c3['sample'] != local['sample'] or c3['chrom'] != local['chrom'] or c3['type'] != local['type']:
        return False
    if local['trid'] not in c3['overlapping_local_loci']:
        return False
    if min(c3['length'], abs(local['delta'])) / max(c3['length'], abs(local['delta'])) < .5:
        return False
    names = set().union(*(set(r['names']) for r in c3['records']))
    # Association is deliberately tagged as compatibility, not sequence identity.
    supported = set().union(*(set(v['novel_read_names']['C']) for v in local['evidence'].values()))
    return len(names & supported) >= 2


def run(manifest_path, output, rates=(.05, .06, .08), workers=16, minimum=50):
    if not rates or any(not 0 < r < 1 for r in rates):
        raise ValueError('Error rates must lie between zero and one')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads(Path(manifest_path).read_text())
    jobs = [(s['sample'], str(p), list(rates), minimum) for s in manifest['samples']
            for p in sorted(Path(s['local_directory']).glob('*.json'))]
    if not jobs:
        raise ValueError('No local evidence files')
    events, hashes, warnings = [], {}, []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(process_locus, jobs, chunksize=8):
            events.extend(result['events'])
            hashes[result['path']] = result['input_hash']
            if result.get('warning'):
                warnings.append(dict(path=result['path'], warning=result['warning']))
    c3 = json.loads(Path(manifest['c3_audit']).read_text())
    by_locus = {}
    for e in events:
        by_locus.setdefault((e['sample'], e['trid']), []).append(e)
    for c in c3:
        hits = [e for trid in c['overlapping_local_loci'] for e in by_locus.get((c['sample'], trid), []) if related(c, e)]
        for e in hits:
            if 'C3_COMPATIBLE' not in e['branches']:
                e['branches'].append('C3_COMPATIBLE')
            e['c3_event_ids'].append(c['event_id'])
            e['c3_record_ids'].extend(c['record_ids'])
        # Preserve C3 representations and their original evidence, even without a local match.
        c['compatible_local_event_ids'] = [e['event_id'] for e in hits]
        c.pop('compatible_truth_ids', None)
    summary = {}
    for rate in rates:
        selected = [e for e in events if e['evidence'][str(rate)]['status'] == 'RULE_PASS']
        summary[str(rate)] = dict(rule_pass_events=len(selected),
            rule_pass_loci=len({(e['sample'], e['trid']) for e in selected}),
            sample_events=dict(Counter(e['sample'] for e in selected)),
            sample_loci={s['sample']: len({e['trid'] for e in selected if e['sample'] == s['sample']}) for s in manifest['samples']},
            local_only_events=sum(not e['c3_event_ids'] for e in selected))
        with (out / f'rule_pass_{round(rate * 100):02d}.jsonl').open('w') as f:
            for e in selected:
                f.write(json.dumps(e) + '\n')
    with (out / 'all_local_events.jsonl').open('w') as f:
        for e in events:
            f.write(json.dumps(e) + '\n')
    (out / 'c3_evidence.json').write_text(json.dumps(c3, indent=2))
    with (out / 'unified_candidates.jsonl').open('w') as f:
        for e in events:
            f.write(json.dumps(dict(kind='LOCAL_ALLELE_EVENT', **e)) + '\n')
        for c in c3:
            f.write(json.dumps(dict(kind='C3_REPRESENTATION', **c)) + '\n')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    (out / 'provenance.json').write_text(json.dumps(dict(manifest=manifest, rates=list(rates), workers=workers,
        local_files=len(jobs), local_events=len(events), warnings=warnings, hashes=hashes,
        c3_sha256=hashlib.sha256(Path(manifest['c3_audit']).read_bytes()).hexdigest(),
        module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        criteria=dict(min_child_specific_reads=3, min_parent_depth=10, parent_target_support=0, edit_margin=10,
                      minimum_variant_length=minimum),
        interpretation='RULE_PASS is an uncalibrated candidate tier, not validated DNM'), indent=2))
    (out / 'status.json').write_text(json.dumps(dict(stage='COMPLETE', local_files=len(jobs))))
    return summary
