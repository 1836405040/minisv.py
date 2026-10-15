"""Orientation-aware, sequence-based structural junction review.

This experimental reviewer never promotes a structural event to high confidence:
independent positive controls and parental haplotype validation are still required.
Coordinates are zero-based interbase positions; L/R names the retained reference side.
"""
import hashlib
import json
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import edlib
import pysam


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def junctions(fields):
    info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
    from .trio_types import normalize_record
    fields, info, _length, _error = normalize_record(fields)
    typ = info.get('SVTYPE', '').upper()
    chrom, pos = fields[0], int(fields[1])
    if typ in ('INV', 'DUP'):
        end = int(info.get('END', pos + abs(int(info.get('SVLEN', '0')))))
        if end <= pos:
            raise ValueError('INVALID_INTERVAL')
        if typ == 'INV':
            return [((chrom, pos, 'L'), (chrom, end, 'L')),
                    ((chrom, pos, 'R'), (chrom, end, 'R'))]
        # Only a tandem-duplication junction is tested. Dispersed copies require
        # a destination breakpoint that a symbolic DUP record does not encode.
        return [((chrom, end, 'L'), (chrom, pos, 'R'))]
    if typ != 'BND':
        raise ValueError('UNSUPPORTED_STRUCTURAL_TYPE')
    alt = fields[4]
    match = re.fullmatch(r'([^\[\]]*)([\[\]])([^:\[\]]+):(\d+)\2([^\[\]]*)', alt)
    if not match:
        raise ValueError('UNRESOLVED_BND_ALLELE')
    before, bracket, remote, remote_pos, after = match.groups()
    if bool(before) == bool(after) or (before or after).upper() != fields[3].upper() or len(fields[3]) != 1:
        raise ValueError('BND_INSERTED_SEQUENCE_OR_COMPLEX_ANCHOR')
    local_side = 'L' if before else 'R'
    remote_side = 'R' if bracket == '[' else 'L'
    return [((chrom, pos if before else pos-1, local_side),
             (remote, int(remote_pos)-(remote_side == 'R'), remote_side))]


def canonical(junction_list):
    return tuple(sorted(tuple(sorted(pair)) for pair in junction_list))


def models(fasta, pair, flank=200):
    arms, refs = [], []
    for chrom, boundary, side in pair:
        if chrom not in fasta.references or boundary < flank or boundary+flank > fasta.get_reference_length(chrom):
            raise ValueError('REFERENCE_INTERVAL_UNAVAILABLE')
        seq = fasta.fetch(chrom, boundary-flank, boundary+flank).upper()
        if set(seq)-set('ACGT'):
            raise ValueError('REFERENCE_SEQUENCE_UNRESOLVED')
        arms.append(seq[:flank] if side == 'L' else rc(seq[flank:]))
        refs.append(seq)
    return arms[0]+rc(arms[1]), refs


def distance(model, sequence):
    return min(edlib.align(model, sequence, mode='HW')['editDistance'],
               edlib.align(rc(model), sequence, mode='HW')['editDistance'])


def molecule_evidence(bam, pair, alt, refs, rate, flank=200):
    molecules = {}
    with pysam.AlignmentFile(bam) as source:
        for chrom, boundary, _side in pair:
            if chrom not in source.references:
                raise ValueError('BAM_CONTIG_UNAVAILABLE')
            for read in source.fetch(chrom, max(0,boundary-flank), boundary+flank):
                if read.is_secondary or read.is_unmapped or read.is_duplicate or read.is_qcfail or read.mapping_quality < 20 or not read.query_sequence:
                    continue
                old = molecules.get(read.query_name)
                # Full primary sequence is preferred over hard-clipped supplementary sequence.
                if old is None or len(read.query_sequence) > len(old[0]):
                    molecules[read.query_name] = (read.get_forward_sequence().upper(), read.is_reverse)
    alt_reads, ref_names = [], [[], []]
    for name, (seq, reverse) in molecules.items():
        da = distance(alt, seq)
        dr = [distance(ref, seq) for ref in refs]
        if da <= rate*len(alt) and min(dr)-da >= 10:
            orientation = edlib.align(rc(alt), seq, mode='HW')['editDistance'] < edlib.align(alt, seq, mode='HW')['editDistance']
            alt_reads.append(dict(name=name, reverse=orientation, alt_distance=da, reference_distances=dr))
        for i, d in enumerate(dr):
            if d <= rate*len(refs[i]) and da-d >= 10:
                ref_names[i].append(name)
    return dict(fetched_molecules=len(molecules), alt_specific=len(alt_reads),
                alt_names=[r['name'] for r in alt_reads],
                alt_strands=dict(Counter('reverse' if r['reverse'] else 'forward' for r in alt_reads)),
                reference_names=ref_names, reference_counts=list(map(len,ref_names)), alt_reads=alt_reads)


def review_structural(job):
    event, reference, bams, rate = job
    event = dict(event, junction_evidence=[], audit_flags=list(event['audit_flags']))
    flags = event['audit_flags']
    try:
        with pysam.FastaFile(reference) as fasta:
            for pair in event['junctions']:
                alt, refs = models(fasta, pair)
                if min(distance(alt, ref) for ref in refs) < 10:
                    flags.append('JUNCTION_NOT_DISCRIMINATING')
                evidence = {role:molecule_evidence(bam,pair,alt,refs,rate) for role,bam in bams.items()}
                event['junction_evidence'].append(dict(endpoints=pair,alt_sequence=alt,support=evidence))
                if evidence['C']['alt_specific'] < 3:
                    flags.append('CHILD_JUNCTION_NOT_REPRODUCED')
                if len(evidence['C']['alt_strands']) < 2:
                    flags.append('SINGLE_STRAND_SUPPORT')
                if any(evidence[r]['alt_specific'] for r in ('F','M')):
                    flags.append('PARENT_TARGET_JUNCTION')
                if any(min(evidence[r]['reference_counts']) < 10 for r in ('F','M')):
                    flags.append('PARENT_BREAKPOINT_CALLABILITY_UNRESOLVED')
    except (ValueError, KeyError) as error:
        flags.append(str(error))
    # Junction evidence is deliberately distinct from independent DNM confirmation.
    flags.append('STRUCTURAL_MODEL_VALIDATION_REQUIRED')
    if event['type'] == 'DUP':
        flags.append('DUP_COPY_NUMBER_AND_DESTINATION_UNRESOLVED')
    if any(endpoint[0].removeprefix('chr') not in {str(i) for i in range(1,23)}
           for pair in event['junctions'] for endpoint in pair):
        flags.append('PLOIDY_NOT_VALIDATED')
    flags = sorted(set(flags))
    event['audit_flags'] = flags
    event['tier'] = ('REVIEW_REQUIRED' if set(flags) <= {
        'STRUCTURAL_MODEL_VALIDATION_REQUIRED','DUP_COPY_NUMBER_AND_DESTINATION_UNRESOLVED'}
        else 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE')
    event['junction_support_reproduced'] = bool(event['junction_evidence']) and all(
        j['support']['C']['alt_specific'] >= 3 for j in event['junction_evidence'])
    return event


def structural_events(records, family_id, minimum=50):
    groups, unresolved = {}, []
    for record in records:
        fields = record['original_record'].split('\t')
        info = dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
        try:
            if fields[6] not in ('PASS','.'):
                raise ValueError('EXCLUDED_INPUT_FILTER')
            js = canonical(junctions(fields))
            from .trio_types import normalize_record
            _fields, normalized_info, _length, _error = normalize_record(fields)
            typ = normalized_info['SVTYPE']
            if typ != 'BND' and abs(js[0][0][1]-js[0][1][1]) < minimum:
                raise ValueError('BELOW_MINIMUM_LENGTH')
        except (ValueError, KeyError) as error:
            unresolved.append(dict(record, status=str(error)))
            continue
        key = (record['sample'],typ,js)
        if key not in groups:
            length = None if typ == 'BND' else abs(js[0][0][1]-js[0][1][1])
            event_id = hashlib.sha256(repr(key).encode()).hexdigest()[:16]
            groups[key] = dict(event_id=f"{record['sample']}.STRUCT.{event_id}",
                sample=record['sample'],family_id=family_id,type=typ,delta=None,length=length,
                variant_class='de_novo_SV',basis='REFERENCE_JUNCTION',junctions=js,
                chrom=js[0][0][0],start=js[0][0][1],end=js[0][1][1],
                branches=['STRUCTURAL_JUNCTION_REVIEW'],records=[],audit_flags=[])
        event = groups[key]
        event['records'].append(record)
        if 'IMPRECISE' in fields[7].split(';'):
            event['audit_flags'].append('IMPRECISE_BREAKPOINT')
        if typ == 'BND' and not info.get('MATEID'):
            event['audit_flags'].append('BND_MATE_LINK_UNRESOLVED')
    # A mate ID alone is not proof of a reciprocal record with matching geometry.
    for event in groups.values():
        if event['type'] == 'BND':
            ids = set()
            for r in event['records']:
                info = dict(x.split('=',1) for x in r['original_record'].split('\t')[7].split(';') if '=' in x)
                ids.add((r.get('source'), info.get('MINISV_ORIGINAL_ID',r['id'])))
            for record in event['records']:
                info = dict(x.split('=',1) for x in record['original_record'].split('\t')[7].split(';') if '=' in x)
                if not {(record.get('source'),mid) for mid in info.get('MATEID','').split(',')} & ids:
                    event['audit_flags'].append('BND_MATE_LINK_UNRESOLVED')
    return list(groups.values()), unresolved


def review_pending(records, config, workers, out, rate=.08):
    family_id = config['father']['sample']+'__'+config['mother']['sample']
    events, unresolved = structural_events(records,family_id,config.get('minimum_variant_length',50))
    intervals = []
    if config.get('structural_exclude_bed'):
        with open(config['structural_exclude_bed']) as source:
            intervals = [(a[0],int(a[1]),int(a[2])) for line in source
                         if line.strip() and not line.startswith('#') for a in [line.split()]]
    retained = []
    for event in events:
        endpoints = [ep for pair in event['junctions'] for ep in pair]
        outside = any(ep[0].removeprefix('chr') not in {str(i) for i in range(1,23)} for ep in endpoints)
        for chrom, left, right in intervals:
            positions = [ep[1] for ep in endpoints if ep[0] == chrom]
            if positions and min(positions) < right and max(positions) >= left:
                outside = True
        if outside:
            unresolved.append(dict(event,status='OUTSIDE_REGULAR_REGION'))
        else:
            retained.append(event)
    events = retained
    children = {c['sample']:c for c in config['children']}
    jobs = [(e,config['reference'],dict(C=children[e['sample']]['bam'],
        F=config['father']['bam'],M=config['mother']['bam']),rate) for e in events]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        reviewed = list(pool.map(review_structural,jobs))
    # Shared junctions in siblings are retained as a transmission/mosaicism flag.
    for event in reviewed:
        if any(other['sample'] != event['sample'] and other['junctions'] == event['junctions'] for other in reviewed):
            event['audit_flags'] = sorted(set(event['audit_flags']+['JUNCTION_SHARED_WITH_SIBLING']))
            event['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    for name, rows in [('structural_reviewed.jsonl',reviewed),('structural_unresolved.jsonl',unresolved)]:
        with (Path(out)/name).open('w') as target:
            for row in rows:
                target.write(json.dumps(row)+'\n')
    return reviewed


def run_manifest(manifest, output, workers=8, rate=.08):
    """Review caller structural candidates without rerunning the indel pipeline."""
    import gzip
    import time
    config = json.loads(Path(manifest).read_text())
    if workers < 1 or not 0 < rate < 1:
        raise ValueError('Invalid workers/rate')
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    (out/'status.json').write_text(json.dumps(dict(stage='RUNNING')))
    try:
        with pysam.FastaFile(config['reference']) as ref:
            lengths = dict(zip(ref.references,ref.lengths))
        for member in [config['father'],config['mother']]+config['children']:
            with pysam.AlignmentFile(member['bam']) as bam:
                if not bam.has_index() or any(lengths.get(c) != n for c,n in zip(bam.references,bam.lengths)):
                    raise ValueError('BAM index/reference mismatch: '+member['sample'])
        records = []
        for child in config['children']:
            paths = list(child.get('candidate_vcfs',[]))
            if child.get('candidate_vcf'):
                paths.append(child['candidate_vcf'])
            if not paths:
                raise ValueError('No candidate VCF: '+child['sample'])
            for path in dict.fromkeys(paths):
                opener = gzip.open if str(path).endswith('.gz') else open
                with opener(path,'rt') as source:
                    for line in source:
                        if line.startswith('#') or not line.strip():
                            continue
                        fields = line.rstrip().split('\t')
                        info = dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
                        if info.get('SVTYPE','').upper() not in ('INS','DEL'):
                            records.append(dict(sample=child['sample'],id=fields[2],source=str(path),
                                type=info.get('SVTYPE'),original_record=line.rstrip()))
        with (out/'structural_input_records.jsonl').open('w') as target:
            for record in records:
                target.write(json.dumps(record)+'\n')
        reviewed = review_pending(records,config,workers,out,rate)
        summary = dict(stage='COMPLETE',input_records=len(records),reviewed_events=len(reviewed),
            tiers=dict(Counter(e['tier'] for e in reviewed)),
            types=dict(Counter(e['type'] for e in reviewed)),
            child_junction_reproduced=dict(Counter(e['type'] for e in reviewed if e['junction_support_reproduced'])),
            flags=dict(Counter(f for e in reviewed for f in e['audit_flags'])))
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        (out/'provenance.json').write_text(json.dumps(dict(config=config,workers=workers,rate=rate,
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),truth_used=False,
            time=time.strftime('%Y-%m-%dT%H:%M:%S%z')),indent=2))
        (out/'status.json').write_text(json.dumps(summary,indent=2))
        return summary
    except Exception as error:
        (out/'status.json').write_text(json.dumps(dict(stage='FAILED',error=str(error))))
        raise
