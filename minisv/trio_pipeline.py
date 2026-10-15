"""One-entry workflow from prepared family inputs to reviewed DNM candidates."""
import concurrent.futures
import gzip
import hashlib
import json
import shutil
import subprocess
import time
from collections import Counter
from pathlib import Path

import pysam

from .trio_finalize import run as finalize
from .trio_local import compare_locus, full_read_evidence
from .trio_locus import candidate_region, filter_locus
from .trio_review import assign_tiers, review
from .trio_c3_only import group_calls, review_deletion
from .trio_insertion import review_insertion
from .trio_prepare import prepare_candidates
from .trio_types import partition_for_review
from .trio_io import compressed_writer, combine_parents, extract_candidate_reads


def state(out, stage, **values):
    payload = dict(stage=stage,time=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**values)
    temp = out/'status.tmp'
    temp.write_text(json.dumps(payload,indent=2))
    temp.replace(out/'status.json')


def read_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as source:
        for line in source:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            info = dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
            yield fields,info


def family_bams(config, child):
    return dict(C=child['bam'], F=config['father']['bam'], M=config['mother']['bam'])


def generate_local(job):
    directory, child, father, mother, trid, chrom, start, end, bams, output = job
    data = compare_locus(directory,child,chrom,start,end,trid,father=father,mother=mother)
    data['full_reads'] = {role:full_read_evidence(path,chrom,start,end,data['models']) for role,path in bams.items()}
    data['source_bams'] = bams
    Path(output).write_text(json.dumps(data))
    return output


def review_job(job):
    path, events, reference, bams, rate = job
    return review(json.loads(Path(path).read_text()),events,reference,bams,rate)


def review_c3_only(c3, config, children, rate, workers, out):
    groups = group_calls(c3, minimum=config.get('minimum_variant_length', 50))
    family_id = config['father']['sample']+'__'+config['mother']['sample']
    jobs = [(calls,config['reference'],family_bams(config,children[sample]),rate,family_id)
            for (sample, _chrom, _pos, _length),calls in groups.items()]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        reviewed = list(pool.map(review_deletion,jobs))
        if config.get('parent_haplotype_alignments'):
            insertions = [r for r in c3 if r['type'] == 'INS' and r['basis'] == 'REFERENCE'
                          and r['length'] >= max(30,config.get('minimum_variant_length',50))
                          and r.get('ref') and r.get('alt')]
            ins_jobs = [(r,config['reference'],family_bams(config,children[r['sample']]),
                         config['parent_haplotype_alignments'],rate,family_id) for r in insertions]
            reviewed.extend(pool.map(review_insertion,ins_jobs))
    with (out/'c3_only_reviewed.jsonl').open('w') as target:
        for event in reviewed:
            target.write(json.dumps(event)+'\n')
    return reviewed


def run_c3(child, config, out, threads):
    """Reuse explicit evidence or align complete candidate molecules to parents."""
    pending = out/'structural_pending.jsonl'
    pending.write_text('')
    normalized = out/'c3_supported_input.vcf'
    source = child.get('c3_vcf') or child['candidate_vcf']
    count = partition_for_review(source, normalized, pending)
    if child.get('c3_vcf') or not count:
        # Cached evidence must pass the same type routing as fresh candidates.
        if child.get('c3_vcf') and child.get('candidate_vcf'):
            partition_for_review(child['candidate_vcf'], out/'candidate_supported_input.vcf', pending)
        return str(normalized)
    from .cli import EvalOpt,opt
    from .filtercaller import othercaller_filterasm
    from .read_parser import load_reads
    vcf = str(normalized)
    minimum = config.get('minimum_variant_length', 50)
    paf, gsv = child.get('parent_paf'), child.get('parent_gsv')
    if not paf:
        haps = config.get('parent_assemblies', [])
        if len(haps) != 4:
            raise ValueError('Four parental haplotype FASTAs or existing C3/PAF evidence required')
        fq = out/'child_reads.fq.gz'
        extracted = extract_candidate_reads(vcf,child['bam'],fq,threads)
        parents = out/'parents.fa'
        combine_parents(haps,config['father']['sample'],config['mother']['sample'],parents)
        paf = str(out/'parents.paf.gz')
        command = [config.get('minimap2','minimap2'),'--ds','-t',str(threads),'-cx','map-hifi','-s50','-I100g','--secondary=no',str(parents),str(fq)]
        (out/'alignment_command.json').write_text(json.dumps(command,indent=2))
        with (out/'minimap2.log').open('w') as log, compressed_writer(paf,'wb',min(16,threads)) as target:
            with subprocess.Popen(command,stdout=subprocess.PIPE,stderr=log) as proc:
                shutil.copyfileobj(proc.stdout,target)
                if proc.wait():
                    raise RuntimeError('minimap2 failed; see minimap2.log')
        with gzip.open(paf,'rt') as source:
            if extracted['requested'] and not any(line.strip() for line in source):
                raise RuntimeError('No parental alignments for requested candidate molecules')
    if not gsv:
        gsv = str(out/'parents.gsv.gz')
        options = opt(cen={},min_mapq=0,min_mapq_end=0,min_len=minimum,min_frac=.7,max_cnt_10k=5)
        with compressed_writer(gsv,'wt',min(16,threads)) as target:
            load_reads(paf,options,target)
    legacy = out/'legacy_C3.vcf'
    options = EvalOpt(min_len=minimum,min_count=3,read_len_ratio=.6,min_len_ratio=.6,print_all=True)
    options.max_diff_len=50
    with legacy.open('w') as target:
        othercaller_filterasm(vcf,options,vcf,gsv,str(out/'legacy_C3.stat'),'',target)
    selected = out/'locus_C3.vcf'
    filter_locus(str(legacy),child['bam'],gsv,str(selected),str(out/'locus_C3.audit.tsv'),min_count=3,paf=paf)
    return str(selected)


def c3_events(sample,vcf,local_paths):
    """Exact representation dedup; local sequence events receive final dedup."""
    catalog = []
    for path in local_paths:
        data = json.loads(Path(path).read_text())
        parts = data['trid'].split('_')
        catalog.append((parts[0],int(parts[1]),int(parts[2]),data['trid']))
    groups = {}
    for a,info in read_records(vcf):
        if info.get('SVTYPE') not in ('INS', 'DEL'):
            raise ValueError('Route structural records through partition_for_review before c3_events')
        names = sorted(set(info.get('RNAMES','').split(','))-{'','.'})
        basis = info.get('DELTA_BASIS',info.get('SVLEN_SOURCE','REFERENCE'))
        key = (a[0],a[1],a[3],a[4],info['SVTYPE'],info['SVLEN'],basis)
        if basis != 'REFERENCE':
            key += (tuple(names),)
        groups.setdefault(key,[]).append(dict(id=a[2],names=names,info=info))
    rows = []
    for n,(key,records) in enumerate(groups.items(),1):
        chrom,pos,ref,alt,typ,length,basis = key[:7]
        trids = {r['info']['TRID'] for r in records if r['info'].get('TRID')}
        trids.update(trid for c,s,e,trid in catalog if c==chrom and s<=int(pos)-1<e)
        rows.append(dict(event_id=f'{sample}.C3.{n}',sample=sample,chrom=chrom,pos0=int(pos)-1,
                         type=typ,length=abs(int(length)),basis=basis,ref=ref,alt=alt,record_ids=[r['id'] for r in records],
                         records=records,overlapping_local_loci=sorted(trids),status='C3_ONLY_NOT_FINAL_CONFIDENCE'))
    return rows


def validate(config):
    minimum = config.get('minimum_variant_length',50)
    if isinstance(minimum,bool) or not isinstance(minimum,int) or minimum < 1:
        raise ValueError('minimum_variant_length must be a positive integer')
    if config.get('confidence_policy','strict') != 'strict':
        raise ValueError('Only strict confidence_policy is supported')
    if config.get('parent_readability_policy','marker_assisted') not in ('strict','marker_assisted'):
        raise ValueError('parent_readability_policy must be strict or marker_assisted')
    if config.get('analysis_scope','regular') != 'regular':
        raise ValueError('Only regular analysis_scope is supported')
    ids = [config['father']['sample'],config['mother']['sample']]+[c['sample'] for c in config['children']]
    if len(ids) != len(set(ids)) or not config['children']:
        raise ValueError('Family sample identifiers must be unique and include children')
    if any('/' in s or '\\' in s or s in ('.','..') for s in ids):
        raise ValueError('Sample identifiers cannot contain path separators')
    with pysam.FastaFile(config['reference']) as ref:
        lengths = dict(zip(ref.references,ref.lengths))
    for member in [config['father'],config['mother']]+config['children']:
        with pysam.AlignmentFile(member['bam']) as bam:
            if not bam.has_index():
                raise ValueError(f'BAM index missing: {member["sample"]}')
            if any(lengths.get(c) != n for c,n in zip(bam.references,bam.lengths)):
                raise ValueError(f'BAM/reference contig mismatch: {member["sample"]}')
    for child in config['children']:
        candidate_input = child.get('candidate_vcf') or child.get('candidate_vcfs') or child.get('trgt_signal')
        if child.get('trgt_signal') and not config.get('trgt_directory'):
            raise ValueError('trgt_signal requires trgt_directory')
        if not child.get('local_directory') and not (candidate_input and config.get('trgt_directory')):
            raise ValueError('Supply local_directory or candidate_vcf and trgt_directory')
        if not child.get('c3_vcf') and not candidate_input:
            raise ValueError('C3 VCF or candidate VCF is required')


def run(manifest,output,workers=16,rate=.08,recompute_local=False):
    if not 0 < rate < 1 or workers < 1:
        raise ValueError('Invalid rate/workers')
    config = json.loads(Path(manifest).read_text())
    validate(config)
    out = Path(output)
    out.mkdir(parents=True,exist_ok=False)
    state(out,'PREPARING')
    try:
        local_dirs, jobs, c3 = {}, [], []
        for child in config['children']:
            sample = child['sample']
            directory = out/sample
            directory.mkdir()
            if child.get('candidate_vcfs') or child.get('trgt_signal'):
                state(out,'PREPARE_CANDIDATES',sample=sample)
                child['candidate_vcf'] = prepare_candidates(child,config,directory)
            if child.get('local_directory') and not recompute_local:
                local_dir = Path(child['local_directory'])
                paths = sorted(local_dir.glob('*.json'))
                if not paths:
                    raise ValueError(f'No cached evidence for {sample}')
                for path in paths:
                    data = json.loads(path.read_text())
                    expected = dict(C=sample,F=config['father']['sample'],M=config['mother']['sample'])
                    if any(data['family'][r]['sample'] != s for r,s in expected.items()):
                        raise ValueError(f'Cached family mismatch: {path}')
            else:
                local_dir = directory/'local'
                local_dir.mkdir()
                loci = {}
                for fields,info in read_records(child['candidate_vcf']):
                    if info.get('TRID'):
                        start,end = candidate_region(fields,info)
                        loci[info['TRID']] = (fields[0],start,end)
                for trid,(chrom,start,end) in loci.items():
                    jobs.append((config['trgt_directory'],sample,config['father']['sample'],config['mother']['sample'],
                                 trid,chrom,start,end,family_bams(config,child),str(local_dir/f'{trid}.json')))
            local_dirs[sample] = str(local_dir)
        state(out,'LOCAL_EVIDENCE',new_loci=len(jobs),workers=workers)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            for _ in pool.map(generate_local,jobs,chunksize=4):
                pass
        for child in config['children']:
            state(out,'C3_FILTER',sample=child['sample'])
            selected = run_c3(child,config,out/child['sample'],workers)
            c3.extend(c3_events(child['sample'],selected,Path(local_dirs[child['sample']]).glob('*.json')))
        pending_records = []
        for child in config['children']:
            path = out/child['sample']/'structural_pending.jsonl'
            pending_records.extend(dict(json.loads(line), sample=child['sample'])
                                   for line in path.read_text().splitlines())
        with (out/'structural_pending.jsonl').open('w') as target:
            for record in pending_records:
                target.write(json.dumps(record)+'\n')
        (out/'c3_audit.json').write_text(json.dumps(c3))
        fm = dict(samples=[dict(sample=s,local_directory=p) for s,p in local_dirs.items()],c3_audit=str(out/'c3_audit.json'))
        (out/'finalize_manifest.json').write_text(json.dumps(fm,indent=2))
        state(out,'FINALIZE')
        finalize(out/'finalize_manifest.json',out/'discovery',rates=(rate,),workers=workers,
                 minimum=config.get('minimum_variant_length',50))
        selected = [json.loads(line) for line in (out/'discovery'/f'rule_pass_{round(rate*100):02d}.jsonl').open()]
        grouped = {}
        family_id = config['father']['sample']+'__'+config['mother']['sample']
        for event in selected:
            event['family_id'] = family_id
            grouped.setdefault((event['sample'],event['trid']),[]).append(event)
        review_jobs = []
        children = {c['sample']:c for c in config['children']}
        for (sample,trid),events in grouped.items():
            review_jobs.append((str(Path(local_dirs[sample])/f'{trid}.json'),events,config['reference'],family_bams(config,children[sample]),rate))
        state(out,'RAW_BAM_REVIEW',events=len(selected),loci=len(grouped))
        reviewed = []
        with (out/'review_read_evidence.jsonl').open('w') as target:
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
                for events,raw in pool.map(review_job,review_jobs):
                    reviewed.extend(events)
                    target.write(json.dumps(dict(sample=events[0]['sample'],trid=events[0]['trid'],reads=raw))+'\n')
        assign_tiers(reviewed)
        state(out,'C3_ONLY_RAW_BAM_REVIEW')
        c3_only = review_c3_only(json.loads((out/'discovery'/'c3_evidence.json').read_text()),
                                 config,children,rate,workers,out)
        reviewed.extend(c3_only)
        for event in reviewed:
            event['variant_class'] = 'de_novo_SV' if abs(event['delta']) >= 50 else 'de_novo_indel'
        with (out/'strict_reviewed.jsonl').open('w') as target:
            for event in reviewed:
                target.write(json.dumps(event)+'\n')
        if config.get('parent_readability_policy','marker_assisted') == 'marker_assisted':
            from .trio_phase_review import apply_phase_fallback
            state(out,'MARKER_ASSISTED_PARENT_REVIEW')
            reviewed=[apply_phase_fallback(e,family_bams(config,children[e['sample']]),
                config.get('parent_haplotype_alignments',{}),rate) for e in reviewed]
            with (out/'marker_assisted_reviewed.jsonl').open('w') as target:
                for e in reviewed:
                    if 'marker_assisted_review' in e:
                        target.write(json.dumps(e)+'\n')
        outputs = dict(HIGH_CONFIDENCE_CANDIDATE='high_confidence_candidates.jsonl',REVIEW_REQUIRED='review_required.jsonl',
                       INSUFFICIENT_OR_CONFLICTING_EVIDENCE='insufficient_or_conflicting.jsonl')
        for tier,name in outputs.items():
            with (out/name).open('w') as target:
                for event in reviewed:
                    if event['tier']==tier:
                        target.write(json.dumps(event)+'\n')
        with (out/'all_reviewed.jsonl').open('w') as target:
            for event in reviewed:
                target.write(json.dumps(event)+'\n')
        summary = dict(structural_pending_records=len(pending_records),
                       structural_pending_types=dict(Counter(e['type'] for e in pending_records)),
                       rule_pass_events=len(selected),c3_only_reviewed=len(c3_only),
                       tiers=dict(Counter(e['tier'] for e in reviewed)),
                       per_sample={s:dict(Counter(e['tier'] for e in reviewed if e['sample']==s)) for s in children},
                       flags=dict(Counter(f for e in reviewed for f in e['audit_flags'])),
                       local_files=sum(len(list(Path(p).glob('*.json'))) for p in local_dirs.values()),
                       c3_representations=len(c3),workers=workers,error_rate=rate,
                       minimum_variant_length=config.get('minimum_variant_length',50),
                       confidence_policy=config.get('confidence_policy','strict'),
                       parent_readability_policy=config.get('parent_readability_policy','marker_assisted'),
                       interpretation='Experimental high-confidence candidates, not validated DNM or calibrated posterior probabilities')
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        report = ['# Trio Pipeline Result', '',
            'High-confidence here denotes an experimental candidate tier, not a validated de novo mutation.', '',
            '|Sample|High-confidence candidates|Review required|Insufficient/conflicting|',
            '|---|---:|---:|---:|']
        for sample,counts in summary['per_sample'].items():
            report.append(f"|{sample}|{counts.get('HIGH_CONFIDENCE_CANDIDATE',0)}|{counts.get('REVIEW_REQUIRED',0)}|{counts.get('INSUFFICIENT_OR_CONFLICTING_EVIDENCE',0)}|")
        report += ['', f'Structural records pending type-specific review: {len(pending_records)} (structural_pending.jsonl; record count, not deduplicated events).', 'All local candidates and C3-only evidence are retained under discovery/.',
                   'Parent-relative events are JSON records with full sequences, not reference-relative VCF.',
                   'Prepared caller/TRGT inputs are required; upstream alignment, assembly and callers are not launched by this entry point.',
                   'Nonautosomal ploidy, mosaic likelihood and independent-family validation remain incomplete.']
        (out/'REPORT.md').write_text('\n'.join(report)+'\n')
        files = {str(p):dict(size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in
                 [Path(config['reference'])]+[Path(m['bam']) for m in [config['father'],config['mother']]+config['children']]}
        (out/'provenance.json').write_text(json.dumps(dict(config=config,input_file_metadata=files,
            source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('trio_*.py')},
            recompute_local=recompute_local,truth_used=False),indent=2))
        state(out,'COMPLETE',**summary)
        return summary
    except Exception as error:
        state(out,'FAILED',error=str(error))
        raise
