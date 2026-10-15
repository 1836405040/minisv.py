"""Conservative candidate-tier gates; not a calibrated DNM probability model."""
import re
from collections import defaultdict
from statistics import median

import pysam


def segment(chrom, pos, strand, cigar, mapq):
    parts = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
    if not parts or int(mapq) < 20:
        return None
    length = sum(n for n, op in parts if op in 'MISH=X')
    leading = 0
    for n, op in parts:
        if op not in 'SH': break
        leading += n
    aligned = sum(n for n, op in parts if op in 'MI=X')
    reference = sum(n for n, op in parts if op in 'MDN=X')
    if sum(n for n,op in parts if op in 'M=X') < 100:
        return None
    qstart, qend = leading, leading+aligned
    if strand == '-':
        qstart, qend = length-qend, length-qstart
    return dict(chrom=chrom,start=int(pos),end=int(pos)+reference,
                strand=strand,qstart=qstart,qend=qend,length=length)


def endpoint_match(s, endpoint):
    chrom, boundary, side = endpoint
    edge = s['end'] if side == 'L' else s['start']
    if s['chrom'] != chrom or abs(edge-boundary) > 50:
        return None
    toward = (side == 'L') == (s['strand'] == '+')
    return (s['qend'] if toward else s['qstart'], toward)


def split_support(bam, pair, allowed_names):
    segments = defaultdict(list)
    allowed_names = set(allowed_names)
    with pysam.AlignmentFile(bam) as source:
        for chrom,boundary,_ in pair:
            for read in source.fetch(chrom,max(0,boundary-250),boundary+250):
                if read.query_name not in allowed_names or read.is_secondary or read.is_duplicate or read.is_qcfail or read.mapping_quality < 20:
                    continue
                s = segment(chrom,read.reference_start,'-' if read.is_reverse else '+',read.cigarstring,read.mapping_quality)
                if s: segments[read.query_name].append(s)
                if read.has_tag('SA'):
                    for item in read.get_tag('SA').rstrip(';').split(';'):
                        a = item.split(',')
                        if len(a) < 6: continue
                        s = segment(a[0],int(a[1])-1,a[2],a[3],a[4])
                        if s: segments[read.query_name].append(s)
    supported = []
    for name, parts in segments.items():
        found = False
        for a in parts:
            for b in parts:
                x,y = endpoint_match(a,pair[0]),endpoint_match(b,pair[1])
                if x and y and a != b and a['length']==b['length'] and x[1]!=y[1] and abs(x[0]-y[0]) <= 50:
                    found = True
        if found: supported.append(name)
    return sorted(supported)


def parent_haplotypes(bam, endpoint, allowed_names):
    chrom,boundary,_ = endpoint
    allowed_names = set(allowed_names)
    groups = defaultdict(lambda: {1:set(),2:set()})
    with pysam.AlignmentFile(bam) as source:
        for read in source.fetch(chrom,max(0,boundary-200),boundary+200):
            if (read.query_name not in allowed_names or read.is_secondary or read.is_supplementary
                or read.is_duplicate or read.is_qcfail or read.mapping_quality < 20
                or not read.has_tag('HP') or not read.has_tag('PS')): continue
            # Require the reference allele to be aligned continuously through
            # the junction; a reference motif elsewhere in the read is insufficient.
            if not any(s <= boundary-200 and e >= boundary+200 for s,e in read.get_blocks()): continue
            hp = read.get_tag('HP')
            if hp in (1,2): groups[str(read.get_tag('PS'))][hp].add(read.query_name)
    phases = {ps:{str(h):sorted(names) for h,names in data.items()} for ps,data in groups.items()}
    passed = any(all(len(data[h]) >= 3 for h in (1,2)) for data in groups.values())
    return dict(passed=passed,phase_sets=phases)


def window_depth(bam, chrom, start, end):
    if start < 0: return None
    total = 0
    with pysam.AlignmentFile(bam) as source:
        if chrom not in source.references or end > source.get_reference_length(chrom): return None
        for r in source.fetch(chrom,start,end):
            if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality < 30: continue
            total += sum(max(0,min(e,end)-max(s,start)) for s,e in r.get_blocks())
    return total/(end-start)


def duplication_dosage(event,bams):
    endpoints = event['junctions'][0]
    chrom = endpoints[0][0]
    start,end = sorted(ep[1] for ep in endpoints)
    if end-start < 1000:
        return dict(passed=False,reason='INTERVAL_TOO_SHORT_FOR_WINDOW_DOSAGE')
    inside = [(start+int((end-start)*f)-100,start+int((end-start)*f)+100) for f in (.1,.3,.5,.7,.9)]
    outside = [(start-d-200,start-d) for d in (500,1000,1500)]+[(end+d,end+d+200) for d in (500,1000,1500)]
    rows = {}
    for role,bam in bams.items():
        a=[window_depth(bam,chrom,s,e) for s,e in inside]
        b=[window_depth(bam,chrom,s,e) for s,e in outside]
        valid=all(v is not None and v >= 10 for v in a+b)
        stable=valid and max(b)/min(b) <= 1.6 and max(a)/min(a) <= 1.6
        ratio=median(a)/median(b) if valid else None
        passed=bool(stable and ((1.25 <= ratio <= 1.8) if role=='C' else (.8 <= ratio <= 1.2)))
        rows[role]=dict(inside=a,flanks=b,ratio=ratio,passed=passed)
    return dict(passed=all(r['passed'] for r in rows.values()),samples=rows,inside_windows=inside,flank_windows=outside)


def reciprocal_mates(event):
    rows = {}
    for r in event['records']:
        fields=r['original_record'].split('\t')
        info=dict(x.split('=',1) for x in fields[7].split(';') if '=' in x)
        key=(r.get('source'),info.get('MINISV_ORIGINAL_ID',r['id']))
        rows[key]=set(info.get('MATEID','').split(','))-{''}
    return bool(rows) and all(any((source,mate) in rows and ident in rows[(source,mate)] and mate!=ident
        for mate in mates) for (source,ident),mates in rows.items())


def assess(event,bams,policy='audit'):
    if policy not in ('audit','strict_v1'): raise ValueError('Unknown structural_confidence_policy')
    flags=set(event['audit_flags'])-{'STRUCTURAL_MODEL_VALIDATION_REQUIRED','DUP_COPY_NUMBER_AND_DESTINATION_UNRESOLVED'}
    audit=dict(policy=policy,junctions=[])
    if len(event['junction_evidence']) != len(event['junctions']): flags.add('INCOMPLETE_JUNCTION_REVIEW')
    for evidence in event['junction_evidence']:
        pair=evidence['endpoints'];support=evidence['support']
        names=split_support(bams['C'],pair,support['C']['alt_names'])
        directions={r['reverse'] for r in support['C']['alt_reads'] if r['name'] in names}
        if len(names)<3 or len(directions)<2: flags.add('CHILD_SPLIT_TOPOLOGY_UNRESOLVED')
        parents={r:[parent_haplotypes(bams[r],ep,support[r]['reference_names'][i])
                   for i,ep in enumerate(pair)] for r in ('F','M')}
        if any(not p['passed'] for data in parents.values() for p in data): flags.add('PARENT_PHASED_CALLABILITY_UNRESOLVED')
        audit['junctions'].append(dict(split_names=names,parent_haplotypes=parents))
    if event['type']=='DUP':
        audit['dosage']=duplication_dosage(event,bams)
        if not audit['dosage']['passed']: flags.add('DUP_DOSAGE_UNRESOLVED')
    if event['type']=='BND':
        audit['reciprocal_mates']=reciprocal_mates(event)
        if not audit['reciprocal_mates']: flags.add('BND_MATE_LINK_UNRESOLVED')
    audit['eligible_high_confidence']=not flags
    result=dict(event,structural_confidence=audit)
    if policy=='audit': flags.add('STRUCTURAL_MODEL_VALIDATION_REQUIRED')
    result['audit_flags']=sorted(flags)
    result['tier']=('HIGH_CONFIDENCE_CANDIDATE' if not flags else
        'REVIEW_REQUIRED' if flags=={'STRUCTURAL_MODEL_VALIDATION_REQUIRED'} else
        'INSUFFICIENT_OR_CONFLICTING_EVIDENCE')
    return result
