"""Audit parental target sequences without treating absent split evidence as negative."""
from collections import defaultdict
import pysam
from .trio_structural_confidence import segment, endpoint_match


def evidence(bam,pair,support):
    names=set(support['alt_names'])
    segments=defaultdict(dict)
    observed=set()
    with pysam.AlignmentFile(bam) as source:
        for chrom,boundary,_ in pair:
            for r in source.fetch(chrom,max(0,boundary-250),boundary+250):
                if r.query_name not in names:continue
                observed.add(r.query_name)
                if r.is_secondary or r.is_duplicate or r.is_qcfail or r.mapping_quality<20:continue
                items=[(chrom,r.reference_start,'-' if r.is_reverse else '+',r.cigarstring,r.mapping_quality,'BAM')]
                if r.has_tag('SA'):
                    for value in r.get_tag('SA').rstrip(';').split(';'):
                        a=value.split(',')
                        if len(a)>=6:items.append((a[0],int(a[1])-1,a[2],a[3],int(a[4]),'SA'))
                for c,pos,strand,cigar,mq,origin in items:
                    s=segment(c,pos,strand,cigar,mq)
                    if s is None:continue
                    key=(c,pos,strand,cigar)
                    if key not in segments[r.query_name]:
                        segments[r.query_name][key]=dict(s,cigar=cigar,mapq=int(mq),origins=[])
                    if origin not in segments[r.query_name][key]['origins']:
                        segments[r.query_name][key]['origins'].append(origin)
    rows=[]
    for name in sorted(names):
        parts=list(segments[name].values());matches=[]
        for a in parts:
            for b in parts:
                x,y=endpoint_match(a,pair[0]),endpoint_match(b,pair[1])
                if x and y and a!=b and a['length']==b['length'] and x[1]!=y[1] and abs(x[0]-y[0])<=50:
                    matches.append(dict(first=a,second=b,query_gap=abs(x[0]-y[0]),
                        both_bam='BAM' in a['origins'] and 'BAM' in b['origins']))
        rows.append(dict(name=name,observed=name in observed,topology_pass=bool(matches),
            both_bam=any(m['both_bam'] for m in matches),segments=parts,matches=matches,
            reason='TOPOLOGY_SUPPORTED' if matches else 'NO_ELIGIBLE_SEGMENTS' if not parts else 'NO_MATCHING_ENDPOINT_PAIR'))
    supported={r['name'] for r in rows if r['topology_pass']}
    # Independently compare the diagnostic extractor with the frozen F09 API.
    directions={r['reverse'] for r in support['alt_reads'] if r['name'] in supported}
    bam_names={r['name'] for r in rows if r['both_bam']}
    bam_directions={r['reverse'] for r in support['alt_reads'] if r['name'] in bam_names}
    return dict(sequence_support=len(names),topology_support=len(supported),
        topology_names=sorted(supported),directions=sorted(directions),
        strong=len(supported)>=3 and len(directions)==2,
        both_bam_support=len(bam_names),both_bam_strong=len(bam_names)>=3 and len(bam_directions)==2,
        molecules=rows)


def review_parent_topology(event,bams):
    reviews=[]
    for junction in event['junction_evidence']:
        reviews.append(dict(endpoints=junction['endpoints'],parents={role:evidence(bams[role],junction['endpoints'],junction['support'][role]) for role in ('F','M')}))
    complete=len(reviews)==len(event['junctions'])
    strong=[role for role in ('F','M') if complete and reviews and all(j['parents'][role]['strong'] for j in reviews)]
    corroborated=[role for role in strong if all(j['parents'][role]['both_bam_strong'] for j in reviews)]
    any_topology=any(j['parents'][role]['topology_support'] for j in reviews for role in ('F','M'))
    has_sequence=any(j['parents'][role]['sequence_support'] for j in reviews for role in ('F','M'))
    classification=('INCOMPLETE_PARENT_TOPOLOGY_REVIEW' if not complete or not reviews else
        'NO_PARENT_TARGET_SEQUENCE_OBSERVED' if not has_sequence else
        'PARENT_STRUCTURAL_JUNCTIONS_SUPPORTED' if strong else
        'LIMITED_OR_PARTIAL_PARENT_TOPOLOGY' if any_topology else 'PARENT_SEQUENCE_ONLY_UNRESOLVED')
    return dict(event_id=event['event_id'],sample=event['sample'],type=event['type'],
        classification=classification,strong_parents=strong,bam_corroborated_parents=corroborated,
        junctions=reviews,previous_tier=event['tier'],previous_flags=event['audit_flags'],
        interpretation='Parent-supported adjacency; full haplotype transmission/copy number not independently certified' if strong else
        'Missing split topology does not establish parental absence or a false sequence match')

