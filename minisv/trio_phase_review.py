"""Experimental marker-assisted review of already discovered INS candidates.

HP labels alone do not establish readability. Both parental groups require a
linked high-quality heterozygous marker and separately reproduced sequences.
No event coordinates or truth labels are encoded in the decision rules.
"""
import re
from collections import Counter, defaultdict

import edlib
import pysam

from .trio_insertion import spanning_reads, assembly_paths


def consensus(sequences):
    """Reference-free medoid followed by majority aligned base/gap voting."""
    if not sequences:
        return None
    guide = min(sorted(set(sequences)),key=lambda s:sum(
        edlib.align(s,t)['editDistance'] for t in sequences))
    bases = [Counter() for _ in guide]
    insertions = [Counter() for _ in range(len(guide)+1)]
    for sequence in sequences:
        alignment = edlib.align(sequence,guide,task='path')
        q,t = 0,0
        ins = defaultdict(str)
        for size,op in re.findall(r'(\d+)([=XID])',alignment['cigar']):
            size = int(size)
            if op in '=X':
                for offset in range(size):bases[t+offset][sequence[q+offset]] += 1
                q += size;t += size
            elif op == 'I':
                ins[t] += sequence[q:q+size];q += size
            else:
                for offset in range(size):bases[t+offset][''] += 1
                t += size
        for boundary in range(len(insertions)):
            insertions[boundary][ins[boundary]] += 1
    output=[];ambiguous=[]
    for i in range(len(guide)+1):
        choice,count=insertions[i].most_common(1)[0]
        if count/len(sequences) < .8:ambiguous.append(['insertion',i])
        output.append(choice)
        if i<len(guide):
            choice,count=bases[i].most_common(1)[0]
            if count/len(sequences)<.8:ambiguous.append(['base',i])
            output.append(choice)
    return dict(sequence=''.join(output),ambiguous_columns=ambiguous,depth=len(sequences))


def phase_groups(bam,chrom,left,right,search=15000):
    raw={r['name']:r for r in spanning_reads(bam,chrom,left,right)}
    groups=defaultdict(list);sites=defaultdict(lambda:defaultdict(list));seen=set()
    with pysam.AlignmentFile(bam) as source:
        for read in source.fetch(chrom,left,right):
            if (read.query_name not in raw or read.query_name in seen or read.is_secondary
                    or read.is_supplementary):continue
            seen.add(read.query_name)
            hp=read.get_tag('HP') if read.has_tag('HP') else None
            ps=read.get_tag('PS') if read.has_tag('PS') else None
            row=dict(raw[read.query_name],hp=hp,ps=ps)
            groups[(ps,hp)].append(row)
            if hp not in (1,2) or ps is None:continue
            for q,r in read.get_aligned_pairs(matches_only=True):
                if r<left-search or r>right+search:continue
                if read.query_qualities is None or read.query_qualities[q]<20:continue
                base=read.query_sequence[q].upper()
                if base in 'ACGT':sites[(ps,r)][hp].append((read.query_name,base))
    markers=[]
    for (ps,pos),hps in sites.items():
        votes={hp:Counter(b for _,b in hps[hp]) for hp in (1,2)}
        if any(sum(v.values())<3 for v in votes.values()):continue
        tops={hp:vs.most_common(1)[0] for hp,vs in votes.items()}
        if tops[1][0]==tops[2][0]:continue
        if any(count/sum(votes[hp].values())<.9 for hp,(_,count) in tops.items()):continue
        markers.append(dict(ps=ps,pos0=pos,alleles={str(hp):dict(base=base,count=count,
            depth=sum(votes[hp].values()),names=[n for n,b in hps[hp] if b==base])
            for hp,(base,count) in tops.items()}))
    return raw,groups,markers


def review_phase(event,bams,rate=.08):
    """Only resolves parent readability; every other prior flag remains active."""
    chrom=event['chrom'];anchor=event['start'];left=max(0,anchor-200);right=anchor+200
    result=dict(event_id=event['event_id'],original_tier=event['tier'],
        policy='STRICT_WITH_MARKER_ASSISTED_PARENT_READABILITY_V1',parents={})
    target=event['child_sequence']
    roles={}
    for role in ('F','M'):
        parent_templates=[s for k,s in event['parent_paths'].items() if k.startswith(role)]
        raw,groups,markers=phase_groups(bams[role],chrom,left,right)
        ps_candidates=sorted({m['ps'] for m in markers})
        attempts=[]
        for ps in ps_candidates:
            reports=[]
            for hp in (1,2):
                rows=groups.get((ps,hp),[])
                # Require direct molecule linkage, rather than HP tags alone.
                linked_markers=[m for m in markers if m['ps']==ps]
                marker=min(linked_markers,key=lambda m:abs(m['pos0']-anchor))
                linked_names=set(marker['alleles'][str(hp)]['names'])
                linked_rows=[r for r in rows if r['name'] in linked_names]
                flags=[]
                if len(rows)<3:flags.append('PHASE_GROUP_DEPTH_LT_3')
                if len({r['reverse'] for r in rows})<2:flags.append('PHASE_GROUP_SINGLE_STRAND')
                if len(linked_rows)<3 or len({r['reverse'] for r in linked_rows})<2:
                    flags.append('MARKER_LINKED_MOLECULES_INSUFFICIENT')
                c=consensus([r['sequence'] for r in rows])
                linked_consensus=consensus([r['sequence'] for r in linked_rows])
                strand_consensus={str(strand):consensus([r['sequence'] for r in rows if r['reverse']==strand])
                                  for strand in (False,True)}
                if c is None:flags.append('NO_CONSENSUS')
                else:
                    if (linked_consensus is None or linked_consensus['ambiguous_columns']
                            or linked_consensus['sequence']!=c['sequence']):
                        flags.append('MARKER_LINKED_CONSENSUS_DISCORDANT')
                    if c['ambiguous_columns']:flags.append('CONSENSUS_AMBIGUOUS')
                    if any(s is None or s['sequence']!=c['sequence'] for s in strand_consensus.values()):
                        flags.append('STRAND_CONSENSUS_DISCORDANT')
                    dp=min(edlib.align(c['sequence'],s)['editDistance'] for s in parent_templates)
                    dt=edlib.align(c['sequence'],target)['editDistance']
                    c.update(parent_distance=dp,target_distance=dt)
                    if dp/max(1,len(c['sequence']))>rate or dt-dp<10:
                        flags.append('PARENT_CONSENSUS_NOT_TARGET_NEGATIVE')
                reports.append(dict(hp=hp,depth=len(rows),names=[r['name'] for r in rows],
                    marker_pos0=marker['pos0'],linked_depth=len(linked_rows),
                    linked_names=[r['name'] for r in linked_rows],
                    linked_consensus=linked_consensus,
                    flags=flags,consensus=c,strand_consensus=strand_consensus))
            attempts.append(dict(ps=ps,groups=reports,passed=all(not r['flags'] for r in reports)))
        passed=any(a['passed'] for a in attempts)
        roles[role]=passed
        result['parents'][role]=dict(depth=len(raw),unphased_or_other=len(raw)-sum(
            len(rs) for (ps,hp),rs in groups.items() if ps is not None and hp in (1,2)),
            markers=markers,attempts=attempts,passed=passed)
    flags=[f for f in event['audit_flags'] if f!='PARENT_ALLELE_COVERAGE_UNRESOLVED']
    if not all(roles.values()):flags.append('PHASE_PARENT_READABILITY_UNRESOLVED')
    result['audit_flags']=flags
    result['tier']='HIGH_CONFIDENCE_CANDIDATE' if not flags else 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    result['original_strict_tier']=event['tier']
    result['interpretation']='Experimental alternative readability audit; not independent biological validation'
    return result


def apply_phase_fallback(event,bams,haplotypes,rate=.08):
    """Resolve only eligible concrete-INS readability failures; retain provenance."""
    if (event.get('type')!='INS' or event.get('basis')!='REFERENCE'
            or not event.get('parent_paths')
            or 'PARENT_ALLELE_COVERAGE_UNRESOLVED' not in event.get('audit_flags',[])):
        return event
    result=dict(event)
    audit=review_phase(event,bams,rate)
    anchor=event['start'];left=max(0,anchor-200);right=anchor+200
    names=set(event['flank_recheck']['C']['specific_names'])
    rows=[r for r in spanning_reads(bams['C'],event['chrom'],left,right) if r['name'] in names]
    c=consensus([r['sequence'] for r in rows])
    strands={str(s):consensus([r['sequence'] for r in rows if r['reverse']==s]) for s in (False,True)}
    paths,_=assembly_paths(haplotypes,event['chrom'],left,right,anchor)
    offset=anchor-left
    inserted=event['child_sequence'][offset:offset+event['delta']]
    models={'reference_background':event['child_sequence']}
    if paths:
        models.update({k:s[:at]+inserted+s[at:] for k,(s,at) in paths.items()})
    matched=[k for k,s in models.items() if c and c['sequence']==s]
    passed=bool(c and not c['ambiguous_columns'] and matched and all(
        s and not s['ambiguous_columns'] and s['sequence']==c['sequence'] for s in strands.values()))
    audit['child_consensus']=dict(passed=passed,consensus=c,strands=strands,
        matched_backgrounds=matched,names=[r['name'] for r in rows])
    result['marker_assisted_review']=audit
    result['original_strict_tier']=event['tier']
    result['original_strict_flags']=list(event['audit_flags'])
    if audit['tier']=='HIGH_CONFIDENCE_CANDIDATE' and passed:
        result['audit_flags']=audit['audit_flags']
        result['tier']=audit['tier']
        result['confidence_basis']=audit['policy']
    return result
