"""Sample-independent raw-BAM checks for sequence-based trio candidates."""
from collections import Counter

import edlib
import pysam

from .trio_locus import query_interval


def parent_allele_coverage(local, rate):
    result = {}
    for role in ('F', 'M'):
        unique = {}
        for key, seq in local['models'].items():
            if key.startswith(role):
                unique.setdefault(seq, []).append(key)
        result[role] = []
        for seq, keys in unique.items():
            key = keys[0]
            other = [k for k, s in local['models'].items() if k.startswith(role) and s != seq]
            reads = [r for r in local['full_reads'][role] if r['distances'][key]/max(1,r['length']) <= rate
                     and (not other or min(r['distances'][k] for k in other)-r['distances'][key] >= 5)]
            result[role].append(dict(models=keys, discriminating_reads=len(reads)))
    return result


def review(local, events, reference, bams, rate=.08):
    event = events[0]
    with pysam.FastaFile(reference) as ref:
        left = ref.fetch(event['chrom'], max(0,event['start']-200), event['start']).upper()
        right = ref.fetch(event['chrom'], event['end'], event['end']+200).upper()
    flank_length = len(left)+len(right)
    models = {k: left+s.upper()+right for k,s in local['models'].items()}
    raw = {}
    for role in ('C','F','M'):
        rows, seen = [], set()
        with pysam.AlignmentFile(bams[role]) as source:
            for read in source.fetch(event['chrom'], max(0,event['start']-250), event['end']+250):
                if read.query_name in seen or read.mapping_quality < 20:
                    continue
                start, end = max(0,event['start']-200), event['end']+200
                if query_interval(read,start,end,50) is None:
                    continue
                iv = query_interval(read,start,end,0)
                sequence = read.get_forward_sequence()[iv[0]:iv[1]]
                if read.is_reverse:
                    sequence = sequence.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
                distances = {k:edlib.align(sequence.upper(),s)['editDistance'] for k,s in models.items()}
                rows.append(dict(name=read.query_name, length=len(sequence), distances=distances,
                                 reverse=read.is_reverse, mapq=read.mapping_quality,
                                 hp=read.get_tag('HP') if read.has_tag('HP') else None,
                                 ps=read.get_tag('PS') if read.has_tag('PS') else None))
                seen.add(read.query_name)
        raw[role] = rows
    coverage = parent_allele_coverage(local,rate)
    output = []
    for event in events:
        key = event['child_model']
        supports = {}
        for role, reads in raw.items():
            target = [r for r in reads if r['distances'][key] == min(r['distances'].values())
                      and r['distances'][key]/max(1,r['length']-flank_length) <= rate]
            specific = [r for r in target if min(v for k,v in r['distances'].items() if k.startswith(('F','M')))-r['distances'][key] >= 10]
            supports[role] = dict(depth=len(reads),target=len(target),specific=len(specific),
                target_names=[r['name'] for r in target], specific_names=[r['name'] for r in specific],
                strands=dict(Counter('reverse' if r['reverse'] else 'forward' for r in specific)))
        flags = []
        if supports['F']['target'] or supports['M']['target']:
            flags.append('PARENT_TARGET_AFTER_FLANK_RECHECK')
        if supports['C']['specific'] < 3:
            flags.append('CHILD_SUPPORT_NOT_REPRODUCED_WITH_FLANKS')
        if min(supports['F']['depth'],supports['M']['depth']) < 10:
            flags.append('PARENT_FLANK_DEPTH_LT_MIN')
        if any(a['discriminating_reads'] < 3 for alleles in coverage.values() for a in alleles):
            flags.append('PARENT_ALLELE_COVERAGE_UNRESOLVED')
        if len(supports['C']['strands']) < 2:
            flags.append('SINGLE_STRAND_SUPPORT')
        if event['origin_ambiguous']:
            flags.append('ORIGIN_AMBIGUOUS')
        if event['chrom'].removeprefix('chr') not in {str(n) for n in range(1,23)}:
            flags.append('PLOIDY_NOT_VALIDATED')
        output.append(dict(event, flank_recheck=supports,parent_allele_coverage=coverage,audit_flags=flags))
    return output,raw


def assign_tiers(events):
    for event in events:
        event['same_sequence_siblings'] = sorted({other['sample'] for other in events
            if other['sample'] != event['sample'] and other['trid'] == event['trid']
            and other['child_sequence'] == event['child_sequence']
            and other.get('family_id') == event.get('family_id')})
        if event['same_sequence_siblings']:
            event['audit_flags'].append('EXACT_TARGET_ALLELE_SHARED_WITH_SIBLING')
        event['audit_flags'] = sorted(set(event['audit_flags']))
        flags = set(event['audit_flags'])
        if not flags:
            tier = 'HIGH_CONFIDENCE_CANDIDATE'
        elif flags & {'PARENT_TARGET_AFTER_FLANK_RECHECK','PARENT_ALLELE_COVERAGE_UNRESOLVED',
                      'PARENT_FLANK_DEPTH_LT_MIN','PLOIDY_NOT_VALIDATED'}:
            tier = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
        else:
            tier = 'REVIEW_REQUIRED'
        event['tier'] = tier
    return events
