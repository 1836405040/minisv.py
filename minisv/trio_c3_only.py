"""Reference-allele review of C3 calls without a compatible catalog event."""

from collections import Counter, defaultdict

import edlib
import pysam

from .trio_locus import query_interval


def group_calls(rows, minimum=50):
    by_chrom = defaultdict(list)
    for row in rows:
        if row['type'] != 'DEL' or row['length'] < minimum or row['basis'] != 'REFERENCE':
            continue
        by_chrom[(row['sample'], row['chrom'])].append(row)
    groups = {}
    for (sample, chrom), calls in by_chrom.items():
        clusters = []
        for call in sorted(calls, key=lambda x: (x['pos0'], -x['length'])):
            start, end = call['pos0'] + 1, call['pos0'] + 1 + call['length']
            match = next((cluster for cluster in clusters
                          if abs(start - cluster[0]['pos0'] - 1) <= 20
                          and abs(end - cluster[0]['pos0'] - 1 - cluster[0]['length']) <= 20
                          and min(call['length'], cluster[0]['length']) /
                          max(call['length'], cluster[0]['length']) >= .8), None)
            if match is None:
                clusters.append([call])
            else:
                match.append(call)
        for cluster in clusters:
            if not any(call.get('compatible_local_event_ids') for call in cluster):
                representative = cluster[0]
                groups[(sample, chrom, representative['pos0'], representative['length'])] = cluster
    return groups


def review_deletion(job):
    calls, reference, bams, rate, family_id = job
    first = calls[0]
    chrom, pos0, length = first['chrom'], first['pos0'], first['length']
    start = pos0 + 1
    end = start + length
    left = max(0, start - 200)
    right = end + 200
    record_ids = sorted({rid for call in calls for rid in call['record_ids']})
    result = dict(event_id=f"{first['sample']}.C3ONLY.{chrom}.{start}.{end}.DEL",
                  sample=first['sample'], family_id=family_id, chrom=chrom,
                  start=start, end=end, type='DEL', delta=-length,
                  basis='REFERENCE', branches=['C3_ONLY_REFERENCE_ALLELE'],
                  c3_event_ids=[call['event_id'] for call in calls],
                  c3_record_ids=record_ids, audit_flags=[])
    with pysam.FastaFile(reference) as fasta:
        if chrom not in fasta.references or right > fasta.get_reference_length(chrom):
            result['audit_flags'].append('REFERENCE_INTERVAL_UNAVAILABLE')
            result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
            return result
        ref = fasta.fetch(chrom, left, right).upper()
    alt = ref[:start-left] + ref[end-left:]
    result['reference_allele'] = ref
    result['child_sequence'] = alt
    result['trid'] = f'C3ONLY_{chrom}_{start}_{end}'
    if len(ref) - len(alt) != length or 'N' in ref:
        result['audit_flags'].append('REFERENCE_ALLELE_UNRESOLVED')
        result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
        return result
    supports = {}
    for role, bam in bams.items():
        seen, reads = set(), []
        with pysam.AlignmentFile(bam) as source:
            for read in source.fetch(chrom, left, right):
                if (read.query_name in seen or read.mapping_quality < 20 or read.is_duplicate
                        or read.is_qcfail):
                    continue
                interval = query_interval(read, left, right, 0)
                if interval is None or read.reference_start > left or read.reference_end < right:
                    continue
                seq = read.get_forward_sequence()[interval[0]:interval[1]].upper()
                if read.is_reverse:
                    seq = seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
                d_alt = edlib.align(seq, alt)['editDistance']
                d_ref = edlib.align(seq, ref)['editDistance']
                reads.append(dict(name=read.query_name, reverse=read.is_reverse,
                                  alt_distance=d_alt, ref_distance=d_ref,
                                  alt_specific=d_ref-d_alt >= 10 and d_alt/max(1,len(seq)) <= rate,
                                  ref_specific=d_alt-d_ref >= 10 and d_ref/max(1,len(seq)) <= rate))
                seen.add(read.query_name)
        alt_reads = [r for r in reads if r['alt_specific']]
        supports[role] = dict(depth=len(reads), alt_specific=len(alt_reads),
                              ref_specific=sum(r['ref_specific'] for r in reads),
                              alt_names=[r['name'] for r in alt_reads],
                              alt_strands=dict(Counter('reverse' if r['reverse'] else 'forward'
                                                       for r in alt_reads)))
    result['flank_recheck'] = supports
    flags = result['audit_flags']
    if supports['C']['alt_specific'] < 3:
        flags.append('CHILD_SUPPORT_NOT_REPRODUCED_WITH_FLANKS')
    if min(supports[r]['depth'] for r in ('F','M')) < 10:
        flags.append('PARENT_FLANK_DEPTH_LT_MIN')
    if any(supports[r]['alt_specific'] for r in ('F','M')):
        flags.append('PARENT_TARGET_AFTER_FLANK_RECHECK')
    if any(supports[r]['ref_specific'] < 3 for r in ('F','M')):
        flags.append('PARENT_REFERENCE_ALLELE_UNRESOLVED')
    if len(supports['C']['alt_strands']) < 2:
        flags.append('SINGLE_STRAND_SUPPORT')
    if chrom.removeprefix('chr') not in {str(i) for i in range(1,23)}:
        flags.append('PLOIDY_NOT_VALIDATED')
    result['audit_flags'] = sorted(set(flags))
    if not flags:
        result['tier'] = 'HIGH_CONFIDENCE_CANDIDATE'
    elif set(flags) & {'PARENT_FLANK_DEPTH_LT_MIN','PARENT_TARGET_AFTER_FLANK_RECHECK',
                       'PARENT_REFERENCE_ALLELE_UNRESOLVED','PLOIDY_NOT_VALIDATED'}:
        result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    else:
        result['tier'] = 'REVIEW_REQUIRED'
    return result
