"""Concrete C3 insertion review against four assembly-derived parent paths.

No truth coordinates are used. Target testing stays local; parent readability
can use larger windows without treating neighboring variants as the target.
"""
import re
from collections import Counter

import edlib
import pysam

from .trio_locus import query_interval


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def project_target(fields, position):
    """Project reference boundary into oriented assembly; reject gap interiors."""
    tags = {t[:2]: t[5:] for t in fields[12:]}
    cs = tags.get('cs', '')
    tokens = re.findall(r':\d+|=[A-Za-z]+|\*[A-Za-z]{2}|[+-][A-Za-z]+', cs)
    if not cs or ''.join(tokens) != cs:
        return None
    target, query = int(fields[7]), 0
    points = []
    for token in tokens:
        op = token[0]
        if op in ':=*':
            length = int(token[1:]) if op == ':' else (1 if op == '*' else len(token)-1)
            if target <= position <= target+length:
                points.append(query+position-target)
            target += length
            query += length
        elif op == '+':
            if target == position:
                points.extend((query, query+len(token)-1))
            query += len(token)-1
        else:
            length = len(token)-1
            if target < position < target+length:
                return None
            if position in (target, target+length):
                points.append(query)
            target += length
    # A gap exactly at a boundary has two valid projections, not one.
    if target != int(fields[8]) or query != int(fields[3])-int(fields[2]):
        return None
    return points[0] if points and len(set(points)) == 1 else None


def assembly_paths(haplotypes, chrom, left, right, anchor):
    paths, provenance = {}, {}
    for key, spec in haplotypes.items():
        matches = []
        with open(spec['paf']) as source:
            for line in source:
                f = line.rstrip().split('\t')
                if (f[5] != chrom or int(f[11]) < 20 or int(f[7]) > left
                        or int(f[8]) < right or 'tp:A:S' in f[12:]):
                    continue
                bounds = [project_target(f, p) for p in (left, right, anchor)]
                if any(p is None for p in bounds):
                    continue
                lo, hi, at = bounds
                with pysam.FastaFile(spec['fasta']) as fasta:
                    if f[4] == '+':
                        sequence = fasta.fetch(f[0], int(f[2])+lo, int(f[2])+hi).upper()
                    else:
                        sequence = reverse_complement(fasta.fetch(
                            f[0], int(f[3])-hi, int(f[3])-lo).upper())
                matches.append((sequence, at-lo, dict(contig=f[0], strand=f[4],
                    mapq=int(f[11]), fasta=spec['fasta'], paf=spec['paf'])))
        # Multiple alignment choices remain unresolved even if one looks good.
        if len(matches) != 1 or 'N' in matches[0][0]:
            return None, dict(unresolved=key, alignment_choices=len(matches))
        sequence, at, origin = matches[0]
        paths[key] = (sequence, at)
        provenance[key] = origin
    if set(paths) != {'F0', 'F1', 'M0', 'M1'}:
        return None, dict(unresolved='FOUR_PARENT_HAPLOTYPES_REQUIRED')
    return paths, provenance


def spanning_reads(bam, chrom, left, right):
    rows, seen = [], set()
    with pysam.AlignmentFile(bam) as source:
        for read in source.fetch(chrom, left, right):
            if (read.query_name in seen or read.mapping_quality < 20
                    or read.is_duplicate or read.is_qcfail):
                continue
            interval = query_interval(read, left, right, 0)
            if interval is None:
                continue
            sequence = read.get_forward_sequence()[interval[0]:interval[1]].upper()
            if read.is_reverse:
                sequence = reverse_complement(sequence)
            rows.append(dict(name=read.query_name, reverse=read.is_reverse,
                             sequence=sequence))
            seen.add(read.query_name)
    return rows


def allele_readability(paths, reads, rate):
    result = {}
    for role in ('F', 'M'):
        unique = {}
        for key, (sequence, _) in paths.items():
            if key.startswith(role):
                unique.setdefault(sequence, []).append(key)
        result[role] = []
        for sequence, keys in unique.items():
            names = []
            for read in reads[role]:
                distance = edlib.align(read['sequence'], sequence)['editDistance']
                others = [edlib.align(read['sequence'], other)['editDistance']
                          for other in unique if other != sequence]
                if (distance/max(1,len(read['sequence'])) <= rate
                        and (not others or min(others)-distance >= 5)):
                    names.append(read['name'])
            result[role].append(dict(models=keys, discriminating_reads=len(names), names=names))
    return result


def review_insertion(job):
    call, reference, bams, haplotypes, rate, family_id = job
    chrom, pos0 = call['chrom'], call['pos0']
    ref, alt = call['ref'].upper(), call['alt'].upper()
    length = len(alt)-len(ref)
    result = dict(event_id=f"{call['sample']}.C3INS.{chrom}.{pos0}.{length}",
        sample=call['sample'], family_id=family_id, chrom=chrom, start=pos0+len(ref),
        end=pos0+len(ref), type='INS', delta=length, basis='REFERENCE',
        variant_class='de_novo_indel' if length < 50 else 'de_novo_SV',
        branches=['C3_ASSEMBLY_BACKGROUND_INS'], c3_event_ids=call.get('c3_event_ids', []),
        c3_record_ids=call.get('record_ids', []), audit_flags=[])
    flags = result['audit_flags']
    if (not set(ref+alt) <= set('ACGT') or not alt.startswith(ref) or length < 30
            or length != call['length']):
        flags.append('INSERTION_REPRESENTATION_UNRESOLVED')
        result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
        return result
    anchor = pos0+len(ref)
    left, right = max(0,anchor-200), anchor+200
    with pysam.FastaFile(reference) as fasta:
        if fasta.fetch(chrom,pos0,anchor).upper() != ref:
            flags.append('REFERENCE_ANCHOR_MISMATCH')
            result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
            return result
        refseq = fasta.fetch(chrom,left,right).upper()
    paths, provenance = assembly_paths(haplotypes,chrom,left,right,anchor)
    result['assembly_projection'] = provenance
    if paths is None:
        flags.append('PARENT_ASSEMBLY_PATH_UNRESOLVED')
        result['tier'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
        return result
    inserted = alt[len(ref):]
    target = refseq[:anchor-left]+inserted+refseq[anchor-left:]
    # Test the target on every parental background as well as the reference.
    alternatives = [target]+[seq[:at]+inserted+seq[at:] for seq,at in paths.values()]
    parental = [seq for seq,_ in paths.values()]
    if any(edlib.align(seq,refseq)['editDistance']-edlib.align(seq,target)['editDistance'] >= 10
           for seq in parental):
        flags.append('PARENT_ASSEMBLY_CONTAINS_TARGET')
    result['reference_allele'], result['child_sequence'] = refseq, target
    result['parent_paths'] = {k:s for k,(s,_) in paths.items()}
    raw = {role:spanning_reads(bam,chrom,left,right) for role,bam in bams.items()}
    support, details = {}, {}
    for role, reads in raw.items():
        details[role] = []
        for read in reads:
            da = min(edlib.align(read['sequence'],s)['editDistance'] for s in alternatives)
            dp = min(edlib.align(read['sequence'],s)['editDistance'] for s in parental)
            details[role].append(dict(name=read['name'],reverse=read['reverse'],
                alt_distance=da,parent_distance=dp,
                specific=dp-da >= 10 and da/max(1,len(read['sequence'])) <= rate,
                target=da <= dp and da/max(1,len(read['sequence'])) <= rate))
        specific = [r for r in details[role] if r['specific']]
        support[role] = dict(depth=len(reads),specific=len(specific),
            target=sum(r['target'] for r in details[role]),
            specific_names=[r['name'] for r in specific],
            strands=dict(Counter('reverse' if r['reverse'] else 'forward' for r in specific)))
    result['flank_recheck'], result['read_evidence'] = support, details
    attempts = []
    for flank in (200,500,1000,2000,5000):
        pl, pr = max(0,anchor-flank), anchor+flank
        ppaths, origin = (paths,provenance) if flank == 200 else assembly_paths(haplotypes,chrom,pl,pr,anchor)
        if ppaths is None:
            attempts.append(dict(flank=flank,unresolved=origin))
            continue
        preads = raw if flank == 200 else {r:spanning_reads(bams[r],chrom,pl,pr) for r in ('F','M')}
        coverage = allele_readability(ppaths,preads,rate)
        attempt = dict(flank=flank,coverage=coverage,depth={r:len(preads[r]) for r in ('F','M')})
        attempts.append(attempt)
        if (min(attempt['depth'].values()) >= 10
                and all(a['discriminating_reads'] >= 3 for rows in coverage.values() for a in rows)):
            result['parent_allele_coverage'] = coverage
            result['parent_readability_flank'] = flank
            break
    result['parent_readability_attempts'] = attempts
    if 'parent_allele_coverage' not in result:
        flags.append('PARENT_ALLELE_COVERAGE_UNRESOLVED')
    if support['C']['specific'] < 3:
        flags.append('CHILD_SUPPORT_NOT_REPRODUCED_WITH_FLANKS')
    if len(support['C']['strands']) < 2:
        flags.append('SINGLE_STRAND_SUPPORT')
    if min(support[r]['depth'] for r in ('F','M')) < 10:
        flags.append('PARENT_FLANK_DEPTH_LT_MIN')
    if any(support[r]['target'] for r in ('F','M')):
        flags.append('PARENT_TARGET_AFTER_FLANK_RECHECK')
    if chrom.removeprefix('chr') not in {str(i) for i in range(1,23)}:
        flags.append('PLOIDY_NOT_VALIDATED')
    result['tier'] = ('HIGH_CONFIDENCE_CANDIDATE' if not flags else
        'INSUFFICIENT_OR_CONFLICTING_EVIDENCE' if any(f.startswith('PARENT_') or f == 'PLOIDY_NOT_VALIDATED'
                                                    for f in flags) else 'REVIEW_REQUIRED')
    return result
