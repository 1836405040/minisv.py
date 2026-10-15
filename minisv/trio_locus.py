"""Associate parental assembly evidence with a candidate's interval on a read."""
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import pysam

from .filtercaller import parse_msvasm_withtype
from .type import simple_type


def query_interval(read, start, end, flank=100):
    """Project a fully spanned reference interval to forward molecule coordinates.

    Hard-clipped alignments are not comparable with FASTQ-derived GSV offsets.
    Insertions at interval boundaries are included by walking the CIGAR path.
    """
    left, right = max(0, start - flank), end + flank
    if (read.is_unmapped or read.is_secondary or read.is_supplementary
            or not read.cigartuples or any(op == 5 for op, _ in read.cigartuples)
            or read.reference_start > left or read.reference_end < right):
        return None
    ref, query, points = read.reference_start, 0, []
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            lo, hi = max(left, ref), min(right, ref + length)
            if lo <= hi:
                points.extend((query + lo - ref, query + hi - ref))
            ref += length
            query += length
        elif op == 1:
            if left <= ref <= right:
                points.extend((query, query + length))
            query += length
        elif op in (2, 3):
            if ref <= right and ref + length >= left:
                if op == 3:
                    return None
                points.append(query)
            ref += length
        elif op == 4:
            query += length
    if not points:
        return None
    lo, hi = min(points), max(points)
    if read.is_reverse:
        lo, hi = read.query_length - hi, read.query_length - lo
    return lo, hi


def matches_interval(evidence, interval, flag, length, max_diff=50, ratio=.6):
    if interval is None or not (flag & evidence.flag):
        return False
    lo, hi = sorted((evidence.st, evidence.en))
    # Require the complete evidence span, not an incidental overlap.
    if lo < interval[0] or hi > interval[1]:
        return False
    a, b = abs(length), abs(evidence.len)
    return abs(a - b) <= max_diff or (min(a, b) >= ratio * max(a, b))


def candidate_region(fields, info):
    """VCF POS is one-based; TRID intervals are zero-based catalog coordinates."""
    trid = info.get('TRID', '')
    match = re.fullmatch(r'(.+)_(\d+)_(\d+)(?:_.*)?', trid)
    if match and match[1] == fields[0]:
        return int(match[2]), int(match[3])
    pos = int(fields[1]) - 1
    return pos, max(pos + 1, int(info.get('END', fields[1])))


def filter_locus(vcf, bam, gsv, output, audit, min_count=3, flank=100,
                 min_mapq=20, max_diff=50, ratio=.6, paf=None):
    inputs = {Path(p).resolve() for p in (vcf, bam, gsv, paf) if p}
    if Path(output).resolve() in inputs or Path(audit).resolve() in inputs or Path(output).resolve() == Path(audit).resolve():
        raise ValueError('output and audit must be distinct from each other and inputs')
    evidence = parse_msvasm_withtype(gsv)
    aligned = set(evidence)
    if paf:
        paf_open = gzip.open if str(paf).endswith('.gz') else open
        with paf_open(paf, 'rt') as source:
            aligned.update(line.split('\t', 1)[0] for line in source if line.strip())
    opener = gzip.open if str(vcf).endswith('.gz') else open
    with opener(vcf, 'rt') as source, pysam.AlignmentFile(bam) as alignments, \
            open(output, 'w') as out, open(audit, 'w', newline='') as stats:
        writer = csv.DictWriter(stats, delimiter='\t', fieldnames=[
            'id', 'status', 'requested_reads', 'projectable_reads',
            'matching_reads', 'matching_read_names'])
        writer.writeheader()
        for line in source:
            if line.startswith('#'):
                out.write(line)
                continue
            fields = line.rstrip().split('\t')
            info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
            names = set(info.get('RNAMES', '').split(',')) - {'', '.'}
            typ = info.get('SVTYPE', '')
            intervals = defaultdict(list)
            supported = set()
            status = ('UNKNOWN_MISSING_CONTIG' if fields[0] not in alignments.references
                      else 'UNKNOWN_UNSUPPORTED_TYPE')
            if typ in ('INS', 'DEL', 'DUP', 'INV') and fields[0] in alignments.references:
                start, end = candidate_region(fields, info)
                for read in alignments.fetch(fields[0], max(0, start - flank), end + flank):
                    if read.query_name not in names or read.mapping_quality < min_mapq:
                        continue
                    region = query_interval(read, start, end, flank)
                    if region is not None:
                        intervals[read.query_name].append(region)
                flag = simple_type(typ, fields[0], fields[0])
                length = int(info['SVLEN'])
                for name, regions in intervals.items():
                    # Multiple primary alignments with the same name are ambiguous.
                    if len(regions) == 1 and any(matches_interval(
                            item, regions[0], flag, length, max_diff, ratio)
                            for item in evidence.get(name, [])):
                        supported.add(name)
                projectable = sum(len(v) == 1 for v in intervals.values())
                status = ('PASS_LOCUS_EVIDENCE' if len(supported) >= min_count else
                          'UNKNOWN_INSUFFICIENT_PROJECTION' if projectable < min_count else
                          'UNKNOWN_MISSING_PARENT_ALIGNMENT' if any(
                              name not in aligned for name in names) else
                          'FAIL_INSUFFICIENT_LOCUS_EVIDENCE')
            if fields[6] not in ('.', 'PASS'):
                status = 'EXCLUDED_INPUT_FILTER'
            writer.writerow(dict(id=fields[2], status=status, requested_reads=len(names),
                                 projectable_reads=sum(len(v) == 1 for v in intervals.values()),
                                 matching_reads=len(supported),
                                 matching_read_names=','.join(sorted(supported))))
            if status == 'PASS_LOCUS_EVIDENCE':
                out.write(line)
