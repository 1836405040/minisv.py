"""Complete candidate molecule extraction and uniquely named parent references."""
import contextlib
import gzip
import io
import json
import shutil
import subprocess

import pysam

from .trio_locus import candidate_region
from .trio_prepare import open_text


@contextlib.contextmanager
def compressed_writer(path, mode='wb', threads=16):
    if mode not in ('wb', 'wt'):
        raise ValueError('compressed_writer requires wb or wt')
    pigz = shutil.which('pigz')
    if not pigz:
        with gzip.open(path, mode, compresslevel=1) as stream:
            yield stream
        return
    with open(path, 'wb') as destination:
        process = subprocess.Popen([pigz, '-1', '-p', str(threads)], stdin=subprocess.PIPE, stdout=destination)
        stream = io.TextIOWrapper(process.stdin) if mode == 'wt' else process.stdin
        try:
            yield stream
        finally:
            try:
                stream.close()
            finally:
                code = process.wait()
            if code:
                raise RuntimeError(f'pigz failed with code {code}')


def extract_candidate_reads(vcf, bam, output, threads):
    names, regions = set(), {}
    with open_text(vcf) as source:
        for line in source:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip().split('\t')
            if fields[6] not in ('PASS', '.'):
                continue
            info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
            names.update(set(info.get('RNAMES', '').split(',')) - {'', '.'})
            start, end = candidate_region(fields, info)
            regions.setdefault(fields[0], []).append((max(0, start - 1000), end + 1000))
    found = set()
    with compressed_writer(output, 'wt', min(16, threads)) as target:
        def write_read(read):
            if read.query_name not in names or read.query_name in found or read.is_secondary or read.is_supplementary:
                return
            seq = read.get_forward_sequence()
            if not seq:
                return
            qualities = read.get_forward_qualities() if read.query_qualities is not None else None
            qual = ''.join(chr(q + 33) for q in qualities) if qualities is not None else 'I' * len(seq)
            target.write(f'@{read.query_name}\n{seq}\n+\n{qual}\n')
            found.add(read.query_name)
        with pysam.AlignmentFile(bam, threads=min(8, threads)) as source:
            for chrom, intervals in regions.items():
                merged = []
                for start, end in sorted(intervals):
                    if merged and start <= merged[-1][1]:
                        merged[-1][1] = max(end, merged[-1][1])
                    else:
                        merged.append([start, end])
                for start, end in merged:
                    for read in source.fetch(chrom, start, end):
                        write_read(read)
        missing_before_scan = len(names - found)
        if missing_before_scan:
            with pysam.AlignmentFile(bam, threads=min(32, threads)) as source:
                for read in source.fetch(until_eof=True):
                    write_read(read)
    result = dict(requested=len(names), found=len(found), missing=len(names - found),
                  recovered_by_full_scan=missing_before_scan - len(names - found))
    output.with_name('extraction.json').write_text(json.dumps(result, indent=2))
    if names - found:
        raise ValueError(f'{len(names - found)} candidate molecules missing from BAM')
    return result


def combine_parents(haps, father, mother, output):
    if len(haps) != 4:
        raise ValueError('Four parental haplotype FASTAs required')
    with output.open('w') as target:
        for path, sample, hap in zip(haps, [father, father, mother, mother], [1, 2, 1, 2]):
            seen = set()
            with open_text(path) as source:
                for line in source:
                    if line.startswith('>'):
                        name = line[1:].split()[0]
                        if name in seen:
                            raise ValueError(f'Duplicate contig {name} within {path}')
                        seen.add(name)
                        target.write(f'>{sample}.hap{hap}.{line[1:].rstrip()}\n')
                    else:
                        target.write(line.rstrip() + '\n')
            if not seen:
                raise ValueError(f'No FASTA sequences in {path}')
