import gzip
import json
import random
import shutil

import pysam
import pytest

from minisv.trio_io import combine_parents, compressed_writer, extract_candidate_reads
from minisv.trio_pipeline import read_records, run_c3
from minisv.trio_prepare import prepare_candidates


def write_vcf(path, rows):
    path.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(rows))


def test_merge_filters_and_keeps_source_identity(tmp_path):
    vcf = tmp_path / 'caller.vcf'
    write_vcf(vcf, [f'chr1\t101\tsame\tN\t<INS>\t.\t{flt}\tSVTYPE=INS;SVLEN={size};RNAMES={names}\n'
        for flt, size, names in [('PASS', 60, 'r1,r2'), ('PASS', 49, 'r1,r2'),
                                 ('LowQual', 60, 'r1,r2'), ('PASS', 60, 'r1,r1')]])
    output = tmp_path / 'out'
    output.mkdir()
    result = prepare_candidates(dict(candidate_vcfs=[str(vcf), str(vcf)]), {}, output)
    records = list(read_records(result))
    assert len(records) == 2
    assert records[0][0][2] != records[1][0][2]
    assert json.loads((output / 'candidate_summary.json').read_text())['total_records'] == 2


def test_complete_molecules_recovered_outside_candidate_window(tmp_path):
    bam = tmp_path / 'child.bam'
    with pysam.AlignmentFile(bam, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 10000}]}) as target:
        for name, start, flag, sequence in [('near', 50, 0, 'A' * 100), ('far', 8000, 16, 'C' * 100)]:
            read = pysam.AlignedSegment()
            read.query_name, read.reference_id, read.reference_start = name, 0, start
            read.flag, read.mapping_quality = flag, 60
            read.cigartuples, read.query_sequence = [(0, 100)], sequence
            target.write(read)
    pysam.index(str(bam))
    vcf = tmp_path / 'input.vcf'
    write_vcf(vcf, ['chr1\t101\tx\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60;RNAMES=near,far\n'])
    fq = tmp_path / 'reads.fq.gz'
    result = extract_candidate_reads(vcf, str(bam), fq, 1)
    assert result == dict(requested=2, found=2, missing=0, recovered_by_full_scan=1)
    with pysam.FastxFile(str(fq)) as source:
        sequences = {r.name: r.sequence for r in source}
    assert sequences == {'near': 'A' * 100, 'far': 'G' * 100}
    write_vcf(vcf, ['chr1\t101\tx\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60;RNAMES=absent\n'])
    with pytest.raises(ValueError, match='missing from BAM'):
        extract_candidate_reads(vcf, str(bam), fq, 1)


def test_compression_fallback_and_parent_names(tmp_path, monkeypatch):
    monkeypatch.setattr('minisv.trio_io.shutil.which', lambda _: None)
    path = tmp_path / 'test.gz'
    with compressed_writer(path, 'wt') as target:
        target.write('test\n')
    with gzip.open(path, 'rt') as source:
        assert source.read() == 'test\n'
    hap = tmp_path / 'hap.fa'
    hap.write_text('>contig\nACGT')
    combined = tmp_path / 'parents.fa'
    combine_parents([hap] * 4, 'dad', 'mom', combined)
    assert [line for line in combined.read_text().splitlines() if line.startswith('>')] == [
        '>dad.hap1.contig', '>dad.hap2.contig', '>mom.hap1.contig', '>mom.hap2.contig']
    hap.write_text('>contig\nACGT\n>contig\nACGT\n')
    with pytest.raises(ValueError, match='Duplicate contig'):
        combine_parents([hap] * 4, 'dad', 'mom', combined)


@pytest.mark.skipif(not shutil.which('minimap2'), reason='minimap2 executable required')
def test_fresh_parent_alignment_to_locus_c3(tmp_path):
    rng = random.Random(31)
    reference = ''.join(rng.choices('ACGT', k=6000))
    insertion = ''.join(rng.choices('ACGT', k=60))
    hap = tmp_path / 'hap.fa'
    hap.write_text('>contig\n' + reference + '\n')
    bam = tmp_path / 'child.bam'
    with pysam.AlignmentFile(bam, 'wb', header={'SQ': [{'SN': 'chr1', 'LN': 6000}]}) as target:
        for i in range(3):
            read = pysam.AlignedSegment()
            read.query_name, read.reference_id, read.reference_start = f'r{i}', 0, 0
            read.mapping_quality = 60
            read.cigartuples = [(0, 3000), (1, 60), (0, 3000)]
            read.query_sequence = reference[:3000] + insertion + reference[3000:]
            target.write(read)
    pysam.index(str(bam))
    vcf = tmp_path / 'candidates.vcf'
    write_vcf(vcf, ['chr1\t3001\tx\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60;SUPPORT=3;RNAMES=r0,r1,r2\n'])
    output = tmp_path / 'run'
    output.mkdir()
    result = run_c3(dict(candidate_vcf=str(vcf), bam=str(bam)),
                   dict(parent_assemblies=[str(hap)] * 4, father=dict(sample='dad'), mother=dict(sample='mom')),
                   output, 2)
    assert len(list(read_records(result))) == 1
    assert json.loads((output / 'extraction.json').read_text())['found'] == 3
    with gzip.open(output / 'parents.paf.gz', 'rt') as source:
        assert all(line.split('\t')[5].startswith(('dad.hap', 'mom.hap')) for line in source)
