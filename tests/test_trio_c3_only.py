import random

import pysam

from minisv.trio_c3_only import group_calls, review_deletion
from minisv.trio_pipeline import review_c3_only


def fixture(tmp_path, parent_alt=False, short_child=False, chrom='chr1'):
    bases = random.Random(17).choices('ACGT', k=1000)
    sequence = ''.join(bases)
    reference = tmp_path/'ref.fa'
    reference.write_text(f'>{chrom}\n{sequence}\n')
    pysam.faidx(str(reference))
    bams = {}
    for role in ('C', 'F', 'M'):
        path = tmp_path/f'{role}.bam'
        with pysam.AlignmentFile(path, 'wb', header={'SQ': [{'SN': chrom, 'LN': 1000}]}) as bam:
            for i in range(12):
                alt = (role == 'C' and i < 5) or (parent_alt and role == 'F' and i == 0)
                read = pysam.AlignedSegment()
                read.query_name = f'{role}{i}'
                read.reference_id = 0
                read.reference_start = 350 if short_child and role == 'C' else 0
                read.mapping_quality = 60
                read.flag = 16 if i % 2 else 0
                if short_child and role == 'C':
                    read.cigartuples = [(0, 200)]
                    aligned = sequence[350:550]
                elif alt:
                    read.cigartuples = [(0, 401), (2, 60), (0, 539)]
                    aligned = sequence[:401] + sequence[461:]
                else:
                    read.cigartuples = [(0, 1000)]
                    aligned = sequence
                read.query_sequence = aligned
                bam.write(read)
        pysam.index(str(path))
        bams[role] = str(path)
    return str(reference), bams


def calls(chrom='chr1'):
    return [dict(event_id=f'child.C3.{i}', sample='child', chrom=chrom,
                 pos0=400, length=60, type='DEL', basis='REFERENCE',
                 record_ids=[f'caller{i}'], compatible_local_event_ids=[])
            for i in (1, 2)]


def test_two_callers_become_one_high_confidence_event(tmp_path):
    groups = group_calls(calls())
    assert len(groups) == 1
    ref, bams = fixture(tmp_path)
    result = review_deletion((next(iter(groups.values())), ref, bams, .08, 'family'))
    assert result['tier'] == 'HIGH_CONFIDENCE_CANDIDATE', result['flank_recheck']
    assert result['c3_record_ids'] == ['caller1', 'caller2']
    assert result['flank_recheck']['C']['alt_specific'] == 5
    assert result['flank_recheck']['F']['ref_specific'] == 12


def test_local_match_suppresses_duplicate_c3_only_event():
    rows = calls()
    rows[0]['compatible_local_event_ids'] = ['local.1']
    assert group_calls(rows) == {}


def test_nearby_caller_representations_share_one_event():
    rows = calls()
    rows[1]['pos0'] += 4
    rows[1]['length'] -= 5
    groups = group_calls(rows)
    assert len(groups) == 1
    assert len(next(iter(groups.values()))) == 2
    rows[1]['compatible_local_event_ids'] = ['local.1']
    assert group_calls(rows) == {}


def test_parent_alt_support_prevents_high_confidence(tmp_path):
    ref, bams = fixture(tmp_path, parent_alt=True)
    result = review_deletion((calls(), ref, bams, .08, 'family'))
    assert result['tier'] == 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    assert 'PARENT_TARGET_AFTER_FLANK_RECHECK' in result['audit_flags']


def test_nonspanning_child_reads_do_not_count(tmp_path):
    ref, bams = fixture(tmp_path, short_child=True)
    result = review_deletion((calls(), ref, bams, .08, 'family'))
    assert result['flank_recheck']['C']['depth'] == 0
    assert result['tier'] != 'HIGH_CONFIDENCE_CANDIDATE'


def test_sex_chromosome_requires_ploidy_review(tmp_path):
    ref, bams = fixture(tmp_path, chrom='chrX')
    result = review_deletion((calls('chrX'), ref, bams, .08, 'family'))
    assert result['tier'] == 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    assert 'PLOIDY_NOT_VALIDATED' in result['audit_flags']


def test_pipeline_review_writes_one_event_for_nearby_callers(tmp_path):
    ref, bams = fixture(tmp_path)
    rows = calls()
    rows[1]['pos0'] += 4
    rows[1]['length'] -= 5
    config = dict(reference=ref, father=dict(sample='dad', bam=bams['F']),
                  mother=dict(sample='mom', bam=bams['M']))
    children = {'child': dict(sample='child', bam=bams['C'])}
    reviewed = review_c3_only(rows, config, children, .08, 1, tmp_path)
    assert len(reviewed) == 1
    assert reviewed[0]['tier'] == 'HIGH_CONFIDENCE_CANDIDATE'
    assert len((tmp_path/'c3_only_reviewed.jsonl').read_text().splitlines()) == 1
