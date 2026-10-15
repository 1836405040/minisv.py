import random

import pysam
import pytest

from minisv.trio_insertion import project_target, review_insertion


def fixture(tmp_path, parent_target=False, father_snp=False, neighbor=False, length=34):
    seq = ''.join(random.Random(4).choices('ACGT',k=1000))
    ins = ('ACGTTGCA'*7)[:length]
    ref = tmp_path/'ref.fa'
    ref.write_text('>chr1\n'+seq+'\n')
    pysam.faidx(str(ref))
    haps = {}
    for key in ('F0','F1','M0','M1'):
        model = seq
        cs = ':1000'
        if key == 'F1' and father_snp:
            model = seq[:300]+('A' if seq[300] != 'A' else 'C')+seq[301:]
            cs = ':300*'+seq[300].lower()+model[300].lower()+':699'
        fa = tmp_path/f'{key}.fa'
        fa.write_text('>ctg\n'+model+'\n')
        pysam.faidx(str(fa))
        paf = tmp_path/f'{key}.paf'
        paf.write_text(f'ctg\t1000\t0\t1000\t+\tchr1\t1000\t0\t1000\t1000\t1000\t60\tcs:Z:{cs}\n')
        haps[key] = dict(fasta=str(fa),paf=str(paf))
    bams = {}
    for role in ('C','F','M'):
        path = tmp_path/f'{role}.bam'
        with pysam.AlignmentFile(path,'wb',header={'SQ':[{'SN':'chr1','LN':1000}]}) as bam:
            for i in range(12):
                read = pysam.AlignedSegment()
                read.query_name = f'{role}{i}'
                read.reference_id = 0
                read.reference_start = 0
                read.mapping_quality = 60
                read.flag = 16 if i % 2 else 0
                target = role == 'C' and i < 6 or parent_target and role == 'F' and i == 0
                if target:
                    read.query_sequence = seq[:401]+ins+seq[401:]
                    read.cigartuples = [(0,401),(1,len(ins)),(0,599)]
                elif neighbor and role == 'F':
                    read.query_sequence = seq[:92]+'A'*26+seq[92:]
                    read.cigartuples = [(0,92),(1,26),(0,908)]
                else:
                    read.query_sequence = seq
                    read.cigartuples = [(0,1000)]
                bam.write(read)
        pysam.index(str(path))
        bams[role] = str(path)
    call = dict(sample='child',chrom='chr1',pos0=400,ref=seq[400],
                alt=seq[400]+ins,length=len(ins),record_ids=['caller'])
    return call,str(ref),bams,haps,.08,'family'


def test_34_bp_ins_requires_molecules_and_four_parent_paths(tmp_path):
    result = review_insertion(fixture(tmp_path))
    assert result['delta'] == 34
    assert result['tier'] == 'HIGH_CONFIDENCE_CANDIDATE'
    assert result['variant_class'] == 'de_novo_indel'
    assert result['flank_recheck']['C']['specific'] == 6
    assert len(result['flank_recheck']['C']['strands']) == 2
    assert result['parent_allele_coverage']['F'][0]['models'] == ['F0','F1']


def test_parent_target_prevents_strict_high_confidence(tmp_path):
    result = review_insertion(fixture(tmp_path,parent_target=True))
    assert 'PARENT_TARGET_AFTER_FLANK_RECHECK' in result['audit_flags']
    assert result['tier'] == 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'


def test_distinct_but_unreadable_parent_paths_remain_unresolved(tmp_path):
    result = review_insertion(fixture(tmp_path,father_snp=True))
    assert 'PARENT_ALLELE_COVERAGE_UNRESOLVED' in result['audit_flags']


def test_parent_neighbor_is_not_the_target(tmp_path):
    result = review_insertion(fixture(tmp_path,neighbor=True))
    assert result['tier'] == 'HIGH_CONFIDENCE_CANDIDATE'
    assert result['flank_recheck']['F']['target'] == 0


def test_symbolic_insertion_and_29_bp_rejected(tmp_path):
    job = list(fixture(tmp_path))
    job[0]['alt'] = '<INS>'
    assert review_insertion(job)['tier'] == 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    job[0]['alt'] = job[0]['ref']+'A'*29
    job[0]['length'] = 29
    assert review_insertion(job)['tier'] == 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'


@pytest.mark.parametrize('length', [30, 33, 49, 50])
def test_indel_sv_boundary_with_actual_bam_molecules(tmp_path, length):
    result=review_insertion(fixture(tmp_path,length=length))
    assert result['tier']=='HIGH_CONFIDENCE_CANDIDATE'
    assert result['delta']==length
    assert result['variant_class']==('de_novo_indel' if length<50 else 'de_novo_SV')


def test_assembly_projection_rejects_gap_boundary_and_invalid_span():
    f = 'ctg 1034 0 1034 + chr1 1000 0 1000 1000 1034 60'.split()
    f.append('cs:Z::401+'+'a'*34+':599')
    assert project_target(f,401) is None
    assert project_target(f,500) == 534
    f[8] = '999'
    assert project_target(f,500) is None
