import json

import pysam
import pytest

from minisv.trio_pipeline import run
from minisv.trio_review import assign_tiers, parent_allele_coverage


def test_tiers_preserve_risks_and_family_boundaries():
    base = dict(sample='child',trid='chr1_300_400_x',child_sequence='AAA',family_id='family1')
    events = [dict(base,audit_flags=[]),dict(base,sample='other',family_id='family2',audit_flags=[]),
              dict(base,sample='uncallable',child_sequence='CCC',audit_flags=['PARENT_ALLELE_COVERAGE_UNRESOLVED']),
              dict(base,sample='strand',child_sequence='GGG',audit_flags=['SINGLE_STRAND_SUPPORT'])]
    assign_tiers(events)
    assert [e['tier'] for e in events] == ['HIGH_CONFIDENCE_CANDIDATE','HIGH_CONFIDENCE_CANDIDATE',
                                          'INSUFFICIENT_OR_CONFLICTING_EVIDENCE','REVIEW_REQUIRED']
    assert not events[0]['same_sequence_siblings']


def test_parent_coverage_deduplicates_identical_alleles():
    data = dict(models={'F0.0':'AAA','F1.0':'AAA','M0.0':'CCC','M1.0':'GGG'},
                full_reads={'F':[dict(length=3,distances={'F0.0':0,'F1.0':0})]*3,
                            'M':[dict(length=3,distances={'M0.0':0,'M1.0':6})]*3})
    result = parent_allele_coverage(data,.08)
    assert len(result['F']) == 1
    assert result['F'][0]['discriminating_reads'] == 3
    assert [r['discriminating_reads'] for r in result['M']] == [3,0]


def fixture_files(root):
    reference = root/'ref.fa'
    reference.write_text('>chr1\n'+'A'*1000+'\n')
    pysam.faidx(str(reference))
    members = {}
    for role,sample in [('C','kid'),('F','dad'),('M','mom')]:
        bam = root/f'{sample}.bam'
        with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':1000}]}) as target:
            for i in range(12):
                read = pysam.AlignedSegment()
                read.query_name=f'{sample}{i}'
                read.reference_id=0
                read.reference_start=0
                read.mapping_quality=60
                read.flag=16 if i%2 else 0
                read.cigartuples=[(0,400),(1,60),(0,600)] if role=='C' and i<4 else [(0,1000)]
                read.query_sequence='A'*1060 if role=='C' and i<4 else 'A'*1000
                target.write(read)
        pysam.index(str(bam))
        members[role]=dict(sample=sample,bam=str(bam))
    local_dir=root/'local'
    local_dir.mkdir()
    models={'C0.0':'A'*160,'C1.0':'A'*100,'F0.0':'A'*100,'F1.0':'A'*100,'M0.0':'A'*100,'M1.0':'A'*100}
    reads={}
    for role,sample in [('C','kid'),('F','dad'),('M','mom')]:
        reads[role]=[]
        for i in range(12):
            n=160 if role=='C' and i<4 else 100
            distances={k:abs(n-len(s)) for k,s in models.items()}
            reads[role].append(dict(name=f'{sample}{i}',length=n,distances=distances))
    local=dict(child='kid',trid='chr1_300_400_x',models=models,full_reads=reads,
               family={role:dict(sample=member['sample']) for role,member in members.items()},
               transmissions=[dict(total_edit_distance=60,paternal_child='C0.0',father='F0.0',maternal_child='C1.0',mother='M0.0')])
    (local_dir/'chr1_300_400_x.json').write_text(json.dumps(local))
    c3=root/'c3.vcf'
    c3.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    config=dict(reference=str(reference),father=members['F'],mother=members['M'],
                children=[dict(members['C'],local_directory=str(local_dir),c3_vcf=str(c3))])
    manifest=root/'manifest.json'
    manifest.write_text(json.dumps(config))
    return manifest,config


def test_complete_pipeline_from_generic_family_cache_and_bams(tmp_path):
    manifest,config=fixture_files(tmp_path)
    summary=run(manifest,tmp_path/'result',workers=1)
    assert summary['rule_pass_events']==1
    assert summary['tiers']=={'HIGH_CONFIDENCE_CANDIDATE':1}
    assert json.loads((tmp_path/'result/status.json').read_text())['stage']=='COMPLETE'
    assert len((tmp_path/'result/high_confidence_candidates.jsonl').read_text().splitlines())==1


def test_pipeline_rejects_cached_family_mismatch(tmp_path):
    manifest,config=fixture_files(tmp_path)
    config['father']['sample']='wrong_dad'
    manifest.write_text(json.dumps(config))
    with pytest.raises(ValueError,match='Cached family mismatch'):
        run(manifest,tmp_path/'result',workers=1)
    assert json.loads((tmp_path/'result/status.json').read_text())['stage']=='FAILED'


@pytest.mark.parametrize('merge_inputs', [False, True])
def test_recompute_local_and_c3_from_prepared_generic_inputs(tmp_path, merge_inputs):
    import gzip
    import shutil
    manifest,config=fixture_files(tmp_path)
    directory=tmp_path/'trgt'
    directory.mkdir()
    for sample in ('kid','dad','mom'):
        vcf=directory/f'{sample}.sorted.vcf'
        ref='A'*101
        alt='A'*161
        gt='0/1:100,160' if sample=='kid' else '0/0:100,100'
        vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=TRID,Number=1,Type=String,Description="TR">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            '##FORMAT=<ID=AL,Number=.,Type=Integer,Description="AL">\n'
            f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n'
            f'chr1\t300\t.\t{ref}\t{alt}\t.\tPASS\tTRID=chr1_300_400_x\tGT:AL\t{gt}\n')
        pysam.tabix_index(str(vcf),preset='vcf',force=True)
        shutil.copyfile(tmp_path/f'{sample}.bam',directory/f'{sample}.spanning.sorted.bam')
        shutil.copyfile(tmp_path/f'{sample}.bam.bai',directory/f'{sample}.spanning.sorted.bam.bai')
    candidates=tmp_path/'candidates.vcf'
    candidates.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t301\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60;SUPPORT=4;TRID=chr1_300_400_x;RNAMES=kid0,kid1,kid2,kid3\n')
    gsv=tmp_path/'parents.gsv.gz'
    with gzip.open(gsv,'wt') as target:
        for i in range(4):
            strand='-' if i%2 else '+'
            start,end=(600,660) if i%2 else (400,460)
            target.write(f'p\t400\t400\tkid{i}\t60\t{strand}\tSVTYPE=INS;SVLEN=60;qoff_l={start};qoff_r={end}\n')
    paf=tmp_path/'parents.paf'
    paf.write_text(''.join(f'kid{i}\t1060\t0\t1060\t+\tp\t1000\t0\t1000\t1000\t1060\t60\n' for i in range(4)))
    child=config['children'][0]
    child.pop('local_directory')
    child.pop('c3_vcf')
    child.update(candidate_vcf=str(candidates),parent_gsv=str(gsv),parent_paf=str(paf))
    if merge_inputs:
        child['candidate_vcfs'] = [child.pop('candidate_vcf')]
    config['trgt_directory']=str(directory)
    manifest.write_text(json.dumps(config))
    summary=run(manifest,tmp_path/'fresh',workers=1,recompute_local=True)
    assert summary['local_files']==1
    assert summary['rule_pass_events']==1
    assert summary['c3_representations']==1
    assert (tmp_path/'fresh/kid/local/chr1_300_400_x.json').exists()
