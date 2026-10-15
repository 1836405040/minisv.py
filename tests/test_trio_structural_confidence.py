import pysam
import pytest
from minisv.trio_structural_confidence import segment, endpoint_match, duplication_dosage, assess
from minisv.trio_structural import review_structural, rc
from tests.test_trio_structural import fixture


def phased_split_fixture(tmp_path,typ='INV',missing_hp=False,bad_sa=False):
    event,ref,members,_=fixture(tmp_path,typ)
    if typ=='BND':
        event['audit_flags']=[]
        event['records']=[dict(id='a',source='x',original_record='chr1\t1000\ta\tN\tN[chr2:2000[\t.\tPASS\tSVTYPE=BND;MATEID=b'),
            dict(id='b',source='x',original_record='chr2\t2000\tb\tN\t]chr1:1000]N\t.\tPASS\tSVTYPE=BND;MATEID=a')]
    for role,m in members.items():
        raw=tmp_path/(role+'.new.bam')
        with pysam.AlignmentFile(m['bam']) as src, pysam.AlignmentFile(raw,'wb',header=src.header) as dst:
            for r in src:
                if role!='C' and not missing_hp:
                    r.set_tag('HP',1+int(r.query_name.split('.')[-1])%2);r.set_tag('PS',10)
                if role=='C' and '.alt.' in r.query_name:
                    j,i=map(int,r.query_name.split('.')[-2:]);a,b=event['junctions'][j]
                    # Build reference-oriented alignments of the two physical read halves.
                    alt=r.query_sequence
                    a_reverse=(a[2]=='R') != bool(i%2)
                    b_reverse=(b[2]=='L') != bool(i%2)
                    physical=rc(alt) if i%2 else alt
                    r.query_sequence=rc(physical) if a_reverse else physical
                    r.flag=16 if a_reverse else 0
                    r.reference_start=a[1]-200 if a[2]=='L' else a[1]
                    r.cigarstring='200M200S' if a[2]=='L' else '200S200M'
                    bpos=(b[1]-200 if b[2]=='L' else b[1])+1
                    bcigar='200M200S' if b[2]=='L' else '200S200M'
                    if not bad_sa:r.set_tag('SA',f'{b[0]},{bpos},{"-" if b_reverse else "+"},{bcigar},60,0;')
                dst.write(r)
        dest=tmp_path/(role+'.sorted.bam');pysam.sort('-o',str(dest),str(raw));pysam.index(str(dest));m['bam']=str(dest)
    return event,ref,{r:m['bam'] for r,m in members.items()}


@pytest.mark.parametrize('typ',['INV','BND'])
def test_strict_positive_requires_sequence_topology_and_phased_parents(tmp_path,typ):
    event,ref,bams=phased_split_fixture(tmp_path,typ)
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    assert result['tier']=='HIGH_CONFIDENCE_CANDIDATE',result['audit_flags']
    assert result['structural_confidence']['eligible_high_confidence']
    audited=assess(result,bams,'audit')
    assert audited['tier']=='REVIEW_REQUIRED'


@pytest.mark.parametrize('kw,flag',[({'missing_hp':True},'PARENT_PHASED_CALLABILITY_UNRESOLVED'),
    ({'bad_sa':True},'CHILD_SPLIT_TOPOLOGY_UNRESOLVED')])
def test_missing_orthogonal_evidence_prevents_promotion(tmp_path,kw,flag):
    event,ref,bams=phased_split_fixture(tmp_path,**kw)
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    assert flag in result['audit_flags']
    assert result['tier']=='INSUFFICIENT_OR_CONFLICTING_EVIDENCE'


def test_existing_parent_conflict_cannot_be_cleared(tmp_path):
    event,ref,bams=phased_split_fixture(tmp_path)
    event['audit_flags'].append('PARENT_TARGET_JUNCTION')
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    assert 'PARENT_TARGET_JUNCTION' in result['audit_flags']
    assert result['tier']!='HIGH_CONFIDENCE_CANDIDATE'


def test_split_cigar_hard_clips_and_orientation():
    s=segment('chr1',800,'-','200M200H',60)
    assert endpoint_match(s,('chr1',1000,'L'))==(200,False)
    assert s['length']==400
    assert segment('chr1',800,'+','90M310S',60) is None
    assert segment('chr1',800,'+','200M200S',10) is None


@pytest.mark.parametrize('extra_child,extra_parent,passed',[(10,0,True),(0,0,False),(10,10,False)])
def test_duplication_dosage_from_aligned_bases(tmp_path,extra_child,extra_parent,passed):
    bams={}
    for role in ('C','F','M'):
        raw=tmp_path/(role+'.bam')
        with pysam.AlignmentFile(raw,'wb',header={'SQ':[{'SN':'chr1','LN':10000}]}) as out:
            for i in range(20):
                r=pysam.AlignedSegment();r.query_name=f'b{i}';r.query_sequence='A'*10000
                r.reference_id=0;r.reference_start=0;r.mapping_quality=60;r.cigarstring='10000M';out.write(r)
            for i in range(extra_child if role=='C' else extra_parent):
                r=pysam.AlignedSegment();r.query_name=f'e{i}';r.query_sequence='A'*2000
                r.reference_id=0;r.reference_start=3000;r.mapping_quality=60;r.cigarstring='2000M';out.write(r)
        pysam.index(str(raw));bams[role]=str(raw)
    result=duplication_dosage(dict(junctions=[(('chr1',3000,'R'),('chr1',5000,'L'))]),bams)
    assert result['passed']==passed
    assert result['samples']['C']['ratio']==(20+extra_child)/20


def test_tandem_duplication_full_strict_path(tmp_path):
    import random
    from minisv.trio_structural import structural_events,models
    rng=random.Random(12);sequence=''.join(rng.choices('ACGT',k=10000))
    ref=tmp_path/'ref.fa';ref.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(ref))
    line='chr1\t3000\tdup\tN\t<DUP>\t60\tPASS\tSVTYPE=DUP;END=5000'
    events,_=structural_events([dict(sample='kid',id='dup',original_record=line)],'dad__mom')
    event=events[0];pair=event['junctions'][0]
    with pysam.FastaFile(ref) as fa:alt,_=models(fa,pair)
    bams={}
    for role in ('C','F','M'):
        raw=tmp_path/(role+'.raw.bam')
        with pysam.AlignmentFile(raw,'wb',header={'SQ':[{'SN':'chr1','LN':10000}]}) as out:
            for i in range(20):
                r=pysam.AlignedSegment();r.query_name=f'{role}.ref{i}';r.query_sequence=sequence
                r.reference_id=0;r.reference_start=0;r.cigarstring='10000M';r.mapping_quality=60
                r.set_tag('HP',1+i%2);r.set_tag('PS',10);out.write(r)
            if role=='C':
                for i in range(10):
                    r=pysam.AlignedSegment();r.query_name=f'extra{i}';r.query_sequence=sequence[3000:5000]
                    r.reference_id=0;r.reference_start=3000;r.cigarstring='2000M';r.mapping_quality=60;out.write(r)
                for i in range(4):
                    r=pysam.AlignedSegment();r.query_name=f'alt{i}';r.query_sequence=rc(alt)
                    r.reference_id=0;r.reference_start=3000;r.cigarstring='200S200M';r.mapping_quality=60;r.flag=16 if i%2==0 else 0
                    r.set_tag('SA',f'chr1,4801,{"-" if i%2==0 else "+"},200M200S,60,0;');out.write(r)
        dest=tmp_path/(role+'.bam');pysam.sort('-o',str(dest),str(raw));pysam.index(str(dest));bams[role]=str(dest)
    result=review_structural((event,str(ref),bams,.08,'strict_v1'))
    assert result['tier']=='HIGH_CONFIDENCE_CANDIDATE',result['audit_flags']
    assert result['structural_confidence']['dosage']['passed']
