import json
import random
from pathlib import Path

import pysam
import pytest

from minisv.trio_structural import canonical, junctions, models, rc, review_structural, structural_events, run_manifest


def fields(typ, alt=None, extra=';END=2000', pos=1000, chrom='chr1', ident='v1'):
    return [chrom,str(pos),ident,'N',alt or '<'+typ+'>','60','PASS','SVTYPE='+typ+extra]


@pytest.mark.parametrize('alt,expected', [
    ('N[chr2:2000[', (('chr1',1000,'L'),('chr2',1999,'R'))),
    ('N]chr2:2000]', (('chr1',1000,'L'),('chr2',2000,'L'))),
    ('[chr2:2000[N', (('chr1',999,'R'),('chr2',1999,'R'))),
    (']chr2:2000]N', (('chr1',999,'R'),('chr2',2000,'L'))),
])
def test_all_four_bnd_orientations(alt,expected):
    assert junctions(fields('BND',alt)) == [expected]


def test_inv_and_dup_coordinates():
    assert junctions(fields('INV')) == [(('chr1',1000,'L'),('chr1',2000,'L')), (('chr1',1000,'R'),('chr1',2000,'R'))]
    assert junctions(fields('DUP')) == [(('chr1',2000,'L'),('chr1',1000,'R'))]


def test_reciprocal_bnd_dedup_and_mate_validation():
    a=fields('BND','N[chr2:2000[',';MATEID=v2')
    b=fields('BND',']chr1:1000]N',';MATEID=v1',2000,'chr2','v2')
    records=[dict(sample='kid',id=f[2],original_record='\t'.join(f)) for f in (a,b)]
    events, unresolved=structural_events(records,'family')
    assert len(events)==1 and not unresolved
    assert not events[0]['audit_flags']
    assert len(events[0]['records'])==2
    events,_=structural_events(records[:1],'family')
    assert 'BND_MATE_LINK_UNRESOLVED' in events[0]['audit_flags']


@pytest.mark.parametrize('alt', ['N[chr2:2000]', 'NN[chr2:2000[', '<BND>'])
def test_complex_bnd_is_explicitly_unresolved(alt):
    f=fields('BND',alt)
    events, unresolved=structural_events([dict(sample='kid',id='v1',original_record='\t'.join(f))],'f')
    assert not events and len(unresolved)==1


def fixture(tmp_path, typ='INV', inherited=False, no_parent=False, one_junction=False):
    rng=random.Random(91)
    sequences={c:''.join(rng.choices('ACGT',k=6000)) for c in ('chr1','chr2')}
    reference=tmp_path/'ref.fa'
    reference.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in sequences.items()))
    pysam.faidx(str(reference))
    f=fields(typ, 'N[chr2:2000[' if typ=='BND' else None)
    events,_=structural_events([dict(sample='kid',id='v1',original_record='\t'.join(f))],'dad__mom')
    event=events[0]
    with pysam.FastaFile(reference) as fa:
        alleles=[models(fa,pair)[0] for pair in event['junctions']]
    members={}
    for role,sample in [('C','kid'),('F','dad'),('M','mom')]:
        raw=tmp_path/(sample+'.raw.bam')
        with pysam.AlignmentFile(raw,'wb',header={'SQ':[{'SN':c,'LN':6000} for c in sequences]}) as bam:
            if not (no_parent and role!='C'):
                endpoints=set(ep for pair in event['junctions'] for ep in pair)
                for chrom,boundary,side in endpoints:
                    for i in range(12):
                        r=pysam.AlignedSegment();r.query_name=f'{sample}.ref.{chrom}.{boundary}.{i}'
                        r.reference_id=list(sequences).index(chrom);r.reference_start=boundary-200
                        r.query_sequence=sequences[chrom][boundary-200:boundary+200]
                        r.cigarstring='400M';r.mapping_quality=60
                        bam.write(r)
            if role=='C' or (inherited and role=='F'):
                for j,alt in enumerate(alleles):
                    if one_junction and j:
                        continue
                    chrom,boundary,side=event['junctions'][j][0]
                    for i in range(4):
                        r=pysam.AlignedSegment();r.query_name=f'{sample}.alt.{j}.{i}'
                        r.reference_id=list(sequences).index(chrom);r.reference_start=boundary-200
                        r.query_sequence=alt;r.flag=16 if i%2 else 0
                        r.cigarstring='200M200S';r.mapping_quality=60
                        bam.write(r)
        path=tmp_path/(sample+'.bam')
        pysam.sort('-o',str(path),str(raw));pysam.index(str(path))
        members[role]=dict(sample=sample,bam=str(path))
    return event,str(reference),members,f


@pytest.mark.parametrize('typ',['INV','DUP','BND'])
def test_positive_junctions_reproduced_but_not_unvalidated_high_confidence(tmp_path,typ):
    event,ref,members,f=fixture(tmp_path,typ)
    result=review_structural((event,ref,{r:m['bam'] for r,m in members.items()},.08))
    assert result['junction_support_reproduced']
    assert all(j['support']['C']['alt_specific']==4 for j in result['junction_evidence'])
    assert all(len(j['support']['C']['alt_strands'])==2 for j in result['junction_evidence'])
    assert 'PARENT_TARGET_JUNCTION' not in result['audit_flags']
    assert 'PARENT_BREAKPOINT_CALLABILITY_UNRESOLVED' not in result['audit_flags']
    assert result['tier']!='HIGH_CONFIDENCE_CANDIDATE'


@pytest.mark.parametrize('kw,flag', [({'inherited':True},'PARENT_TARGET_JUNCTION'),
    ({'no_parent':True},'PARENT_BREAKPOINT_CALLABILITY_UNRESOLVED'),
    ({'one_junction':True},'CHILD_JUNCTION_NOT_REPRODUCED')])
def test_negative_and_incomplete_controls(tmp_path,kw,flag):
    event,ref,members,_=fixture(tmp_path,**kw)
    result=review_structural((event,ref,{r:m['bam'] for r,m in members.items()},.08))
    assert flag in result['audit_flags']
    assert result['tier']=='INSUFFICIENT_OR_CONFLICTING_EVIDENCE'


def test_standalone_manifest_preserves_provenance_and_masks(tmp_path):
    event,ref,members,f=fixture(tmp_path)
    source=tmp_path/'input.vcf';source.write_text('\t'.join(f)+'\n')
    mask=tmp_path/'excluded.bed';mask.write_text('chr1\t1500\t1600\n')
    config=dict(reference=ref,father=members['F'],mother=members['M'],
        children=[dict(members['C'],candidate_vcf=str(source))],structural_exclude_bed=str(mask))
    manifest=tmp_path/'manifest.json';manifest.write_text(json.dumps(config))
    summary=run_manifest(manifest,tmp_path/'result',1)
    assert summary['reviewed_events']==0
    assert json.loads((tmp_path/'result/structural_unresolved.jsonl').read_text())['status']=='OUTSIDE_REGULAR_REGION'
    assert not json.loads((tmp_path/'result/provenance.json').read_text())['truth_used']
    with pytest.raises(FileExistsError):
        run_manifest(manifest,tmp_path/'result',1)


def test_preparation_keeps_reciprocal_mate_identity(tmp_path):
    from minisv.trio_prepare import prepare_candidates
    a=fields('BND','N[chr2:2000[',';MATEID=v2;RNAMES=r1,r2')
    b=fields('BND',']chr1:1000]N',';MATEID=v1;RNAMES=r1,r2',2000,'chr2','v2')
    source=tmp_path/'pairs.vcf'
    source.write_text('\n'.join('\t'.join(f) for f in (a,b))+'\n')
    combined=prepare_candidates(dict(candidate_vcfs=[str(source)]),{},tmp_path)
    records=[dict(sample='kid',source=combined,id=line.split('\t')[2],original_record=line)
             for line in Path(combined).read_text().splitlines() if not line.startswith('#')]
    events,unresolved=structural_events(records,'family')
    assert len(events)==1 and not unresolved and not events[0]['audit_flags']


def test_bnd_type_inferred_from_brackets():
    f=fields('BND','N[chr2:2000[')
    f[7]='MATEID=v2'
    assert junctions(f)==[(('chr1',1000,'L'),('chr2',1999,'R'))]
