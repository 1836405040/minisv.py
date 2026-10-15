import pytest
from minisv.trio_structural import review_structural
from minisv.trio_parent_topology import review_parent_topology
from tests.test_trio_structural_confidence import phased_split_fixture


@pytest.mark.parametrize('missing_sa,expected',[(False,'PARENT_STRUCTURAL_JUNCTIONS_SUPPORTED'),(True,'PARENT_SEQUENCE_ONLY_UNRESOLVED')])
def test_parent_sequence_and_topology_are_separate(tmp_path,missing_sa,expected):
    event,ref,bams=phased_split_fixture(tmp_path,bad_sa=missing_sa)
    # Use the same synthetic variant-bearing molecules as a transmitting parent.
    bams['F']=bams['C']
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    audit=result['structural_confidence']['parent_topology']
    assert audit['classification']==expected
    assert 'PARENT_TARGET_JUNCTION' in result['audit_flags']
    assert result['tier']=='INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    assert not result['structural_confidence']['eligible_high_confidence']
    assert audit['junctions'][0]['parents']['F']['molecules']


def test_no_parent_target_does_not_replace_callability(tmp_path):
    event,ref,bams=phased_split_fixture(tmp_path,missing_hp=True)
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    assert result['structural_confidence']['parent_topology']['classification']=='NO_PARENT_TARGET_SEQUENCE_OBSERVED'
    assert 'PARENT_PHASED_CALLABILITY_UNRESOLVED' in result['audit_flags']
    assert result['tier']!='HIGH_CONFIDENCE_CANDIDATE'


def test_partial_inversion_and_missing_review_are_not_negative(tmp_path):
    event,ref,bams=phased_split_fixture(tmp_path)
    bams['F']=bams['C']
    result=review_structural((event,ref,bams,.08,'strict_v1'))
    for role in ('F','M'):
        support=result['junction_evidence'][1]['support'][role]
        support['alt_names']=[];support['alt_reads']=[]
    audit=review_parent_topology(result,bams)
    assert audit['classification']=='LIMITED_OR_PARTIAL_PARENT_TOPOLOGY'
    result['junction_evidence']=[]
    assert review_parent_topology(result,bams)['classification']=='INCOMPLETE_PARENT_TOPOLOGY_REVIEW'
