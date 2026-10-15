import json

from tests.test_trio_finalize import fixture
from minisv.trio_finalize import classify, process_locus
from minisv.trio_prepare import prepare_candidates
from minisv.trio_pipeline import read_records
from minisv.trio_c3_only import group_calls


def test_33_boundary_is_inclusive_and_default_stays_50(tmp_path):
    data = fixture()
    assert classify(data, 'C0.0', -33, .08, minimum=33)['status'] == 'RULE_PASS'
    assert 'BELOW_33_BP' in classify(data, 'C0.0', -32, .08, minimum=33)['reasons']
    assert 'BELOW_50_BP' in classify(data, 'C0.0', -33, .08)['reasons']
    data['models']['C0.0'] = 'A' * 127
    data['transmissions'][0]['total_edit_distance'] = 33
    path = tmp_path / 'evidence.json'
    path.write_text(json.dumps(data))
    assert len(process_locus(('child', str(path), [.08], 33))['events']) == 1
    assert process_locus(('child', str(path), [.08]))['events'] == []


def test_34_setting_excludes_33_without_changing_evidence_rules():
    data = fixture()
    assert classify(data,'C0.0',-34,.08,minimum=34)['status'] == 'RULE_PASS'
    assert 'BELOW_34_BP' in classify(data,'C0.0',-33,.08,minimum=34)['reasons']


def test_prepare_and_c3_only_keep_33_but_exclude_32(tmp_path):
    vcf = tmp_path / 'calls.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' +
        ''.join(f'chr1\t{n*10}\tins{n}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={n};RNAMES=r1,r2\n' for n in (32,33,49,50)))
    out = tmp_path / 'out'
    out.mkdir()
    combined = prepare_candidates(dict(candidate_vcfs=[str(vcf)]), dict(minimum_variant_length=33), out)
    assert [int(info['SVLEN']) for _, info in read_records(combined)] == [33,49,50]
    calls = [dict(type='DEL', length=n, basis='REFERENCE', sample='child', chrom='chr1',
                  pos0=n*100, compatible_local_event_ids=[]) for n in (32,33,49,50)]
    assert len(group_calls(calls, minimum=33)) == 3
    assert len(group_calls(calls)) == 1
