import gzip
import json
from pathlib import Path

import pytest

from minisv.eval import gc_parse_sv
from minisv.trio_prepare import prepare_candidates
from minisv.trio_pipeline import run_c3, c3_events, run
from minisv.trio_types import normalize_record

HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def record(typ, extra='', alt=None):
    return f'chr1\t1000\t{typ}\tN\t{alt or "<"+typ+">"}\t60\tPASS\tSVTYPE={typ};SUPPORT=3;RNAMES=r1,r2,r3{extra}\n'


@pytest.mark.parametrize('typ', ['INV', 'DUP', 'DEL'])
@pytest.mark.parametrize('length', ['', ';SVLEN=0', ';SVLEN=.'])
def test_end_only_intervals_survive_core_parser(tmp_path, typ, length):
    path = tmp_path/'input.vcf'
    path.write_text(HEADER+record(typ, ';END=1050'+length))
    calls = gc_parse_sv(str(path), 50, 3, False, False)
    assert len(calls) == 1
    assert calls[0].SVTYPE == typ
    assert calls[0].SVLEN == (-50 if typ == 'DEL' else 50)
    assert not gc_parse_sv(str(path), 51, 3, False, False)


@pytest.mark.parametrize('extra,reason', [('', 'UNRESOLVED_SVLEN'), (';END=900','UNRESOLVED_SVLEN'), (';SVLEN=bad','INVALID_SVLEN_OR_END')])
def test_unresolved_interval_is_explicit(extra, reason):
    assert normalize_record(record('INV', extra).strip().split('\t'))[3] == reason


def test_prepare_preserves_bnd_without_fake_length_and_audits_tra(tmp_path):
    source = tmp_path/'input.vcf'
    source.write_text(HEADER+record('BND', alt='N[chr2:2000[')+record('INV',';END=1050')+record('TRA'))
    result = prepare_candidates(dict(candidate_vcfs=[str(source)]), {}, tmp_path)
    text = Path(result).read_text()
    assert 'SVTYPE=BND' in text and 'SVTYPE=INV' in text
    assert 'SVLEN=50' in text
    audit = [json.loads(x) for x in (tmp_path/'candidate_preparation_audit.jsonl').read_text().splitlines()]
    assert [x['status'] for x in audit] == ['RETAINED','RETAINED','UNSUPPORTED_TRA_REPRESENTATION']


@pytest.mark.parametrize('cached', [False, True])
def test_structural_only_gzip_routes_without_alignment(tmp_path, cached):
    source = tmp_path/'input.vcf.gz'
    with gzip.open(source, 'wt') as f:
        f.write(HEADER+record('INV',';END=2000')+record('DUP',';END=3000')+record('BND',alt='N[chr2:2000['))
    result = run_c3({('c3_vcf' if cached else 'candidate_vcf'):str(source)}, {}, tmp_path, 1)
    assert c3_events('child', result, []) == []
    pending = [json.loads(x) for x in (tmp_path/'structural_pending.jsonl').read_text().splitlines()]
    assert [x['type'] for x in pending] == ['INV','DUP','BND']
    assert all(x['status']=='PENDING_TYPE_SPECIFIC_REVIEW' for x in pending)
    assert not (tmp_path/'child_reads.fq.gz').exists()


def test_mixed_cached_pipeline_reports_pending_without_promoting(tmp_path):
    from tests.test_trio_pipeline import fixture_files
    manifest, config = fixture_files(tmp_path)
    with open(config['children'][0]['c3_vcf'], 'a') as target:
        target.write(record('INV',';END=2000')+record('BND',alt='N[chr2:2000['))
    summary = run(manifest, tmp_path/'result', workers=1)
    assert summary['structural_pending_records'] == 2
    assert summary['structural_pending_types'] == {'INV':1,'BND':1}
    assert summary['tiers'] == {'HIGH_CONFIDENCE_CANDIDATE':1, 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE':2}
    assert len((tmp_path/'result/structural_pending.jsonl').read_text().splitlines()) == 2


def test_c3_events_requires_explicit_type_routing(tmp_path):
    path = tmp_path/'input.vcf'
    path.write_text(HEADER+record('INV',';END=2000'))
    with pytest.raises(ValueError, match='partition_for_review'):
        c3_events('child', path, [])
