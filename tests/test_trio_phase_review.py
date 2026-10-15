from minisv.trio_phase_review import consensus,review_phase


def test_consensus_keeps_true_insertion_and_removes_single_read_error():
    seq='ACGT'*20
    rows=[seq,seq,seq,seq,seq[:30]+'TT'+seq[30:]]
    c=consensus(rows)
    assert c['sequence']==seq
    assert c['ambiguous_columns']==[]
    inserted=seq[:30]+'A'*34+seq[30:]
    assert consensus([inserted]*4)['sequence']==inserted


def test_consensus_exposes_mixed_alleles():
    c=consensus(['AACCGG']*2+['AATCGG']*2)
    assert c['ambiguous_columns']


def test_marker_review_does_not_erase_parent_target(monkeypatch):
    sequence='ACGT'*100
    target=sequence[:200]+'T'*34+sequence[200:]
    def phased(*args):
        groups={(1,hp):[dict(name=f'{hp}.{i}',sequence=sequence,reverse=bool(i%2))
                        for i in range(4)] for hp in (1,2)}
        marker=dict(ps=1,pos0=501,alleles={str(hp):dict(names=[r['name'] for r in groups[(1,hp)]]) for hp in (1,2)})
        return {r['name']:r for rs in groups.values() for r in rs},groups,[marker]
    monkeypatch.setattr('minisv.trio_phase_review.phase_groups',phased)
    e=dict(event_id='e',tier='INSUFFICIENT_OR_CONFLICTING_EVIDENCE',chrom='chr1',start=500,
        parent_paths=dict(F0=sequence,F1=sequence,M0=sequence,M1=sequence),
        child_sequence=target,audit_flags=['PARENT_TARGET_AFTER_FLANK_RECHECK',
        'PARENT_ALLELE_COVERAGE_UNRESOLVED'])
    r=review_phase(e,dict(F='f',M='m'))
    assert all(v['passed'] for v in r['parents'].values())
    assert r['tier']!='HIGH_CONFIDENCE_CANDIDATE'
    assert 'PARENT_TARGET_AFTER_FLANK_RECHECK' in r['audit_flags']


def test_missing_phase_markers_never_pass(monkeypatch):
    monkeypatch.setattr('minisv.trio_phase_review.phase_groups',lambda *a:({}, {}, []))
    e=dict(event_id='e',tier='INSUFFICIENT_OR_CONFLICTING_EVIDENCE',chrom='chr1',start=500,
        parent_paths=dict(F0='ACGT'),child_sequence='ACGT'+'T'*34,
        audit_flags=['PARENT_ALLELE_COVERAGE_UNRESOLVED'])
    assert review_phase(e,dict(F='f',M='m'))['tier']!='HIGH_CONFIDENCE_CANDIDATE'


def test_hp_tags_without_direct_marker_link_never_pass(monkeypatch):
    sequence='ACGT'*100
    groups={(1,hp):[dict(name=f'{hp}.{i}',sequence=sequence,reverse=bool(i%2))
                    for i in range(4)] for hp in (1,2)}
    marker=dict(ps=1,pos0=501,alleles={str(hp):dict(names=['other']) for hp in (1,2)})
    monkeypatch.setattr('minisv.trio_phase_review.phase_groups',lambda *a:(
        {r['name']:r for rs in groups.values() for r in rs},groups,[marker]))
    e=dict(event_id='e',tier='INSUFFICIENT_OR_CONFLICTING_EVIDENCE',chrom='chr1',start=500,
        parent_paths=dict(F0=sequence,F1=sequence,M0=sequence,M1=sequence),
        child_sequence=sequence[:200]+'T'*34+sequence[200:],
        audit_flags=['PARENT_ALLELE_COVERAGE_UNRESOLVED'])
    r=review_phase(e,dict(F='f',M='m'))
    assert r['tier']!='HIGH_CONFIDENCE_CANDIDATE'
    assert all('MARKER_LINKED_MOLECULES_INSUFFICIENT' in g['flags']
        for p in r['parents'].values() for a in p['attempts'] for g in a['groups'])
