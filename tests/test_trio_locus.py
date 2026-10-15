from types import SimpleNamespace

import pysam

from minisv.trio_locus import query_interval, matches_interval, candidate_region
from minisv.trio_alleles import allele_hypotheses


def alignment(cigar, reverse=False):
    read = pysam.AlignedSegment()
    read.query_name = 'r'
    read.reference_start = 100
    read.cigartuples = cigar
    read.query_sequence = 'A' * sum(n for op,n in cigar if op in (0,1,4,7,8))
    read.flag = 16 if reverse else 0
    return read


def test_projection_includes_insertion_and_reverse():
    forward = alignment([(0,100),(1,50),(0,200)])
    reverse = alignment([(0,100),(1,50),(0,200)], True)
    assert query_interval(forward,200,210,0) == (100,160)
    assert query_interval(reverse,200,210,0) == (190,250)


def test_uncallable_projection():
    assert query_interval(alignment([(0,100)]),150,250,0) is None
    assert query_interval(alignment([(5,10),(0,300)]),150,250,0) is None
    assert query_interval(alignment([(0,100),(3,50),(0,100)]),150,300,0) is None


def test_wrong_read_position_does_not_support():
    item = SimpleNamespace(flag=1,len=100,st=9000,en=9100)
    assert not matches_interval(item,(100,300),1,100)
    item.st,item.en=100,200
    assert matches_interval(item,(100,300),1,100)
    assert not matches_interval(item,(100,300),2,100)


def test_catalog_and_vcf_coordinate_conventions():
    fields=['chr1','101']
    assert candidate_region(fields,{'END':'150'}) == (100,150)
    assert candidate_region(fields,{'TRID':'chr1_100_200_trsolve'}) == (100,200)


def test_parent_delta_not_mean_difference():
    row=dict(denovo_status='Y:+',index='0',child_AL='771,1885',
             father_AL='699,1723',mother_AL='1885,3057',allele_origin='F:1',
             mean_diff_father='49')
    assert [h['delta'] for h in allele_hypotheses(row)] == [72]


def test_unknown_origin_keeps_hypotheses_not_verdict():
    row=dict(denovo_status='Y:?',index='1',child_AL='32,224',
             father_AL='32,296',mother_AL='32,32',allele_origin='?')
    h=allele_hypotheses(row)
    assert -72 in [a['delta'] for a in h]
    assert all(not a['origin_resolved'] for a in h)


def test_end_to_end_locus_filter(tmp_path):
    import gzip
    import csv
    from minisv.trio_locus import filter_locus
    bam = tmp_path/'reads.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':10000}]}) as out:
        for i in range(3):
            read=alignment([(0,100),(1,100),(0,300)])
            read.reference_id=0
            read.mapping_quality=60
            read.query_name=f'r{i}'
            out.write(read)
    pysam.index(str(bam))
    vcf=tmp_path/'in.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   'chr1\t200\tgood\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=100;RNAMES=r0,r1,r2\n'
                   'chr1\t400\twrong\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=100;RNAMES=r0,r1,r2\n')
    gsv=tmp_path/'e.gsv.gz'
    with gzip.open(gsv,'wt') as out:
        for i in range(3):out.write(f'p\t500\t500\tr{i}\t60\t+\tSVTYPE=INS;SVLEN=100;qoff_l=100;qoff_r=200\n')
    output=tmp_path/'out.vcf';audit=tmp_path/'audit.tsv'
    filter_locus(str(vcf),str(bam),str(gsv),str(output),str(audit),flank=10)
    assert '\tgood\t' in output.read_text()
    assert '\twrong\t' not in output.read_text()
    with audit.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert [r['matching_reads'] for r in rows]==['3','0']


def test_candidate_builder_preserves_reads_and_parent_delta(tmp_path):
    import csv
    from minisv.trio_alleles import build_candidates
    row=dict(chrom='chr1',start='200',end='250',trid='chr1_200_250_trsolve',
             denovo_status='Y:+',index='0',child_AL='771',father_AL='699',
             mother_AL='771',allele_origin='F:1',mean_diff_father='49')
    signal=tmp_path/'signal.tsv'
    with signal.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
    bam=tmp_path/'span.bam'
    with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':10000}]}) as f:
        read=alignment([(0,300)]);read.reference_id=0
        read.set_tag('TR',row['trid']);read.set_tag('AL',0);f.write(read)
    out=tmp_path/'out.vcf';audit=tmp_path/'audit.jsonl'
    build_candidates(str(signal),str(bam),str(out),str(audit))
    assert 'SVLEN=72;' in out.read_text()
    assert 'RNAMES=r;' in out.read_text()
    assert 'chr1\t201\t' in out.read_text()
