import pysam
from minisv.trio_local import edit_summary, load_alleles


def test_edit_direction_and_total():
    result=edit_summary('ACGGGT','ACT')
    assert result['delta']==3
    assert result['distance']==3
    assert edit_summary('ACT','ACGGGT')['delta']==-3


def test_phase_order_does_not_change_length_mapping(tmp_path):
    vcf=tmp_path/'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
      '##INFO=<ID=TRID,Number=1,Type=String,Description="locus">\n'
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
      '##FORMAT=<ID=AL,Number=.,Type=Integer,Description="lengths">\n'
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
      'chr1\t100\t.\tTAAA\tTAAAAA\t.\tPASS\tTRID=chr1_100_103_x\tGT:AL\t1|0:3,5\n')
    zipped=pysam.tabix_index(str(vcf),preset='vcf',force=True)
    a=load_alleles(zipped,'S','chr1',100,103,'chr1_100_103_x')
    assert a['gt']==[1,0]
    assert a['alleles'][0]['sequences']==['AAA']
    assert a['alleles'][1]['sequences']==['AAAAA']
