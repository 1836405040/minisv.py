import json
from click.testing import CliRunner

from minisv.cli import cli as minisv_cli
cli = minisv_cli.commands['denovo']
from tests.test_trio_pipeline import fixture_files


def test_public_commands_and_identity():
    assert set(cli.commands) == {'init','check','run','export','candidates','locus-filter',
                                 'local-compare','finalize','structural'}
    result=CliRunner().invoke(cli,['--version'])
    assert result.exit_code == 0
    assert 'minisv denovo' in result.output
    assert {'extract','merge','sv-trio-filter','trio-run','denovo'} <= set(minisv_cli.commands)


def test_init_multiple_children_and_no_overwrite(tmp_path):
    runner=CliRunner()
    output=tmp_path/'family.json'
    args=['init','-o',str(output),'--father','dad','--mother','mom','--child','one','--child','two']
    assert runner.invoke(cli,args).exit_code == 0
    data=json.loads(output.read_text())
    assert [c['sample'] for c in data['children']] == ['one','two']
    assert runner.invoke(cli,args).exit_code != 0
    assert runner.invoke(cli,['init','-o',str(tmp_path/'bad.json'),'--child','father']).exit_code != 0


def test_cli_check_run_and_export(tmp_path):
    manifest,config=fixture_files(tmp_path)
    runner=CliRunner()
    assert runner.invoke(cli,['check','-m',str(manifest)]).exit_code == 0
    out=tmp_path/'result'
    result=runner.invoke(cli,['run','--manifest',str(manifest),'--output',str(out),'--workers','1'])
    assert result.exit_code == 0, result.output
    source=out/'all_reviewed.jsonl'
    export=tmp_path/'events.tsv'
    assert runner.invoke(cli,['export','-i',str(source),'-o',str(export)]).exit_code == 0
    lines=export.read_text().splitlines()
    assert len(lines)==2
    assert 'PARENT_ALLELE' in lines[1]
    assert 'HIGH_CONFIDENCE_CANDIDATE' in lines[1]
    assert runner.invoke(cli,['export','-i',str(source),'-o',str(source)]).exit_code != 0
    config['children'][0]['c3_vcf']=str(tmp_path/'missing.vcf')
    manifest.write_text(json.dumps(config))
    assert runner.invoke(cli,['check','-m',str(manifest)]).exit_code != 0
