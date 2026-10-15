"""Public interface for the de novo specialization of Minisv."""
import csv
import json
import shutil
from pathlib import Path
from importlib.metadata import PackageNotFoundError, version

import rich_click as click


try:
    VERSION = version('minisv')
except PackageNotFoundError:
    VERSION = '0.2.0'


@click.group('denovo', help='Family de novo SV discovery and sequence review.')
@click.version_option(VERSION, prog_name='minisv denovo')
def cli():
    pass


@cli.command('init')
@click.option('--output', '-o', required=True, type=click.Path(path_type=Path))
@click.option('--father', default='father', show_default=True)
@click.option('--mother', default='mother', show_default=True)
@click.option('--child', multiple=True, default=('child',), show_default=True)
def init_manifest(output, father, mother, child):
    """Write a family manifest template; repeat --child for siblings."""
    names = [father, mother, *child]
    if len(names) != len(set(names)) or any(not s or any(c in s for c in '/\\\t\n ') or s in ('.','..') for s in names):
        raise click.ClickException('Use unique sample identifiers without whitespace or path separators')
    config = dict(reference='/data/reference/ref.fa',
        father=dict(sample=father, bam=f'/data/bam/{father}.bam'),
        mother=dict(sample=mother, bam=f'/data/bam/{mother}.bam'),
        trgt_directory='/data/trgt', minimap2='minimap2',
        parent_assemblies=[f'/data/assembly/{s}.hap{h}.fa' for s in (father,mother) for h in (1,2)],
        children=[dict(sample=s,bam=f'/data/bam/{s}.bam',
            candidate_vcfs=[f'/data/candidates/{s}.{caller}.minisv.vcf' for caller in ('sniffles2','longcallD')],
            trgt_signal=f'/data/trgt/{s}.trgt_denovo.tsv') for s in child])
    try:
        with output.open('x') as target:
            json.dump(config,target,indent=2)
            target.write('\n')
    except OSError as error:
        raise click.ClickException(str(error)) from error
    click.echo(str(output))


@cli.command('check')
@click.option('--manifest', '-m', required=True, type=click.Path(exists=True,path_type=Path))
def check_manifest(manifest):
    """Check family inputs, BAM/reference agreement and alignment tools."""
    from .trio_pipeline import validate
    try:
        config = json.loads(manifest.read_text())
        validate(config)
        paths = []
        for child in config['children']:
            for key in ('candidate_vcf','trgt_signal','parent_paf','parent_gsv','c3_vcf'):
                if child.get(key):
                    paths.append(child[key])
            paths.extend(child.get('candidate_vcfs',[]))
            if child.get('local_directory'):
                if not list(Path(child['local_directory']).glob('*.json')):
                    raise ValueError(f'No local evidence: {child["local_directory"]}')
            else:
                directory = Path(config['trgt_directory'])
                for member in (config['father'], config['mother'], child):
                    paths.extend([str(directory/f'{member["sample"]}.sorted.vcf.gz'),
                                  str(directory/f'{member["sample"]}.spanning.sorted.bam')])
            if not child.get('c3_vcf') and not child.get('parent_paf'):
                haps = config.get('parent_assemblies',[])
                if len(haps) != 4:
                    raise ValueError('Four parental haplotype FASTAs required')
                paths.extend(haps)
                if not shutil.which(config.get('minimap2','minimap2')):
                    raise ValueError('minimap2 executable not found')
        missing = sorted({p for p in paths if not Path(p).is_file()})
        if missing:
            raise ValueError('Missing input files: ' + ', '.join(missing))
    except (OSError, ValueError, KeyError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(dict(status='OK',children=len(config['children']),
                              compression='pigz' if shutil.which('pigz') else 'gzip')))


@cli.command('export')
@click.option('--input', '-i', 'source', required=True, type=click.Path(exists=True,path_type=Path))
@click.option('--output', '-o', required=True, type=click.Path(path_type=Path))
def export_events(source, output):
    """Export reviewed parent-relative JSONL events to a compact TSV."""
    fields = ['event_id','sample','chrom','start','end','type','delta','basis','tier',
              'child_specific','father_target','mother_target','c3_associated','audit_flags']
    try:
        # Exclusive creation protects input files and existing exports.
        with source.open() as stream, output.open('x',newline='') as target:
            writer=csv.DictWriter(target,fieldnames=fields,delimiter='\t')
            writer.writeheader()
            for line in stream:
                if not line.strip():
                    continue
                event=json.loads(line)
                row={key:event.get(key,'') for key in fields[:9]}
                support=event.get('flank_recheck',{})
                row.update(child_specific=support.get('C',{}).get('specific',''),
                    father_target=support.get('F',{}).get('target',''),
                    mother_target=support.get('M',{}).get('target',''),
                    c3_associated=bool(event.get('c3_event_ids')),
                    audit_flags=','.join(event.get('audit_flags',[])))
                writer.writerow(row)
    except (OSError,ValueError) as error:
        raise click.ClickException(str(error)) from error


