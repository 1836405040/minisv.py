"""Sequence comparisons for catalog-anchored trio loci; no truth-directed scoring."""
import itertools
import json
import re
from collections import Counter
from pathlib import Path

import edlib
import pysam

from .trio_locus import query_interval


def full_read_evidence(bam, chrom, start, end, models, min_mapq=20):
    """Reference-CIGAR projection is diagnostic and can be ambiguous in repeats."""
    rows=[];seen=set()
    with pysam.AlignmentFile(str(bam)) as source:
        for read in source.fetch(chrom,max(0,start-100),end+100):
            if read.query_name in seen or read.mapping_quality<min_mapq:continue
            if query_interval(read,start,end,50) is None:continue
            interval=query_interval(read,start,end,0)
            if interval is None:continue
            seq=read.get_forward_sequence()[interval[0]:interval[1]]
            if read.is_reverse:
                seq=seq.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
            if not seq:continue
            distances={key:edlib.align(seq,value,mode='NW',task='distance')['editDistance']
                       for key,value in models.items()}
            best=min(distances.values())
            rows.append(dict(name=read.query_name,length=len(seq),mapq=read.mapping_quality,
                             best_distance=best,best_models=[k for k,v in distances.items() if v==best],
                             distances=distances))
            seen.add(read.query_name)
    return rows


def edit_summary(child, parent):
    result = edlib.align(child, parent, mode='NW', task='path')
    counts = Counter()
    for length, op in re.findall(r'(\d+)([=XID])', result['cigar'] or ''):
        counts[op] += int(length)
    return dict(distance=result['editDistance'], delta=len(child)-len(parent),
                cigar=result['cigar'], operations=dict(counts))


def load_alleles(vcf, sample, chrom, start, end, trid):
    with pysam.VariantFile(str(vcf)) as source:
        records = [r for r in source.fetch(chrom, max(0,start-1),end)
                   if r.info.get('TRID') == trid]
        if len(records) != 1:
            raise ValueError(f'{sample}: expected one record for {trid}')
        r=records[0]; data=r.samples[sample]
        # AL/spanning-BAM labels can be in length order while GT is phased.
        # Keep same-length sequence alternatives rather than guessing a mapping.
        called=[(i,r.alleles[i][1:]) for i in data['GT'] if i is not None]
        alleles=[]
        for index,length in enumerate(data['AL']):
            choices=sorted({seq for _,seq in called if len(seq)==length})
            if not choices:
                raise ValueError(f'{sample}: AL length {length} not represented by GT sequences')
            alleles.append(dict(index=index,length=length,sequences=choices,
                                sequence_mapping_unique=len(choices)==1))
        return dict(sample=sample,gt=list(data['GT']),alleles=alleles)


def read_evidence(bam, chrom, start, end, trid, models):
    rows=[]; seen=set()
    with pysam.AlignmentFile(str(bam)) as source:
        for read in source.fetch(chrom,max(0,start-1),end):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or not read.has_tag('TR') or read.get_tag('TR')!=trid
                    or not read.has_tag('FL') or read.query_name in seen):
                continue
            left,right=map(int,read.get_tag('FL'))
            seq=read.query_sequence
            if not seq or left+right>=len(seq):continue
            # TRGT's trimmed spanning reads carry reference-oriented sequence.
            # Reverse records need an explicit tag convention before use.
            if read.is_reverse:continue
            seq=seq[left:len(seq)-right if right else len(seq)]
            distances={key:edlib.align(seq,value,mode='NW',task='distance')['editDistance']
                       for key,value in models.items()}
            best=min(distances.values())
            rows.append(dict(name=read.query_name,allele_label=read.get_tag('AL') if read.has_tag('AL') else None,
                             mapq=read.mapping_quality,length=len(seq),best_distance=best,
                             best_models=[k for k,v in distances.items() if v==best],distances=distances))
            seen.add(read.query_name)
    return rows


def compare_locus(directory, child, chrom, start, end, trid, father='NA12877', mother='NA12878'):
    directory=Path(directory)
    family={role:load_alleles(directory/f'{s}.sorted.vcf.gz',s,chrom,start,end,trid)
            for role,s in [('F',father),('M',mother),('C',child)]}
    models={}
    for role,sample in family.items():
        for allele in sample['alleles']:
            for version,seq in enumerate(allele['sequences']):
                models[f'{role}{allele["index"]}.{version}']=seq
    pairs=[]
    for c,cs in models.items():
        if not c.startswith('C'):continue
        for p,ps in models.items():
            if p.startswith('C'):continue
            pairs.append(dict(child=c,parent=p,**edit_summary(cs,ps)))
    transmissions=[]
    child_groups=[[k for k in models if k.startswith(f'C{i}.')] for i in range(2)]
    if len(family['C']['alleles'])==2:
        for c0,c1 in itertools.product(*child_groups):
            for f in [k for k in models if k.startswith('F')]:
                for m in [k for k in models if k.startswith('M')]:
                    for paternal,maternal in [(c0,c1),(c1,c0)]:
                        fp=next(p for p in pairs if p['child']==paternal and p['parent']==f)
                        mp=next(p for p in pairs if p['child']==maternal and p['parent']==m)
                        transmissions.append(dict(paternal_child=paternal,father=f,
                            maternal_child=maternal,mother=m,total_edit_distance=fp['distance']+mp['distance'],
                            paternal_delta=fp['delta'],maternal_delta=mp['delta']))
    transmissions.sort(key=lambda x:x['total_edit_distance'])
    reads={role:read_evidence(directory/f'{s["sample"]}.spanning.sorted.bam',chrom,start,end,trid,models)
           for role,s in family.items()}
    return dict(child=child,trid=trid,family=family,models=models,pairs=pairs,
                transmissions=transmissions,reads=reads,
                interpretation='edit-distance ranking is a diagnostic, not a phased inheritance verdict')


def run_loci(directory, loci, output, bam_manifest=None):
    import csv
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    bams={}
    if bam_manifest:
        with open(bam_manifest) as source:
            bams={row['sample']:row['bam'] for row in csv.DictReader(source,delimiter='\t')}
    with open(loci) as source:
        for row in csv.DictReader(source,delimiter='\t'):
            result=compare_locus(directory,row['child'],row['chrom'],int(row['start']),int(row['end']),row['trid'],
                                 father=row.get('father') or 'NA12877', mother=row.get('mother') or 'NA12878')
            result['source_directory']=str(directory)
            if bam_manifest:
                result['full_reads']={}
                result['source_bams']={}
                for role,sample in result['family'].items():
                    bam=bams[sample['sample']]
                    result['source_bams'][role]=bam
                    result['full_reads'][role]=full_read_evidence(
                        bam,row['chrom'],int(row['start']),int(row['end']),result['models'])
            (out/f'{row["child"]}.{row["trid"]}.json').write_text(json.dumps(result,indent=2))
