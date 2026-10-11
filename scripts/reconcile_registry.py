#!/usr/bin/env python3
"""Map NCBI BioSample, genome-info and *explicitly validated* manuscript cohort.

Never equate all registered BioSamples (167) to the study denominator (160).
The script fails when a provided cohort has duplicates or non-registry IDs.
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
from common import table, save_tsv, unique_index

HEAD=['isolate_id','sample_name','biosample_accession','bioproject','wgs_root_accession','ncbi_status',
      'genome_filename','registered_organism','manuscript_cohort','taxonomy_review_needed','source_warning']


def build(attrs, genomes, cohort=None, expected_cohort=160):
    bio=unique_index(attrs,'accession',source='BioSample attributes')
    names={}
    for r in attrs:
        nm=r.get('sample_name','')
        if nm in names:raise ValueError('Repeated BioSample sample_name '+nm)
        names[nm]=r
    gc={}
    for r in genomes:
        k=r.get('biosample_accession','')
        if not k or k in gc:raise ValueError('Duplicated/blank genome BioSample '+k)
        gc[k]=r
    cohort_set=set(cohort) if cohort is not None else None
    if cohort_set is not None:
        if len(cohort_set)!=len(cohort):raise ValueError('Duplicated IDs in cohort whitelist')
        if len(cohort_set)!=expected_cohort:raise ValueError(f'Expected {expected_cohort} verified manuscript IDs, found {len(cohort_set)}')
    out=[];represented=set();warnings=[]
    for b in attrs:
        sample=b['sample_name']
        mt=re.fullmatch(r'Ecoli_Turkey_(\d+)',sample)
        isolate=(mt.group(1) if mt else (b.get('strain') or '').strip())
        if not isolate:raise ValueError(f'Cannot infer isolate code from {sample}')
        if isolate in represented:raise ValueError(f'Duplicated isolate code {isolate}')
        represented.add(isolate)
        linked=gc.get(b['accession'])
        if linked is None:warnings.append(f'No genome record for {isolate} {b["accession"]}')
        published=(b.get('organism') or '')
        if published!='Escherichia coli':warnings.append(f'Organism taxon requires review: {isolate}, {published}')
        row={'isolate_id':isolate,'sample_name':sample,'biosample_accession':b['accession'],
          'bioproject':b.get('bioproject_accession',''),
          'wgs_root_accession':linked.get('genome_acc','') if linked else '',
          'ncbi_status':linked.get('status','') if linked else '',
          'genome_filename':linked.get('filename','') if linked else '',
          'registered_organism':published,
          'manuscript_cohort':('Yes' if isolate in cohort_set else 'No') if cohort_set is not None else 'Not verified',
          'taxonomy_review_needed':'Yes' if published!='Escherichia coli' else 'Verify against current NCBI record',
          'source_warning':'BioSample and WGS registration do not verify species or original MIC assignment'}
        out.append(row)
    if cohort_set is not None and cohort_set-represented:
        raise ValueError('Whitelisted isolate IDs absent from NCBI registry: '+','.join(sorted(cohort_set-represented)))
    return out,warnings


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--biosample',required=True,help='NCBI exported attributes TSV')
    ap.add_argument('--genome-info',required=True,help='NCBI genome-info-*.tsv')
    ap.add_argument('--cohort',help='Verified isolate whitelist; one column isolate_id')
    ap.add_argument('--expected',type=int,default=160)
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    b=table(args.biosample,['accession','sample_name','bioproject_accession','organism'])
    g=table(args.genome_info,['genome_acc','biosample_accession','status','filename'])
    ids=[z['isolate_id'] for z in table(args.cohort,['isolate_id'])] if args.cohort else None
    rows,warnings=build(b,g,ids,args.expected)
    save_tsv(args.out,HEAD,rows)
    for w in warnings:print('REVIEW: '+w)
    print(f'BioSamples={len(rows)}; explicit manuscript cohort size={len(ids) if ids else "NOT VERIFIED"}')

if __name__=='__main__':main()
