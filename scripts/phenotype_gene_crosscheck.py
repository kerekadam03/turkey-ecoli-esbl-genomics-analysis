#!/usr/bin/env python3
"""Genotype/phenotype contingency counts ONLY for annotation-complete isolates.

Genuine status from laboratory, not inferred from CTX-M presence. No 'strong
concordance' interpretation is embedded, because conditional directions differ.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from common import table,save_tsv,unique_index,tri_bool

OUT=['phenotypic_group','blaCTXM_carriage','isolates_n','group_denominator_n','fraction_percent',
     'interpretation']

def crosscheck(cohort, phenotypes, calls, completeness):
    ids=list(cohort)
    phen=unique_index(phenotypes,'isolate_id');complete=unique_index(completeness,'isolate_id')
    positives={r['isolate_id'] for r in calls if r['gene'].lower().startswith('blactx-m') and r['status'].lower()=='verified'}
    counts=Counter();den=Counter()
    for code in ids:
        status=tri_bool(phen.get(code,{}).get('laboratory_esbl_phenotype'))
        if status is None:continue
        if tri_bool(complete.get(code,{}).get('annotation_complete')) is not True:continue
        p='ESBL-positive' if status else 'ESBL-negative'
        c='CTX-M-present' if code in positives else 'CTX-M-not-detected'
        counts[(p,c)]+=1;den[p]+=1
    for p in ('ESBL-positive','ESBL-negative'):
        for c in ('CTX-M-present','CTX-M-not-detected'):
            n=counts[(p,c)];N=den[p]
            yield {'phenotypic_group':p,'blaCTXM_carriage':c,'isolates_n':n,
                'group_denominator_n':N,'fraction_percent':f'{n/N*100:.2f}' if N else '',
                'interpretation':'Descriptive joint distribution; complete annotation required to infer not-detected'}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',required=True)
    p.add_argument('--phenotypes',required=True)
    p.add_argument('--calls',required=True)
    p.add_argument('--complete',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args(argv)
    cohort=[r['isolate_id'] for r in table(a.cohort,['isolate_id'])]
    if len(cohort)!=len(set(cohort)):raise ValueError('Repeated cohort isolate ID')
    rows=list(crosscheck(cohort,table(a.phenotypes,['isolate_id','laboratory_esbl_phenotype']),
       table(a.calls,['isolate_id','gene','status']),table(a.complete,['isolate_id','annotation_complete'])))
    save_tsv(a.out,OUT,rows)
    print(f'Genotype-phenotype contingency with annotation completeness: {a.out}')

if __name__=='__main__':main()
