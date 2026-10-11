#!/usr/bin/env python3
"""Build Figure 3 and Figure 6 source TSVs only from verified gene calls.

Do not use the absolute number of 160 until 160 valid, completely annotated
isolates have been confirmed.  Missing annotation is not gene absence.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from common import table,save_tsv,unique_index,tri_bool


def build(cohort,genes,completeness,phenotypes):
    all_ids=list(cohort);eligible={i for i in all_ids if tri_bool(completeness.get(i,{}).get('annotation_complete')) is True}
    by_gene=defaultdict(set);classes=defaultdict(set);per_isolate=defaultdict(set)
    for call in genes:
        who=call['isolate_id']
        if who not in all_ids:raise ValueError(f'Unregistered gene call {who}')
        if who not in eligible or call.get('status','').lower()!='verified':continue
        gene=call['gene'];by_gene[gene].add(who);per_isolate[who].add(gene);classes[gene].add(call['call_class'])
    fig3=[{'gene':gene,'carriers_n':len(ids),'denominator_n':len(eligible),'call_class':';'.join(sorted(classes[gene]))}
          for gene,ids in sorted(by_gene.items())]
    fig6=[]
    for who in all_ids:
        if who not in eligible:continue
        esbl=phenotypes.get(who,{}).get('laboratory_esbl_phenotype','Unknown')
        if esbl not in ('Positive','Negative'):continue
        fig6.append({'isolate_id':who,'laboratory_esbl_phenotype':esbl,
          'n_distinct_verified_gene_labels':len(per_isolate[who])})
    return fig3,fig6


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('cohort','calls','complete','phenotypes','figure3','figure6'):
        ap.add_argument('--'+key,required=True)
    a=ap.parse_args(argv)
    all_ids=[r['isolate_id'] for r in table(a.cohort,['isolate_id'])]
    if len(all_ids)!=len(set(all_ids)):raise ValueError('Duplicate cohort ID')
    calls=table(a.calls,['isolate_id','gene','call_class','status'])
    comp=unique_index(table(a.complete,['isolate_id','annotation_complete']),'isolate_id')
    phen=unique_index(table(a.phenotypes,['isolate_id','laboratory_esbl_phenotype']),'isolate_id')
    fg3,fg6=build(all_ids,calls,comp,phen)
    save_tsv(a.figure3,['gene','carriers_n','denominator_n','call_class'],fg3)
    save_tsv(a.figure6,['isolate_id','laboratory_esbl_phenotype','n_distinct_verified_gene_labels'],fg6)

if __name__=='__main__':main()
