#!/usr/bin/env python3
"""Descriptive annotation summaries from manually reconciled isolate-level gene calls.

Only calls explicitly given a verified category contribute to acquired-ARG totals.
Missing calls are never imputed as verified absence.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from common import table,save_tsv,unique_index,tri_bool

REQUIRED=['isolate_id','gene','call_class','status','plasmid_predicted','mge_proximal']
CATEGORIES={'acquired','intrinsic','efflux','regulatory','other','uncertain'}
HEAD=['isolate_id','n_distinct_verified_gene_labels','n_acquired','n_intrinsic','n_efflux',
      'n_regulatory','n_other','n_uncertain','n_plasmid_predicted','n_mge_proximal',
      'n_harmonized_gene_calls','interpretation']


def summarize(calls, cohort):
    by_id=defaultdict(list)
    eligible=set(cohort)
    for r in calls:
        code=r['isolate_id']
        if code not in eligible:raise ValueError(f'Gene call isolate {code} not in declared cohort')
        if r['call_class'].strip().lower() not in CATEGORIES:
            raise ValueError(f'Unreviewed gene class {r["call_class"]!r}; reconcile RGI/ABRicate first')
        if r['status'].strip().lower()!='verified':
            raise ValueError(f'Not a verified gene call for {code}/{r["gene"]}')
        by_id[code].append(r)
    out=[]
    for code in cohort:
        items=by_id[code]
        # Collapse multiple hits to unique gene labels at isolate level;
        # allele-level vs family-level IDs depend on the input schema.
        unique={z['gene'] for z in items}
        categories={cat:{z['gene'] for z in items if z['call_class'].strip().lower()==cat} for cat in CATEGORIES}
        plasmid={z['gene'] for z in items if tri_bool(z['plasmid_predicted']) is True}
        mge={z['gene'] for z in items if tri_bool(z['mge_proximal']) is True}
        out.append({'isolate_id':code,'n_distinct_verified_gene_labels':len(unique),
          **{f'n_{cat}':len(categories[cat]) for cat in CATEGORIES},
          'n_plasmid_predicted':len(plasmid),'n_mge_proximal':len(mge),
          'n_harmonized_gene_calls':len(items),
          'interpretation':'No genes listed does NOT imply gene-free; verify input exhaustiveness'})
    return out


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--calls',required=True)
    ap.add_argument('--cohort',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    c=table(args.calls,REQUIRED)
    cohort=[r['isolate_id'] for r in table(args.cohort,['isolate_id'])]
    if len(cohort)!=len(set(cohort)):raise ValueError('Cohort contains duplicate isolate IDs')
    save_tsv(args.out,HEAD,summarize(c,cohort))

if __name__=='__main__':main()
