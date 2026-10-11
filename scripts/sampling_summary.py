#!/usr/bin/env python3
"""Transparent sample-level summaries with Wilson intervals (not farm-adjusted).

No differences are tested; no generalization from selected isolates to population
prevalence. Farm/flock observations may be correlated or repeated.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from common import table,save_tsv,unique_index,tri_bool,wilson_interval

GROUPS=['sampling_site','farm_id','flock_id','region','age_group','production_type']
REQUIRED=['isolate_id','sample_id','bird_id']+GROUPS
OUT=['stratum','category','all_isolates_n','known_esbl_status_n','esbl_positive_n','esbl_positive_pct',
     'wilson_low_pct','wilson_high_pct','missing_esbl_status_n','warning']


def compute(meta,phenotypes):
    cohort=unique_index(meta,'isolate_id',source='validated cohort')
    if len({r['sample_id'] for r in meta if r['sample_id']})!=len([r for r in meta if r['sample_id']]):
        raise ValueError('Duplicate sample_id: clarify whether same bird/sample was sampled multiple times')
    result=[]
    for grouping in GROUPS:
        groups=defaultdict(list)
        for r in meta:groups[(r.get(grouping) or 'Not reported').strip()].append(r)
        for level, rows in sorted(groups.items()):
            if level=='Not reported': pass
            stats=[phenotypes.get(r['isolate_id']) for r in rows]
            known=[x for x in stats if x is not None]
            pos=sum(bool(x) for x in known)
            lo,hi=wilson_interval(pos,len(known))
            result.append({'stratum':grouping,'category':level,'all_isolates_n':len(rows),
              'known_esbl_status_n':len(known),'esbl_positive_n':pos,
              'esbl_positive_pct':f'{pos/len(known)*100:.2f}' if known else '',
              'wilson_low_pct':f'{lo*100:.2f}' if known else '',
              'wilson_high_pct':f'{hi*100:.2f}' if known else '',
              'missing_esbl_status_n':len(rows)-len(known),
              'warning':'Unadjusted descriptive proportion; farm/flock clustering not accounted for'})
    return result


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata',required=True)
    ap.add_argument('--phenotypes',required=True,help='One row per isolate_id, laboratory_esbl_phenotype=Positive/Negative/Unknown')
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    meta=table(args.metadata,REQUIRED)
    statuses=unique_index(table(args.phenotypes,['isolate_id','laboratory_esbl_phenotype']),'isolate_id')
    if set(statuses)-{r['isolate_id'] for r in meta}:raise ValueError('Phenotype IDs not in cohort')
    phenotypes={ident:tri_bool(row['laboratory_esbl_phenotype']) for ident,row in statuses.items()}
    save_tsv(args.out,OUT,compute(meta,phenotypes))

if __name__=='__main__':main()
