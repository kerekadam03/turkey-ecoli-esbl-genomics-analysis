#!/usr/bin/env python3
"""Gene co-occurrence n/N (%) by independently assigned ESBL phenotype.

No imputation of gene absence when an isolate has incomplete annotation.
No p-values or co-selection/physical linkage claims are generated.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from common import table,save_tsv,unique_index,tri_bool

OUT=['gene','call_class','ESBL_positive_n','ESBL_positive_N','ESBL_positive_pct',
     'ESBL_negative_n','ESBL_negative_N','ESBL_negative_pct','interpretation']

def summarize(cohort, phenotypes, calls, completeness):
    ids=list(cohort)
    phen=unique_index(phenotypes,'isolate_id');compl=unique_index(completeness,'isolate_id')
    good={code for code in ids if tri_bool(compl.get(code,{}).get('annotation_complete')) is True}
    who={code:tri_bool(phen.get(code,{}).get('laboratory_esbl_phenotype')) for code in good}
    pos={code for code in who if who[code] is True}; neg={code for code in who if who[code] is False}
    genes=defaultdict(set);classes=defaultdict(set)
    for c in calls:
        if c['isolate_id'] not in ids:raise ValueError('Call outside cohort')
        if c['status'].strip().lower()!='verified':continue
        if c['isolate_id'] in good:
            genes[c['gene']].add(c['isolate_id']); classes[c['gene']].add(c['call_class'])
    for gene,carriers in sorted(genes.items()):
        p=len(pos&carriers);n=len(neg&carriers)
        yield {'gene':gene,'call_class':';'.join(sorted(classes[gene])),
          'ESBL_positive_n':p,'ESBL_positive_N':len(pos),'ESBL_positive_pct':f'{100*p/len(pos):.2f}' if pos else '',
          'ESBL_negative_n':n,'ESBL_negative_N':len(neg),'ESBL_negative_pct':f'{100*n/len(neg):.2f}' if neg else '',
          'interpretation':'Isolate-level co-occurrence only; not physical co-localization or co-selection'}

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('cohort','phenotypes','calls','complete','out'):p.add_argument('--'+arg,required=True)
    a=p.parse_args(argv)
    c=[x['isolate_id'] for x in table(a.cohort,['isolate_id'])]
    if len(c)!=len(set(c)):raise ValueError('Repeated cohort ID')
    rows=list(summarize(c,table(a.phenotypes,['isolate_id','laboratory_esbl_phenotype']),
       table(a.calls,['isolate_id','gene','status','call_class']),
       table(a.complete,['isolate_id','annotation_complete'])))
    save_tsv(a.out,OUT,rows)
    print(f'Summarized {len(rows)} verified gene labels; denominator excludes unverified annotations')

if __name__=='__main__':main()
