#!/usr/bin/env python3
"""Audit blaCTX-M-carrier predictions and MGE hits; NEVER claim direct transfer."""
from __future__ import annotations
import argparse
from collections import defaultdict
from common import table,save_tsv,unique_index,tri_bool

REQ=['isolate_id','gene','contig_id','plasmid_predicted','mge_element','mge_hit_id','distance_bp']
OUT=['isolate_id','blaCTXM_gene_labels','blactxm_on_plasmid_predicted_contig','has_mge_proximity',
     'mge_hits_n','mge_element_names','comment']


def ctxm(s):
    s=(s or '').lower().replace(' ','').replace('_','-')
    return s.startswith('blactx-m') or s.startswith('ctx-m')


def summarize(rows):
    d=defaultdict(list)
    for r in rows:
        if ctxm(r['gene']):d[r['isolate_id']].append(r)
    ret=[]
    for isolate,items in sorted(d.items(),key=lambda t:t[0]):
        genes=sorted({x['gene'] for x in items})
        p=[tri_bool(x['plasmid_predicted']) for x in items]
        m=[x for x in items if x['mge_element'].strip()]
        hits={(x['contig_id'],x['mge_hit_id']) for x in m}
        if any(not x[1] for x in hits):raise ValueError(f'MGE hit ID missing for {isolate}; cannot count hits reproducibly')
        ret.append({'isolate_id':isolate,'blaCTXM_gene_labels':'; '.join(genes),
          'blactxm_on_plasmid_predicted_contig':'Yes' if any(z is True for z in p) else 'No' if all(z is False for z in p) else 'Unknown',
          'has_mge_proximity':'Yes' if m else 'No',
          'mge_hits_n':len(hits),'mge_element_names':'; '.join(sorted({x['mge_element'] for x in m})),
          'comment':'Predicted contig origin/proximity; not plasmid reconstruction nor evidence of horizontal transfer'})
    return ret


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    rows=summarize(table(args.input,REQ))
    save_tsv(args.out,OUT,rows)
    print(f'CTX-M carriers={len(rows)}; MGE-near carriers={sum(x["has_mge_proximity"]=="Yes" for x in rows)}; MGE hits={sum(int(x["mge_hits_n"]) for x in rows)}')

if __name__=='__main__':main()
