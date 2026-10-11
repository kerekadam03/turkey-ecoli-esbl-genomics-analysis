#!/usr/bin/env python3
"""Convert ABRicate TSV into a normalized per-hit ledger WITHOUT claiming RGI agreement.

ABRicate nucleotide identity and coverage are separate from RGI CARD Strict model
classification. The two must not be merged by taking an undocumented union.
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path
from common import table,save_tsv

COLS=['isolate_id','gene','contig_id','start','end','identity_pct','coverage_pct','database','method',
      'source_accession','call_class','plasmid_predicted','mge_proximal','status','method_note']

def normalize_file(path, min_identity=90., min_cov=90.):
    rows=table(path, ['#FILE','SEQUENCE','GENE','%IDENTITY','%COVERAGE'])
    result=[]
    for r in rows:
        pi=float(r['%IDENTITY']);pc=float(r['%COVERAGE'])
        if pi<min_identity or pc<min_cov:continue
        file=Path(r['#FILE']).name
        matched=re.match(r'^(\d+)(?:_.*)?\.(?:fa|fasta|fna)',file,re.I)
        if not matched: raise ValueError(f'No isolate ID parsed from filename: {file}')
        result.append(dict(isolate_id=matched.group(1),gene=r['GENE'],contig_id=r['SEQUENCE'],
             start=r.get('START',''),end=r.get('END',''),identity_pct=pi,coverage_pct=pc,
             database=r.get('DATABASE','CARD'),method='ABRicate nucleotide alignment',
             source_accession=r.get('ACCESSION',''),call_class='UNREVIEWED',plasmid_predicted='Unknown',
             mge_proximal='Unknown',status='Candidate; manually reconcile with RGI',
             method_note='ABRicate >=90% identity and coverage; RGI Strict not inferred'))
    return result

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--min-id',type=float,default=90.)
    ap.add_argument('--min-coverage',type=float,default=90.)
    args=ap.parse_args(argv)
    rows=normalize_file(args.input,args.min_id,args.min_coverage)
    save_tsv(args.out,COLS,rows)
    print(f'Candidate calls retained: {len(rows)}; none classified as CARD/RGI Strict automatically')

if __name__=='__main__':main()
