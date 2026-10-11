#!/usr/bin/env python3
"""Check Panaroo gene-category thresholds and concatenated alignment length.

Panaroo's --core_threshold for alignments is not the same as a separately
reported descriptive >=99% strict-core category. Warn when they are confused.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from common import table


def alignment_dimensions(fasta):
    records={};title=None;partial=[]
    with Path(fasta).open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if title is not None: records[title]=len(''.join(partial))
                title=line[1:].strip().split()[0];partial=[]
            else:partial.append(line.strip())
    if title is not None:records[title]=len(''.join(partial))
    return {'alignment_taxa':len(records), 'alignment_bp':next(iter(set(records.values()))) if len(set(records.values()))==1 else None,
            'lengths_equal':len(set(records.values()))==1,'alignment_ids':sorted(records)}


def gene_categories(records, genome_count):
    categories={'strict_core_99_to_100':0,'soft_core_95_to_lt99':0,'shell_15_to_lt95':0,'cloud_lt15':0}
    for r in records:
        n=float(r['No. isolates'])
        if n>genome_count: raise ValueError('Panaroo row has more isolates than genome count')
        fraction=n/genome_count
        if fraction>=0.99:categories['strict_core_99_to_100']+=1
        elif fraction>=0.95:categories['soft_core_95_to_lt99']+=1
        elif fraction>=0.15:categories['shell_15_to_lt95']+=1
        else:categories['cloud_lt15']+=1
    categories['all_gene_clusters']=len(records)
    return categories


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--presence-absence',required=True,help='Panaroo gene_presence_absence.csv')
    ap.add_argument('--genomes',type=int,required=True,help='Verified exact number of genomes included in Panaroo')
    ap.add_argument('--alignment',help='core_gene_alignment.aln FASTA')
    ap.add_argument('--panaroo-core-threshold',type=float,default=None,help='Read from actual Panaroo command/log')
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    if args.genomes<1:raise ValueError('Genome count must be >0')
    with open(args.presence_absence,encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    if not rows or 'No. isolates' not in rows[0]:raise ValueError('Missing Panaroo No. isolates field')
    summary=gene_categories(rows,args.genomes)
    if args.alignment:summary.update(alignment_dimensions(args.alignment))
    summary['input_genomes']=args.genomes
    summary['actual_panaroo_core_threshold_supplied']=args.panaroo_core_threshold
    summary['interpretation']='Gene frequencies and concatenated alignment are distinct objects; confirm exact loci in Panaroo logs.'
    issues=[]
    if summary['strict_core_99_to_100']==0 and args.alignment:
        issues.append('No >=99% gene clusters but an alignment was provided; verify --core_threshold, aligned loci and missingness')
    if args.alignment and summary['alignment_bp'] is None:
        issues.append('Alignment taxa have unequal sequence lengths; cannot infer a valid concatenated alignment')
    if args.panaroo_core_threshold is None:issues.append('Original Panaroo --core_threshold not documented')
    summary['validation_issues']=issues
    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    Path(args.out).write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f"Clusters={summary['all_gene_clusters']}; flags={len(issues)}")

if __name__=='__main__':main()
