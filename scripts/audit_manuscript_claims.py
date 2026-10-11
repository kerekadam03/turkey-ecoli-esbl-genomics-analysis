#!/usr/bin/env python3
"""Compare descriptive observations against *separately* documented manuscript claims.

A value appearing in a manuscript cannot stand in for the underlying isolate-level
records. This tool returns REVIEW, not invented validation, if inputs are absent.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from pathlib import Path
from common import table,save_tsv,tri_bool

CHECKS=['study_isolate_n','phenotypic_esbl_positive_n','blactxm_positive_n',
        'external_reference_genome_n','farm_n','plasmid_predicted_CTXM_carriers_n','MGE_CTXM_carriers_n']

def audit(config, metadata=None, phenotype=None, ctxm=None, reference=None):
    expected=config.get('MANUSCRIPT_CLAIMS_TO_VERIFY',{})
    measured={}
    if metadata is not None:
        measured['study_isolate_n']=len(metadata)
        farms={r.get('farm_id') for r in metadata if r.get('farm_id')}
        if all(r.get('farm_id') for r in metadata):measured['farm_n']=len(farms)
    if phenotype is not None:
        measured['phenotypic_esbl_positive_n']=sum(tri_bool(x.get('laboratory_esbl_phenotype')) is True for x in phenotype)
    if ctxm is not None:
        measured['blactxm_positive_n']=len(ctxm)
        if all(x.get('blactxm_on_plasmid_predicted_contig') in ('Yes','No') for x in ctxm):
            measured['plasmid_predicted_CTXM_carriers_n']=sum(x['blactxm_on_plasmid_predicted_contig']=='Yes' for x in ctxm)
        measured['MGE_CTXM_carriers_n']=sum(x.get('has_mge_proximity')=='Yes' for x in ctxm)
    if reference is not None:
        measured['external_reference_genome_n']=len(reference)
    return [{'check':key,'manuscript_claim':expected.get(key,''),
      'measured_from_supplied_input':measured.get(key,''),
      'assessment':'NO EVIDENCE' if key not in measured else 'MATCH' if measured[key]==expected.get(key) else 'MISMATCH',
      'important':'Manuscript claims alone never establish actual source values'} for key in CHECKS]


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--claims-json',required=True)
    ap.add_argument('--cohort')
    ap.add_argument('--phenotypes')
    ap.add_argument('--ctxm-summary')
    ap.add_argument('--references')
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    config=json.loads(Path(args.claims_json).read_text(encoding='utf-8'))
    meta=table(args.cohort,['isolate_id']) if args.cohort else None
    phen=table(args.phenotypes,['isolate_id','laboratory_esbl_phenotype']) if args.phenotypes else None
    ctxm=table(args.ctxm_summary,['isolate_id']) if args.ctxm_summary else None
    refs=table(args.references,['accession']) if args.references else None
    save_tsv(args.out,['check','manuscript_claim','measured_from_supplied_input','assessment','important'],audit(config,meta,phen,ctxm,refs))
    print('Audited values against claims; review all NO EVIDENCE and MISMATCH rows before submission.')

if __name__=='__main__':main()
