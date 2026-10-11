#!/usr/bin/env python3
"""Independent cefotaxime/clavulanate MIC synergy calculation.

This is NOT equivalent to documenting the full CLSI ESBL screen and dual-agent
confirmatory method. It preserves the laboratory's independent phenotype call.
"""
from __future__ import annotations
import argparse
import math
from common import table, save_tsv, unique_index, tri_bool

REQ = ['isolate_id','cefotaxime_mic_ug_ml','cefotaxime_clavulanate_mic_ug_ml',
       'initial_esbl_screen_positive','laboratory_esbl_phenotype']
OUT = REQ+['fold_reduction','log2_fold_reduction','cefotaxime_screen_ge2','cefotaxime_synergy_ge8',
           'screen_plus_synergy','full_CLSI_two_agent_confirmation','interpretation_note']


def positive_mic(raw):
    raw=str(raw or '').strip().replace(',','.')
    if not raw: return None
    if raw.startswith(('<','>','≤','≥')): return None  # cannot assign a precise ratio from censored values
    try:
        x=float(raw)
        if x<=0 or not math.isfinite(x): raise ValueError('MIC must be finite and >0')
        return x
    except ValueError as e: raise ValueError(f'Invalid exact MIC {raw!r}: {e}')


def calculate(row):
    ctx=positive_mic(row['cefotaxime_mic_ug_ml']); cla=positive_mic(row['cefotaxime_clavulanate_mic_ug_ml'])
    screen=tri_bool(row.get('initial_esbl_screen_positive'))
    if ctx is None or cla is None:
        ratio=None;synergy='Not evaluable';log=''
    else:
        ratio=ctx/cla;synergy='Yes' if ratio>=8 else 'No';log=f'{math.log2(ratio):.4f}'
    if screen is None or ratio is None:
        combination='Not established'
    elif screen and ratio>=8:combination='Screen-positive and CTX synergy-positive'
    else:combination='Not confirmed using supplied CTX screen and synergy data'
    z=dict(row)
    z.update({'fold_reduction':f'{ratio:.5f}' if ratio is not None else '',
     'log2_fold_reduction':log,'cefotaxime_screen_ge2':('Yes' if ctx>=2 else 'No') if ctx is not None else 'Unknown',
     'cefotaxime_synergy_ge8':synergy,'screen_plus_synergy':combination,
     'full_CLSI_two_agent_confirmation':'Not evaluated: ceftazidime-based confirmation not supplied',
     'interpretation_note':'Ratio-only cefotaxime+clavulanate evidence; phenotype never inferred from genes'})
    return z


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mic',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args(argv)
    rows=table(args.mic,REQ);unique_index(rows,'isolate_id',source=args.mic)
    save_tsv(args.out,OUT,[calculate(r) for r in rows])

if __name__=='__main__':main()
