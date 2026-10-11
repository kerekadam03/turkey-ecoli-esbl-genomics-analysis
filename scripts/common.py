"""Shared, dependency-free TSV helpers for the turkey ESBL reproducibility package.

New reference implementation (October 2026), not an archived command from the study.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path


def table(path, required=()):
    """Read UTF-8/UTF-8-SIG TSV; require exact column headings and unique heading names."""
    p = Path(path)
    with p.open('r', encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        cols = reader.fieldnames or []
        if len(cols) != len(set(cols)):
            raise ValueError(f'{p}: duplicated column headings')
        missing = set(required) - set(cols)
        if missing:
            raise ValueError(f'{p}: missing columns {sorted(missing)}')
        records = list(reader)
    return records


def save_tsv(path, columns, rows):
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=list(columns), extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def unique_index(rows, key, *, source='input'):
    out = {}
    for r in rows:
        value = (r.get(key) or '').strip()
        if not value:
            raise ValueError(f'{source}: row without {key}')
        if value in out:
            raise ValueError(f'{source}: duplicated {key}={value}')
        out[value] = r
    return out


def tri_bool(v):
    v = str(v or '').strip().lower()
    if v in ('yes', 'y', 'true', '1', 'positive', 'present'): return True
    if v in ('no', 'n', 'false', '0', 'negative', 'absent'): return False
    if v in ('', 'unknown', 'na', 'not tested', 'not collected', 'missing', 'uncertain'): return None
    raise ValueError(f'Unrecognized yes/no/unknown: {v!r}')


def wilson_interval(successes, trials, z=1.959963984540054):
    """Two-sided Wilson score confidence interval; descriptive, unadjusted for clustering."""
    if trials <= 0:
        return (math.nan, math.nan)
    p = successes / trials
    denominator = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / denominator
    half = z/denominator * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials))
    return (max(0.0, center-half), min(1.0, center+half))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''):
            h.update(b)
    return h.hexdigest()
