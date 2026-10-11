#!/usr/bin/env python3
"""Generate iTOL annotation datasets from a **verified** tree-tip metadata TSV.

Do not substitute old 166-tip metadata for a new 160-isolate analysis.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from common import table,unique_index

FIELDS={
 'origin':{'Study isolate':'#406B91','Reference genome':'#8B9299'},
 'sampling_site':{'Cloacal':'#0076A9','Tracheal':'#2C944E','Unknown':'#D9DEE5'},
 'esbl_status':{'Positive':'#CD402E','Negative':'#BAC7D1','Unknown':'#D9DEE5'},
 'ctxm_type':{'CTX-M-1':'#813AAA','CTX-M-27':'#E28D30','None':'#E7E7E7','Unknown':'#D9DEE5'},
 'plasmid_prediction':{'Yes':'#35B3D0','No':'#D2DEE2','Unknown':'#D9DEE5'},
 'mge_proximity':{'Yes':'#4D9F49','No':'#CFDCCF','Unknown':'#D9DEE5'},
}

def make(rows,outdir):
    path=Path(outdir);path.mkdir(parents=True,exist_ok=True)
    unique_index(rows,'tree_tip')
    for field,palette in FIELDS.items():
        if any((r.get(field) or 'Unknown') not in palette for r in rows):
            raise ValueError(f'Unexpected {field}; update the palette explicitly')
        lines=['DATASET_COLORSTRIP','SEPARATOR TAB',f'DATASET_LABEL\t{field}',
               'COLOR\t#596373','STRIP_WIDTH\t20','SHOW_INTERNAL\t0','MARGIN\t0','BORDER_WIDTH\t0',
               'DATA']
        for r in rows:
            value=r.get(field) or 'Unknown'
            lines.append(f'{r["tree_tip"]}\t{palette[value]}\t{value}')
        (path/f'itol_{field}.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return sorted(path.glob('itol_*.txt'))

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tree-metadata',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args(argv)
    rows=table(a.tree_metadata,['tree_tip']+list(FIELDS))
    print('Generated',len(make(rows,a.outdir)),'annotation files for',len(rows),'tree tips')

if __name__=='__main__':main()
