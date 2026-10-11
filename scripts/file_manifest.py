#!/usr/bin/env python3
"""SHA-256 provenance manifest of input files (optionally outside GitHub)."""
from pathlib import Path
from datetime import datetime,timezone
from common import sha256,save_tsv
import argparse

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',nargs='+',required=True)
    p.add_argument('--out',required=True)
    x=p.parse_args(argv)
    rows=[];time=datetime.now(timezone.utc).isoformat()
    for part in x.inputs:
        files=sorted(Path(part).rglob('*')) if Path(part).is_dir() else [Path(part)]
        for f in files:
            if f.is_file():rows.append({'file_path':str(f),'size_bytes':f.stat().st_size,'sha256':sha256(f),'manifest_created_utc':time})
    save_tsv(x.out,['file_path','size_bytes','sha256','manifest_created_utc'],rows)

if __name__=='__main__':main()
