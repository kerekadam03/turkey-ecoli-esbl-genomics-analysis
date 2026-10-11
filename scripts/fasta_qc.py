#!/usr/bin/env python3
"""Non-destructive assembly QC for .fasta/.fna/.fa (and ZIPs containing them).

The outputs are measurements of the sequences supplied to this command; they do
NOT demonstrate that an assembly is uncontaminated, taxonomically valid, or part
of the 160-genome manuscript cohort.
"""
from __future__ import annotations
import argparse
import gzip
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from common import save_tsv

SUFFIXES = ('.fasta', '.fa', '.fna', '.fas', '.fasta.gz', '.fa.gz', '.fna.gz')
HEAD = ['isolate_id','filename','source_archive','assembly_size_bp','contig_count','n50_bp','gc_percent_acgt',
        'ambiguous_N_bp','shortest_contig_bp','longest_contig_bp','contigs_lt_200','contigs_lt_1000',
        'duplicated_ids','invalid_characters','ncbi_size_screen_only','comments']
ALLOWED = set('ACGTRYSWKMBDHVN')


def read_fasta(handle, label):
    """Yield (header id, sequence); reject malformed files and empty records."""
    header = None
    buf = []
    for binary in handle:
        line = binary.decode('ascii').strip() if isinstance(binary, bytes) else binary.strip()
        if not line: continue
        if line.startswith('>'):
            if header is not None:
                seq = ''.join(buf).upper()
                if not seq: raise ValueError(f'{label}: empty sequence for {header}')
                yield header, seq
            header = line[1:].split()[0]
            if not header: raise ValueError(f'{label}: empty FASTA identifier')
            buf=[]
        else:
            if header is None: raise ValueError(f'{label}: nucleotides before first FASTA header')
            buf.append(line)
    if header is None: raise ValueError(f'{label}: no FASTA records')
    seq=''.join(buf).upper()
    if not seq: raise ValueError(f'{label}: empty final sequence')
    yield header, seq


def summarize(records, name, source=''):
    sizes=[]; chars=Counter(); identifiers=Counter()
    for identifier, sequence in records:
        identifiers[identifier]+=1
        sizes.append(len(sequence))
        chars.update(sequence)
    if not sizes: raise ValueError(f'{name}: no contigs')
    total=sum(sizes)
    halfway=(total+1)//2; partial=0; n50=0
    for length in sorted(sizes,reverse=True):
        partial+=length
        if partial>=halfway: n50=length;break
    gc=chars['G']+chars['C']
    acgt=sum(chars[c] for c in 'ACGT')
    invalid=''.join(sorted(set(chars).difference(ALLOWED)))
    specimen=Path(name).name
    match=re.match(r'^(\d+)(?:_spades_megahit[._]gam)?\.(?:fasta|fa|fna|fas)(?:\.gz)?$',specimen,re.I)
    isolate=match.group(1) if match else specimen.rsplit('.',1)[0]
    size_screen=('below 4 Mb' if total<4_000_000 else 'above 6.3 Mb' if total>6_300_000 else 'inside size range')
    return {'isolate_id':isolate,'filename':specimen,'source_archive':source,
      'assembly_size_bp':total,'contig_count':len(sizes),'n50_bp':n50,
      'gc_percent_acgt':f'{100*gc/max(acgt,1):.4f}', 'ambiguous_N_bp':chars['N'],
      'shortest_contig_bp':min(sizes),'longest_contig_bp':max(sizes),
      'contigs_lt_200':sum(x<200 for x in sizes),'contigs_lt_1000':sum(x<1000 for x in sizes),
      'duplicated_ids':sum(v-1 for v in identifiers.values() if v>1),
      'invalid_characters':invalid,'ncbi_size_screen_only':size_screen,
      'comments':'NOT a taxonomic, contamination, or cohort-membership validation'}


def scan_paths(paths):
    for p in (Path(x) for x in paths):
        if p.is_dir():
            yield from scan_paths(sorted(q for q in p.rglob('*') if q.is_file()))
        elif p.suffix.lower()=='.zip':
            with zipfile.ZipFile(p) as z:
                for info in z.infolist():
                    if info.is_dir() or not info.filename.lower().endswith(SUFFIXES): continue
                    if info.filename.lower().endswith('.gz'):
                        with z.open(info) as raw, gzip.GzipFile(fileobj=raw) as source:
                            yield summarize(read_fasta(source, info.filename), info.filename,p.name)
                    else:
                        with z.open(info) as raw:
                            yield summarize(read_fasta(raw,info.filename),info.filename,p.name)
        elif p.name.lower().endswith(SUFFIXES):
            if p.name.lower().endswith('.gz'):
                with gzip.open(p,'rb') as f: yield summarize(read_fasta(f,p.name),p.name)
            else:
                with p.open('rb') as f: yield summarize(read_fasta(f,p.name),p.name)
        else:
            raise ValueError(f'Unsupported input: {p}')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inputs',nargs='+',help='FASTA files, input directories, or ZIPs')
    ap.add_argument('--out',required=True,help='Output TSV path')
    args=ap.parse_args(argv)
    rows=list(scan_paths(args.inputs))
    seen=Counter(x['isolate_id'] for x in rows)
    if any(n>1 for n in seen.values()):
        raise SystemExit('WARNING: duplicate isolate IDs across input FASTAs; refusing misleading QC table.')
    save_tsv(args.out,HEAD,rows)
    print(f'QC inspected {len(rows)} assemblies; under4Mb={sum(int(x["assembly_size_bp"])<4_000_000 for x in rows)}; over6.3Mb={sum(int(x["assembly_size_bp"])>6_300_000 for x in rows)}',file=sys.stderr)

if __name__=='__main__':main()
