#!/usr/bin/env bash
# Policy check; never changes a FASTA. Not the actual original NCBI correction code.
set -euo pipefail
if [[ $# -lt 2 ]]; then echo "Usage: $0 OUT_QC_TSV ASSEMBLY.zip [MORE_ASSEMBLIES...]" >&2;exit 2;fi
OUT="$1";shift
python scripts/fasta_qc.py "$@" --out "$OUT"
echo 'Never drop contigs solely to force assembly size into a target range.' >&2
echo 'Investigate taxonomic mixtures, contamination, true read depth, plasmid content, and reassembly.' >&2
