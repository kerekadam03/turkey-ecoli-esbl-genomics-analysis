#!/usr/bin/env bash
# RECONSTRUCTED reference commands, NOT verified historical runs.
# Example: bash workflow/02_resistome_and_mobility_reference.sh SAMPLE assembly.fasta
# Requires explicitly installed tools and local databases; inspect versions first.
set -euo pipefail
if [[ $# -ne 2 ]]; then echo "Usage: $0 SAMPLE assembly.fasta" >&2; exit 2; fi
SAMPLE="$1"; ASSEMBLY="$2"; : "${THREADS:=8}"; : "${PLASFLOW_THRESHOLD:=0.7}"
for tool in rgi abricate;do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 3; };done
mkdir -p "results/$SAMPLE/resistome" "logs"
# CARD RGI prediction. RGI Strict is a CARD-specific model label, NOT a 90% BLAST filter.
rgi main --input_sequence "$ASSEMBLY" \
  --output_file "results/$SAMPLE/resistome/rgi_card" \
  --input_type contig --alignment_tool DIAMOND --threads "$THREADS" --local \
  2>&1 | tee "logs/${SAMPLE}.rgi.log"
# Independent nucleotide alignment screen. Do not silently union it with RGI calls.
abricate --db card --minid 90 --mincov 90 "$ASSEMBLY" \
  > "results/$SAMPLE/resistome/abricate_CARD.tsv"
# PlasFlow classification is a contig-level prediction. 0.7 is the tool default,
# NOT demonstrated to be the cutoff originally used in this study.
if command -v PlasFlow.py >/dev/null 2>&1; then
  PlasFlow.py --input "$ASSEMBLY" \
    --output "results/$SAMPLE/resistome/plasflow.tsv" \
    --threshold "$PLASFLOW_THRESHOLD" \
    2>&1 | tee "logs/${SAMPLE}.plasflow.log"
else
  echo 'PlasFlow not installed; contig classification not performed.' >&2
fi
# MobileElementFinder local database must be installed separately.
if command -v mefinder >/dev/null 2>&1; then
  mefinder find --contig "$ASSEMBLY" "results/$SAMPLE/resistome/mobileelementfinder" \
    2>&1 | tee "logs/${SAMPLE}.mefinder.log"
else
  echo 'MobileElementFinder not installed; skipping. Manual MGE distances require actual output.' >&2
fi
echo 'Annotation outputs generated. Confirm database versions, hit reconciliation and genomic coordinates before interpretation.'
