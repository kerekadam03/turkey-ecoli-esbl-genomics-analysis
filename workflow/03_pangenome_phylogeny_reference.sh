#!/usr/bin/env bash
# RECONSTRUCTED reference, NOT the original Panaroo/IQ-TREE 3.0.1 invocation.
# Must be run only after 160-isolate cohort + external reference accessions
# and sequence annotations are independently verified.
# Usage: bash workflow/03_pangenome_phylogeny_reference.sh input-gff-dir
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: $0 DIRECTORY_OF_GFF3" >&2; exit 2; fi
for tool in panaroo;do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 3; };done
if command -v iqtree3 >/dev/null;then IQTREE=iqtree3
elif command -v iqtree >/dev/null;then IQTREE=iqtree
else echo 'Missing IQ-TREE executable' >&2;exit 3;fi
: "${THREADS:=16}"; : "${CORE_THRESHOLD:=0.95}"
mkdir -p results/panaroo logs results/phylogeny
shopt -s nullglob
GFFS=("$1"/*.gff "$1"/*.gff3)
if [[ ${#GFFS[@]} -lt 2 ]];then echo 'At least two verified GFF3 inputs required' >&2; exit 2;fi
# The choice of 95% core threshold is a NEW REANALYSIS CHOICE unless
# an original Panaroo command/log confirms it. Historical text only reports
# 99% 'strict core' and 95-99% 'soft core' counts; not the alignment threshold.
panaroo -i "${GFFS[@]}" -o results/panaroo/ --clean-mode strict \
  -a core --aligner mafft --core_threshold "$CORE_THRESHOLD" \
  -t "$THREADS" --remove-invalid-genes \
  2>&1 | tee logs/panaroo_reanalysis.log
ALN=results/panaroo/core_gene_alignment.aln
[[ -s "$ALN" ]] || { echo 'Panaroo core alignment absent or empty; inspect the logs.' >&2;exit 4;}
# NEW SUPPORT-BASED ANALYSIS: the manuscript's original -fast command did not
# explicitly request UFBoot or SH-aLRT; do not say support values were obtained
# previously. This command is designed to address reviewer feedback.
"$IQTREE" -s "$ALN" -m MFP -B 1000 -alrt 1000 -T "$THREADS" \
  --prefix results/phylogeny/core_tree_with_support \
  2>&1 | tee logs/iqtree_reanalysis_with_support.log
echo 'NEW phylogenetic analysis generated; confirm alignment, recombination, cluster definitions and manuscript cohort.' >&2
