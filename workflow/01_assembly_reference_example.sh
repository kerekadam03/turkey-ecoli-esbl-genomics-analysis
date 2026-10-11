#!/usr/bin/env bash
# RECONSTRUCTED WORKFLOW EXAMPLE; NOT AN ARCHIVE OF HISTORICAL STUDY COMMANDS.
# Usage: bash workflow/01_assembly_reference_example.sh SAMPLE R1.fastq.gz R2.fastq.gz
# The original GAM-NGS merge commands and parameters are not available.
set -euo pipefail
if [[ $# -ne 3 ]]; then echo "Usage: $0 SAMPLE R1.fastq.gz R2.fastq.gz" >&2; exit 2; fi
SAMPLE="$1"; R1="$2"; R2="$3"
: "${THREADS:=16}"; : "${MEM_GB:=64}"
mkdir -p "results/${SAMPLE}/" "logs/" 
for command in fastqc fastp spades.py megahit quast.py; do
  command -v "$command" >/dev/null || { echo "MISSING tool: $command" >&2; exit 3; }
done
FASTQ_OUT="results/${SAMPLE}/clean"
mkdir -p "$FASTQ_OUT"
fastqc --threads "$THREADS" --outdir "$FASTQ_OUT" "$R1" "$R2" \
  2>&1 | tee "logs/${SAMPLE}.fastqc.log"
fastp -i "$R1" -I "$R2" \
  -o "$FASTQ_OUT/${SAMPLE}_R1.clean.fastq.gz" \
  -O "$FASTQ_OUT/${SAMPLE}_R2.clean.fastq.gz" \
  -w "$THREADS" --json "$FASTQ_OUT/fastp.json" \
  --html "$FASTQ_OUT/fastp.html" 2>&1 | tee "logs/${SAMPLE}.fastp.log"
R1C="$FASTQ_OUT/${SAMPLE}_R1.clean.fastq.gz"
R2C="$FASTQ_OUT/${SAMPLE}_R2.clean.fastq.gz"
spades.py -1 "$R1C" -2 "$R2C" -t "$THREADS" -m "$MEM_GB" \
  -o "results/${SAMPLE}/spades" 2>&1 | tee "logs/${SAMPLE}.spades.log"
megahit -1 "$R1C" -2 "$R2C" -t "$THREADS" \
  -o "results/${SAMPLE}/megahit" 2>&1 | tee "logs/${SAMPLE}.megahit.log"
quast.py "results/${SAMPLE}/spades/contigs.fasta" \
  "results/${SAMPLE}/megahit/final.contigs.fa" \
  -o "results/${SAMPLE}/quast_compare" 2>&1 | tee "logs/${SAMPLE}.quast.log"
echo 'Candidate assemblies and comparative QUAST results generated. The selection or merging of final assemblies is outside this reference workflow.'
