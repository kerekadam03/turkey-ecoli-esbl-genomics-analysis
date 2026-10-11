# Code provenance and research data

## Scope

This repository contains a documented Python toolkit for processing study inputs, checking data consistency, and generating descriptive plots. The Python tools and their example inputs were developed as research-support software in October 2026 and have automated tests.

The `workflow/` shell scripts document **reference commands for independent reanalysis** with standard bioinformatics tools. They are not original laboratory machine logs and should not be cited as proof of the precise commands historically used for the published assembly and annotation datasets.

## Research inputs and provenance

| Component | Sources needed for independent reanalysis | Repository role |
| --- | --- | --- |
| Laboratory MICs | Paired cefotaxime and cefotaxime–clavulanate values plus independent laboratory classifications | `scripts/analysis_esbl.py` computes ratios and reports screening status |
| Resistance gene identification | Original per-contig gene annotations and database identifiers | Normalization and descriptive summaries; phenotype and genotype kept separate |
| Genomic mobility context | Named contigs, gene positions, prediction scores, MGE annotations and coordinates | Reports plasmid-predicted origin and proximity statistics |
| Assembly quality | FASTA assemblies or archive files | Contig count, total assembly size, N50, GC% and sequence-format checks |
| Comparative phylogenomics | Genome set, reference accession list, Panaroo alignment and IQ-TREE outputs | Input-quality and alignment checks plus annotation layers |
| Cohort and metadata | Curated isolate inclusion table, sample identifiers and sampling metadata | Exact-match reconciliation and descriptive strata |

The original wet-laboratory data and sequencing records are distinct from source code. Genome/BioSample accessions should be checked against the publicly released NCBI records before citing a specific version. Gene-context inferences from short-read sequencing do not establish horizontal transfer.

## Versions and environment

The software versions described in the manuscript and the independent example configurations are listed in [`SOFTWARE_REFERENCE.md`](SOFTWARE_REFERENCE.md). Exact original execution flags and database snapshots should be documented from instrument and pipeline records when the complete historical analysis is archived.

## Demonstration data

Files in `data/example/` are synthetic fixtures used to verify the command-line interface. No biological conclusions or manuscript tables should be drawn from them. The repository does not include restricted raw reads or identifiable farm/flock source records.
