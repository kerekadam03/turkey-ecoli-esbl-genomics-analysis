# Turkey-derived *Escherichia coli*: ESBL phenotypes and genomic context analysis

Research code, quality-control utilities, and documented analysis workflows associated with the *One Health* manuscript on extended-spectrum β-lactamase (ESBL)-associated resistance in turkey-derived *Escherichia coli* from Hungary.

**Study:** Phenotypic ESBL profiles and genomic context of *bla*CTX-M in turkey-derived *Escherichia coli* (manuscript in revision)  
**Journal:** *One Health*  
**Manuscript:** ONEHLT-D-26-01460  
**NCBI BioProject:** [PRJNA1475617](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1475617)

## Overview

This repository provides computational tools for processing and checking phenotypic antimicrobial susceptibility data, genomic annotation summaries, predicted mobile-element contexts, sequencing-assembly quality metrics, and phylogenetic metadata. Its purpose is to facilitate reproducible data processing, independent validation, and transparent reporting of analysis settings.

The code distinguishes between directly recorded observations and computational inferences. In particular, plasmid-origin prediction and proximity to mobile genetic elements are reported as genomic context annotations; they are not presented as experimental evidence of horizontal gene transfer.

## Repository contents

| Location | Purpose |
| --- | --- |
| [`scripts/`](scripts/) | Python utilities for isolate matching, MIC ratios, ARG summaries, MGE-context summaries, pangenome/phylogeny checks, and figures |
| [`workflow/`](workflow/) | Documented example commands for sequence assembly, resistance annotation, mobility analysis, and phylogenetic reconstruction |
| [`config/`](config/) | Input settings and comparison targets |
| [`data/example/`](data/example/) | Small demonstration datasets for running and testing the tools; these are **not** experimental observations |
| [`reference_metadata/`](reference_metadata/) | Reference-genome identifiers retained from comparative tree metadata |
| [`tests/`](tests/) | Automated Python tests |
| [`docs/`](docs/) | Input formats, software and workflow provenance, and code quality-control guidance |

## Getting started

**Requirements:** Python 3.10 or newer. Core analysis utilities use Python's standard library. For optional plots, install the packages in `requirements-optional.txt`.

Run the included automated tests from the repository's root directory:

```bash
python -m unittest discover -s tests -v
```

Run a demonstration of the paired cefotaxime/cefotaxime–clavulanate MIC analysis:

```bash
python scripts/analysis_esbl.py \
  --mic data/example/mic_pairs.tsv \
  --out results/example_esbl.tsv
```

Summarize *bla*CTX-M sequence-context annotations:

```bash
python scripts/mge_context_summary.py \
  --input data/example/ctxm_context.tsv \
  --out results/example_ctxm.tsv
```

Generate a demonstration susceptibility plot (requires `matplotlib`):

```bash
python scripts/plot_reanalysis.py \
  --figure 1 \
  --input data/example/figure1_categories.tsv \
  --output results/example_figure1.png
```

All files under `data/example/` are illustrative. Figures and numerical summaries intended for a manuscript must be regenerated using the appropriate checked experimental data and documented interpretation criteria.

## Input data and reproducibility

Required input columns are documented in [`docs/INPUT_TSV_SCHEMA.md`](docs/INPUT_TSV_SCHEMA.md). The scripts accept explicitly identified study isolates and recorded data rather than generating missing measurements from manuscript-level percentages.

A complete re-execution of the genome-analysis workflow additionally requires the relevant sequence inputs, software and database versions, assembly configurations, annotation parameters, and phylogenetic input files. The shell scripts under `workflow/` are **documented reference implementations** for independent reanalysis, not archived execution logs of the original sequencing runs. The distinction between available utilities and original laboratory or bioinformatics records is described in [`docs/PROVENANCE.md`](docs/PROVENANCE.md).

Sequence assemblies and associated sample records should be retrieved from their verified NCBI accessions where publicly available. Isolate inclusion and metadata must be matched to the analysis cohort before combining phenotype and genotype data. Do not upload identifiable or restricted sampling information to this repository.

## Analysis principles

- **Susceptibility interpretation:** clinical breakpoints and epidemiological cutoffs are handled as distinct interpretive frameworks.
- **Phenotypic ESBL testing:** paired MIC data and the laboratory's independently determined phenotype are retained separately from gene detection.
- **Resistance genes:** gene annotations distinguish acquired determinants from intrinsic, efflux-associated, and regulatory categories when the curated source calls support this distinction.
- **Plasmid and mobile-element context:** predicted contig origin, gene proximity, and observed sequence context are distinguished from experimentally confirmed transfer.
- **Phylogenetic analysis:** tree interpretation requires documented genome membership, reference sequences, alignments, and support measures.

## Citation and use

For use of the computational resources, cite the version or commit of this repository and the associated research article after its publication. Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

This repository provides research software and analysis documentation; the accompanying experimental datasets and formal study conclusions are described in the manuscript and corresponding verified data deposits. No software-use license is implied unless one is added by the repository owner.
