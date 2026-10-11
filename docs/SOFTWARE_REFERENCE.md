# Software environment and reference workflows

The following names and versions reflect the methods stated in the manuscript. The commands in `workflow/` illustrate independent analysis stages; they are separate from historical run logs.

| Analysis stage | Software stated in manuscript | Repository resource |
| --- | --- | --- |
| Read QC | FastQC 0.11.9, fastp 0.23.2, Bloocoo 1.0.7 | `workflow/01_assembly_reference_example.sh` (FastQC/fastp stages) |
| Genome assembly and assessment | SPAdes 4.0.0; MEGAHIT 1.2.9; GAM-NGS 1.1b; QUAST 5.2 | `workflow/01_assembly_reference_example.sh` (candidate assembly comparison; does not merge assemblies) |
| Resistance gene annotation | RGI 5.1.0; ABRicate; CARD | `workflow/02_resistome_and_mobility_reference.sh`; `scripts/normalize_abricate.py` |
| Mobile sequence context | MobileElementFinder 1.0.3; PlasFlow 1.1 | `workflow/02_resistome_and_mobility_reference.sh`; `scripts/mge_context_summary.py` |
| Comparative genomics | Panaroo 1.6.0 | `workflow/03_pangenome_phylogeny_reference.sh`; `scripts/audit_pangenome_phylogeny.py` |
| Phylogenetics | IQ-TREE 3.0.1, ModelFinder, iTOL v6 | `workflow/03_pangenome_phylogeny_reference.sh` and `scripts/make_itol_annotation_tracks.py` |
| Statistical processing | R 4.5.3 (reported methods); Python 3.10+ (repository) | `scripts/` for independent descriptive input checking |

For reproducible bioinformatics work, save the exact command lines, software versions, reference database version/date, checksums of inputs, per-step logs, and any imposed thresholds. A value in an illustrative `.env` file is a **configuration example**, not evidence that the same value was used for the historical study.

Do not substitute genome-size filtering for taxonomic or assembly-quality assessment. In a contamination investigation, investigate read-level evidence and assembly identity.
