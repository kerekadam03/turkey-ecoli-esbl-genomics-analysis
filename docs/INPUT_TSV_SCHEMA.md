# Minimal normalized input TSV schemas

All files use UTF-8 TSV headers, one row per record and stable **isolate_id**.
Missing values must be blank or explicitly `Unknown` — never invented.

## `data/private/verified_study_cohort.tsv`

`isolate_id` (160 distinct verified IDs; requires independent proof of inclusion)

## `data/private/mic_pairs.tsv`

`isolate_id`, `cefotaxime_mic_ug_ml`, `cefotaxime_clavulanate_mic_ug_ml`,
`initial_esbl_screen_positive`, `laboratory_esbl_phenotype`.

For the MIC ratio script, non-censored positive MICs are required for exact
ratios. `<=`, `>=` or values outside assay dilution ranges are marked unevaluable;
**do not convert without modelling interval censoring**.
The `initial_esbl_screen_positive` field is independently measured and not
inferred from WGS. Document whether an initial cefotaxime screen at >=2 μg/mL
was performed and which standard governed it.

## `data/private/gene_calls_curated.tsv`

`isolate_id`, `gene`, `call_class` (`acquired`,`intrinsic`,`efflux`,`regulatory`,`other`,`uncertain`),
`status` (`verified`), `plasmid_predicted`, `mge_proximal` (`Yes`,`No`,`Unknown`).

This requires an independently checked reconciliation of the original CARD
RGI and ABRicate outputs. Do not fill `call_class` by assuming all CARD genes
are clinically acquired ARGs.

## `data/private/annotation_completeness.tsv`

`isolate_id`, `annotation_complete` (`Yes`, `No`, `Unknown`).
Only verified complete annotation supports absence/not detected claims.

## `data/private/ctxm_context.tsv`

`isolate_id`, `gene`, `contig_id`, `plasmid_predicted`, `mge_element`,
`mge_hit_id`, `distance_bp`.

For several elements near the same *bla* gene, repeat rows with distinct
`mge_hit_id`. The MGE `distance_bp` should be calculated from documented,
contig-specific coordinates using the actual database rule. A blank MGE name
means **no described hit**, not a negative experimental transfer assay.

## `data/private/sampling_metadata.tsv`

`isolate_id`, `sample_id`, `bird_id`, `sampling_site`, `farm_id`, `flock_id`,
`region`, `age_group`, `production_type`.

If a person-invariant or farm-coded link exists, anonymize it before release.
Unknown farm/flock IDs must be clearly marked as missing, not replaced by
settlement names or flock-size bands.

## `data/private/figure1_categories.tsv`

`drug`, `category` (Susceptible / Intermediate / Resistant), `n`, `N`,
`interpretation_type` (`clinical` only is plotted in Figure 1).

**No cefotaxime+clavulanate column as a therapeutic S/I/R result.**
Include ECOFF WT/NWT results in a different output, not overlaid as S/I/R.

## `data/private/reference_genomes.tsv`

`accession`, `source`, `host`, `selection_rationale`, `download_date`, `version`.
The historical GCA list alone is not a complete accession provenance record.

## `data/private/panaroo/gene_presence_absence.csv`

Official Panaroo output with `No. isolates` and a validated number of input
genomes; the `--core_threshold` of the *original analysis* must be independently established.
