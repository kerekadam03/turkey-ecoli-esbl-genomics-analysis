.PHONY: test demo clean
PYTHON ?= python3

test:
	$(PYTHON) -m compileall -q scripts tests
	$(PYTHON) -m unittest discover -s tests -v

demo:
	mkdir -p results
	$(PYTHON) scripts/analysis_esbl.py --mic data/example/mic_pairs.tsv --out results/demo_esbl.tsv
	$(PYTHON) scripts/mge_context_summary.py --input data/example/ctxm_context.tsv --out results/demo_ctxm.tsv
	$(PYTHON) scripts/build_figure_inputs.py --cohort data/example/sampling_metadata.tsv --phenotypes data/example/mic_pairs.tsv --calls data/example/reconciled_gene_calls.tsv --complete data/example/annotation_completeness.tsv --figure3 results/demo_fig3_input.tsv --figure6 results/demo_fig6_input.tsv
	$(PYTHON) scripts/gene_frequency_by_esbl.py --cohort data/example/sampling_metadata.tsv --phenotypes data/example/mic_pairs.tsv --calls data/example/reconciled_gene_calls.tsv --complete data/example/annotation_completeness.tsv --out results/demo_gene_by_esbl.tsv
	$(PYTHON) scripts/phenotype_gene_crosscheck.py --cohort data/example/sampling_metadata.tsv --phenotypes data/example/mic_pairs.tsv --calls data/example/reconciled_gene_calls.tsv --complete data/example/annotation_completeness.tsv --out results/demo_genotype_phenotype.tsv
	$(PYTHON) scripts/sampling_summary.py --metadata data/example/sampling_metadata.tsv --phenotypes data/example/mic_pairs.tsv --out results/demo_sampling.tsv
	$(PYTHON) scripts/audit_manuscript_claims.py --claims-json config/manuscript_claims_to_verify.json --cohort data/example/sampling_metadata.tsv --phenotypes data/example/mic_pairs.tsv --ctxm-summary results/demo_ctxm.tsv --references reference_metadata/legacy_reference_tips_normalized.tsv --out results/demo_claim_check.tsv
	@echo 'Synthetic demo completed; these are NOT study results.'

clean:
	rm -rf results/*
