# Software verification

Checked on 2026-10-11 before repository packaging.

## Automated checks

- `python -m compileall -q scripts tests`: Python syntax compilation.
- `python -m unittest discover -s tests -v`: nine automated tests covering parsers, sample matching, MIC calculations, MGE counts, gene annotation summaries, and FASTA QC.
- `make demo`: small independent example-input workflow.
- `bash -n workflow/*.sh` (each file): shell script syntax.
- Example plotting with `matplotlib` for Figures 1, 2, 3, 4, 6, and 7.

The test fixtures evaluate code behavior. To independently reproduce reported study figures, use the actual validated data and the originally documented analytical settings as described in the methods and accompanying data records.
