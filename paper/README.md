# LocusPick manuscript analyses

Scripts that produce Figs. 2–7 of the LocusPick methods paper. They are not
part of the R package (`paper/` is in `.Rbuildignore`).

Run everything from this folder: `cd paper`.

## Test mode first

```bash
LOCUSPICK_PAPER_TEST=1 bash run_all.sh      # ~5 min, no downloads
```

Test mode builds a small simulated genome with every input laid out like the
real files (ENCODE manifest, nf-core/chipseq output folder, Hnf4a tables,
Seurat/Signac objects) and runs all six figure scripts on it. The simulation
contains each pitfall the scripts look for, so every code path runs. Results go
to `results_test/`, and every figure is stamped *TEST MODE*. **Do not report
those numbers.**

## Real data: order of steps

| Step | Command | Needs | Produces |
|---|---|---|---|
| 1 | `Rscript 00_select_encode.R` | internet (encodeproject.org) | `data/encode/manifest.tsv` |
| 2 | `bash 00_download.sh` | curl | references, ENCODE files, hg38+dm6 genome |
| 3 | fill `templates/gse60104_design.csv`, then `bash 00_process_gse60104.sh` | nextflow, singularity, samtools | `data/gse60104/chipseq/` |
| 4 | download GSE215915, fill `templates/gse215915_samples.csv`, then `Rscript 00_process_pbmc.R` | Azimuth, EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38 | `data/gse215915/pbmc_multiome.rds` |
| 5 | score papers into `templates/audit_rubric.csv` (see `templates/AUDIT_GUIDE.md`) | two raters | — |
| 6 | `bash run_all.sh` (or one script at a time) | | `results/fig2` … `results/fig7` |

Fig. 7 (and the Hnf4a parts of Figs. 4–5) use the existing files under
`/vast/igc/analyses/garrett/P01/kat` (paths in `config.R`) and need the
chromDicts package to read the `.rds` tracks.

## Scripts

| Script | Figure | Question |
|---|---|---|
| `01_fig2_audit.R` | 2 | How do published papers choose and report example loci? (two-rater audit, Cohen's κ) |
| `02_fig3_shortcuts.R` | 3 | How often do shortcut rules (strongest peak, largest fold change, largest signal change, random) pick a locus that fails a basic check? ENCODE HepG2 vs K562 |
| `03_fig4_global_shift.R` | 4 | Global shifts: ChIP-Rx H3K79me2 ± DOT1L inhibitor (CPM vs spike-in) and Hnf4a CUT&RUN young vs old |
| `04_fig5_robustness.R` | 5 | Are the picks stable to replicates, weights, distance cap and normalization? How typical are they? |
| `05_fig6_singlecell.R` | 6 | Single-cell multiome (PBMC): links vs nearest gene, pseudobulk vs cell-level DE under a null, paired vs unpaired, donor numbers |
| `06_fig7_case_study.R` | 7 | Hnf4a aging liver: bulk four-group selection, the original hand-picked genes, single-nucleus per cell type |

`config.R` holds every path and setting; `R/common.R` the shared helpers;
`R/simulate.R` the test data.

## Placeholders and things to fill in

* **Hnf4a CUT&RUN replicate bigWigs** — `HNF4A$cnr_reps` in `config.R` points
  at `PLACEHOLDER_*` files. Until they exist, Fig. 5 skips the Hnf4a
  leave-one-replicate-out analysis (with a message) and uses the chromDict
  tracks for everything else.
* **`HNF4A$manual_picks`** — genes shown in the original hand-picked panels,
  named by group, e.g. `c(repressed_age_up = "Saa1")`. Fig. 7 reports their
  rank, percentile and QC flags.
* **sn column names** — `HNF4A$sn_cols` (same as the earlier test script).
* **`templates/gse60104_design.csv`** — one row per library: `gsm`, `srr`
  (several runs separated by `;`), `sample` (nf-core sample name, e.g.
  `DMSO`; the pipeline appends `_REP1`), `replicate`, `condition` (must match
  `CHIPRX$test` / `CHIPRX$ref` in `config.R`), `antibody`, `is_input`,
  `control`, `control_replicate`. nf-core/fetchngs no longer accepts GEO IDs,
  so take SRR IDs from the SRA Run Selector.
* **`templates/gse215915_samples.csv`** — `donor`, `severity`
  (`decompensated`/`stable`), `h5`, `fragments`.

## Not tested here

The `00_*` steps could not be run where the scripts were written (no access to
ENCODE, GEO or SRA). They check their inputs and stop with a message rather
than continuing on bad data, but expect to adjust them on first run:

* ENCODE portal field names in `00_select_encode.R`;
* the hg38 Umap bigWig URL (fallback: build it as for mm10);
* the nf-core/chipseq output paths used by Fig. 4 (`bwa/mergedLibrary/…`);
* whether GSE215915 includes fragment files for every sample.

## What each figure is expected to show, and what would change the story

* **Fig. 3.** LocusPick's own checks are applied to every rule, so LocusPick
  is ~0% on those four by construction. The informative numbers are the
  shortcut rates and the replicate-reproducibility column, which LocusPick
  does not filter on.
* **Fig. 4.** A per-condition scale factor shifts every log ratio by a
  constant, so the comparison against all peaks (`shift_vs_background`) and
  ratio-based rankings should barely move between CPM and spike-in, while
  "is there a loss?" and difference-based rankings can. If the real data
  disagree, the per-replicate spike-in factors are the first thing to check.
* **Fig. 6B.** Under permuted donor labels pseudobulk DE should find ~0
  genes and cell-level tests many; if pseudobulk also finds many, the
  permutation is not breaking the design (check donors vs libraries).
