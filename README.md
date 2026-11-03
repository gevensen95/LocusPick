# LocusPick

Transparent selection of representative loci for genome browser figures.

Example tracks in papers are usually picked by eye. LocusPick makes that
choice reproducible and auditable:

1. **Link** binding peaks (ChIP-seq, CUT&RUN, ATAC-seq) to genes, or check an
   existing peak-to-gene annotation against your gene models.
2. **Group** genes by two differential expression contrasts (e.g. a factor
   knockout × aging), giving activated/repressed × up/down groups.
3. **Check** each candidate for common pitfalls: distant or inconsistent gene
   assignments, low mappability (paralog clusters), gaps in signal coverage and
   outdated gene symbols.
4. **Rank** candidates with explicit weights you can state in a legend.
5. **Test** whether a change in binding is specific to the group or genome-wide.
6. **Plot** multi-track loci with gene models and expression insets.

Every candidate and flag is kept in a table for supplementary material.

## Installation

```r
# install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
remotes::install_github("gevensen95/LocusPick", build_vignettes = TRUE)
```

## Quick start

```r
library(LocusPick)
ex <- locusPickExample()

cand   <- buildCandidates(ex$peaks, ex$genes, ex$groups,
                          signal = ex$signal, mappability = ex$mappability)
qcSummary(cand)
ranked <- rankLoci(cand, c(A_mean = 1, abs_lfc_perturbation = 1,
                           abs_lfc_condition = 1))
top    <- topLoci(ranked, n = 3)
plotLocus(top[1, ], ex$signal, ex$genes,
          trackGroups = list(Binding = c("A", "B")),
          expression = ex$expression)
```

See `vignette("LocusPick")` for the full workflow and
`vignette("preparing-inputs")` for ways to produce each input from BAM files,
peak callers and DE tools.

## Inputs

| What | Accepted |
|---|---|
| Peaks | BED / narrowPeak path, `data.frame`, `GRanges` |
| Gene models | GTF/GFF path or imported `GRanges` |
| Signal | bigWig (path or URL), in-memory `GRanges`, chromDict, or any function |
| Mappability | Umap bigWig (`inst/scripts/make_umap_bw.sh` builds one) |
| Expression | DE tables (any column names, via `deTable()`); long table for insets |
