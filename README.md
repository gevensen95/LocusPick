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

Pass your tools' output as it is:

```r
library(LocusPick)
res <- locusPick(
    peaks        = c("rep1_peaks.narrowPeak", "rep2_peaks.narrowPeak"),  # consensus
    perturbation = DESeq2::results(dds_ko, contrast = c("genotype", "KO", "WT")),
    condition    = "old_vs_young_edgeR.csv", conditionDirection = "old/young",
    genes        = "gencode.vM25.annotation.gtf",
    orgDb        = org.Mm.eg.db::org.Mm.eg.db,
    signal       = list(CNR_young = c("y1.bw", "y2.bw"), CNR_old = c("o1.bw", "o2.bw")),
    mappability  = "k50.Umap.MultiTrackMappability.bw")
res                          # checks, QC, signal-change specificity, top loci
plot(res, "Rad51b")
reportLoci(res, "review.pdf")
```

Try it on the bundled example data:

```r
ex <- locusPickExample()
cand   <- buildCandidates(ex$peaks, ex$genes, ex$groups,
                          signal = ex$signal, mappability = ex$mappability)
ranked <- rankLoci(cand, c(A_mean = 1, abs_lfc_perturbation = 1, abs_lfc_condition = 1))
plotLocus(topLoci(ranked, n = 1)[1, ], ex$signal, ex$genes,
          trackGroups = list(Binding = c("A", "B")), expression = ex$expression)
```

Single-cell multiome or scRNA + scATAC data:

```r
res <- scLocusPick(seurat_obj, cellType = "cell_type", sample = "sample",
                   perturbation = c("genotype", "KO", "WT"),
                   condition = c("age", "old", "young"),
                   genes = "gencode.vM25.annotation.gtf")
```

Either contrast can be omitted, or both (rank on signal alone).

See `vignette("LocusPick")` for the full workflow, `vignette("preparing-inputs")`
for every supported input format, and `vignette("single-cell")` for single-cell
data.

## Documentation

Rendered vignettes, readable without installing the package:

| Vignette | Markdown | PDF |
|---|---|---|
| Choosing representative loci | [LocusPick.md](docs/vignettes/LocusPick.md) | [PDF](docs/vignettes/LocusPick.pdf) |
| Preparing input data | [preparing-inputs.md](docs/vignettes/preparing-inputs.md) | [PDF](docs/vignettes/preparing-inputs.pdf) |
| Single-cell data | [single-cell.md](docs/vignettes/single-cell.md) | [PDF](docs/vignettes/single-cell.pdf) |

## Inputs

| What | Accepted directly |
|---|---|
| DE results | DESeq2, edgeR, limma objects or tables; any CSV/TSV; Ensembl/Entrez IDs (with an OrgDb) |
| Peaks | MACS2/MACS3, Genrich, SEACR, HOMER, DiffBind, nf-core consensus, BED; replicates merged to consensus |
| Signal | bigWig or BAM, one or several per condition (averaged); nf-core results folders; chromDicts; any function |
| Gene models | GTF/GFF path or `GRanges` |
| Expression insets | `DESeqDataSet`, `DGEList`, `SummarizedExperiment`, or a matrix |
| Mappability | Umap bigWig (`inst/scripts/make_umap_bw.sh` builds one) |
| Single cell | Seurat + Signac, `SingleCellExperiment`, Cell Ranger / ATAC / ARC folders, matrices; fragments files; peak-gene links |
