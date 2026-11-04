sc <- locusPickSingleCellExample()

test_that("pseudobulk sums counts per sample and cell type", {
    pb <- pseudobulk(sc$rna, sample = "sample", cellType = "cell_type", meta = sc$meta)
    expect_s3_class(pb, "LocusPickPseudobulk")
    expect_equal(ncol(pb$counts), 16)
    expect_equal(sum(pb$counts), sum(sc$rna))
    cells <- rownames(sc$meta)[sc$meta$sample == "WT_young_1" & sc$meta$cell_type == "Kupffer"]
    expect_equal(unname(pb$counts[, "Kupffer|WT_young_1"]),
                 unname(Matrix::rowSums(sc$rna[, cells])))
    expect_true(all(c("genotype", "age") %in% names(pb$colData)))   # sample-level columns kept
    expect_message(pseudobulk(sc$rna, "sample", "cell_type", meta = sc$meta, minCells = 30),
                   "dropped 16")
})

test_that("pseudobulkDE recovers the simulated effects and requires replicates", {
    skip_if_not_installed("DESeq2")
    de <- suppressMessages(pseudobulkDE(sc$rna, c("genotype", "KO", "WT"), sample = "sample",
                                        cellType = "cell_type", meta = sc$meta,
                                        covariates = "age"))
    expect_s3_class(de, "LocusPickDEList")
    h <- de$Hepatocyte
    expect_lt(h$lfc[h$gene == "geneA"], -1)
    expect_lt(h$padj[h$gene == "geneA"], 0.01)
    expect_gt(h$lfc[h$gene == "geneB"], 0.8)
    expect_equal(attr(h, "contrast"), "genotype KO vs WT")
    expect_false(any(de$Kupffer$padj[de$Kupffer$gene %in% c("geneA", "geneB")] < 0.05))
    expect_warning(expect_message(
        skip <- pseudobulkDE(sc$rna, c("genotype", "KO", "WT"), sample = "sample",
                             cellType = "cell_type", meta = sc$meta, minSamples = 5),
        "skipped"), "enough samples")
    expect_length(skip, 0)
})

test_that("pseudobulkDE works with edgeR", {
    skip_if_not_installed("edgeR")
    de <- suppressMessages(pseudobulkDE(sc$rna, c("age", "old", "young"), sample = "sample",
                                        cellType = "cell_type", meta = sc$meta,
                                        method = "edgeR", cellTypes = "Hepatocyte"))
    expect_named(de, "Hepatocyte")
    expect_lt(de$Hepatocyte$lfc[de$Hepatocyte$gene == "geneA"], -0.5)
})

test_that("fragmentSource reads the selected cells and normalizes", {
    skip_if_not_installed("Rsamtools")
    hep <- rownames(sc$meta)[sc$meta$cell_type == "Hepatocyte"]
    raw <- fragmentSource(sc$fragments, cells = hep, normalize = "none")
    whole <- GenomicRanges::GRanges("chr1:1-60000")
    s <- fetchSignal(raw, whole)
    fr <- utils::read.delim(sc$fragments, header = FALSE, comment.char = "#")
    fr <- fr[fr$V4 %in% hep, ]
    expect_equal(sum(GenomicRanges::width(s) * s$score), sum(fr$V3 - fr$V2))
    expect_message(cpm <- fragmentSource(sc$fragments, cells = hep), "Counting")
    s2 <- fetchSignal(cpm, whole)
    expect_equal(sum(GenomicRanges::width(s2) * s2$score),
                 sum(fr$V3 - fr$V2) * 1e6 / nrow(fr))
    given <- fragmentSource(sc$fragments, cells = hep, total = nrow(fr))
    expect_equal(signalAtPeaks(sc$peaks, list(a = cpm))$a_mean,
                 signalAtPeaks(sc$peaks, list(a = given))$a_mean)
    expect_error(fragmentSource(f <- tempfile()), "not found")
})

test_that("scSignal builds one track per group", {
    skip_if_not_installed("Rsamtools")
    tr <- suppressMessages(scSignal(NULL, c("cell_type", "age"), fragments = sc$fragments,
                                    meta = sc$meta))
    expect_setequal(names(tr), c("Hepatocyte_young", "Hepatocyte_old",
                                 "Kupffer_young", "Kupffer_old"))
    s <- signalAtPeaks(sc$peaks, tr)
    expect_gt(s$Hepatocyte_young_mean[1], 1.5 * s$Hepatocyte_old_mean[1])    # simulated loss
    expect_gt(s$Kupffer_young_mean[7], 5 * s$Hepatocyte_young_mean[7])     # Kupffer-specific
    expect_equal(attr(tr, "groups")$n_cells, rep(100L, 4))
})

test_that("read10x reads a Cell Ranger ARC outs folder", {
    d <- file.path(tempdir(), "arc", "outs")
    m <- file.path(d, "filtered_feature_bc_matrix")
    dir.create(m, recursive = TRUE, showWarnings = FALSE)
    peakNames <- paste0("chr1:", c(4800, 19800), "-", c(5300, 20300))
    cnt <- rbind(as.matrix(sc$rna[1:3, 1:5]), matrix(1:10, 2))
    Matrix::writeMM(Matrix::Matrix(cnt, sparse = TRUE), file.path(m, "matrix.mtx"))
    writeLines(colnames(sc$rna)[1:5], file.path(m, "barcodes.tsv"))
    utils::write.table(data.frame(c(paste0("ENSG", 1:3), peakNames),
                                  c(rownames(sc$rna)[1:3], peakNames),
                                  c(rep("Gene Expression", 3), rep("Peaks", 2))),
                       file.path(m, "features.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = FALSE)
    file.copy(sc$fragments, file.path(d, "atac_fragments.tsv.gz"), overwrite = TRUE)
    file.copy(paste0(sc$fragments, ".tbi"), file.path(d, "atac_fragments.tsv.gz.tbi"),
              overwrite = TRUE)
    x <- suppressMessages(read10x(dirname(d)))
    expect_s3_class(x, "LocusPick10x")
    expect_equal(dim(x$rna), c(3L, 5L))
    expect_equal(rownames(x$atac), peakNames)
    expect_equal(length(x$peaks), 2)
    expect_match(x$fragments, "atac_fragments")
    pb <- suppressMessages(pseudobulk(x, sample = "sample", meta = sc$meta, minCells = 1))
    expect_equal(sum(pb$counts), sum(cnt[1:3, ]))
})

test_that("readDE warns about cell-level FindMarkers output", {
    fm <- data.frame(p_val = 1e-10, avg_log2FC = 1.2, pct.1 = 0.5, pct.2 = 0.2,
                     p_val_adj = 1e-8, row.names = "geneA")
    expect_warning(readDE(fm, direction = "KO/WT", verbose = FALSE), "inflated")
})

test_that("scLocusPick runs on matrices (multiome-style, with links)", {
    skip_if_not_installed("DESeq2"); skip_if_not_installed("Rsamtools")
    res <- suppressMessages(scLocusPick(
        sc$rna, atac = NULL, cellType = "cell_type", sample = "sample",
        perturbation = c("genotype", "KO", "WT"), condition = c("age", "old", "young"),
        genes = sc$genes, peaks = sc$peaks, links = sc$links,
        fragments = sc$fragments, meta = sc$meta, atacMeta = sc$meta))
    expect_s3_class(res, "LocusPickResult")
    expect_equal(res$cellTypes, "Hepatocyte")       # no DE genes in Kupffer cells
    top <- res$top[res$top$group == "Hepatocyte | activated_age_down", ]
    expect_equal(top$gene[1], "geneA")
    expect_true("link_score" %in% names(res$weights))
    expect_setequal(unique(res$expression$group), c("WT young", "WT old", "KO young", "KO old"))
    p <- plot(res, "geneA")
    expect_s3_class(p, "patchwork")
    pdf <- tempfile(fileext = ".pdf")
    expect_true(all(suppressMessages(reportLoci(res, pdf))))
})

test_that("scLocusPick runs on unpaired data and with one or no contrast", {
    skip_if_not_installed("DESeq2"); skip_if_not_installed("Rsamtools")
    # 'unpaired': ATAC cells are a different set of barcodes with their own metadata
    atacMeta <- sc$meta
    res <- suppressMessages(scLocusPick(
        sc$rna, atac = NULL, cellType = "cell_type", sample = "sample",
        condition = c("age", "old", "young"), covariates = "genotype", genes = sc$genes,
        peaks = sc$peaks, links = NULL, fragments = sc$fragments, meta = sc$meta,
        atacMeta = atacMeta))
    expect_true(all(grepl("^Hepatocyte \\| age_", res$top$group)))
    none <- suppressMessages(scLocusPick(
        sc$rna, atac = NULL, cellType = "cell_type", sample = "sample", genes = sc$genes,
        peaks = sc$peaks, links = NULL, fragments = sc$fragments, meta = sc$meta,
        atacMeta = atacMeta, weights = c(accessibility_mean = 1)))
    expect_setequal(none$cellTypes, c("Hepatocyte", "Kupffer"))
    k <- none$top[none$top$cell_type == "Kupffer", ]
    expect_equal(k$gene[1], "geneF")                # the Kupffer-specific open peak
})

test_that("scLocusPick reads a Seurat + Signac object", {
    skip_if_not_installed("Signac"); skip_if_not_installed("SeuratObject")
    skip_if_not_installed("DESeq2")
    so <- suppressWarnings(suppressMessages(locusPickSingleCellExample(seurat = TRUE)))$object
    expect_s4_class(so, "Seurat")
    res <- suppressMessages(scLocusPick(so, cellType = "cell_type", sample = "sample",
                                        perturbation = c("genotype", "KO", "WT"),
                                        genes = sc$genes))
    expect_true("link_score" %in% names(res$candidates))    # links taken from the object
    expect_equal(res$top$gene[res$top$group == "Hepatocyte | activated"][1], "geneA")
    tr <- suppressMessages(scSignal(so, c("cell_type", "age")))
    expect_length(tr, 4)
})
