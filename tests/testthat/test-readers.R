ex <- locusPickExample()
f <- function(x) system.file("extdata", x, package = "LocusPick")
tmp <- function(lines, name) {
    p <- file.path(tempdir(), name); writeLines(lines, p); p
}

# ---- peaks -------------------------------------------------------------------

test_that("narrowPeak is read 0-based with summit converted", {
    p <- tmp(c("chr1\t4799\t5300\tp1\t100\t.\t5.2\t12.1\t9.8\t251",
               "chr1\t19799\t20300\tp2\t80\t.\t4.1\t10.0\t8.0\t200"), "s1_peaks.narrowPeak")
    gr <- readPeaks(p, verbose = FALSE)
    expect_equal(GenomicRanges::start(gr)[1], 4800)
    expect_equal(gr$peak_format[1], "narrowPeak")
    expect_equal(gr$summit[1], 4800 + 251)
    expect_equal(gr$qvalue[2], 8)
})

test_that("MACS2 xls is read 1-based", {
    p <- tmp(c("# Command line: macs2 callpeak -t a.bam", "# name = s1",
               "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname",
               "chr1\t4800\t5300\t501\t5050\t40\t12.1\t6.5\t9.8\ts1_peak_1"), "s1_peaks.xls")
    gr <- readPeaks(p, verbose = FALSE)
    expect_equal(GenomicRanges::start(gr), 4800)
    expect_equal(gr$fold_enrichment, 6.5)
    expect_equal(gr$peak_format, "macs2_xls")
})

test_that("SEACR output is detected", {
    p <- tmp(c("chr1\t4799\t5300\t120.5\t3.2\tchr1:5000-5100",
               "chr1\t26999\t27400\t80.1\t2.1\tchr1:27100-27200"), "s1.stringent.bed")
    gr <- readPeaks(p, verbose = FALSE)
    expect_equal(gr$peak_format[1], "seacr")
    expect_equal(gr$max_region[1], "chr1:5000-5100")
    expect_equal(GenomicRanges::start(gr)[1], 4800)
})

test_that("HOMER peaks and annotatePeaks output are read", {
    p <- tmp(c("# HOMER Peaks", "# total peaks = 1",
               "#PeakID\tchr\tstart\tend\tstrand\tNormalized Tag Count\tfocus ratio\tfindPeaks Score",
               "peak1\tchr1\t4800\t5300\t+\t40.0\t0.9\t55.2"), "peaks.txt")
    gr <- readPeaks(p, verbose = FALSE)
    expect_equal(gr$peak_format, "homer")
    expect_equal(GenomicRanges::start(gr), 4800)
    expect_equal(gr$homer_score, 55.2)
    a <- tmp(c(paste0("PeakID (cmd=annotatePeaks.pl peaks.txt mm10)\tChr\tStart\tEnd\tStrand\t",
                      "Peak Score\tFocus Ratio/Region Size\tAnnotation\tDetailed Annotation\t",
                      "Distance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\t",
                      "Nearest Refseq\tNearest Ensembl\tGene Name"),
               paste0("peak1\tchr1\t4800\t5300\t+\t55\tNA\tpromoter-TSS (NM_1)\tx\t50\tNM_1\t1\t",
                      "Mm.1\tNM_1\tENSMUSG1\tgeneA")), "annotated.txt")
    ga <- readPeaks(a, verbose = FALSE)
    expect_equal(ga$peak_format, "homer_annotated")
    expect_equal(ga$homer_gene, "geneA")
    expect_equal(ga$homer_annotation, "promoter-TSS")
    expect_equal(ga$homer_distance_to_tss, 50)
})

test_that("DiffBind reports are recognized (GRanges and saved table)", {
    db <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(4800, 19800), c(5300, 20300)))
    db$Conc <- c(6.1, 5.0); db$Conc_Young <- c(6.8, 5.2); db$Conc_Old <- c(5.1, 4.8)
    db$Fold <- c(1.7, 0.4); db$p.value <- c(1e-6, 0.2); db$FDR <- c(1e-4, 0.5)
    gr <- readPeaks(db, verbose = FALSE)
    expect_equal(gr$peak_format[1], "diffbind")
    expect_equal(gr$diff_fold, c(1.7, 0.4))
    expect_true("diff_conc_young" %in% names(S4Vectors::mcols(gr)))
    p <- file.path(tempdir(), "db.tsv")
    utils::write.table(as.data.frame(db), p, sep = "\t", quote = FALSE, row.names = FALSE)
    g2 <- readPeaks(p, verbose = FALSE)
    expect_equal(g2$peak_format[1], "diffbind")
    expect_equal(GenomicRanges::start(g2)[1], 4800)
})

test_that("several files give consensus peaks; blacklist removes peaks", {
    a <- tmp(c("chr1\t4799\t5300", "chr1\t19799\t20300", "chr1\t30000\t30200"), "r1.bed")
    b <- tmp(c("chr1\t4850\t5250", "chr1\t40000\t40200"), "r2.bed")
    cons <- readPeaks(c(a, b), verbose = FALSE)
    expect_equal(length(cons), 1)
    expect_equal(cons$n_reps, 2L)
    any1 <- readPeaks(c(a, b), minReps = 1, verbose = FALSE)
    expect_equal(length(any1), 4)
    bl <- GenomicRanges::GRanges("chr1:30000-30100")
    expect_equal(length(readPeaks(a, blacklist = bl, verbose = FALSE)), 2)
})

test_that("findPipelineOutputs parses nf-core-like layouts", {
    d <- file.path(tempdir(), "nfcore")
    unlink(d, recursive = TRUE)
    dirs <- c("03_peak_calling/03_bed_to_bigwig", "03_peak_calling/04_called_peaks/seacr",
              "03_peak_calling/04_called_peaks/macs2", "03_peak_calling/05_consensus_peaks",
              "work/ab/cd")
    for (x in dirs) dir.create(file.path(d, x), recursive = TRUE)
    file.create(file.path(d, c(
        "03_peak_calling/03_bed_to_bigwig/young_R1.bigWig",
        "03_peak_calling/03_bed_to_bigwig/igg_R1.bigWig",
        "03_peak_calling/04_called_peaks/seacr/young_R1.seacr.peaks.stringent.bed",
        "03_peak_calling/04_called_peaks/macs2/young_R1_peaks.narrowPeak",
        "03_peak_calling/05_consensus_peaks/young.consensus.peaks.bed",
        "work/ab/cd/young_R1.bigWig")))
    o <- findPipelineOutputs(d)
    expect_equal(sum(o$type == "bigwig"), 2)           # work/ excluded
    expect_true(o$control[o$sample == "igg_R1"])
    expect_equal(o$caller[grepl("seacr", o$file)], "seacr")
    expect_equal(o$caller[grepl("narrowPeak", o$file)], "macs2")
    expect_equal(o$sample[o$type == "consensus"], "young")
    expect_setequal(o$sample[o$type == "bigwig"], c("young_R1", "igg_R1"))
    expect_equal(o$sample[grepl("seacr", o$file)], "young_R1")
})

# ---- signal ------------------------------------------------------------------

test_that("meanSource averages replicates; vectors of files work as one track", {
    w <- GenomicRanges::GRanges("chr1:4800-5300")
    s <- signalAtPeaks(w, list(A = ex$signal$A, B = ex$signal$B,
                               AB = c(f("toy_A.bw"), f("toy_B.bw"))))
    expect_equal(s$AB_mean, (s$A_mean + s$B_mean) / 2, tolerance = 1e-6)
    tr <- signalFromFiles(f(c("toy_A.bw", "toy_B.bw")), c(A = "toy_A", both = "toy_"))
    expect_named(tr, c("A", "both"))
})

test_that("bamSource gives CPM coverage, drops duplicates and low MAPQ", {
    skip_if_not_installed("Rsamtools"); skip_if_not_installed("GenomicAlignments")
    bam <- f("toy.bam")
    expect_message(src <- bamSource(bam), "paired")
    s <- signalAtPeaks(asPeaks(f("toy_peaks.bed")), list(bam = src))
    expect_true(all(s$bam_cov == 1))                     # BAM covers every base
    expect_gt(s$bam_mean[1], s$bam_mean[3])              # pile-up at peak1
    raw <- suppressMessages(bamSource(bam, normalize = "none", removeDuplicates = FALSE,
                                      minMapq = 0))
    s2 <- signalAtPeaks(asPeaks(f("toy_peaks.bed")), list(raw = raw))
    expect_gt(s2$raw_mean[1], 0)
    # every fragment covers 150 bp: total coverage = 150 x fragments kept
    full <- GenomicRanges::GRanges("chr1:1-60000")
    tot <- sum(GenomicRanges::width(fetchSignal(raw, full)) * fetchSignal(raw, full)$score)
    expect_equal(tot, 150 * 650)
})

# ---- DE readers --------------------------------------------------------------

make_counts <- function() {
    set.seed(3)
    m <- matrix(rnbinom(600, mu = 100, size = 10), nrow = 100,
                dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
    m[1:10, 4:6] <- m[1:10, 4:6] * 6L      # up in KO
    m[11:20, 4:6] <- m[11:20, 4:6] %/% 6L  # down in KO
    m
}

test_that("readDE reads DESeq2 objects and reports the contrast", {
    skip_if_not_installed("DESeq2")
    m <- make_counts()
    cd <- data.frame(genotype = factor(rep(c("WT", "KO"), each = 3), levels = c("WT", "KO")))
    dds <- suppressMessages(DESeq2::DESeq(
        DESeq2::DESeqDataSetFromMatrix(m, cd, ~ genotype), quiet = TRUE))
    res <- DESeq2::results(dds, contrast = c("genotype", "KO", "WT"))
    expect_message(d <- readDE(res), "genotype KO vs WT")
    expect_equal(attr(d, "contrast"), "genotype KO vs WT")
    expect_gt(d$lfc[d$gene == "g1"], 1)
    ex <- expressionFrom(dds, group = "genotype", panel = "KO")
    expect_equal(levels(factor(ex$group, unique(ex$group))), c("WT", "KO"))
    expect_equal(nrow(ex), 600)
})

test_that("readDE reads edgeR and limma objects", {
    skip_if_not_installed("edgeR"); skip_if_not_installed("limma")
    m <- make_counts()
    g <- factor(rep(c("WT", "KO"), each = 3), levels = c("WT", "KO"))
    y <- edgeR::DGEList(m, group = g)
    y <- edgeR::normLibSizes(y)
    des <- stats::model.matrix(~ g)
    fit <- edgeR::glmQLFit(edgeR::estimateDisp(y, des), des)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    expect_message(d1 <- readDE(qlf), "edgeR")
    expect_message(d2 <- readDE(edgeR::topTags(qlf, n = Inf)), "contrast")
    expect_equal(d1$lfc[d1$gene == "g1"], d2$lfc[d2$gene == "g1"])
    v <- limma::voom(y, des)
    lf <- limma::eBayes(limma::lmFit(v, des))
    expect_error(readDE(lf), "coef")
    d3 <- suppressWarnings(readDE(lf, coef = 2, verbose = FALSE))
    expect_gt(d3$lfc[d3$gene == "g1"], 1)
    expect_equal(attr(d3, "contrast"), "gKO")
    e <- expressionFrom(y, group = "group", panel = "KO")
    expect_equal(unique(e$group), c("WT", "KO"))
})

test_that("readDE detects table formats and warns about unknown direction", {
    tab <- data.frame(Gene = c("a", "b"), logFC = c(1, -1), AveExpr = 5,
                      t = 3, P.Value = 0.01, adj.P.Val = c(0.02, 0.03))
    expect_warning(d <- readDE(tab, verbose = FALSE), "direction")
    expect_equal(attr(d, "source"), "limma-style table")
    expect_silent(readDE(tab, direction = "KO/WT", verbose = FALSE))
    p <- file.path(tempdir(), "deseq.csv")
    utils::write.csv(data.frame(baseMean = 1, log2FoldChange = c(2, -1), lfcSE = 1,
                                stat = 1, pvalue = 0.01, padj = c(0.01, 0.2),
                                row.names = c("x", "y")), p)
    d2 <- readDE(p, direction = "KO/WT", verbose = FALSE)
    expect_equal(d2$gene, c("x", "y"))
    expect_error(readDE(data.frame(a = 1, b = 2)), "fold-change")
    expect_error(readDE(data.frame(gene = "ENSMUSG00000017950.16", logFC = 1, FDR = 0.01),
                        direction = "x"), "orgDb")
})

test_that("readDE converts Ensembl IDs with an OrgDb", {
    skip_if_not_installed("org.Hs.eg.db")
    tab <- data.frame(gene = c("ENSG00000101076.17", "ENSG00000141510"),
                      log2FoldChange = c(-2, 1), padj = c(1e-5, 1e-3))
    d <- suppressMessages(readDE(tab, orgDb = org.Hs.eg.db::org.Hs.eg.db, direction = "KO/WT"))
    expect_setequal(d$gene, c("HNF4A", "TP53"))
})

test_that("defineGroups accepts raw tool output directly", {
    de <- read.delim(f("toy_de.tsv"))
    ko <- data.frame(gene = de$gene, log2FoldChange = de$ko_lfc, padj = de$ko_padj)
    age <- data.frame(gene = de$gene, logFC = de$cond_lfc, FDR = de$cond_padj)
    g <- suppressWarnings(suppressMessages(defineGroups(ko, age)))
    expect_equal(g$group, ex$groups$group)
})

# ---- wrapper -----------------------------------------------------------------

test_that("locusPick runs end to end and its outputs plot", {
    de <- read.delim(f("toy_de.tsv"))
    res <- suppressMessages(locusPick(
        peaks = f("toy_peaks.bed"),
        perturbation = data.frame(gene = de$gene, log2FoldChange = de$ko_lfc, padj = de$ko_padj),
        condition = data.frame(gene = de$gene, logFC = de$cond_lfc, FDR = de$cond_padj),
        perturbationDirection = "KO/WT", conditionDirection = "old/young",
        genes = f("toy.gtf"),
        signal = list(CNR_young = f("toy_A.bw"), CNR_old = f("toy_B.bw")),
        mappability = f("toy_map.bw"), expression = ex$expression))
    expect_s3_class(res, "LocusPickResult")
    expect_equal(res$trackGroups, list(CNR = c("CNR_young", "CNR_old")))
    expect_true("background (all peaks)" %in% res$shift$group)
    expect_equal(res$top$gene[res$top$group == "activated_condition_down"][1], "geneA")
    expect_output(print(res), "Top loci")
    expect_s3_class(plot(res, "geneB"), "patchwork")
    pdf <- tempfile(fileext = ".pdf")
    ok <- suppressMessages(reportLoci(res, pdf))
    expect_true(all(ok))
    expect_true(file.exists(sub("\\.pdf$", "_candidates.tsv", pdf)))
})

test_that("scaleSource multiplies scores and rejects bad factors", {
    bw <- system.file("extdata", "toy_A.bw", package = "LocusPick")
    w <- GenomicRanges::GRanges("chr1:5000-5040")
    a <- fetchSignal(bw, w); b <- fetchSignal(scaleSource(bw, 2.5), w)
    expect_equal(b$score, a$score * 2.5)
    expect_error(scaleSource(bw, 0), "positive")
    expect_error(scaleSource(bw, c(1, 2)), "single")
})
