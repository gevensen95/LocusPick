ex <- locusPickExample()

test_that("gene models are read with spans, TSS and exons", {
    gm <- ex$genes
    expect_s3_class(gm, "locusGeneModels")
    expect_setequal(gm$genes$gene_name, paste0("gene", LETTERS[1:6]))
    b <- gm$genes[gm$genes$gene_name == "geneB"]
    expect_equal(as.character(GenomicRanges::strand(b)), "-")
    expect_equal(b$tss, GenomicRanges::end(b))  # minus strand TSS = end
    expect_true(all(gm$exons$gene_name %in% gm$genes$gene_name))
})

test_that("asPeaks reads BED as 1-based and adds peak ids", {
    expect_equal(GenomicRanges::start(ex$peaks)[1], 4800)
    expect_equal(ex$peaks$peak_id[1], "chr1:4800-5300")
    df <- asPeaks(data.frame(chr = "chr1", start = 99, end = 200), zeroBased = TRUE)
    expect_equal(GenomicRanges::start(df), 100)
})

test_that("assignPeaks classifies features and respects maxDist", {
    a <- assignPeaks(ex$peaks, ex$genes, maxDist = 5000)
    expect_equal(a$feature[a$peak_id == "chr1:4800-5300" & a$gene == "geneA"], "promoter")
    expect_equal(a$feature[a$peak_id == "chr1:27000-27400" & a$gene == "geneC"], "gene body")
    expect_true(all(a$distance <= 5000))
    near <- assignPeaks(ex$peaks, ex$genes, mode = "nearest")
    expect_false(any(duplicated(near$peak_id)) &&
                 any(near$distance > stats::ave(near$distance, near$peak_id, FUN = min)))
    expect_equal(near$gene[near$peak_id == "chr1:19800-20300"], "geneB")
    none <- assignPeaks(ex$peaks, ex$genes, maxDist = 0)
    expect_false("chr1:58500-58900" %in% none$peak_id)
})

test_that("checkAssignment flags distant and unknown genes", {
    pk <- ex$peaks
    pk$gene <- c("geneA", "geneB", "geneA", "geneD", "geneE", "geneE", "geneF", "notAGene")
    chk <- checkAssignment(pk, ex$genes, geneCol = "gene")
    expect_true(chk$near_gene[1])
    expect_false(chk$near_gene[3])            # peak at 27 kb is 15 kb past geneA
    expect_false(chk$gene_in_models[8])
    expect_true(is.na(chk$distance[8]))
})

test_that("defineGroups assigns regulation and direction", {
    g <- ex$groups
    expect_equal(as.character(g$group[g$gene == "geneA"]), "activated_condition_down")
    expect_equal(as.character(g$group[g$gene == "geneB"]), "repressed_condition_up")
    expect_false("geneF" %in% g$gene)  # not significant
    gain <- defineGroups(
        deTable(data.frame(gene = "x", l = 2, p = 0.01), lfc = "l", padj = "p"),
        deTable(data.frame(gene = "x", l = 1, p = 0.01), lfc = "l", padj = "p"),
        perturbationType = "gain")
    expect_equal(gain$regulation, "activated")
})

test_that("deTable keeps the most significant duplicate", {
    d <- data.frame(gene = c("a", "a"), lfc = c(1, 2), padj = c(0.5, 0.01))
    expect_message(out <- deTable(d, lfc = "lfc"), "duplicated")
    expect_equal(out$lfc, 2)
})

test_that("signal summaries handle coverage gaps and zeros", {
    s <- signalAtPeaks(ex$peaks, ex$signal)
    hole <- s$peak_id == "chr1:44000-44400"
    expect_true(is.na(s$A_mean[hole]))
    expect_equal(s$A_cov[hole], 0)
    expect_true(all(s$A_cov[!hole] > 0.99))
    expect_gt(s$A_mean[1], s$B_mean[1])
    m <- mappabilityAtPeaks(ex$peaks, ex$mappability)
    expect_equal(m[ex$peaks$peak_id == "chr1:34600-35100"], 0.3, tolerance = 1e-6)
    expect_equal(m[1], 1)
})

test_that("rangesSource and functionSource agree with bigWig", {
    w <- GenomicRanges::GRanges("chr1:4800-5300")
    bw <- fetchSignal(ex$signal$A, w)
    rs <- rangesSource(bw)
    fs <- functionSource(function(chr, start, end) {
        x <- bw[GenomicRanges::end(bw) >= start & GenomicRanges::start(bw) <= end]
        data.frame(start = GenomicRanges::start(x), end = GenomicRanges::end(x), score = x$score)
    })
    a <- signalAtPeaks(w, list(bw = ex$signal$A, gr = rs, fn = fs))
    expect_equal(a$gr_mean, a$bw_mean)
    expect_equal(a$fn_mean, a$bw_mean)
})

test_that("buildCandidates flags failing loci without dropping them", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal,
                            mappability = ex$mappability)
    expect_true(all(c("eligible", "mappable", "has_signal", "near_gene") %in% names(cand)))
    expect_false(cand$mappable[cand$gene == "geneD" & cand$peak_id == "chr1:34600-35100"])
    expect_false(any(cand$has_signal[cand$peak_id == "chr1:44000-44400"]))
    expect_true(cand$eligible[cand$gene == "geneA"])
    # nearest mode: geneB's promoter peak is not offered as a geneC example
    expect_false(any(cand$gene == "geneC" & cand$peak_id == "chr1:19800-20300"))
    all <- buildCandidates(ex$peaks, ex$genes, ex$groups, assignMode = "all")
    expect_true(any(all$gene == "geneC" & all$peak_id == "chr1:19800-20300"))
    q <- qcSummary(cand)
    expect_equal(sum(q$peaks), nrow(cand))
})

test_that("buildCandidates checks a user-supplied assignment", {
    pk <- ex$peaks
    pk$gene <- c("geneA", "geneB", "geneA", "geneD", "geneE", "geneE", "geneF", "geneF")
    cand <- buildCandidates(pk, ex$genes, ex$groups, geneCol = "gene")
    bad <- cand$peak_id == "chr1:27000-27400"
    expect_false(cand$near_gene[bad])
    expect_false(cand$eligible[bad])
})

test_that("rankLoci never scores on partial data and respects direction", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal,
                            mappability = ex$mappability)
    r <- rankLoci(cand, c(A_mean = 1))
    expect_true(all(is.na(r$score[!r$eligible])))
    g <- r[r$group == "activated_condition_down" & r$eligible, ]
    expect_equal(g$gene[which.max(g$score)], "geneA")
    rl <- rankLoci(cand, c(A_mean = -1))
    gl <- rl[rl$group == "activated_condition_down" & rl$eligible, ]
    expect_false(gl$gene[which.max(gl$score)] == "geneA" &&
                 gl$peak_id[which.max(gl$score)] == "chr1:4800-5300")
    cand$A_mean[cand$gene == "geneA"] <- NA
    cand$eligible[cand$gene == "geneA"] <- TRUE
    r2 <- rankLoci(cand, c(A_mean = 1))
    expect_true(is.na(r2$score[r2$gene == "geneA"]))
    expect_error(rankLoci(cand, c(nope = 1)), "missing column")
    perGroup <- stats::setNames(rep(list(c(A_mean = 1)), 4), levels(ex$groups$group))
    expect_silent(rankLoci(cand, perGroup))
})

test_that("topLoci keeps one peak per gene", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    t1 <- topLoci(rankLoci(cand, c(A_mean = 1)), n = 5)
    expect_false(any(duplicated(t1[, c("group", "gene")])))
})

test_that("signalShift summarizes per group and background", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    bg <- signalAtPeaks(ex$peaks, ex$signal)
    s <- signalShift(cand, "A_mean", "B_mean", background = bg)
    expect_true("background (all peaks)" %in% s$group)
    expect_true(all(s$frac_a_gt_b >= 0 & s$frac_a_gt_b <= 1))
})

test_that("plotLocus and reportLoci produce output", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    p <- plotLocus(cand[1, ], ex$signal, ex$genes,
                   trackGroups = list(Signal = c("A", "B")), expression = ex$expression)
    expect_s3_class(p, "patchwork")
    f <- tempfile(fileext = ".png")
    ggplot2::ggsave(f, p, width = 8, height = 5, dpi = 50)
    expect_gt(file.size(f), 0)
    # a locus with no signal draws a placeholder rather than failing
    hole <- cand[cand$peak_id == "chr1:44000-44400", ][1, ]
    expect_s3_class(plotLocus(hole, ex$signal, ex$genes), "patchwork")
    pdf <- tempfile(fileext = ".pdf")
    ok <- reportLoci(cand[1:2, ], pdf, signal = ex$signal, genes = ex$genes,
                     candidates = cand)
    expect_true(all(ok))
    expect_true(file.exists(sub("\\.pdf$", "_candidates.tsv", pdf)))
    expect_error(plotLocus(cand[1:2, ], ex$signal, ex$genes), "one row")
})

test_that("harmonizeSymbols maps outdated symbols", {
    skip_if_not_installed("org.Mm.eg.db")
    x <- harmonizeSymbols(c("Hnf4a", "1600002H07Rik", "notAGene", NA),
                          org.Mm.eg.db::org.Mm.eg.db)
    expect_equal(x[c(1, 3)], c("Hnf4a", "notAGene"))
    expect_equal(x[2], "Tedc2")
    expect_true(is.na(x[4]))
    expect_true("1600002H07Rik" %in% attr(x, "changes")$original)
})

test_that("checkInputs passes good inputs and flags a build mismatch", {
    out <- NULL
    capture.output(out <- checkInputs(ex$peaks, ex$genes, ex$signal, ex$groups))
    expect_true(all(out$status != "fail"))
    bad <- ex$genes
    bad$genes <- GenomicRanges::shift(bad$genes, 1e6)
    capture.output(out2 <- checkInputs(ex$peaks, bad))
    expect_equal(out2$status[grepl("within", out2$check)], "fail")
    ens <- GenomeInfoDb::renameSeqlevels(ex$peaks, c(chr1 = "1"))
    capture.output(out3 <- checkInputs(ens, ex$genes))
    expect_equal(out3$status[1], "warn")
    expect_equal(out3$status[2], "ok")   # naming converted automatically
})
