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

test_that("harmonizeSymbols works with a human OrgDb", {
    skip_if_not_installed("org.Hs.eg.db")
    x <- harmonizeSymbols(c("TCF14", "HNF4A", "MLL", "notAGene"),
                          org.Hs.eg.db::org.Hs.eg.db)
    expect_equal(as.vector(x), c("HNF4A", "HNF4A", "KMT2A", "notAGene"))
    expect_equal(nrow(attr(x, "changes")), 2)
})

test_that("one contrast or none is supported", {
    de <- read.delim(system.file("extdata", "toy_de.tsv", package = "LocusPick"))
    ko  <- deTable(de, lfc = "ko_lfc", padj = "ko_padj")
    age <- deTable(de, lfc = "cond_lfc", padj = "cond_padj")
    gp <- defineGroups(ko)
    expect_equal(levels(gp$group), c("activated", "repressed"))
    expect_false("lfc_condition" %in% names(gp))
    gc <- defineGroups(condition = age, conditionName = "age")
    expect_equal(levels(gc$group), c("age_down", "age_up"))
    expect_error(defineGroups(), "at least one contrast")

    c1 <- buildCandidates(ex$peaks, ex$genes, gp, signal = ex$signal)
    expect_true("abs_lfc_perturbation" %in% names(c1))
    expect_false("abs_lfc_condition" %in% names(c1))
    c0 <- buildCandidates(ex$peaks, ex$genes, groups = NULL, signal = ex$signal)
    expect_equal(as.character(unique(c0$group)), "all")
    expect_true("geneF" %in% c0$gene)      # no DE filter

    f <- function(x) system.file("extdata", x, package = "LocusPick")
    r1 <- suppressMessages(locusPick(f("toy_peaks.bed"), perturbation = ko, genes = f("toy.gtf"),
                                     signal = list(young = f("toy_A.bw"), old = f("toy_B.bw"))))
    expect_equal(names(r1$weights), c("young_mean", "abs_lfc_perturbation"))
    expect_equal(r1$contrasts[["condition"]], "none")
    r0 <- suppressMessages(locusPick(f("toy_peaks.bed"), genes = f("toy.gtf"),
                                     signal = list(young = f("toy_A.bw"), old = f("toy_B.bw")),
                                     weights = c(signal_change = 1)))
    expect_equal(r0$top$gene[1], "geneA")   # largest young - old difference
    expect_error(suppressMessages(locusPick(f("toy_peaks.bed"), genes = f("toy.gtf"))),
                 "Nothing to rank")
})

test_that("peak-gene links replace distance assignment", {
    links <- data.frame(peak = c("chr1-4800-5300", "chr1-27000-27400", "chr1-27000-27400"),
                        gene = c("geneA", "geneA", "geneE"),
                        score = c(0.6, 0.4, 0.2), pvalue = c(1e-5, 1e-3, 0.04))
    cand <- buildCandidates(ex$peaks, ex$genes, groups = NULL, links = links,
                            maxLinkDist = 25000)
    expect_equal(nrow(cand), 3)
    expect_true(all(c("link_score", "link_pvalue") %in% names(cand)))
    far <- cand$gene == "geneA" & cand$peak_id == "chr1:27000-27400"
    expect_equal(cand$distance[far], 14999)
    expect_true(cand$near_gene[far])
    strict <- buildCandidates(ex$peaks, ex$genes, groups = NULL, links = links,
                              maxLinkDist = 10000)
    expect_false(strict$near_gene[strict$gene == "geneE"])   # 20.4 kb away
})

test_that("rankLoci reports percentiles and signalShift tests specificity", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    r <- rankLoci(cand, c(A_mean = 1, abs_lfc_perturbation = -1))
    expect_true(all(c("pct_A_mean", "pct_abs_lfc_perturbation") %in% names(r)))
    expect_true(all(r$pct_A_mean[r$eligible] > 0 & r$pct_A_mean[r$eligible] <= 1))
    g <- r[r$group == "activated_condition_down" & r$eligible, ]
    expect_equal(g$pct_A_mean[which.max(g$A_mean)], 1)
    # a global shift: every background peak 2-fold lower in B
    set.seed(1)
    bg <- data.frame(A = rexp(200) + 1); bg$B <- bg$A / 2
    d <- data.frame(group = "g", peak_id = paste0("p", 1:10), A = rexp(10) + 1)
    d$B <- d$A / 2
    s <- signalShift(d, "A", "B", background = bg)
    expect_lt(s$p_vs_zero[1], 0.01)              # changed ...
    expect_gt(s$p_vs_background[1], 0.05)        # ... but not more than everywhere
    d$B <- d$A / 8
    s2 <- signalShift(d, "A", "B", background = bg)
    expect_lt(s2$p_vs_background[1], 0.01)       # specific
    expect_equal(s2$shift_vs_background[1], 2, tolerance = 0.05)
})

test_that("weightStability and topOverlap measure stability", {
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    w <- c(A_mean = 1, abs_lfc_perturbation = 1, abs_lfc_condition = 1)
    ws <- weightStability(cand, w, n = 1, B = 30)
    expect_true(all(ws$summary$mean_jaccard >= 0 & ws$summary$mean_jaccard <= 1))
    expect_true(all(ws$loci$top_freq >= 0 & ws$loci$top_freq <= 1))
    single <- ws$summary$mean_jaccard[ws$summary$group == "repressed_condition_up"]
    expect_equal(single, 1)                       # one candidate: always on top
    r1 <- rankLoci(cand, w); r2 <- rankLoci(cand, c(B_mean = 1))
    ov <- topOverlap(list(full = r1, alt = r2), n = 2)
    expect_true(all(ov$jaccard >= 0 & ov$jaccard <= 1))
    expect_equal(unique(ov$ranking), "alt")
})

test_that("weightStability accepts per-group weight lists", {
    ex <- locusPickExample()
    cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
    w <- stats::setNames(rep(list(c(A_mean = 1, abs_lfc_perturbation = 1)), length(unique(cand$group))),
                         unique(as.character(cand$group)))
    ws <- weightStability(cand, w, n = 1, B = 10)
    expect_true(all(ws$summary$mean_jaccard >= 0 & ws$summary$mean_jaccard <= 1))
})
