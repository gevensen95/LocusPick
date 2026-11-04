## =============================================================================
## Fig. 4. Global shifts: a locus that "loses signal" is only an example of a
## group-specific change if the rest of the genome does not lose it too.
##
## A. ChIP-Rx (GSE60104, Orlando et al. 2014): H3K79me2 in Jurkat cells after
##    DOT1L inhibition (EPZ5676) vs DMSO, with a Drosophila spike-in. The
##    inhibitor lowers H3K79me2 genome-wide. CPM tracks (nf-core/chipseq
##    default, scaled to total mapped reads) partly hide that; spike-in tracks
##    (reads per million fly reads) show it. We compare:
##      - per-peak log2(EPZ5676 / DMSO) under both normalizations
##      - signalShift() by baseline-signal quartile, against all peaks
##      - the 10 "largest loss" loci picked under each normalization
##    What to expect: a per-condition scale factor shifts every log ratio by
##    the same constant. So (i) "is there a loss?" (p_vs_zero, median vs 0)
##    depends on the normalization, (ii) the comparison against all peaks
##    (shift_vs_background) and picks ranked on the log ratio barely change,
##    (iii) picks ranked on a signal DIFFERENCE (LocusPick's signal_change
##    column) can change. Panel B and the Jaccard table show which holds.
## B. Hnf4a CUT&RUN young vs old (chromDicts): signalShift() per gene group
##    against all consensus peaks, i.e. is the age loss at G1/G2 examples
##    larger than the genome-wide loss?
##
## Inputs: nf-core/chipseq output for GSE60104 (00_process_gse60104.sh) and the
## design table templates/gse60104_design.csv; Hnf4a paths in config.R.
## Outputs (results/fig4/): spikein_factors.tsv, shift_by_quartile.tsv,
## largest_loss_picks.tsv, hnf4a_shift_by_group.tsv, fig4.pdf/png
## =============================================================================

source("config.R"); source("R/common.R")
out <- figDir("fig4")
set.seed(SEED)

## ---- A. ChIP-Rx ------------------------------------------------------------------------
needFiles(CHIPRX$design, "GSE60104 design table")
design <- utils::read.csv(CHIPRX$design, stringsAsFactors = FALSE)
design <- design[!as.logical(design$is_input) & design$condition %in% c(CHIPRX$ref, CHIPRX$test), ]
if (!nrow(design) || any(design$sample == ""))
    stop("Fill in templates/gse60104_design.csv (sample, condition, ...) first.", call. = FALSE)
design$sample <- paste0(design$sample, "_REP", design$replicate)    # nf-core/chipseq naming

fo <- findPipelineOutputs(CHIPRX$nfcore)
bw <- fo[fo$type == "bigwig" & fo$sample %in% design$sample, ]
if (!all(design$sample %in% bw$sample))
    stop("No bigWig for: ", paste(setdiff(design$sample, bw$sample), collapse = ", "), call. = FALSE)

## Spike-in factors from samtools idxstats of the filtered BAMs:
## factor = all mapped reads / fly reads, so CPM x factor = reads per million fly reads.
ix <- list.files(CHIPRX$nfcore, "\\.idxstats$", recursive = TRUE, full.names = TRUE)
ix <- ix[grepl("mLb\\.clN", basename(ix))]
sf <- do.call(rbind, lapply(ix, function(f) {
    d <- utils::read.delim(f, header = FALSE, col.names = c("chr", "len", "mapped", "unmapped"))
    d <- d[d$chr != "*", ]
    fly <- sum(d$mapped[startsWith(d$chr, CHIPRX$spikePrefix)])
    data.frame(sample = sub("\\.mLb\\.clN.*$", "", basename(f)), total = sum(d$mapped),
               fly = fly, human = sum(d$mapped) - fly)
}))
sf <- merge(design[, c("sample", "condition", "replicate")], sf, by = "sample")
if (nrow(sf) != nrow(design) || any(sf$fly == 0))
    stop("Spike-in read counts missing for some samples; check the idxstats files and spikePrefix.",
         call. = FALSE)
sf$fly_fraction <- sf$fly / sf$total
sf$factor <- sf$total / sf$fly
writeTSV(sf, file.path(out, "spikein_factors.tsv"))
print(sf)

trackOf <- function(cond, spike) {
    s <- sf[sf$condition == cond, ]
    src <- lapply(seq_len(nrow(s)), function(i) {
        f <- bw$file[match(s$sample[i], bw$sample)]
        if (spike) scaleSource(f, s$factor[i]) else bigWigSource(f)
    })
    if (length(src) == 1) src[[1]] else meanSource(src)
}
SIG <- list(CPM     = list(ref = trackOf(CHIPRX$ref, FALSE), test = trackOf(CHIPRX$test, FALSE)),
            spikein = list(ref = trackOf(CHIPRX$ref, TRUE),  test = trackOf(CHIPRX$test, TRUE)))

pkFiles <- fo$file[fo$type == "peaks" & fo$sample %in% design$sample]
pk <- readPeaks(pkFiles, minReps = 2, blacklist = readBlacklist(REF$hg38$blacklist), verbose = FALSE)
pk <- pk[!startsWith(as.character(seqnames(pk)), CHIPRX$spikePrefix)]
say(length(pk), " consensus H3K79me2 peaks (>= 2 samples)")
gm <- readGeneModels(needFiles(REF$hg38$gtf))

perPeak <- list(); shifts <- list(); tops <- list()
for (nm in names(SIG)) {
    s <- SIG[[nm]]
    cand <- buildCandidates(pk, gm, groups = NULL, signal = s, mappability = REF$hg38$umap)
    cand <- cand[cand$eligible, ]
    cand$log2_test_ref <- log2((cand$test_mean + 0.01) / (cand$ref_mean + 0.01))
    q <- stats::quantile(cand$ref_mean, 0:4 / 4)
    cand$group <- paste0("baseline Q", findInterval(cand$ref_mean, q, rightmost.closed = TRUE, all.inside = TRUE))
    bg <- signalAtPeaks(pk, s)
    sh <- signalShift(cand, "test_mean", "ref_mean", background = bg)
    sh$normalization <- nm
    shifts[[nm]] <- sh
    perPeak[[nm]] <- data.frame(peak_id = cand$peak_id, group = cand$group,
                                log2_test_ref = cand$log2_test_ref, normalization = nm)
    ## "largest loss" picks, one group: ranked on the log ratio or the difference
    cand$loss <- -cand$log2_test_ref
    cand$loss_diff <- cand$ref_mean - cand$test_mean
    cand$group <- "all"
    for (rk in c("loss", "loss_diff")) {
        tp <- topLoci(rankLoci(cand, stats::setNames(1, rk)), n = N_TOP)
        tp <- tp[, !startsWith(names(tp), "pct_")]
        tp$normalization <- nm; tp$ranked_on <- rk
        tops[[paste(nm, rk)]] <- tp
    }
}
shift <- do.call(rbind, shifts)
writeTSV(shift, file.path(out, "shift_by_quartile.tsv"))
tops <- do.call(rbind, tops)
writeTSV(tops[, c("normalization", "ranked_on", "gene", "peak_id", "loss", "loss_diff",
                 "ref_mean", "test_mean", "rank")], file.path(out, "largest_loss_picks.tsv"))
jac <- do.call(rbind, lapply(c("loss", "loss_diff"), function(rk) data.frame(
    ranked_on = rk, jaccard_CPM_vs_spikein = jaccard(
        tops$peak_id[tops$normalization == "CPM" & tops$ranked_on == rk],
        tops$peak_id[tops$normalization == "spikein" & tops$ranked_on == rk]))))
writeTSV(jac, file.path(out, "largest_loss_jaccard.tsv"))
print(jac)

## How specific do the CPM picks look under each normalization?
cpmTop <- tops$peak_id[tops$normalization == "CPM" & tops$ranked_on == "loss_diff"]
spec <- do.call(rbind, lapply(names(SIG), function(nm) {
    d <- perPeak[[nm]]; d$group <- ifelse(d$peak_id %in% cpmTop, "CPM top picks", "other")
    bg <- signalAtPeaks(pk, SIG[[nm]])
    sp <- signalAtPeaks(pk[pk$peak_id %in% cpmTop], SIG[[nm]]); sp$group <- "CPM top picks"
    r <- signalShift(sp, "test_mean", "ref_mean", background = bg); r$normalization <- nm; r
}))
writeTSV(spec, file.path(out, "cpm_picks_vs_background.tsv"))

pp <- do.call(rbind, perPeak)
pp$normalization <- factor(pp$normalization, c("CPM", "spikein"), c("CPM", "Spike-in"))
med <- stats::aggregate(log2_test_ref ~ normalization, pp, stats::median)
pA <- ggplot(pp, aes(log2_test_ref, colour = normalization)) +
    geom_density() + geom_vline(xintercept = 0, linetype = 2) +
    geom_vline(data = med, aes(xintercept = log2_test_ref, colour = normalization)) +
    scale_colour_manual(values = c(CPM = PAL[["CPM"]], `Spike-in` = PAL[["spikein"]])) +
    labs(x = sprintf("log2(%s / %s) per peak", CHIPRX$test, CHIPRX$ref), y = "Density",
         colour = NULL, title = "A  H3K79me2 ChIP-Rx") + theme_paper()
sh <- shift[shift$group != "background (all peaks)", ]
sh$normalization <- factor(sh$normalization, c("CPM", "spikein"), c("CPM", "Spike-in"))
pB <- ggplot(sh, aes(group, median_log2_a_over_b, colour = normalization)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = position_dodge(0.4)) +
    geom_hline(data = transform(shift[shift$group == "background (all peaks)", ],
                                normalization = factor(normalization, c("CPM", "spikein"), c("CPM", "Spike-in"))),
               aes(yintercept = median_log2_a_over_b, colour = normalization), linetype = 3) +
    scale_colour_manual(values = c(CPM = PAL[["CPM"]], `Spike-in` = PAL[["spikein"]])) +
    labs(x = NULL, y = "Median log2 ratio\n(dotted: all peaks)", colour = NULL,
         title = "B  By baseline signal") + theme_paper()

## ---- B. Hnf4a CUT&RUN --------------------------------------------------------------------
orgDb <- if (!TEST && requireNamespace("org.Mm.eg.db", quietly = TRUE)) org.Mm.eg.db::org.Mm.eg.db
de <- hnf4aDE()
groups <- defineGroups(de$ko, de$age, perturbationType = "loss", conditionName = "age")
hk <- hnf4aPeaks()
hg <- readGeneModels(needFiles(HNF4A$gtf),
                     harmonize = if (!is.null(orgDb)) function(s) harmonizeSymbols(s, orgDb))
hs <- list(Yng = hnf4aSource(HNF4A$chromDicts[["CNR_Yng"]]), Old = hnf4aSource(HNF4A$chromDicts[["CNR_Old"]]))
hc <- buildCandidates(hk, hg, groups, signal = hs, mappability = HNF4A$umap)
hc <- hc[hc$eligible, ]
hbg <- signalAtPeaks(hk, hs)
hsh <- signalShift(hc, "Yng_mean", "Old_mean", background = hbg)
writeTSV(hsh, file.path(out, "hnf4a_shift_by_group.tsv"))
print(hsh)

hp <- rbind(data.frame(group = as.character(hc$group), r = log2((hc$Yng_mean + 0.01) / (hc$Old_mean + 0.01))),
            data.frame(group = "all peaks", r = log2((hbg$Yng_mean + 0.01) / (hbg$Old_mean + 0.01))))
hp <- hp[is.finite(hp$r), ]
hp$group <- factor(hp$group, c("all peaks", levels(groups$group)))
pC <- ggplot(hp, aes(group, r)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = stats::median(hp$r[hp$group == "all peaks"]), linetype = 3, colour = "grey40") +
    geom_boxplot(outlier.size = 0.3, width = 0.6, fill = "grey90") +
    labs(x = NULL, y = "log2(young / old) CUT&RUN\n(dotted: median of all peaks)",
         title = "C  Hnf4a CUT&RUN, by gene group") +
    theme_paper() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

## ---- Example locus under both normalizations -----------------------------------------------
ex <- tops[tops$normalization == "CPM" & tops$ranked_on == "loss_diff", ][1, ]
tr <- list(DMSO_CPM = SIG$CPM$ref, EPZ_CPM = SIG$CPM$test,
           DMSO_spikein = SIG$spikein$ref, EPZ_spikein = SIG$spikein$test)
names(tr) <- sub("^DMSO", CHIPRX$ref, sub("^EPZ", CHIPRX$test, names(tr)))
pD <- plotLocus(ex, tr, gm, trackGroups = list(CPM = names(tr)[1:2], `Spike-in` = names(tr)[3:4]),
                annotate = FALSE, maxWindow = 6e4)

p <- (pA | pB) / pC + plot_annotation(caption = if (TEST) "TEST MODE: simulated data, numbers are meaningless")
saveFig(p, file.path(out, "fig4"), width = 9, height = 7)
pdf(file.path(out, "fig4_example_locus.pdf"), width = 7, height = 5)
print(pD + plot_annotation(title = sprintf("D  Top CPM 'loss' locus (%s): CPM vs spike-in", ex$gene)))
dev.off()
say("Fig. 4 done: ", out)
