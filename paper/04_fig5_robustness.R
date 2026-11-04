## =============================================================================
## Fig. 5. Are LocusPick's picks robust to reasonable analysis choices?
##
## For every ENCODE factor (HepG2 vs K562) and for Hnf4a CUT&RUN (young vs
## old), the top N_TOP loci per group are recomputed after changing one choice
## at a time and compared with the reference picks (Jaccard of top sets):
##   replicates : leave one replicate out (ENCODE: 2 reps; Hnf4a: placeholder
##                per-replicate bigWigs, skipped until they exist)
##   weights    : weightStability(), every weight drawn log-uniform in [0.5, 2]
##   maxDist    : peak-to-gene distance cap 1, 5 (reference), 10, 50 kb
##   normalization (ChIP-Rx, from Fig. 4): CPM vs spike-in, if fig4 ran
## Typicality: percentile of each top pick within its group for every ranking
## criterion (pct_<column> from rankLoci), showing whether picks are typical
## or extreme members of their group.
##
## Outputs (results/fig5/): overlap.tsv, weight_stability.tsv, typicality.tsv,
## fig5.pdf/png
## =============================================================================

source("config.R"); source("R/common.R")
out <- figDir("fig5")
set.seed(SEED)
B_WEIGHTS <- if (TEST) 50 else 200
DISTS <- c(1e3, 5e3, 1e4, 5e4)

overlap <- list(); stab <- list(); typ <- list()
addOverlap <- function(dataset, analysis, rankings) {
    o <- topOverlap(rankings, n = N_TOP, reference = 1)
    o <- o[o$ranking != names(rankings)[1], ]
    if (nrow(o)) overlap[[length(overlap) + 1]] <<- cbind(dataset = dataset, analysis = analysis, o)
}

## Shared routine: candidates -> reference ranking -> variations.
## `repSignals` is a list of signal lists, one per left-out replicate.
## `post` adds derived ranking columns to every candidate table.
robustness <- function(dataset, peaks, genes, groups, signal, repSignals, mappability, weights,
                       post = identity) {
    say("== ", dataset)
    build <- function(sig, maxDist = 5000)
        post(buildCandidates(peaks, genes, groups, signal = sig, mappability = mappability,
                             maxDist = maxDist))
    cand <- build(signal)
    ref <- rankLoci(cand, weights)

    if (length(repSignals))
        addOverlap(dataset, "replicates", c(list(reference = ref),
                                            lapply(repSignals, function(s) rankLoci(build(s), weights))))
    ws <- weightStability(cand, weights, n = N_TOP, B = B_WEIGHTS, seed = SEED)
    stab[[dataset]] <<- cbind(dataset = dataset, ws$summary)
    rk <- lapply(DISTS, function(d) if (d == 5000) ref else rankLoci(build(signal, d), weights))
    names(rk) <- paste0("maxDist_", DISTS / 1000, "kb")
    rk <- c(rk["maxDist_5kb"], rk[names(rk) != "maxDist_5kb"])
    addOverlap(dataset, "maxDist", rk)

    tp <- topLoci(ref, n = N_TOP)
    pc <- grep("^pct_", names(tp), value = TRUE)
    typ[[dataset]] <<- do.call(rbind, lapply(pc, function(p)
        data.frame(dataset = dataset, group = tp$group, gene = tp$gene,
                   criterion = sub("^pct_", "", p), percentile = tp[[p]])))
    invisible(ref)
}

## ---- ENCODE ---------------------------------------------------------------------------------
m  <- readManifest()
gm <- readGeneModels(needFiles(REF$hg38$gtf))
bl <- readBlacklist(REF$hg38$blacklist)
de <- encodeRNADE(m, REF$hg38$gtf)
cond <- readDE(de, gene = "gene", lfc = "log2FoldChange", padj = "padj",
               direction = paste(ENCODE$cells, collapse = "/"), verbose = FALSE)
cells <- unname(ENCODE$cells)
groups <- defineGroups(condition = cond, conditionName = cells[1])
for (tf in setdiff(unique(m$target), "RNA")) {
    pk <- readPeaks(c(encodeFile(m, tf, cells[1], "peaks"), encodeFile(m, tf, cells[2], "peaks")),
                    minReps = 1, blacklist = bl, verbose = FALSE)
    reps <- lapply(cells, function(cl) encodeFile(m, tf, cl, "fc_rep", 1:2))
    sig  <- stats::setNames(lapply(reps, function(r) meanSource(as.list(r))), cells)
    loro <- lapply(1:2, function(drop) stats::setNames(lapply(reps, function(r) r[-drop]), cells))
    names(loro) <- paste0("without_rep", 1:2)
    w <- stats::setNames(c(1, 1), c(paste0(cells[1], "_mean"), "abs_lfc_condition"))
    robustness(tf, pk, gm, groups, sig, loro, REF$hg38$umap, w)
}

## ---- Hnf4a ------------------------------------------------------------------------------------
orgDb <- if (!TEST && requireNamespace("org.Mm.eg.db", quietly = TRUE)) org.Mm.eg.db::org.Mm.eg.db
hde <- hnf4aDE()
hgroups <- defineGroups(hde$ko, hde$age, perturbationType = "loss", conditionName = "age")
hg <- readGeneModels(needFiles(HNF4A$gtf),
                     harmonize = if (!is.null(orgDb)) function(s) harmonizeSymbols(s, orgDb))
hk <- hnf4aPeaks()
repsOK <- all(file.exists(unlist(HNF4A$cnr_reps)))
if (repsOK) {
    hsig <- lapply(HNF4A$cnr_reps, function(r) meanSource(as.list(r)))
    nRep <- min(lengths(HNF4A$cnr_reps))
    hloro <- lapply(seq_len(nRep), function(d) lapply(HNF4A$cnr_reps, function(r) r[-d]))
    names(hloro) <- paste0("without_rep", seq_len(nRep))
} else {
    message("Hnf4a per-replicate bigWigs not found (placeholders in config.R): ",
            "leave-one-replicate-out skipped; using the chromDict tracks.")
    hsig <- list(Yng = hnf4aSource(HNF4A$chromDicts[["CNR_Yng"]]),
                 Old = hnf4aSource(HNF4A$chromDicts[["CNR_Old"]]))
    hloro <- list()
}
## G1/G2 additionally reward a loss with age, as in the original selection.
base <- c(Yng_mean = 1, abs_lfc_perturbation = 1, abs_lfc_condition = 1)
hw <- list(activated_age_down = c(base, signal_change = 1), repressed_age_up = c(base, signal_change = 1),
           activated_age_up = base, repressed_age_down = base)
hw <- hw[names(hw) %in% unique(as.character(hgroups$group))]
robustness("Hnf4a CUT&RUN", hk, hg, hgroups, hsig, hloro, HNF4A$umap, hw,
           post = function(c) { c$signal_change <- c$Yng_mean - c$Old_mean; c })

## ---- Normalization (from Fig. 4) ----------------------------------------------------------------
f4 <- file.path(RESULTS, "fig4", "largest_loss_jaccard.tsv")
if (file.exists(f4)) {
    j <- utils::read.delim(f4)
    overlap[[length(overlap) + 1]] <- data.frame(dataset = "ChIP-Rx H3K79me2", analysis = "normalization",
                                                 group = "all", ranking = paste0("spikein_", j$ranked_on),
                                                 jaccard = j$jaccard_CPM_vs_spikein)
} else message("Run 03_fig4_global_shift.R first to include the normalization comparison.")

ov <- do.call(rbind, overlap)
ws <- do.call(rbind, stab)
ty <- do.call(rbind, typ)
ty <- ty[!is.na(ty$percentile), ]   # criteria not used in a group
writeTSV(ov, file.path(out, "overlap.tsv"))
writeTSV(ws, file.path(out, "weight_stability.tsv"))
writeTSV(ty, file.path(out, "typicality.tsv"))

## ---- Figure -----------------------------------------------------------------------------------------
pd <- rbind(ov[, c("dataset", "analysis", "ranking", "jaccard")],
            data.frame(dataset = ws$dataset, analysis = "weights", ranking = "random weights",
                       jaccard = ws$mean_jaccard))
pd$analysis <- factor(pd$analysis, c("replicates", "weights", "maxDist", "normalization"))
pA <- ggplot(pd, aes(analysis, jaccard)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey92") +
    geom_jitter(aes(colour = dataset), width = 0.12, height = 0, size = 1.2) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = sprintf("Jaccard, top %d vs reference", N_TOP), colour = NULL,
         title = "A  Stability of the top loci") + theme_paper()
pB <- ggplot(ty, aes(criterion, percentile)) +
    geom_hline(yintercept = c(0.5), linetype = 2, colour = "grey50") +
    geom_violin(fill = "grey92", scale = "width") +
    geom_jitter(width = 0.1, height = 0, size = 0.6, alpha = 0.6) +
    facet_wrap(~dataset, scales = "free_x") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Percentile within group", title = "B  How typical are the picks?") +
    theme_paper() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p <- pA / pB + plot_layout(heights = c(1, 1.2)) +
    plot_annotation(caption = if (TEST) "TEST MODE: simulated data, numbers are meaningless")
saveFig(p, file.path(out, "fig5"), width = 9, height = 7)
say("Fig. 5 done: ", out)
