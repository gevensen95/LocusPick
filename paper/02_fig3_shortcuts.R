## =============================================================================
## Fig. 3. How often do common shortcut rules pick a misleading example locus?
##
## Data: ENCODE TF ChIP-seq and ATAC-seq in HepG2 and K562 (00_select_encode.R),
## with HepG2 vs K562 RNA-seq as the expression contrast (one contrast, so
## groups are "HepG2_up" and "HepG2_down").
##
## For every factor, peaks from both cell lines are pooled and assigned to
## genes in two ways:
##   shortcuts : the nearest TSS at any distance (what annotation tools report),
##               then the top N_RULE (peak, gene) pairs per group by
##                 top_signal - strongest peak (max over the two cell lines)
##                 top_lfc    - largest |log2FC| of the gene (its strongest peak)
##                 top_change - largest signal change in the group's direction
##                 random     - random pairs (200 draws)
##   LocusPick : locusPick() with blacklist removal, three configurations
##     LocusPick                - default weights: HepG2 signal + |log2FC|
##     LocusPick +change        - plus the signal change in the group's
##                                direction (HepG2 - K562 for HepG2_up, the
##                                reverse for HepG2_down)
##     LocusPick +change +reps  - as above, keeping only peaks whose change
##                                agrees in every HepG2 x K562 replicate pair
##                                (buildCandidates(replicates =))
##   The +change configurations were added after the first run showed that
##   the default weights, which ignore the change, pick loci whose change does
##   not reproduce. Report all three.
## Every pick is then judged with the same checks (pairFlags()):
##   peak > 5 kb from the gene body, mean Umap k50 < 0.8, signal track
##   coverage < 50 % of the peak in either line, overlaps the ENCODE blacklist.
## LocusPick filters on these by design, so its failure rate is ~0 by
## construction; the figure's point is the size of the rates for the shortcuts.
## One check is NOT used by LocusPick and so is a fair comparison for all rules:
##   not_reproduced - the sign of the HepG2 - K562 signal difference differs
##                    between replicate 1 and replicate 2 (per-replicate bigWigs).
##                    Fair for every rule except "+change +reps", which filters
##                    on a stricter version of it (circular there, shown for
##                    completeness).
## The signal-gap check cannot fail on ENCODE fold-change bigWigs (they cover
## the genome); it stays in the tables but is left out of the figure.
##
## Outputs (results/fig3/): picks.tsv, failure_rates.tsv, failure_rates_by_factor.tsv,
## fig3_failure_rates.pdf/png, fig3_examples_<factor>.pdf
## =============================================================================

source("config.R"); source("R/common.R")
out <- figDir("fig3")
set.seed(SEED)
RULES <- c("LocusPick", "LocusPick +change", "LocusPick +change +reps",
           "top_signal", "top_lfc", "top_change", "random")
FLAGS <- c(far_from_gene = "Peak > 5 kb from gene", low_mappability = "Low mappability",
           no_signal = "Gap in signal track", blacklisted = "ENCODE blacklist",
           not_reproduced = "Change not reproduced\nacross replicates**")

m  <- readManifest()
needFiles(c(REF$hg38$gtf, REF$hg38$umap, REF$hg38$blacklist), "hg38 reference")
gm <- readGeneModels(REF$hg38$gtf)
bl <- readBlacklist(REF$hg38$blacklist)
de <- encodeRNADE(m, REF$hg38$gtf)
cond <- readDE(de, gene = "gene", lfc = "log2FoldChange", padj = "padj",
               direction = paste(ENCODE$cells, collapse = "/"), verbose = FALSE)
cells <- unname(ENCODE$cells)
factors <- setdiff(unique(m$target), "RNA")
say("Factors: ", paste(factors, collapse = ", "))

## Pick N per group by a score (higher = picked first), one row per gene.
pickBy <- function(tab, score, n = N_RULE) {
    tab$.s <- score
    tab <- tab[!is.na(tab$.s) & !is.na(tab$group), ]
    tab <- tab[order(tab$group, -tab$.s), ]
    tab <- tab[!duplicated(tab[, c("group", "gene")]), ]
    do.call(rbind, lapply(split(tab, tab$group, drop = TRUE), utils::head, n))
}

allPicks <- list()
for (tf in factors) {
    say("== ", tf)
    pk <- readPeaks(c(encodeFile(m, tf, cells[1], "peaks"), encodeFile(m, tf, cells[2], "peaks")),
                    minReps = 1, verbose = FALSE)
    signal <- stats::setNames(list(encodeFile(m, tf, cells[1], "fc_pooled"),
                                   encodeFile(m, tf, cells[2], "fc_pooled")), cells)

    ## ---- LocusPick ------------------------------------------------------------
    res <- suppressMessages(locusPick(pk[!overlapsAny(pk, bl)], condition = cond, genes = gm,
                                      signal = signal, mappability = REF$hg38$umap,
                                      conditionName = cells[1], n = N_RULE))
    lp <- res$top[, c("group", "gene", "peak_id")]
    lp$rule <- "LocusPick"; lp$draw <- 0L

    ## +change: reward the signal change in the group's direction
    grpNames <- unique(as.character(res$candidates$group))
    base <- stats::setNames(c(1, 1), c(paste0(cells[1], "_mean"), "abs_lfc_condition"))
    wChange <- stats::setNames(lapply(grpNames, function(g)
        c(base, signal_change = if (grepl("_up$", g)) 1 else -1)), grpNames)
    reps <- stats::setNames(lapply(cells, function(cl) encodeFile(m, tf, cl, "fc_rep", 1:2)), cells)
    lpVariant <- function(label, replicates = NULL) {
        r <- suppressMessages(locusPick(pk[!overlapsAny(pk, bl)], condition = cond, genes = gm,
                                        signal = signal, mappability = REF$hg38$umap,
                                        conditionName = cells[1], weights = wChange,
                                        replicates = replicates, n = N_RULE))
        o <- r$top[, c("group", "gene", "peak_id")]
        o$rule <- label; o$draw <- 0L
        o
    }
    lpC  <- lpVariant("LocusPick +change")
    lpCR <- lpVariant("LocusPick +change +reps", replicates = reps)

    ## ---- Shortcuts: nearest TSS, no checks -------------------------------------
    nt <- nearestTSS(pk, gm)
    sig <- signalAtPeaks(pk, signal)
    nt <- merge(nt, sig, by = "peak_id")
    grp <- res$groups[, c("gene", "group", "lfc_condition")]
    nt <- merge(nt, grp, by = "gene")
    nt$group <- as.character(nt$group)
    a <- nt[[paste0(cells[1], "_mean")]]; b <- nt[[paste0(cells[2], "_mean")]]
    dirn <- ifelse(grepl("_up$", nt$group), 1, -1)
    sc <- list(top_signal = pmax(a, b, na.rm = TRUE),
               top_lfc    = abs(nt$lfc_condition) + 1e-6 * pmax(a, b, na.rm = TRUE),
               top_change = dirn * (a - b))
    cut <- lapply(names(sc), function(r) {
        p <- pickBy(nt, sc[[r]])[, c("group", "gene", "peak_id")]
        p$rule <- r; p$draw <- 0L; p
    })
    rnd <- lapply(seq_len(if (TEST) 20 else 200), function(i) {
        p <- pickBy(nt, stats::runif(nrow(nt)))[, c("group", "gene", "peak_id")]
        p$rule <- "random"; p$draw <- i; p
    })
    picks <- do.call(rbind, c(list(lp, lpC, lpCR), cut, rnd))

    ## ---- Judge every pick the same way -------------------------------------------
    u <- unique(picks[, c("peak_id", "gene")])
    fl <- pairFlags(u, pk, gm, signal, mappability = REF$hg38$umap, blacklist = bl)
    picks <- merge(picks, fl, by = c("peak_id", "gene"), sort = FALSE)
    picks$far_from_gene   <- !picks$near_gene
    picks$low_mappability <- !picks$mappable
    picks$no_signal       <- !picks$has_signal
    picks$any_flag <- with(picks, far_from_gene | low_mappability | no_signal | blacklisted)
    rr <- stats::setNames(as.list(c(reps[[1]], reps[[2]])), c("a1", "a2", "b1", "b2"))
    rs <- signalAtPeaks(pk[pk$peak_id %in% u$peak_id], rr)
    rs$not_reproduced <- sign(rs$a1_mean - rs$b1_mean) != sign(rs$a2_mean - rs$b2_mean)
    picks$not_reproduced <- rs$not_reproduced[match(picks$peak_id, rs$peak_id)]
    picks$factor <- tf
    allPicks[[tf]] <- picks

    ## ---- Example panel: best shortcut pick vs LocusPick pick, first group -----------
    g1 <- intersect(levels(factor(res$groups$group)), unique(picks$group))[1]
    ex <- rbind(picks[picks$rule == "top_signal" & picks$group == g1 & picks$any_flag, ][1, ],
                picks[picks$rule == "LocusPick" & picks$group == g1, ][1, ])
    ex <- ex[!is.na(ex$peak_id), ]
    if (nrow(ex)) {
        pdf(file.path(out, paste0("fig3_examples_", tf, ".pdf")), width = 8, height = 5)
        for (i in seq_len(nrow(ex))) {
            p <- pk[pk$peak_id == ex$peak_id[i]]
            loc <- data.frame(seqnames = as.character(seqnames(p)), start = start(p), end = end(p),
                              gene = ex$gene[i], peak_id = ex$peak_id[i])
            print(plotLocus(loc, signal, gm, trackGroups = list(`Fold change` = cells),
                            maxWindow = 6e4) +
                      plot_annotation(title = sprintf("%s | %s | rule: %s | flags: %s", tf, g1, ex$rule[i],
                                                      paste(names(FLAGS)[unlist(ex[i, names(FLAGS)])],
                                                            collapse = ", ")),
                                      caption = if (TEST) "TEST MODE"))
        }
        dev.off()
    }
}
picks <- do.call(rbind, allPicks)
writeTSV(picks, file.path(out, "picks.tsv"))

## ---- Failure rates with Wilson 95 % CIs --------------------------------------------------
rates <- function(d, by) {
    long <- do.call(rbind, lapply(c(names(FLAGS), "any_flag"), function(f)
        data.frame(d[, by, drop = FALSE], flag = f, fail = d[[f]])))
    long <- long[!is.na(long$fail), ]
    agg <- stats::aggregate(fail ~ ., data = long, FUN = function(x) c(x = sum(x), n = length(x)))
    agg <- cbind(agg[, setdiff(names(agg), "fail")], x = agg$fail[, "x"], n = agg$fail[, "n"])
    cbind(agg, wilson(agg$x, agg$n))
}
## Random draws are pooled; each draw contributes N_RULE picks per group.
r_all <- rates(picks, "rule")
r_tf  <- rates(picks, c("factor", "rule"))
writeTSV(r_all, file.path(out, "failure_rates.tsv"))
writeTSV(r_tf, file.path(out, "failure_rates_by_factor.tsv"))
print(r_all[r_all$flag == "any_flag", ])

lab <- c(FLAGS[c(1, 2, 4)], any_flag = "Any of these three*", FLAGS[5])
r_all <- r_all[r_all$flag %in% names(lab), ]          # no_signal cannot fail here (see header)
r_all$flag <- factor(lab[r_all$flag], levels = lab)
r_all$rule <- factor(r_all$rule, levels = RULES)
p <- ggplot(r_all, aes(rule, rate, colour = rule)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) +
    facet_wrap(~flag, nrow = 1) +
    scale_colour_manual(values = PAL, guide = "none") +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(x = NULL, y = "Picks failing the check (95% CI)",
         subtitle = sprintf("%d factors, top %d per group; random = %s draws",
                            length(factors), N_RULE, if (TEST) 20 else 200),
         caption = paste("The first three are LocusPick filters, so all LocusPick rules are ~0 there by design.",
                         "*Any flag also includes signal-track gaps (none on ENCODE data).",
                         "**Fair comparison except for '+change +reps', which filters on it.",
                         if (TEST) "TEST MODE: simulated data." else "", sep = "\n")) +
    theme_paper() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
saveFig(p, file.path(out, "fig3_failure_rates"), width = 10.5, height = 3.8)
say("Fig. 3 done: ", out)
