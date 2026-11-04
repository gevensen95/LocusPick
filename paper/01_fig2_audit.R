## =============================================================================
## Fig. 2. How are representative genome-browser loci chosen in published papers?
##
## Input: templates/audit_rubric.csv, filled in by two raters (see
## templates/AUDIT_GUIDE.md). Rows with rater = "consensus" override the two
## raters for that panel.
## Outputs (results/fig2/): agreement.tsv (Cohen's kappa per item),
## disagreements.tsv (to resolve), item_rates.tsv (share "yes", Wilson 95% CI),
## fig2.pdf/png
## =============================================================================

source("config.R"); source("R/common.R")
out <- figDir("fig2")
set.seed(SEED)
ITEMS <- c(selection_rule_stated = "Selection rule stated", quantification_shown = "Quantification shown",
           normalization_stated = "Normalization stated", global_change_addressed = "Global change addressed",
           peak_gene_rule_stated = "Peak-gene rule stated", replicates_shown = "Replicates shown",
           qc_mentioned = "QC mentioned")

rubric <- if (TEST) {
    n <- 60
    truth <- as.data.frame(lapply(ITEMS, function(i) sample(c("yes", "no"), n, TRUE, prob = c(0.3, 0.7))))
    truth$selection_rule_stated <- sample(c("yes", "partial", "no"), n, TRUE, prob = c(0.1, 0.2, 0.7))
    rate <- function(r) { x <- truth; flip <- matrix(stats::runif(n * length(ITEMS)) < 0.1, n)
        for (j in seq_along(ITEMS)) x[flip[, j], j] <- sample(c("yes", "no"), sum(flip[, j]), TRUE)
        cbind(panel_id = sprintf("P%03d", seq_len(n)), assay = sample(c("ChIP", "CUT&RUN", "ATAC"), n, TRUE),
              n_loci_shown = sample(1:4, n, TRUE), rater = r, x) }
    set.seed(SEED); a <- rate("A"); set.seed(SEED + 1); b <- rate("B"); b$assay <- a$assay; b$n_loci_shown <- a$n_loci_shown
    rbind(a, b)
} else {
    r <- utils::read.csv(needFiles(file.path(PAPER, "templates/audit_rubric.csv")), stringsAsFactors = FALSE)
    if (!nrow(r)) stop("templates/audit_rubric.csv has no rows yet.", call. = FALSE)
    r
}
rubric[names(ITEMS)] <- lapply(rubric[names(ITEMS)], function(x) tolower(trimws(as.character(x))))
raters <- setdiff(unique(rubric$rater), "consensus")
if (length(raters) != 2) stop("Expected exactly two raters (plus optional 'consensus'); found: ",
                              paste(raters, collapse = ", "), call. = FALSE)
A <- rubric[rubric$rater == raters[1], ]; B <- rubric[rubric$rater == raters[2], ]
B <- B[match(A$panel_id, B$panel_id), ]
if (any(is.na(B$panel_id))) stop("Some panels were scored by only one rater.", call. = FALSE)

## Agreement per item
agree <- do.call(rbind, lapply(names(ITEMS), function(i) {
    a <- A[[i]]; b <- B[[i]]; a[a == "na"] <- NA; b[b == "na"] <- NA
    cbind(item = i, cohenKappa(a, b))
}))
writeTSV(agree, file.path(out, "agreement.tsv"))
print(agree)

## Final scores: agreed values, consensus rows where given
final <- A
for (i in names(ITEMS)) final[[i]][A[[i]] != B[[i]]] <- NA
cons <- rubric[rubric$rater == "consensus", ]
for (k in seq_len(nrow(cons))) {
    j <- match(cons$panel_id[k], final$panel_id)
    for (i in names(ITEMS)) if (!is.na(cons[[i]][k]) && cons[[i]][k] != "") final[[i]][j] <- cons[[i]][k]
}
dis <- do.call(rbind, lapply(names(ITEMS), function(i) {
    j <- which(is.na(final[[i]]) & A[[i]] != B[[i]])
    if (length(j)) data.frame(panel_id = final$panel_id[j], item = i, rater_1 = A[[i]][j], rater_2 = B[[i]][j])
}))
if (!is.null(dis)) { writeTSV(dis, file.path(out, "disagreements.tsv"))
    message(nrow(dis), " unresolved disagreements: add 'consensus' rows to the rubric.") }

rates <- do.call(rbind, lapply(names(ITEMS), function(i) {
    v <- final[[i]]; v <- v[!is.na(v) & v != "na"]
    cbind(item = i, label = ITEMS[[i]], x = sum(v == "yes"), n = length(v), wilson(sum(v == "yes"), length(v)))
}))
writeTSV(rates, file.path(out, "item_rates.tsv"))

rates$label <- factor(rates$label, rev(ITEMS))
pA <- ggplot(rates, aes(rate, label)) +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = PAL[["LocusPick"]], size = 0.3) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Panels (95% CI)", y = NULL,
         title = sprintf("A  Reporting in %d figure panels", nrow(final))) + theme_paper()
nl <- as.data.frame(table(pmin(as.integer(final$n_loci_shown), 4)))
nl$Var1 <- factor(nl$Var1, labels = sub("^4$", "4+", levels(nl$Var1)))
pB <- ggplot(nl, aes(Var1, Freq)) + geom_col(fill = "grey60", width = 0.6) +
    labs(x = "Example loci per panel", y = "Panels", title = "B  Loci shown") + theme_paper()
p <- (pA | pB) + plot_layout(widths = c(2, 1)) +
    plot_annotation(caption = if (TEST) "TEST MODE: simulated data, numbers are meaningless")
saveFig(p, file.path(out, "fig2"), width = 8, height = 3)
say("Fig. 2 done: ", out)
