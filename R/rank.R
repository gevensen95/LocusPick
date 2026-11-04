#' Rank candidate loci
#'
#' Scores eligible candidates within each group as a weighted mean of
#' percentile ranks. Weights are explicit, so the selection rule can be stated
#' in a figure legend or methods section.
#'
#' A candidate missing any weighted value receives no score: it is never
#' ranked on a subset of the criteria.
#'
#' @param cand Output of [buildCandidates()].
#' @param weights Named numeric vector of column weights. Positive weights
#'   reward high values; negative weights reward low values. Alternatively a
#'   named list of such vectors, one per group (names matching `cand$group`),
#'   to use different criteria for different groups.
#' @param by Column defining groups.
#' @return `cand` with `score` and `rank` (within group; 1 is best), sorted by
#'   group and rank, plus one `pct_<column>` per weighted column: the percentile
#'   (0-1, higher value = higher percentile, whatever the weight's sign) of the
#'   candidate's value among its group's eligible candidates. These show how
#'   typical or extreme a pick is and can be reported with a figure.
#'   Ineligible candidates get `NA`.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal,
#'                         mappability = ex$mappability)
#' ranked <- rankLoci(cand, c(A_mean = 1, abs_lfc_perturbation = 1,
#'                            abs_lfc_condition = 1))
#' head(ranked[, c("group", "gene", "score", "rank")])
#' @export
rankLoci <- function(cand, weights, by = "group") {
    .stopIfMissingCols(cand, c(by, "eligible"), "candidates")
    grp <- as.character(cand[[by]])
    getW <- function(g) {
        if (is.list(weights)) {
            if (!g %in% names(weights)) stop("No weights for group '", g, "'.", call. = FALSE)
            weights[[g]]
        } else weights
    }
    allW <- if (is.list(weights)) unlist(unname(weights)) else weights
    if (is.null(names(allW)) || any(names(allW) == ""))
        stop("weights must be named by column.", call. = FALSE)
    .stopIfMissingCols(cand, unique(names(allW)), "candidates")

    cand$score <- NA_real_
    for (col in unique(names(allW))) cand[[paste0("pct_", col)]] <- NA_real_
    for (g in unique(grp[!is.na(grp)])) {
        idx <- which(grp == g & cand$eligible)
        if (!length(idx)) next
        w <- getW(g)
        for (col in names(w))
            cand[[paste0("pct_", col)]][idx] <- .pctRank(as.numeric(cand[[col]][idx]))
        terms <- vapply(names(w), function(col)
            .pctRank(sign(w[[col]]) * as.numeric(cand[[col]][idx])),
            numeric(length(idx)))
        terms <- matrix(terms, nrow = length(idx))
        cand$score[idx] <- as.numeric(terms %*% abs(w)) / sum(abs(w))
    }
    cand$rank <- stats::ave(-cand$score, grp, FUN = function(s)
        rank(s, na.last = "keep", ties.method = "min"))
    cand <- cand[order(cand[[by]], cand$rank, na.last = TRUE), ]
    rownames(cand) <- NULL
    attr(cand, "weights") <- weights
    cand
}

#' Select the top loci per group
#'
#' @param ranked Output of [rankLoci()].
#' @param n Number of loci per group.
#' @param onePerGene Keep only the best peak per gene.
#' @param by Column defining groups.
#' @return The top rows of `ranked`.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' topLoci(rankLoci(cand, c(A_mean = 1)), n = 2)
#' @export
topLoci <- function(ranked, n = 10, onePerGene = TRUE, by = "group") {
    .stopIfMissingCols(ranked, c(by, "gene", "score"), "ranked candidates")
    r <- ranked[!is.na(ranked$score), ]
    r <- r[order(r[[by]], -r$score), ]
    if (onePerGene) r <- r[!duplicated(r[, c(by, "gene")]), ]
    r <- do.call(rbind, lapply(split(r, r[[by]], drop = TRUE), utils::head, n))
    rownames(r) <- NULL
    r
}

#' Is a signal change specific to a group?
#'
#' A locus showing reduced binding in one condition is only informative if
#' the reduction is not genome-wide. For each group, this summarizes the
#' per-peak log2 ratio of two signal columns, tests it against zero, and,
#' with a background set of peaks, tests whether it differs from the
#' background.
#'
#' Tests are two-sided Wilcoxon tests: signed-rank against zero
#' (`p_vs_zero`) and rank-sum against the background ratios
#' (`p_vs_background`). Groups with fewer than 3 peaks get `NA`. A group whose
#' shift is significant against zero but not against the background shows a
#' genome-wide change, not a group-specific one.
#'
#' @param cand A candidate table (or any `data.frame` with the signal columns).
#' @param a,b Names of the two signal columns to compare (e.g. young and old).
#' @param by Grouping column.
#' @param background Optional `data.frame` with columns `a` and `b` for all
#'   peaks, e.g. from [signalAtPeaks()] on the full peak set.
#' @param pseudocount Added before taking log2 ratios.
#' @return `data.frame` with, per group: `n`, `frac_a_gt_b`,
#'   `median_log2_a_over_b`, `p_vs_zero`, and with a background
#'   `shift_vs_background` (difference of medians) and `p_vs_background`.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' signalShift(cand, "A_mean", "B_mean",
#'             background = signalAtPeaks(ex$peaks, ex$signal))
#' @export
signalShift <- function(cand, a, b, by = "group", background = NULL,
                        pseudocount = 0.01) {
    .stopIfMissingCols(cand, c(a, b, by), "candidates")
    lr <- function(d) {
        x <- d[[a]]; y <- d[[b]]
        ok <- !is.na(x) & !is.na(y)
        list(r = log2((x[ok] + pseudocount) / (y[ok] + pseudocount)),
             gt = mean(x[ok] > y[ok]))
    }
    wt <- function(...) tryCatch(suppressWarnings(stats::wilcox.test(...)$p.value),
                                 error = function(e) NA_real_)
    bg <- NULL
    if (!is.null(background)) {
        .stopIfMissingCols(background, c(a, b), "background")
        bg <- lr(background)
    }
    summ <- function(L, label, isBg = FALSE) {
        out <- data.frame(group = label, n = length(L$r), frac_a_gt_b = L$gt,
                          median_log2_a_over_b = stats::median(L$r),
                          p_vs_zero = if (length(L$r) >= 3) wt(L$r, mu = 0) else NA_real_)
        if (!is.null(bg)) {
            out$shift_vs_background <- if (isBg) NA_real_ else
                stats::median(L$r) - stats::median(bg$r)
            out$p_vs_background <- if (isBg || length(L$r) < 3) NA_real_ else wt(L$r, bg$r)
        }
        out
    }
    u <- cand[!duplicated(cand[, c("peak_id", by)]), ]
    parts <- lapply(split(u, u[[by]], drop = TRUE), function(d)
        summ(lr(d), as.character(d[[by]][1])))
    if (!is.null(bg)) parts$background <- summ(bg, "background (all peaks)", isBg = TRUE)
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}
