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
#'   group and rank. Ineligible candidates get `NA`.
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
    for (g in unique(grp[!is.na(grp)])) {
        idx <- which(grp == g & cand$eligible)
        if (!length(idx)) next
        w <- getW(g)
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
#' the reduction is not genome-wide. This compares the change between two
#' signal columns within each group and, optionally, across a background set of
#' peaks.
#'
#' @param cand A candidate table (or any `data.frame` with the signal columns).
#' @param a,b Names of the two signal columns to compare (e.g. young and old).
#' @param by Grouping column.
#' @param background Optional `data.frame` with columns `a` and `b` for all
#'   peaks, e.g. from [signalAtPeaks()] on the full peak set.
#' @param pseudocount Added before taking log2 ratios.
#' @return `data.frame` with, per group: `n`, `frac_a_gt_b`, and
#'   `median_log2_a_over_b`.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' signalShift(cand, "A_mean", "B_mean")
#' @export
signalShift <- function(cand, a, b, by = "group", background = NULL,
                        pseudocount = 0.01) {
    .stopIfMissingCols(cand, c(a, b, by), "candidates")
    summ <- function(d, label) {
        x <- d[[a]]; y <- d[[b]]
        ok <- !is.na(x) & !is.na(y)
        data.frame(group = label, n = sum(ok),
                   frac_a_gt_b = mean(x[ok] > y[ok]),
                   median_log2_a_over_b = stats::median(
                       log2((x[ok] + pseudocount) / (y[ok] + pseudocount))))
    }
    u <- cand[!duplicated(cand[, c("peak_id", by)]), ]
    parts <- lapply(split(u, u[[by]], drop = TRUE), function(d) summ(d, as.character(d[[by]][1])))
    if (!is.null(background)) {
        .stopIfMissingCols(background, c(a, b), "background")
        parts$background <- summ(background, "background (all peaks)")
    }
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
}
