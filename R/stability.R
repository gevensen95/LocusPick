#' How stable are the top loci?
#'
#' `weightStability()` re-ranks the candidates many times with every weight
#' multiplied by a random factor (log-uniform within `range`) and reports how
#' often each locus stays in the top `n` of its group, and the mean Jaccard
#' overlap of the perturbed top `n` with the original.
#'
#' `topOverlap()` compares rankings made from different inputs, for example
#' leaving out one replicate at a time or using another normalization, by the
#' Jaccard overlap of their top `n` loci per group with a reference ranking.
#'
#' @param cand Output of [buildCandidates()] (or [rankLoci()]).
#' @param weights Weights as in [rankLoci()]: a named vector, or a named list
#'   of vectors per group (each weight in each group is perturbed independently).
#' @param n Size of the top set per group.
#' @param B Number of perturbations.
#' @param range Range of the multiplicative factor applied to each weight.
#' @param by Grouping column.
#' @param seed Random seed.
#' @return `weightStability()`: a list with `loci` (`group`, `gene`,
#'   `peak_id`, `rank`, `top_freq`) and `summary` (`group`, `mean_jaccard`).
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal)
#' ws <- weightStability(cand, c(A_mean = 1, abs_lfc_perturbation = 1), n = 1, B = 50)
#' ws$summary
#' @export
weightStability <- function(cand, weights, n = 10, B = 200, range = c(0.5, 2),
                            by = "group", seed = 1) {
    perturb <- function() {
        draw <- function(w) w * exp(stats::runif(length(w), log(range[1]), log(range[2])))
        if (is.list(weights)) lapply(weights, draw) else draw(weights)
    }
    base <- topLoci(rankLoci(cand, weights, by = by), n = n, by = by)
    key <- function(t) paste(t[[by]], t$peak_id, t$gene, sep = "\r")
    baseKey <- key(base)
    set.seed(seed)
    hits <- stats::setNames(numeric(0), character(0))
    jac <- matrix(NA_real_, B, length(unique(base[[by]])),
                  dimnames = list(NULL, as.character(unique(base[[by]]))))
    for (i in seq_len(B)) {
        t <- topLoci(rankLoci(cand, perturb(), by = by), n = n, by = by)
        k <- key(t)
        hits[k] <- ifelse(is.na(hits[k]), 0, hits[k]) + 1
        for (g in colnames(jac)) {
            A <- baseKey[base[[by]] == g]; Bk <- k[t[[by]] == g]
            jac[i, g] <- length(intersect(A, Bk)) / length(union(A, Bk))
        }
    }
    allR <- rankLoci(cand, weights, by = by)
    allR <- allR[!is.na(allR$score), ]
    allR <- allR[!duplicated(allR[, c(by, "gene")]), ]
    k <- key(allR)
    loci <- data.frame(group = allR[[by]], gene = allR$gene, peak_id = allR$peak_id,
                       rank = allR$rank,
                       top_freq = unname(ifelse(is.na(hits[k]), 0, hits[k])) / B,
                       stringsAsFactors = FALSE)
    loci <- loci[order(loci$group, -loci$top_freq, loci$rank), ]
    rownames(loci) <- NULL
    list(loci = loci,
         summary = data.frame(group = colnames(jac), mean_jaccard = colMeans(jac),
                              row.names = NULL))
}

#' @rdname weightStability
#' @param rankings A named list of ranked candidate tables (from [rankLoci()]
#'   or `res$candidates`), e.g. one per left-out replicate.
#' @param reference Name or index of the reference ranking in `rankings`.
#' @return `topOverlap()`: a `data.frame` with `group`, `ranking` and
#'   `jaccard` (overlap of top-`n` loci with the reference).
#' @export
topOverlap <- function(rankings, n = 10, reference = 1, by = "group") {
    if (is.null(names(rankings))) names(rankings) <- paste0("ranking", seq_along(rankings))
    tops <- lapply(rankings, topLoci, n = n, by = by)
    key <- function(t) paste(t$peak_id, t$gene, sep = "\r")
    ref <- tops[[reference]]
    out <- list()
    for (nm in setdiff(names(tops), names(tops)[reference])) {
        t <- tops[[nm]]
        for (g in unique(as.character(ref[[by]]))) {
            A <- key(ref[ref[[by]] == g, ]); Bk <- key(t[t[[by]] == g, ])
            out[[length(out) + 1]] <- data.frame(
                group = g, ranking = nm,
                jaccard = if (length(union(A, Bk))) length(intersect(A, Bk)) /
                    length(union(A, Bk)) else NA_real_)
        }
    }
    do.call(rbind, out)
}
