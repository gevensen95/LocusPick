#' Build the candidate locus table
#'
#' Joins peaks, gene groups, signal and mappability into one row per
#' peak-gene candidate, and adds quality-control flags. Nothing is dropped:
#' failing candidates stay in the table with `eligible = FALSE`, so the
#' selection can be audited.
#'
#' @param peaks Peaks (anything accepted by [asPeaks()]). Extra columns (for
#'   example a peak-set label) are carried into the output.
#' @param genes Gene models from [readGeneModels()].
#' @param groups Gene groups from [defineGroups()] (needs `gene` and `group`).
#' @param signal Optional named list of signal sources (see [signalSources]).
#' @param mappability Optional mappability source (see [mappabilityAtPeaks()]).
#' @param geneCol If the peaks already carry a gene assignment (for example
#'   from HOMER), the column holding it. That assignment is then checked
#'   against `genes` with [checkAssignment()]. If `NULL` (default), peaks are
#'   assigned with [assignPeaks()].
#' @param maxDist Maximum peak-to-gene-body distance (bp).
#' @param minMappability Minimum mean mappability for `mappable`.
#' @param minCoverage Minimum fraction of a peak covered by signal for
#'   `has_signal`.
#' @param requireSignal Names of `signal` tracks that must pass
#'   `minCoverage`. Defaults to all tracks.
#' @param assignMode When peaks are assigned by LocusPick: `"nearest"`
#'   (default) links each peak only to its closest gene, which avoids showing
#'   another gene's promoter peak as an example; `"all"` links it to every
#'   gene within `maxDist`. `n_nearby_genes` records how many genes were in
#'   range, so crowded loci can be spotted either way.
#' @return A `data.frame` with peak coordinates, `gene`, group columns, signal
#'   columns (`<track>_mean`, `<track>_cov`), `mappability`, absolute fold
#'   changes, and the flags `gene_in_models`, `near_gene`, `mappable`,
#'   `has_signal` and `eligible`.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal,
#'                         mappability = ex$mappability)
#' head(cand)
#' @export
buildCandidates <- function(peaks, genes, groups, signal = NULL, mappability = NULL,
                            geneCol = NULL, maxDist = 5000, minMappability = 0.8,
                            minCoverage = 0.5, requireSignal = names(signal),
                            assignMode = c("nearest", "all")) {
    assignMode <- match.arg(assignMode)
    pk <- asPeaks(peaks)
    gm <- .asGeneModels(genes)
    .stopIfMissingCols(groups, c("gene", "group"), "groups")

    if (is.null(geneCol)) {
        a <- assignPeaks(pk, gm, maxDist = maxDist, mode = "all")
        nNear <- table(a$peak_id)
        if (assignMode == "nearest") {
            best <- stats::ave(a$distance, a$peak_id, FUN = min)
            a <- a[a$distance == best, ]
        }
        a$n_nearby_genes <- as.integer(nNear[a$peak_id])
        a$gene_in_models <- TRUE
        a$near_gene <- TRUE
    } else {
        chk <- checkAssignment(pk, gm, geneCol = geneCol, maxDist = maxDist)
        a <- data.frame(peak_id = pk$peak_id,
                        seqnames = as.character(GenomicRanges::seqnames(pk)),
                        start = GenomicRanges::start(pk), end = GenomicRanges::end(pk),
                        gene = chk$gene, distance = chk$distance,
                        gene_in_models = chk$gene_in_models,
                        near_gene = chk$near_gene, stringsAsFactors = FALSE)
    }

    cand <- merge(a, groups, by = "gene")
    if (!nrow(cand)) {
        warning("No peaks link to any grouped gene.", call. = FALSE)
        return(cand)
    }

    # carry extra peak columns (skip ones already present)
    extra <- as.data.frame(S4Vectors::mcols(pk))
    extra <- extra[!duplicated(extra$peak_id), setdiff(names(extra), setdiff(names(cand), "peak_id")),
                   drop = FALSE]
    if (ncol(extra) > 1) cand <- merge(cand, extra, by = "peak_id", all.x = TRUE)

    upk <- pk[match(unique(cand$peak_id), pk$peak_id)]
    if (length(signal)) {
        sig <- signalAtPeaks(upk, signal)
        cand <- merge(cand, sig, by = "peak_id", all.x = TRUE)
        covCols <- paste0(requireSignal, "_cov")
        .stopIfMissingCols(cand, covCols, "signal summary")
        cand$has_signal <- Reduce(`&`, lapply(covCols, function(cc)
            !is.na(cand[[cc]]) & cand[[cc]] >= minCoverage))
    } else {
        cand$has_signal <- TRUE
    }
    if (!is.null(mappability)) {
        mp <- data.frame(peak_id = upk$peak_id,
                         mappability = mappabilityAtPeaks(upk, mappability))
        cand <- merge(cand, mp, by = "peak_id", all.x = TRUE)
        cand$mappable <- !is.na(cand$mappability) & cand$mappability >= minMappability
    } else {
        cand$mappability <- NA_real_
        cand$mappable <- TRUE
    }
    if ("lfc_perturbation" %in% names(cand))
        cand$abs_lfc_perturbation <- abs(cand$lfc_perturbation)
    if ("lfc_condition" %in% names(cand))
        cand$abs_lfc_condition <- abs(cand$lfc_condition)

    cand$eligible <- cand$gene_in_models & cand$near_gene & cand$mappable & cand$has_signal
    cand <- cand[order(cand$group, cand$gene, cand$peak_id), ]
    rownames(cand) <- NULL
    attr(cand, "settings") <- list(maxDist = maxDist, minMappability = minMappability,
                                   minCoverage = minCoverage,
                                   requireSignal = requireSignal)
    cand
}

#' Summarize quality-control flags
#'
#' @param cand Output of [buildCandidates()].
#' @return A `data.frame` with, per group, the number of candidate peaks and
#'   how many fail each check.
#' @examples
#' ex <- locusPickExample()
#' cand <- buildCandidates(ex$peaks, ex$genes, ex$groups, signal = ex$signal,
#'                         mappability = ex$mappability)
#' qcSummary(cand)
#' @export
qcSummary <- function(cand) {
    flags <- c("gene_in_models", "near_gene", "mappable", "has_signal")
    .stopIfMissingCols(cand, c("group", "eligible", flags), "candidates")
    g <- split(cand, cand$group, drop = TRUE)
    out <- do.call(rbind, lapply(names(g), function(k) {
        d <- g[[k]]
        data.frame(group = k, peaks = nrow(d), genes = length(unique(d$gene)),
                   fail_gene_in_models = sum(!d$gene_in_models),
                   fail_near_gene = sum(d$gene_in_models & !d$near_gene),
                   fail_mappable = sum(!d$mappable),
                   fail_has_signal = sum(!d$has_signal),
                   eligible = sum(d$eligible),
                   eligible_genes = length(unique(d$gene[d$eligible])))
    }))
    rownames(out) <- NULL
    out
}
