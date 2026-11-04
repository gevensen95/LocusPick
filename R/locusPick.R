#' Run the full LocusPick workflow in one call
#'
#' Reads tool outputs directly, checks the inputs, builds and ranks candidate
#' loci, and tests whether signal changes are specific. Each step is also
#' available as its own function for finer control.
#'
#' @param peaks Peak file(s) or object (see [readPeaks()]). Several files are
#'   merged into consensus peaks.
#' @param perturbation,condition DE results (see [readDE()]): DESeq2, edgeR or
#'   limma objects or tables. Fold changes must be test / reference. Either
#'   or both may be omitted: with one contrast, groups follow its direction
#'   (see [defineGroups()]); with none, all linked genes form one group and
#'   loci are ranked on signal alone.
#' @param genes Gene models: GTF path, `GRanges`, or [readGeneModels()] output.
#' @param signal Optional named list of tracks. Each element can be a bigWig or
#'   BAM path, several paths (averaged), or a signal source.
#' @param mappability Optional mappability bigWig (see [mappabilityAtPeaks()]).
#' @param orgDb Optional `OrgDb` for ID conversion and symbol harmonization,
#'   e.g. `org.Mm.eg.db::org.Mm.eg.db`. Recommended.
#' @param expression Optional expression for plot insets: a long table from
#'   [expressionFrom()], or a list of them.
#' @param weights Ranking weights (see [rankLoci()]). Default: mean signal of
#'   the first track plus the absolute fold changes of both contrasts.
#' @param trackGroups Panels for plots (see [plotLocus()]). Default: tracks
#'   sharing a prefix before `_` (e.g. `CNR_young`, `CNR_old`) share a panel.
#' @param perturbationDirection,conditionDirection Optional labels confirming
#'   the contrast directions of plain tables, e.g. `"KO/WT"`, `"old/young"`.
#' @param alpha,minLFC,perturbationType,conditionName Passed to
#'   [defineGroups()].
#' @param links Optional peak-gene links (e.g. `Signac::Links()`), used
#'   instead of distance-based assignment; see [buildCandidates()].
#' @param maxDist,maxLinkDist,minMappability,minCoverage,geneCol,assignMode
#'   Passed to [buildCandidates()].
#' @param minReps,blacklist Passed to [readPeaks()].
#' @param n Number of top loci per group.
#' @return An object of class `LocusPickResult`: a list with `candidates`
#'   (ranked, all flags), `top`, `qc`, `shift`, `checks`, `groups`, and the
#'   inputs needed for plotting. Use `plot(res, i)` for one locus and
#'   `reportLoci(res, "file.pdf")` for a review report.
#' @examples
#' f <- function(x) system.file("extdata", x, package = "LocusPick")
#' de <- read.delim(f("toy_de.tsv"))
#' res <- locusPick(
#'     peaks = f("toy_peaks.bed"),
#'     perturbation = data.frame(gene = de$gene, log2FoldChange = de$ko_lfc,
#'                               padj = de$ko_padj),
#'     condition = data.frame(gene = de$gene, logFC = de$cond_lfc,
#'                            FDR = de$cond_padj),
#'     perturbationDirection = "KO/WT", conditionDirection = "old/young",
#'     genes = f("toy.gtf"),
#'     signal = list(young = f("toy_A.bw"), old = f("toy_B.bw")),
#'     mappability = f("toy_map.bw"))
#' res
#' @export
locusPick <- function(peaks, perturbation = NULL, condition = NULL, genes, signal = NULL,
                      mappability = NULL, orgDb = NULL, expression = NULL,
                      weights = NULL, trackGroups = NULL,
                      perturbationDirection = NULL, conditionDirection = NULL,
                      alpha = 0.05, minLFC = 0, perturbationType = c("loss", "gain"),
                      conditionName = "condition", links = NULL, maxDist = 5000,
                      maxLinkDist = 5e5,
                      minMappability = 0.8, minCoverage = 0.5, geneCol = NULL,
                      assignMode = c("nearest", "all"), minReps = NULL,
                      blacklist = NULL, n = 10) {
    perturbationType <- match.arg(perturbationType)
    assignMode <- match.arg(assignMode)
    harm <- if (!is.null(orgDb)) function(s) harmonizeSymbols(s, orgDb) else NULL

    message("== Reading inputs")
    pk <- readPeaks(peaks, minReps = minReps, blacklist = blacklist)
    if (is.null(geneCol) && "homer_gene" %in% colnames(S4Vectors::mcols(pk)))
        message("  Peaks carry a HOMER gene annotation; LocusPick assigns genes itself. ",
                "Use geneCol = \"homer_gene\" to check HOMER's assignment instead.")
    gm <- if (inherits(genes, "locusGeneModels")) genes else
        readGeneModels(genes, harmonize = harm)
    pert <- if (!is.null(perturbation)) .readDEfor(perturbation, orgDb, perturbationDirection)
    cond <- if (!is.null(condition)) .readDEfor(condition, orgDb, conditionDirection)
    groups <- if (is.null(pert) && is.null(cond)) {
        message("  No DE results given: all linked genes form one group, ranked on signal.")
        NULL
    } else {
        defineGroups(pert, cond, alpha = alpha, minLFC = minLFC,
                     perturbationType = perturbationType, conditionName = conditionName)
    }
    if (!is.null(signal)) {
        if (is.null(names(signal)) || any(names(signal) == ""))
            stop("signal must be a named list, e.g. list(young = ..., old = ...).", call. = FALSE)
        signal <- lapply(signal, .asSource)
    }
    if (is.list(expression) && !is.data.frame(expression))
        expression <- do.call(rbind, expression)

    message("== Checking inputs")
    checks <- checkInputs(pk, gm, signal, groups, maxDist = maxDist)
    if (any(checks$status == "fail"))
        warning("Some input checks failed; see res$checks before trusting the results.",
                call. = FALSE)

    message("== Building candidates")
    cand <- buildCandidates(pk, gm, groups, signal = signal, mappability = mappability,
                            geneCol = geneCol, links = links, maxDist = maxDist,
                            maxLinkDist = maxLinkDist,
                            minMappability = minMappability, minCoverage = minCoverage,
                            assignMode = assignMode)
    tr <- names(signal)
    if (length(tr) >= 2)
        cand$signal_change <- cand[[paste0(tr[1], "_mean")]] - cand[[paste0(tr[2], "_mean")]]
    if (is.null(weights)) {
        weights <- c(abs_lfc_perturbation = 1, abs_lfc_condition = 1)
        weights <- weights[names(weights) %in% names(cand)]
        if (length(tr)) weights <- c(stats::setNames(1, paste0(tr[1], "_mean")), weights)
        if (!length(weights))
            stop("Nothing to rank on: give signal tracks, DE results, or `weights`.",
                 call. = FALSE)
    }
    message(paste(strwrap(paste0("  ranking weights: ", paste(names(weights), unlist(weights),
                                                           sep = " = ", collapse = ", ")),
                          width = 78, exdent = 4), collapse = "\n"))
    ranked <- rankLoci(cand, weights)
    top <- topLoci(ranked, n = n)
    qc <- qcSummary(cand)

    shift <- NULL
    if (length(tr) >= 2) {
        bg <- signalAtPeaks(pk, signal[tr[1:2]])
        shift <- signalShift(cand, paste0(tr[1], "_mean"), paste0(tr[2], "_mean"),
                             background = bg)
    }
    if (is.null(trackGroups) && length(tr)) trackGroups <- .defaultTrackGroups(tr)

    structure(list(candidates = ranked, top = top, qc = qc, shift = shift,
                   checks = checks, groups = groups, peaks = pk, genes = gm,
                   signal = signal, expression = expression, trackGroups = trackGroups,
                   weights = weights,
                   contrasts = c(perturbation = if (is.null(pert)) "none" else .contrastLabel(pert),
                                 condition = if (is.null(cond)) "none" else .contrastLabel(cond))),
              class = "LocusPickResult")
}

.readDEfor <- function(x, orgDb, direction) {
    std <- is.data.frame(x) && all(c("gene", "lfc", "padj") %in% names(x))
    if (std) {
        if (!is.null(orgDb)) x$gene <- harmonizeSymbols(x$gene, orgDb)
        if (!is.null(direction)) attr(x, "contrast") <- direction
        return(x)
    }
    readDE(x, orgDb = orgDb, direction = direction)
}

.contrastLabel <- function(d) {
    a <- attr(d, "contrast")
    if (is.null(a)) "not stated" else a
}

.defaultTrackGroups <- function(tracks) {
    pre <- ifelse(grepl("_", tracks), sub("_[^_]*$", "", tracks), "Signal")
    split(tracks, factor(pre, levels = unique(pre)))
}

#' @export
print.LocusPickResult <- function(x, ...) {
    cat("LocusPick result\n")
    cat("  perturbation:", x$contrasts[["perturbation"]], "\n")
    cat("  condition:   ", x$contrasts[["condition"]], "\n")
    if (!is.null(x$cellTypes)) cat("  cell types:", paste(x$cellTypes, collapse = ", "), "\n")
    cat("  peaks:", length(x$peaks), "| candidate peak-gene pairs:",
        nrow(x$candidates), "| eligible:", sum(x$candidates$eligible), "\n")
    cat(strwrap(paste("  weights:", paste(names(x$weights), unlist(x$weights), sep = "=",
                                          collapse = ", ")), width = 78, exdent = 4),
        sep = "\n")
    cat("\n")
    print(x$qc, row.names = FALSE)
    if (!is.null(x$shift)) {
        cat("\nSignal change (first vs second track):\n")
        print(x$shift, row.names = FALSE, digits = 3)
    }
    cat("\nTop loci:\n")
    show <- intersect(c("group", "gene", "peak_id", "feature", "score"), names(x$top))
    print(utils::head(x$top[, show], 12), row.names = FALSE, digits = 3)
    if (any(x$checks$status != "ok"))
        cat("\nNote: some input checks were not ok; see $checks.\n")
    invisible(x)
}

#' Plot a locus from a LocusPick result
#'
#' @param x A `LocusPickResult`.
#' @param which Row of `x$top` to plot, or a gene symbol.
#' @param ... Passed to [plotLocus()], e.g. `annotate = FALSE`, `colors`,
#'   `yLimits`.
#' @return A `patchwork` object.
#' @examples
#' f <- function(x) system.file("extdata", x, package = "LocusPick")
#' de <- read.delim(f("toy_de.tsv"))
#' res <- suppressMessages(locusPick(
#'     peaks = f("toy_peaks.bed"),
#'     perturbation = deTable(de, lfc = "ko_lfc", padj = "ko_padj"),
#'     condition = deTable(de, lfc = "cond_lfc", padj = "cond_padj"),
#'     genes = f("toy.gtf"), signal = list(young = f("toy_A.bw"), old = f("toy_B.bw"))))
#' plot(res, "geneA")
#' @export
plot.LocusPickResult <- function(x, which = 1, ...) {
    row <- if (is.character(which)) {
        i <- match(which, x$top$gene)
        if (is.na(i)) {
            r <- x$candidates[x$candidates$gene == which & !is.na(x$candidates$score), ]
            if (!nrow(r)) stop("No eligible locus for gene '", which, "'.", call. = FALSE)
            r[1, , drop = FALSE]
        } else x$top[i, , drop = FALSE]
    } else x$top[which, , drop = FALSE]
    if (is.null(x$signal)) stop("No signal tracks in this result.", call. = FALSE)
    expr <- x$expression
    if (!is.null(expr) && "cell_type" %in% names(row) && row$cell_type %in% expr$panel)
        expr <- expr[expr$panel == row$cell_type, , drop = FALSE]
    plotLocus(row, x$signal, x$genes, trackGroups = x$trackGroups,
              expression = expr, ...)
}
