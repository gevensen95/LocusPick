#' Prepare expression data for plot insets
#'
#' Converts expression objects into the long table used by [plotLocus()]
#' (`gene`, `panel`, `group`, `value`).
#'
#' * `DESeqDataSet`: normalized counts (size factors are estimated if needed).
#' * `DGEList`: counts per million.
#' * `SummarizedExperiment`: the first assay, or the one named by `assay`.
#' * matrix or `data.frame` of genes × samples: used as is.
#'
#' @param x An expression object or matrix (genes in rows, samples in columns).
#' @param group Sample groups: a column name of the sample annotation
#'   (`colData` / `samples`), or a vector with one entry per sample.
#' @param panel Label for this experiment (one facet in the inset).
#' @param levels Optional order of groups (control first).
#' @param orgDb Optional `OrgDb` to convert Ensembl/Entrez IDs to symbols and
#'   harmonize outdated symbols.
#' @param assay For `SummarizedExperiment`: which assay to use.
#' @return A long `data.frame`. Combine several experiments with `rbind()`.
#' @examples
#' m <- matrix(c(10, 12, 30, 32, 5, 6, 2, 1), nrow = 2,
#'             dimnames = list(c("geneA", "geneB"), paste0("s", 1:4)))
#' expressionFrom(m, group = c("WT", "WT", "KO", "KO"), panel = "Knockout")
#' @export
expressionFrom <- function(x, group, panel = "Expression", levels = NULL,
                           orgDb = NULL, assay = NULL) {
    sampleInfo <- NULL
    if (is(x, "DESeqDataSet")) {
        sf <- .requireCall("DESeq2", "sizeFactors")(x)
        if (is.null(sf)) x <- .requireCall("DESeq2", "estimateSizeFactors")(x)
        sampleInfo <- as.data.frame(.requireCall("SummarizedExperiment", "colData")(x))
        mat <- .requireCall("DESeq2", "counts")(x, normalized = TRUE)
    } else if (inherits(x, "DGEList")) {
        sampleInfo <- x$samples
        mat <- .requireCall("edgeR", "cpm")(x)
    } else if (is(x, "SummarizedExperiment")) {
        se <- "SummarizedExperiment"
        sampleInfo <- as.data.frame(.requireCall(se, "colData")(x))
        mat <- if (is.null(assay)) .requireCall(se, "assay")(x) else
            .requireCall(se, "assay")(x, assay)
    } else if (is.matrix(x) || is.data.frame(x)) {
        mat <- as.matrix(x)
    } else {
        stop("Unsupported expression input of class '", class(x)[1], "'.", call. = FALSE)
    }
    mat <- as.matrix(mat)
    if (is.null(rownames(mat))) stop("Expression matrix needs gene row names.", call. = FALSE)

    if (length(group) == 1 && ncol(mat) > 1) {
        if (is.null(sampleInfo) || !group %in% names(sampleInfo))
            stop("Sample annotation has no column '", group, "'.", call. = FALSE)
        grp <- sampleInfo[[group]]
        if (is.factor(grp) && is.null(levels)) levels <- levels(grp)
        grp <- as.character(grp)
    } else {
        if (length(group) != ncol(mat))
            stop("group has ", length(group), " entries but there are ", ncol(mat),
                 " samples.", call. = FALSE)
        grp <- as.character(group)
    }
    if (!is.null(levels)) {
        o <- order(match(grp, levels))
        mat <- mat[, o, drop = FALSE]; grp <- grp[o]
    }

    genes <- rownames(mat)
    idType <- .idType(genes)
    if (idType != "symbol") {
        if (is.null(orgDb)) stop("Row names look like ", idType, " IDs; pass `orgDb`.",
                                 call. = FALSE)
        genes <- .idsToSymbols(genes, orgDb, idType)
    } else if (!is.null(orgDb)) {
        genes <- harmonizeSymbols(genes, orgDb)
    }
    keep <- !is.na(genes) & !duplicated(genes)
    mat <- mat[keep, , drop = FALSE]; genes <- genes[keep]
    data.frame(gene = rep(genes, times = ncol(mat)), panel = panel,
               group = rep(grp, each = nrow(mat)), value = as.vector(mat),
               stringsAsFactors = FALSE)
}
