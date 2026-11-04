#' Aggregate single-cell counts into pseudobulk samples
#'
#' Sums counts over the cells of each sample (and cell type), giving one
#' pseudobulk profile per biological sample. Pseudobulk profiles are the
#' recommended input for differential expression between conditions: tests
#' that treat each cell as a replicate overstate significance.
#'
#' @param x A Seurat object, `SingleCellExperiment`, [read10x()] result, or a
#'   genes x cells count matrix.
#' @param sample Metadata column identifying biological samples.
#' @param cellType Optional metadata column of cell types; pseudobulks are made
#'   per cell type and sample.
#' @param meta Optional cell metadata (`data.frame` with barcodes as row names
#'   or a `barcode` column); required for a matrix or 10x input.
#' @param assay Assay to use. Default: the object's default assay (Seurat's
#'   `DefaultAssay()`, the `counts` assay of a `SingleCellExperiment`, the RNA
#'   matrix of a 10x folder).
#' @param minCells Pseudobulks from fewer cells are dropped.
#' @return An object of class `LocusPickPseudobulk`: a list with `counts`
#'   (features x pseudobulks) and `colData` (`sample`, `cell_type`, `n_cells`,
#'   and every metadata column that is constant within samples, such as
#'   condition labels).
#' @examples
#' sc <- locusPickSingleCellExample()
#' pb <- pseudobulk(sc$rna, sample = "sample", cellType = "cell_type", meta = sc$meta)
#' head(pb$colData)
#' @export
pseudobulk <- function(x, sample, cellType = NULL, meta = NULL, assay = NULL,
                       minCells = 10) {
    counts <- .scCounts(x, assay)
    md <- .scMeta(x, meta, cells = colnames(counts))
    .stopIfMissingCols(md, c(sample, cellType), "cell metadata")
    ct <- if (is.null(cellType)) rep("all", nrow(md)) else as.character(md[[cellType]])
    smp <- as.character(md[[sample]])
    ok <- !is.na(ct) & !is.na(smp)
    key <- factor(paste(ct, smp, sep = "\r")[ok])
    ind <- Matrix::sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1,
                                dims = c(length(key), nlevels(key)))
    pb <- as.matrix(counts[, ok, drop = FALSE] %*% ind)
    lv <- strsplit(levels(key), "\r")
    cd <- data.frame(cell_type = vapply(lv, `[`, "", 1), sample = vapply(lv, `[`, "", 2),
                     n_cells = as.integer(tabulate(as.integer(key), nlevels(key))),
                     stringsAsFactors = FALSE)
    # carry sample-level annotations (constant within each sample)
    mdo <- md[ok, , drop = FALSE]
    for (col in setdiff(names(mdo), c(sample, cellType))) {
        v <- mdo[[col]]
        if (!(is.character(v) || is.factor(v) || is.logical(v))) next
        per <- tapply(as.character(v), smp[ok], function(z) length(unique(z)) == 1)
        if (all(per)) cd[[col]] <- as.character(v)[match(cd$sample, smp[ok])]
    }
    colnames(pb) <- rownames(cd) <- paste(cd$cell_type, cd$sample, sep = "|")
    keep <- cd$n_cells >= minCells
    if (any(!keep))
        message("pseudobulk: dropped ", sum(!keep), " of ", length(keep),
                " pseudobulks with < ", minCells, " cells")
    structure(list(counts = pb[, keep, drop = FALSE], colData = cd[keep, , drop = FALSE]),
              class = "LocusPickPseudobulk")
}

#' @export
print.LocusPickPseudobulk <- function(x, ...) {
    cat("LocusPickPseudobulk:", nrow(x$counts), "features x", ncol(x$counts),
        "pseudobulks\n")
    print(table(cell_type = x$colData$cell_type))
    invisible(x)
}

#' Pseudobulk differential expression per cell type
#'
#' Aggregates cells into pseudobulk samples (see [pseudobulk()]) and tests a
#' contrast within each cell type with DESeq2 or edgeR. Each result is a
#' standard LocusPick DE table, ready for [defineGroups()].
#'
#' At least `minSamples` biological samples per condition are required in a
#' cell type; cell types without enough replicates are skipped with a
#' message. (Cell-level tests such as Seurat's `FindMarkers()` treat cells as
#' replicates and give inflated p-values; [readDE()] accepts their output but
#' warns.)
#'
#' @param x A Seurat object, `SingleCellExperiment`, [read10x()] result, count
#'   matrix, or a [pseudobulk()] result.
#' @param contrast Character vector `c(column, test, reference)`, e.g.
#'   `c("genotype", "KO", "WT")`. Fold changes are test / reference.
#' @param sample,cellType,meta,assay,minCells See [pseudobulk()]. The default
#'   assay is used unless `assay` is given.
#' @param method `"DESeq2"` or `"edgeR"` (quasi-likelihood).
#' @param covariates Optional sample-level columns to adjust for (e.g. the
#'   other contrast in a two-factor design).
#' @param minSamples Minimum samples per condition within a cell type.
#' @param cellTypes Optional subset of cell types to test.
#' @return A named list (one DE table per cell type) of class
#'   `LocusPickDEList`, with the pseudobulk object as attribute `pseudobulk`.
#' @examples
#' sc <- locusPickSingleCellExample()
#' de <- pseudobulkDE(sc$rna, contrast = c("genotype", "KO", "WT"), sample = "sample",
#'                    cellType = "cell_type", meta = sc$meta, covariates = "age")
#' head(de$Hepatocyte)
#' @export
pseudobulkDE <- function(x, contrast, sample = "sample", cellType = NULL, meta = NULL,
                         assay = NULL, method = c("DESeq2", "edgeR"), covariates = NULL,
                         minCells = 10, minSamples = 2, cellTypes = NULL) {
    method <- match.arg(method)
    .requirePkg(method, "for pseudobulk differential expression")
    if (!is.character(contrast) || length(contrast) != 3)
        stop("contrast must be c(column, test, reference), e.g. c(\"genotype\", \"KO\", \"WT\").",
             call. = FALSE)
    pb <- if (inherits(x, "LocusPickPseudobulk")) x else
        pseudobulk(x, sample = sample, cellType = cellType, meta = meta, assay = assay,
                   minCells = minCells)
    col <- contrast[1]; test <- contrast[2]; ref <- contrast[3]
    .stopIfMissingCols(pb$colData, c(col, covariates), "pseudobulk sample annotation")
    cts <- unique(pb$colData$cell_type)
    if (!is.null(cellTypes)) cts <- intersect(cts, cellTypes)

    out <- list()
    for (ct in cts) {
        sel <- pb$colData$cell_type == ct & pb$colData[[col]] %in% c(test, ref)
        cd <- pb$colData[sel, , drop = FALSE]
        cnt <- round(pb$counts[, sel, drop = FALSE])
        nT <- sum(cd[[col]] == test); nR <- sum(cd[[col]] == ref)
        if (nT < minSamples || nR < minSamples) {
            message(sprintf("%s: skipped (%d %s and %d %s samples; need >= %d each)",
                            ct, nT, test, nR, ref, minSamples))
            next
        }
        for (v in setdiff(names(cd), c(col, "n_cells")))
            if (is.character(cd[[v]])) cd[[v]] <- factor(cd[[v]])
        cd[[col]] <- factor(cd[[col]], levels = c(ref, test))
        cov <- covariates[vapply(covariates, function(v)
            length(unique(cd[[v]])) > 1, logical(1))]
        cov <- setdiff(cov, col)
        form <- stats::as.formula(paste("~", paste(c(cov, col), collapse = " + ")))
        cnt <- cnt[rowSums(cnt) >= 10, , drop = FALSE]
        de <- if (method == "DESeq2") {
            dds <- suppressMessages(.requireCall("DESeq2", "DESeqDataSetFromMatrix")(cnt, cd, form))
            dds <- suppressMessages(.requireCall("DESeq2", "DESeq")(dds, quiet = TRUE))
            res <- .requireCall("DESeq2", "results")(dds, contrast = c(col, test, ref))
            suppressWarnings(readDE(res, verbose = FALSE))
        } else {
            y <- .requireCall("edgeR", "DGEList")(cnt)
            design <- stats::model.matrix(form, cd)
            keep <- .requireCall("edgeR", "filterByExpr")(y, design = design)
            y <- .requireCall("edgeR", "normLibSizes")(y[keep, , keep.lib.sizes = FALSE])
            y <- .requireCall("edgeR", "estimateDisp")(y, design)
            fit <- .requireCall("edgeR", "glmQLFit")(y, design)
            qlf <- .requireCall("edgeR", "glmQLFTest")(fit, coef = paste0(col, test))
            suppressWarnings(readDE(qlf, verbose = FALSE))
        }
        attr(de, "contrast") <- paste(col, test, "vs", ref)
        attr(de, "source") <- paste("pseudobulk", method)
        message(sprintf("%s: %d vs %d samples, %s; %d of %d genes padj < 0.05",
                        ct, nT, nR, deparse(form), sum(de$padj < 0.05, na.rm = TRUE),
                        nrow(de)))
        out[[ct]] <- de
    }
    if (!length(out)) warning("No cell type had enough samples for DE.", call. = FALSE)
    structure(out, class = c("LocusPickDEList", "list"), pseudobulk = pb,
              contrast = contrast)
}

#' @export
print.LocusPickDEList <- function(x, ...) {
    cat("Pseudobulk DE (", paste(attr(x, "contrast")[2:3], collapse = " vs "),
        ") for ", length(x), " cell type(s):\n", sep = "")
    for (ct in names(x))
        cat("  ", ct, ": ", nrow(x[[ct]]), " genes, ",
            sum(x[[ct]]$padj < 0.05, na.rm = TRUE), " with padj < 0.05\n", sep = "")
    invisible(x)
}
