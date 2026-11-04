# Single-cell input helpers: Seurat (+ Signac), SingleCellExperiment, and
# Cell Ranger / Cell Ranger ARC output folders.

#' Read a Cell Ranger or Cell Ranger ARC output folder
#'
#' Reads the filtered feature-barcode matrix (MTX folder, or `.h5` when
#' DropletUtils is installed), the ATAC fragments file, the peak BED file and
#' per-barcode metrics from a Cell Ranger (`count`), Cell Ranger ATAC or
#' Cell Ranger ARC (multiome) `outs/` folder.
#'
#' Cell Ranger does not annotate cell types or samples, so downstream
#' functions also need a `meta` table with one row per barcode (row names or a
#' `barcode` column) holding your cell type, sample and condition labels.
#'
#' @param dir The `outs/` folder (or a folder containing it).
#' @return A list of class `LocusPick10x` with `rna` (genes x cells counts,
#'   or `NULL`), `atac` (peaks x cells counts, or `NULL`), `peaks` (`GRanges`),
#'   `fragments` (path, or `NULL`) and `metrics` (`data.frame`, or `NULL`).
#' @examples
#' # see vignette("single-cell", package = "LocusPick")
#' @export
read10x <- function(dir) {
    if (!dir.exists(dir)) stop("Folder not found: ", dir, call. = FALSE)
    f <- list.files(dir, recursive = TRUE, full.names = TRUE)
    pick <- function(pattern) {
        h <- f[grepl(pattern, f)]
        if (length(h)) h[order(nchar(h))][1] else NULL
    }
    mtx <- pick("filtered_(feature|peak)_bc_matrix/matrix\\.mtx(\\.gz)?$")
    h5 <- pick("filtered_(feature|peak)_bc_matrix\\.h5$")
    frag <- pick("(atac_)?fragments\\.tsv\\.gz$")
    bed <- pick("(atac_peaks|/peaks)\\.bed$")
    met <- pick("(per_barcode_metrics|singlecell)\\.csv$")

    if (!is.null(mtx)) {
        d <- dirname(mtx)
        m <- methods::as(Matrix::readMM(mtx), "CsparseMatrix")
        feat <- utils::read.delim(list.files(d, "^(features|genes|peaks)\\.tsv", full.names = TRUE)[1],
                                  header = FALSE, stringsAsFactors = FALSE)
        bc <- readLines(list.files(d, "^barcodes\\.tsv", full.names = TRUE)[1])
        type <- if (ncol(feat) >= 3) feat[[3]] else rep("Peaks", nrow(feat))
        name <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
        colnames(m) <- bc
        src <- basename(d)
    } else if (!is.null(h5)) {
        .requirePkg("DropletUtils", "to read Cell Ranger .h5 files")
        sce <- DropletUtils::read10xCounts(h5)
        m <- methods::as(.requireCall("SummarizedExperiment", "assay")(sce, "counts"),
                         "CsparseMatrix")
        rd <- as.data.frame(.requireCall("SummarizedExperiment", "rowData")(sce))
        name <- rd$Symbol; type <- if ("Type" %in% names(rd)) rd$Type else "Gene Expression"
        colnames(m) <- as.character(.requireCall("SummarizedExperiment", "colData")(sce)$Barcode)
        feat <- data.frame(V1 = rd$ID); src <- basename(h5)
    } else {
        stop("No filtered feature-barcode matrix found in ", dir, call. = FALSE)
    }
    isRNA <- type == "Gene Expression"
    rna <- atac <- NULL
    if (any(isRNA)) {
        rna <- m[isRNA, , drop = FALSE]
        rownames(rna) <- make.unique(name[isRNA])
    }
    if (any(type == "Peaks")) {
        atac <- m[type == "Peaks", , drop = FALSE]
        rownames(atac) <- as.character(feat[[1]][type == "Peaks"])
    }
    peaks <- if (!is.null(bed)) readPeaks(bed, format = "bed", verbose = FALSE) else
        if (!is.null(atac)) .stringToRanges(rownames(atac)) else NULL
    metrics <- if (!is.null(met)) utils::read.csv(met, stringsAsFactors = FALSE) else NULL
    if (!is.null(frag) && !file.exists(paste0(frag, ".tbi")))
        warning("Fragments index (.tbi) not found for ", frag, call. = FALSE)
    message("read10x: ", src, " | ", ncol(m), " cells | RNA: ",
            if (is.null(rna)) "none" else paste(nrow(rna), "genes"), " | ATAC: ",
            if (is.null(atac)) "none" else paste(nrow(atac), "peaks"),
            " | fragments: ", if (is.null(frag)) "none" else basename(frag))
    structure(list(rna = rna, atac = atac, peaks = peaks, fragments = frag,
                   metrics = metrics), class = "LocusPick10x")
}

#' @export
print.LocusPick10x <- function(x, ...) {
    cat("LocusPick10x:", if (!is.null(x$rna)) ncol(x$rna) else ncol(x$atac), "cells\n")
    cat("  RNA: ", if (is.null(x$rna)) "none" else paste(nrow(x$rna), "genes"), "\n", sep = "")
    cat("  ATAC: ", if (is.null(x$atac)) "none" else paste(nrow(x$atac), "peaks"), "\n", sep = "")
    cat("  fragments: ", if (is.null(x$fragments)) "none" else x$fragments, "\n", sep = "")
    invisible(x)
}

.stringToRanges <- function(s) {
    p <- do.call(rbind, strsplit(sub("^(.*)[:-]([0-9]+)-([0-9]+)$", "\\1\t\\2\t\\3", s), "\t"))
    gr <- GenomicRanges::GRanges(p[, 1], IRanges::IRanges(as.integer(p[, 2]), as.integer(p[, 3])))
    gr$peak_id <- .peakId(gr)
    gr
}

.isSeurat <- function(x) inherits(x, "Seurat")
.isSCE <- function(x) is(x, "SummarizedExperiment")

.isChromatinAssay <- function(x, assay) {
    .isSeurat(x) && inherits(x[[assay]], "ChromatinAssay")
}

# Default assay: Seurat's DefaultAssay, SCE "counts", 10x "rna".
.scAssay <- function(x, assay = NULL, kind = c("rna", "atac")) {
    kind <- match.arg(kind)
    if (!is.null(assay)) return(assay)
    if (.isSeurat(x)) {
        def <- .requireCall("SeuratObject", "DefaultAssay")(x)
        if (kind == "atac" && !.isChromatinAssay(x, def)) {
            ch <- names(x@assays)[vapply(names(x@assays), function(a)
                .isChromatinAssay(x, a), logical(1))]
            if (length(ch)) return(ch[1])
        }
        return(def)
    }
    if (inherits(x, "LocusPick10x")) return(kind)
    if (.isSCE(x)) {
        an <- .requireCall("SummarizedExperiment", "assayNames")(x)
        return(if ("counts" %in% an) "counts" else an[1])
    }
    NULL
}

.scCounts <- function(x, assay = NULL, kind = "rna") {
    assay <- .scAssay(x, assay, kind)
    if (.isSeurat(x)) {
        m <- tryCatch(.requireCall("SeuratObject", "LayerData")(x, assay = assay, layer = "counts"),
                      error = function(e) .requireCall("SeuratObject", "GetAssayData")(
                          x, assay = assay, slot = "counts"))
        if (identical(assay, "SCT"))
            message("  Using SCT corrected counts; raw counts (the RNA assay) are ",
                    "usually preferred for pseudobulk DE.")
        return(m)
    }
    if (inherits(x, "LocusPick10x")) {
        m <- x[[assay]]
        if (is.null(m)) stop("This 10x folder has no '", assay, "' matrix.", call. = FALSE)
        return(m)
    }
    if (.isSCE(x)) return(.requireCall("SummarizedExperiment", "assay")(x, assay))
    if (is.matrix(x) || is(x, "Matrix")) return(x)
    stop("Unsupported single-cell object of class '", class(x)[1], "'.", call. = FALSE)
}

# Cell metadata aligned to the cells (columns) of x. `meta` (a data.frame with
# barcodes as row names or in a `barcode` column) adds or overrides columns.
.scMeta <- function(x, meta = NULL, cells = NULL) {
    md <- if (.isSeurat(x)) x@meta.data else
        if (.isSCE(x)) as.data.frame(.requireCall("SummarizedExperiment", "colData")(x)) else
        NULL
    if (is.null(cells)) cells <- if (.isSeurat(x) || .isSCE(x)) colnames(x) else
        if (inherits(x, "LocusPick10x")) colnames(if (!is.null(x$rna)) x$rna else x$atac) else
        colnames(x)
    if (is.null(md)) md <- data.frame(row.names = cells)
    if (inherits(x, "LocusPick10x") && !is.null(x$metrics) && "barcode" %in% names(x$metrics)) {
        mt <- x$metrics[match(cells, x$metrics$barcode), setdiff(names(x$metrics), "barcode"),
                        drop = FALSE]
        rownames(mt) <- cells
        md <- cbind(md[cells, , drop = FALSE], mt)
    }
    if (!is.null(meta)) {
        meta <- as.data.frame(meta)
        if ("barcode" %in% names(meta)) rownames(meta) <- meta$barcode
        hit <- match(cells, rownames(meta))
        if (all(is.na(hit))) stop("No cells of the object are found in `meta` ",
                                  "(match row names or a `barcode` column).", call. = FALSE)
        if (anyNA(hit)) message("  ", sum(is.na(hit)), " cells have no `meta` row and are ignored.")
        add <- meta[hit, , drop = FALSE]
        md <- md[cells, setdiff(names(md), names(add)), drop = FALSE]
        md <- cbind(md, add)
        rownames(md) <- cells
    }
    md[cells, , drop = FALSE]
}
