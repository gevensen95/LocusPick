#' Signal from a single-cell ATAC fragments file
#'
#' Computes pseudobulk accessibility for a set of cells directly from a
#' tabix-indexed fragments file (`fragments.tsv.gz` + `.tbi`, as written by
#' Cell Ranger ATAC/ARC and used by Signac), only over the regions requested.
#'
#' Each fragment counts once (the duplicate-count column is ignored, as in
#' Signac). With `normalize = "CPM"`, coverage is scaled per million fragments
#' of the selected cells; pass their total in `total` (e.g. the sum of a
#' per-cell fragment or in-peak count column), otherwise it is counted by
#' reading the whole file once.
#'
#' @param fragments Path(s) to fragments files. Several files (e.g. one per
#'   sample) are combined.
#' @param cells Barcodes to include, as written in the fragments file(s): a
#'   character vector, or a list with one vector per file. `NULL` uses all.
#' @param total Total fragments of the selected cells, for CPM.
#' @param normalize `"CPM"` or `"none"`.
#' @return A `locusSignal` source.
#' @examples
#' sc <- locusPickSingleCellExample()
#' hep <- rownames(sc$meta)[sc$meta$cell_type == "Hepatocyte"]
#' src <- fragmentSource(sc$fragments, cells = hep)
#' signalAtPeaks(sc$peaks[1:3], list(hepatocyte = src))
#' @export
fragmentSource <- function(fragments, cells = NULL, total = NULL,
                           normalize = c("CPM", "none")) {
    normalize <- match.arg(normalize)
    .requirePkg("Rsamtools", "to read fragments files")
    for (f in fragments) {
        if (!file.exists(f)) stop("Fragments file not found: ", f, call. = FALSE)
        if (!file.exists(paste0(f, ".tbi")))
            stop("Tabix index not found for ", f, ". Create it with ",
                 "Rsamtools::indexTabix(f, format = \"bed\") or tabix -p bed.", call. = FALSE)
    }
    cellList <- if (is.list(cells)) cells else rep(list(cells), length(fragments))
    if (length(cellList) != length(fragments))
        stop("`cells` must have one entry per fragments file.", call. = FALSE)
    if (normalize == "CPM" && is.null(total)) {
        message("Counting fragments of the selected cells in ",
                paste(basename(fragments), collapse = ", "),
                " (pass `total` to skip this step)...")
        total <- sum(mapply(.countFragments, fragments, cellList))
    }
    scale <- if (normalize == "CPM") 1e6 / total else 1
    lv <- unique(unlist(lapply(fragments, Rsamtools::seqnamesTabix)))

    .newSource(function(which) {
        w <- GenomicRanges::reduce(which, ignore.strand = TRUE)
        w <- .matchSeqStyle(w, GenomicRanges::GRanges(lv, IRanges::IRanges(1, 1)))
        w <- w[as.character(GenomicRanges::seqnames(w)) %in% lv]
        if (!length(w)) return(.emptySignal())
        fr <- do.call(c, lapply(seq_along(fragments), function(i) {
            tf <- Rsamtools::TabixFile(fragments[i])
            ok <- as.character(GenomicRanges::seqnames(w)) %in% Rsamtools::seqnamesTabix(tf)
            if (!any(ok)) return(GenomicRanges::GRanges())
            ln <- unlist(Rsamtools::scanTabix(tf, param = w[ok]), use.names = FALSE)
            .parseFragments(ln, cellList[[i]], tag = i)
        }))
        .coverageRuns(fr, w, scale, lv)
    }, paste0("fragments (", if (is.null(cells)) "all" else length(unlist(cells)), " cells)"))
}

.parseFragments <- function(lines, cells, tag) {
    if (!length(lines)) return(GenomicRanges::GRanges())
    p <- strsplit(lines, "\t", fixed = TRUE)
    bc <- vapply(p, `[`, "", 4)
    keep <- if (is.null(cells)) rep(TRUE, length(p)) else bc %in% cells
    if (!any(keep)) return(GenomicRanges::GRanges())
    p <- p[keep]
    st <- as.integer(vapply(p, `[`, "", 2)) + 1L
    en <- as.integer(vapply(p, `[`, "", 3))
    chr <- vapply(p, `[`, "", 1)
    key <- paste(tag, chr, st, en, bc[keep])
    u <- !duplicated(key)        # a fragment spanning two requested regions is returned twice
    GenomicRanges::GRanges(chr[u], IRanges::IRanges(st[u], en[u]))
}

.countFragments <- function(path, cells, chunk = 1e6) {
    con <- gzfile(path, "r"); on.exit(close(con))
    n <- 0
    repeat {
        x <- scan(con, what = list(NULL, NULL, NULL, "", NULL), sep = "\t",
                  nmax = chunk, comment.char = "#", quiet = TRUE, flush = TRUE,
                  fill = TRUE)[[4]]
        if (!length(x)) break
        n <- n + if (is.null(cells)) length(x) else sum(x %in% cells)
    }
    n
}

# Coverage of fragments as runs clipped to the requested regions, with zero
# runs where no fragment falls (the data exist; there is simply no signal).
.coverageRuns <- function(fr, w, scale, lv) {
    GenomicRanges::strand(fr) <- "*"
    GenomeInfoDb::seqlevels(fr, pruning.mode = "coarse") <-
        union(lv, GenomeInfoDb::seqlevels(fr))
    GenomeInfoDb::seqlevels(w) <- GenomeInfoDb::seqlevels(fr)
    cov <- GenomicRanges::coverage(fr) * scale
    runs <- methods::as(cov, "GRanges")
    GenomeInfoDb::seqlengths(runs) <- NA      # requested regions may extend past the last fragment
    h <- GenomicRanges::findOverlaps(runs, w)
    out <- GenomicRanges::pintersect(runs[S4Vectors::queryHits(h)], w[S4Vectors::subjectHits(h)])
    out$score <- runs$score[S4Vectors::queryHits(h)]
    S4Vectors::mcols(out) <- S4Vectors::mcols(out)["score"]
    miss <- GenomicRanges::setdiff(w, out)
    if (length(miss)) { miss$score <- 0; out <- c(out, miss) }
    sort(out)
}

#' Pseudobulk accessibility tracks per cell group
#'
#' Builds one fragment-based signal source per combination of metadata
#' columns (e.g. cell type x condition) from a Signac object, or from
#' fragments files plus cell metadata (`SingleCellExperiment`, [read10x()]
#' output, or a plain metadata table).
#'
#' For Signac objects, fragment files and the mapping between cell names and
#' file barcodes are taken from the object. CPM totals use, when available, a
#' per-cell fragment count column (`atac_fragments`, `passed_filters`,
#' `fragments`) or else in-peak counts (`nCount_<assay>`), which is how Signac
#' normalizes coverage plots.
#'
#' @param x A Seurat/Signac object, `SingleCellExperiment`, [read10x()]
#'   result, or `NULL` when `meta` and `fragments` are given.
#' @param groupBy Metadata column(s) defining groups, e.g.
#'   `c("cell_type", "age")`.
#' @param fragments Fragments file(s); taken from the object if `NULL`.
#' @param meta Optional cell metadata (barcodes as row names or a `barcode`
#'   column).
#' @param assay Chromatin assay of a Signac object (default: the first
#'   `ChromatinAssay`).
#' @param minCells Groups with fewer cells are dropped.
#' @param sep Separator used to build track names from group levels.
#' @return A named list of `locusSignal` sources (e.g. `Hepatocyte_young`),
#'   with attribute `groups` describing each track.
#' @examples
#' sc <- locusPickSingleCellExample()
#' tr <- scSignal(NULL, c("cell_type", "age"), fragments = sc$fragments, meta = sc$meta)
#' names(tr)
#' @export
scSignal <- function(x, groupBy, fragments = NULL, meta = NULL, assay = NULL,
                     minCells = 10, sep = "_") {
    cellMap <- NULL
    if (!is.null(x) && .isSeurat(x)) {
        assay <- .scAssay(x, assay, "atac")
        if (!.isChromatinAssay(x, assay))
            stop("Assay '", assay, "' is not a Signac ChromatinAssay.", call. = FALSE)
        if (is.null(fragments)) {
            frs <- .requireCall("Signac", "Fragments")(x[[assay]])
            if (!length(frs)) stop("The ChromatinAssay has no fragment files.", call. = FALSE)
            gfd <- .requireCall("Signac", "GetFragmentData")
            fragments <- vapply(frs, gfd, "", slot = "path")
            cellMap <- lapply(frs, gfd, slot = "cells")   # names: cell names; values: barcodes
        }
    } else if (!is.null(x) && inherits(x, "LocusPick10x") && is.null(fragments)) {
        fragments <- x$fragments
    }
    if (is.null(fragments)) stop("No fragments file given or found in the object.", call. = FALSE)
    md <- if (is.null(x)) {
        if (is.null(meta)) stop("Give a single-cell object or `meta`.", call. = FALSE)
        m <- as.data.frame(meta)
        if ("barcode" %in% names(m)) rownames(m) <- m$barcode
        m
    } else .scMeta(x, meta)
    .stopIfMissingCols(md, groupBy, "cell metadata")
    if (is.null(cellMap)) {
        bcs <- if ("barcode" %in% names(md)) as.character(md$barcode) else rownames(md)
        cellMap <- rep(list(stats::setNames(bcs, rownames(md))), length(fragments))
    }

    totCol <- intersect(c("atac_fragments", "passed_filters", "fragments",
                          if (!is.null(assay)) paste0("nCount_", assay), "nCount_ATAC",
                          "nCount_peaks"), names(md))[1]
    if (!is.na(totCol)) message("scSignal: CPM using per-cell '", totCol, "'")

    lev <- lapply(groupBy, function(g) {
        v <- md[[g]]; if (is.factor(v)) levels(v) else unique(as.character(v[!is.na(v)]))
    })
    combos <- expand.grid(lev, stringsAsFactors = FALSE)
    names(combos) <- groupBy
    key <- do.call(paste, c(lapply(groupBy, function(g) as.character(md[[g]])), sep = sep))
    out <- list(); info <- list()
    for (i in seq_len(nrow(combos))) {
        nm <- do.call(paste, c(as.list(combos[i, , drop = FALSE]), sep = sep))
        cells <- rownames(md)[key == nm & !is.na(key)]
        if (length(cells) < minCells) next
        perFile <- lapply(cellMap, function(cm) unname(cm[intersect(cells, names(cm))]))
        tot <- if (!is.na(totCol)) sum(as.numeric(md[cells, totCol]), na.rm = TRUE) else NULL
        out[[nm]] <- suppressMessages(fragmentSource(fragments, perFile, total = tot))
        info[[nm]] <- data.frame(track = nm, combos[i, , drop = FALSE], n_cells = length(cells),
                                 row.names = NULL, stringsAsFactors = FALSE)
    }
    if (!length(out)) stop("No group has >= ", minCells, " cells.", call. = FALSE)
    message(paste(strwrap(paste0("scSignal: ", length(out), " tracks (",
                                 paste(sprintf("%s: %d cells", names(out),
                                               vapply(info, function(d) d$n_cells, 1L)),
                                       collapse = "; "), ")"),
                          width = 78, exdent = 2), collapse = "\n"))
    structure(out, groups = do.call(rbind, info))
}
