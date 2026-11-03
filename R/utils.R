#' @importFrom methods is
#' @importFrom stats median setNames
#' @importFrom utils write.table
#' @importFrom ggplot2 .data
NULL

# ---- small internal helpers ---------------------------------------------------

.stopIfMissingCols <- function(df, cols, what = "table") {
    miss <- setdiff(cols, colnames(df))
    if (length(miss))
        stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
             call. = FALSE)
    invisible(TRUE)
}

.requirePkg <- function(pkg, why) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Package '", pkg, "' is required ", why,
             ". Install it with BiocManager::install('", pkg, "').",
             call. = FALSE)
}

.peakId <- function(gr) {
    paste0(as.character(GenomicRanges::seqnames(gr)), ":",
           GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
}

# percentile rank in (0, 1]; NA stays NA
.pctRank <- function(x) {
    r <- rank(x, na.last = "keep", ties.method = "average")
    r / max(r, na.rm = TRUE)
}

# Make seqlevels of `x` follow the style of `ref` (e.g. "1" -> "chr1").
# Silently returns `x` unchanged when the style cannot be determined.
.matchSeqStyle <- function(x, ref) {
    refStyle <- tryCatch(GenomeInfoDb::seqlevelsStyle(ref)[1],
                         error = function(e) NA)
    if (is.na(refStyle)) return(x)
    curStyle <- tryCatch(GenomeInfoDb::seqlevelsStyle(x)[1],
                         error = function(e) NA)
    if (is.na(curStyle) || identical(curStyle, refStyle)) return(x)
    tryCatch({
        GenomeInfoDb::seqlevelsStyle(x) <- refStyle
        x
    }, error = function(e) x)
}

#' Convert peaks to a GRanges
#'
#' Accepts a path to a BED/narrowPeak file, a `data.frame`, or a `GRanges`,
#' and returns a `GRanges` with a `peak_id` column (`chr:start-end`, 1-based)
#' unless one is already present.
#'
#' @param x A file path, `data.frame` or `GRanges`.
#' @param zeroBased For `data.frame` input: are the start coordinates 0-based
#'   (as in a BED file read with `read.delim`)? Files read by path through
#'   `rtracklayer` are converted automatically.
#' @return A `GRanges`.
#' @examples
#' asPeaks(data.frame(chr = "chr1", start = c(100, 500), end = c(200, 650)))
#' @export
asPeaks <- function(x, zeroBased = FALSE) {
    if (is.character(x) && length(x) == 1) {
        if (!file.exists(x)) stop("Peak file not found: ", x, call. = FALSE)
        gr <- tryCatch(rtracklayer::import(x),
                       error = function(e) .readBedLike(x))
    } else if (is.data.frame(x)) {
        gr <- GenomicRanges::makeGRangesFromDataFrame(
            x, keep.extra.columns = TRUE, starts.in.df.are.0based = zeroBased)
    } else if (is(x, "GRanges")) {
        gr <- x
    } else {
        stop("Peaks must be a file path, data.frame or GRanges.", call. = FALSE)
    }
    if (!"peak_id" %in% colnames(S4Vectors::mcols(gr)))
        gr$peak_id <- .peakId(gr)
    gr
}

# Fallback for BED-like files that rtracklayer refuses (extra columns, headers)
.readBedLike <- function(path) {
    df <- utils::read.delim(path, header = FALSE, comment.char = "#",
                            stringsAsFactors = FALSE)
    if (!is.numeric(df[[2]])) df <- df[-1, , drop = FALSE]
    names(df)[1:3] <- c("chr", "start", "end")
    df$start <- as.integer(df$start); df$end <- as.integer(df$end)
    GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE,
                                            starts.in.df.are.0based = TRUE)
}
