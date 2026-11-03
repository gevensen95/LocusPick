#' Signal sources
#'
#' A signal source tells LocusPick how to fetch a genome-wide signal (coverage,
#' mappability, methylation, ...) for a set of regions. All sources return a
#' `GRanges` of intervals with a numeric `score` column.
#'
#' * `bigWigSource()`: a bigWig file (read lazily with `rtracklayer`).
#' * `granges()`-style data: `rangesSource()` wraps an in-memory `GRanges`
#'   with a `score` column (for example a bedGraph you imported).
#' * `functionSource()`: any function `f(chr, start, end)` returning a
#'   `data.frame` with `start`, `end` (or `pos`) and `score` (or `value`).
#' * `chromDictSource()`: a chromDict object. Requires the chromDicts package,
#'   which is not on CRAN or Bioconductor and must be installed separately;
#'   `functionSource()` covers any other storage format.
#'
#' @param path Path or URL of a bigWig file.
#' @param x A `GRanges` with a `score` column.
#' @param f A function of `(chr, start, end)`.
#' @param chromDict A chromDict object.
#' @param binSize Bin width (bp) represented by each value of the chromDict.
#' @return An object of class `locusSignal`.
#' @name signalSources
#' @examples
#' bw <- system.file("extdata", "toy_A.bw", package = "LocusPick")
#' src <- bigWigSource(bw)
#' fetchSignal(src, GenomicRanges::GRanges("chr1:1000-1100"))
NULL

.newSource <- function(fetch, label) {
    structure(list(fetch = fetch, label = label), class = "locusSignal")
}

#' @rdname signalSources
#' @export
bigWigSource <- function(path) {
    if (!grepl("^(https?|ftp)://", path) && !file.exists(path))
        stop("bigWig not found: ", path, call. = FALSE)
    bwf <- rtracklayer::BigWigFile(path)
    lv <- tryCatch(GenomeInfoDb::seqlevels(bwf),
                   error = function(e) stop("Not a valid bigWig: ", path,
                                            " (", conditionMessage(e), ")",
                                            call. = FALSE))
    .newSource(function(which) {
        which <- .matchSeqStyle(which, GenomicRanges::GRanges(lv, IRanges::IRanges(1, 1)))
        which <- which[as.character(GenomicRanges::seqnames(which)) %in% lv]
        if (!length(which)) return(.emptySignal())
        GenomeInfoDb::seqlevels(which, pruning.mode = "coarse") <-
            GenomeInfoDb::seqlevelsInUse(which)
        rtracklayer::import(bwf, which = GenomicRanges::reduce(which), as = "GRanges")
    }, paste("bigWig:", basename(path)))
}

#' @rdname signalSources
#' @export
rangesSource <- function(x) {
    if (!is(x, "GRanges") || !"score" %in% colnames(S4Vectors::mcols(x)))
        stop("x must be a GRanges with a 'score' column.", call. = FALSE)
    .newSource(function(which) {
        which <- .matchSeqStyle(which, x)
        h <- GenomicRanges::findOverlaps(x, which, ignore.strand = TRUE)
        unique(x[S4Vectors::queryHits(h)])
    }, "GRanges")
}

#' @rdname signalSources
#' @export
functionSource <- function(f) {
    stopifnot(is.function(f))
    .newSource(function(which) {
        w <- GenomicRanges::reduce(which)
        res <- lapply(seq_along(w), function(i) {
            d <- f(as.character(GenomicRanges::seqnames(w)[i]),
                   GenomicRanges::start(w)[i], GenomicRanges::end(w)[i])
            .dfToSignal(d, as.character(GenomicRanges::seqnames(w)[i]), binSize = 1)
        })
        res <- res[lengths(res) > 0]
        if (!length(res)) return(.emptySignal())
        do.call(c, res)
    }, "function")
}

#' @rdname signalSources
#' @export
chromDictSource <- function(chromDict, binSize = 20) {
    pkg <- "chromDicts"  # optional, not on CRAN/Bioconductor: not in Suggests
    .requirePkg(pkg, "to read chromDict objects")
    mgd <- getExportedValue(pkg, "MakeGeneData")
    .newSource(function(which) {
        w <- GenomicRanges::reduce(which)
        res <- lapply(seq_along(w), function(i) {
            chr <- as.character(GenomicRanges::seqnames(w)[i])
            d <- tryCatch(mgd(list(x = chromDict), chr, GenomicRanges::start(w)[i],
                              GenomicRanges::end(w)[i], force = TRUE),
                          error = function(e) NULL)
            .dfToSignal(d, chr, binSize = binSize)
        })
        res <- res[lengths(res) > 0]
        if (!length(res)) return(.emptySignal())
        do.call(c, res)
    }, "chromDict")
}

.emptySignal <- function() {
    gr <- GenomicRanges::GRanges()
    gr$score <- numeric()
    gr
}

# Convert a data.frame from a function/chromDict source into signal GRanges.
.dfToSignal <- function(d, chr, binSize) {
    if (is.null(d) || !NROW(d)) return(.emptySignal())
    if (is(d, "GRanges")) return(d)
    sc <- if ("score" %in% names(d)) d$score else d$value
    if (is.null(sc)) stop("Signal data.frame needs a 'score' or 'value' column.",
                          call. = FALSE)
    if (all(c("start", "end") %in% names(d))) {
        st <- d$start; en <- d$end
    } else if ("pos" %in% names(d)) {
        st <- d$pos; en <- d$pos + binSize - 1
    } else {
        stop("Signal data.frame needs 'start'/'end' or 'pos'.", call. = FALSE)
    }
    keep <- !is.na(sc)
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(st[keep], en[keep]))
    gr$score <- as.numeric(sc[keep])
    gr
}

#' @export
print.locusSignal <- function(x, ...) {
    cat("locusSignal source (", x$label, ")\n", sep = "")
    invisible(x)
}

#' Fetch signal for regions
#'
#' @param source A `locusSignal` source (see [signalSources]).
#' @param which A `GRanges` of regions.
#' @return A `GRanges` of signal intervals with a `score` column.
#' @examples
#' src <- bigWigSource(system.file("extdata", "toy_A.bw", package = "LocusPick"))
#' fetchSignal(src, GenomicRanges::GRanges("chr1:1000-1100"))
#' @export
fetchSignal <- function(source, which) {
    if (!inherits(source, "locusSignal")) source <- .asSource(source)
    source$fetch(which)
}

.asSource <- function(x) {
    if (inherits(x, "locusSignal")) return(x)
    if (is.character(x) && length(x) == 1) return(bigWigSource(x))
    if (is(x, "GRanges")) return(rangesSource(x))
    if (is.function(x)) return(functionSource(x))
    stop("Cannot use object of class '", class(x)[1], "' as a signal source.",
         call. = FALSE)
}

#' Summarize signal over peaks
#'
#' For each peak and each source, computes the base-pair-weighted mean signal
#' and the fraction of the peak covered by signal intervals.
#'
#' @param peaks Peaks (anything accepted by [asPeaks()]).
#' @param sources A named list of signal sources, bigWig paths, or `GRanges`.
#' @param missingAsZero Treat bases without a signal interval as 0 rather than
#'   missing. Use `TRUE` for tracks that omit zeros (such as Umap mappability)
#'   and `FALSE` (default) when gaps mean "no data".
#' @return A `data.frame` with `peak_id` plus, for each source `s`, columns
#'   `s_mean` and `s_cov`.
#' @examples
#' pk <- asPeaks(system.file("extdata", "toy_peaks.bed", package = "LocusPick"))
#' bw <- system.file("extdata", c("toy_A.bw", "toy_B.bw"), package = "LocusPick")
#' signalAtPeaks(pk, list(A = bw[1], B = bw[2]))
#' @export
signalAtPeaks <- function(peaks, sources, missingAsZero = FALSE) {
    pk <- asPeaks(peaks)
    if (!is.list(sources) || is.null(names(sources)) || any(names(sources) == ""))
        stop("sources must be a named list.", call. = FALSE)
    out <- data.frame(peak_id = pk$peak_id, stringsAsFactors = FALSE)
    for (nm in names(sources)) {
        sig <- fetchSignal(.asSource(sources[[nm]]), pk)
        s <- .summarizeOver(pk, sig, missingAsZero)
        out[[paste0(nm, "_mean")]] <- s$mean
        out[[paste0(nm, "_cov")]] <- s$cov
    }
    out
}

.summarizeOver <- function(regions, sig, missingAsZero) {
    n <- length(regions)
    w <- GenomicRanges::width(regions)
    num <- numeric(n); covd <- numeric(n)
    if (length(sig)) {
        sig <- .matchSeqStyle(sig, regions)
        h <- GenomicRanges::findOverlaps(regions, sig, ignore.strand = TRUE)
        if (length(h)) {
            q <- S4Vectors::queryHits(h); s <- S4Vectors::subjectHits(h)
            ov <- GenomicRanges::width(GenomicRanges::pintersect(
                regions[q], sig[s], ignore.strand = TRUE))
            num  <- as.numeric(tapply(ov * sig$score[s], factor(q, seq_len(n)), sum))
            covd <- as.numeric(tapply(ov, factor(q, seq_len(n)), sum))
            num[is.na(num)] <- 0; covd[is.na(covd)] <- 0
        }
    }
    covd <- pmin(covd, w)  # overlapping signal intervals should not exceed width
    mean <- if (missingAsZero) num / w else ifelse(covd > 0, num / covd, NA_real_)
    list(mean = mean, cov = covd / w)
}

#' Mean mappability over peaks
#'
#' Convenience wrapper around [signalAtPeaks()] for a mappability track such as
#' Umap multi-read mappability (Karimzadeh et al. 2018). Bases absent from the
#' track are treated as 0 (no uniquely mappable k-mer).
#'
#' @param peaks Peaks (anything accepted by [asPeaks()]).
#' @param mappability A signal source or bigWig path.
#' @return Numeric vector of mean mappability (0-1), one value per peak.
#' @seealso `system.file("scripts", "make_umap_bw.sh", package = "LocusPick")`
#'   builds a Umap bigWig from the bedGraph distributed by the Umap authors.
#' @examples
#' pk <- asPeaks(system.file("extdata", "toy_peaks.bed", package = "LocusPick"))
#' mappabilityAtPeaks(pk, system.file("extdata", "toy_map.bw", package = "LocusPick"))
#' @export
mappabilityAtPeaks <- function(peaks, mappability) {
    signalAtPeaks(peaks, list(map = mappability), missingAsZero = TRUE)$map_mean
}
