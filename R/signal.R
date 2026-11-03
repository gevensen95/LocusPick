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
    if (is.character(x) && length(x) > 1) return(.fileSource(x))
    if (is.list(x) && !is(x, "GRanges") && length(x) > 1 &&
        all(vapply(x, function(e) inherits(e, "locusSignal") || is.character(e), logical(1))))
        return(meanSource(lapply(x, .asSource)))
    if (is.character(x) && length(x) == 1)
        return(if (grepl("\\.bam$", x, ignore.case = TRUE)) bamSource(x) else bigWigSource(x))
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

#' Average several signal sources
#'
#' Combines replicate tracks into one condition track by averaging their
#' scores base by base. Bases missing from a replicate count as 0 for that
#' replicate; bases missing from all replicates stay missing.
#'
#' A character vector of several files, or a list of sources, passed as one
#' element of a `signal` list is averaged automatically, so
#' `signal = list(young = c("y1.bw", "y2.bw"))` works directly.
#'
#' @param sources A list of signal sources or file paths.
#' @return A `locusSignal` source.
#' @examples
#' bw <- system.file("extdata", c("toy_A.bw", "toy_B.bw"), package = "LocusPick")
#' m <- meanSource(as.list(bw))
#' fetchSignal(m, GenomicRanges::GRanges("chr1:5000-5040"))
#' @export
meanSource <- function(sources) {
    srcs <- lapply(sources, .asSource)
    n <- length(srcs)
    .newSource(function(which) {
        w <- GenomicRanges::reduce(which, ignore.strand = TRUE)
        sigs <- lapply(srcs, function(s) {
            g <- s$fetch(w)
            g <- .matchSeqStyle(g, w)
            GenomicRanges::strand(g) <- "*"
            g
        })
        sigs <- sigs[lengths(sigs) > 0]
        if (!length(sigs)) return(.emptySignal())
        d <- GenomicRanges::disjoin(do.call(c, lapply(sigs, GenomicRanges::granges)))
        tot <- numeric(length(d))
        for (s in sigs) {
            h <- GenomicRanges::findOverlaps(d, s)
            add <- tapply(s$score[S4Vectors::subjectHits(h)], S4Vectors::queryHits(h), mean)
            idx <- as.integer(names(add))
            tot[idx] <- tot[idx] + add
        }
        d$score <- tot / n
        d
    }, paste0("mean of ", n, " tracks"))
}

#' Signal computed from a BAM file
#'
#' Computes CPM-normalized fragment coverage directly from an indexed BAM, only
#' over the regions requested. Useful for exploring data before bigWigs exist.
#' For figures, bigWigs made with deepTools (`bamCoverage`) give more control
#' (blacklists, extension, normalization) and are much faster.
#'
#' @param path Indexed BAM file (`.bai` alongside).
#' @param paired `TRUE`/`FALSE`, or `NULL` to detect from the file.
#' @param fragmentLength Single-end only: extend reads to this length.
#' @param minMapq Minimum mapping quality.
#' @param removeDuplicates Drop reads flagged as duplicates.
#' @param normalize `"CPM"` (per million mapped fragments) or `"none"`.
#' @return A `locusSignal` source.
#' @examples
#' bam <- system.file("extdata", "toy.bam", package = "LocusPick")
#' if (requireNamespace("Rsamtools", quietly = TRUE) &&
#'     requireNamespace("GenomicAlignments", quietly = TRUE)) {
#'     src <- bamSource(bam)
#'     fetchSignal(src, GenomicRanges::GRanges("chr1:4900-5000"))
#' }
#' @export
bamSource <- function(path, paired = NULL, fragmentLength = NULL, minMapq = 10,
                      removeDuplicates = TRUE, normalize = c("CPM", "none")) {
    normalize <- match.arg(normalize)
    .requirePkg("Rsamtools", "to read BAM files")
    .requirePkg("GenomicAlignments", "to read BAM files")
    if (!file.exists(path)) stop("BAM not found: ", path, call. = FALSE)
    if (!file.exists(paste0(path, ".bai")) && !file.exists(sub("\\.bam$", ".bai", path)))
        stop("BAM index (.bai) not found for ", path,
             ". Create it with Rsamtools::indexBam() or samtools index.", call. = FALSE)
    if (is.null(paired)) paired <- Rsamtools::testPairedEndBam(path)
    mapped <- sum(Rsamtools::idxstatsBam(path)$mapped)
    nFrag <- if (paired) mapped / 2 else mapped
    scale <- if (normalize == "CPM") 1e6 / nFrag else 1
    message("BAM signal for ", basename(path), ": ", if (paired) "paired" else "single",
            "-end, ", normalize, ", MAPQ >= ", minMapq,
            if (removeDuplicates) ", duplicates removed" else "",
            ". For figures, bigWigs from deepTools bamCoverage are recommended.")
    lv <- GenomeInfoDb::seqlevels(Rsamtools::BamFile(path))

    .newSource(function(which) {
        w <- GenomicRanges::reduce(which, ignore.strand = TRUE)
        w <- .matchSeqStyle(w, GenomicRanges::GRanges(lv, IRanges::IRanges(1, 1)))
        w <- w[as.character(GenomicRanges::seqnames(w)) %in% lv]
        if (!length(w)) return(.emptySignal())
        flag <- Rsamtools::scanBamFlag(
            isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
            isNotPassingQualityControls = FALSE,
            isDuplicate = if (removeDuplicates) FALSE else NA,
            isProperPair = if (paired) TRUE else NA)
        prm <- Rsamtools::ScanBamParam(which = w, flag = flag, mapqFilter = minMapq)
        fr <- if (paired) {
            GenomicRanges::granges(GenomicAlignments::readGAlignmentPairs(
                path, param = prm, use.names = TRUE), on.discordant.seqnames = "drop")
        } else {
            g <- GenomicRanges::granges(GenomicAlignments::readGAlignments(
                path, param = prm, use.names = TRUE))
            if (!is.null(fragmentLength)) g <- GenomicRanges::resize(g, fragmentLength, fix = "start")
            g
        }
        # a read spanning two requested regions is returned once per region
        fr <- fr[!duplicated(paste(names(fr), GenomicRanges::start(fr), GenomicRanges::end(fr)))]
        GenomicRanges::strand(fr) <- "*"
        GenomeInfoDb::seqlevels(fr, pruning.mode = "coarse") <- lv
        GenomeInfoDb::seqlevels(w) <- lv
        cov <- GenomicRanges::coverage(fr) * scale
        runs <- methods::as(cov, "GRanges")
        h <- GenomicRanges::findOverlaps(runs, w)
        out <- GenomicRanges::pintersect(runs[S4Vectors::queryHits(h)], w[S4Vectors::subjectHits(h)])
        out$score <- runs$score[S4Vectors::queryHits(h)]
        S4Vectors::mcols(out) <- S4Vectors::mcols(out)["score"]
        # regions beyond the last read on a chromosome: zero coverage
        miss <- GenomicRanges::setdiff(w, out)
        if (length(miss)) { miss$score <- 0; out <- c(out, miss) }
        sort(out)
    }, paste("BAM:", basename(path)))
}
