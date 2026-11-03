#' Read peaks from common peak callers
#'
#' Reads the direct output of peak callers and differential binding tools,
#' detecting the format automatically, and optionally merges replicates into
#' consensus peaks and removes blacklisted regions.
#'
#' Recognized formats:
#'
#' | Format | Detected by | Extra columns kept |
#' |---|---|---|
#' | MACS2/MACS3, Genrich `narrowPeak`, `broadPeak`, `gappedPeak` | extension | `signal`, `pvalue`, `qvalue`, `summit` |
#' | MACS2 `_peaks.xls` | extension / header | `fold_enrichment`, `pvalue`, `qvalue`, `summit` |
#' | SEACR | 6 columns, last is `chr:start-end` | `total_signal`, `max_signal`, `max_region` |
#' | HOMER `peaks.txt` | `PeakID` header | `homer_score` |
#' | HOMER `annotatePeaks.pl` | `PeakID (cmd=...)` header | `homer_gene`, `homer_annotation`, `homer_distance_to_tss` |
#' | DiffBind `dba.report()` | `Fold` and `FDR` columns | `diff_fold`, `diff_fdr`, `diff_conc` |
#' | nf-core consensus / plain BED | 3+ columns | `name` |
#'
#' Coordinates are converted to 1-based `GRanges` (BED-like files are 0-based;
#' MACS2 `.xls` and HOMER files are 1-based).
#'
#' @param x One or more peak files, or a `GRanges`/`data.frame` (including a
#'   DiffBind report). Several files are merged into consensus peaks.
#' @param format `"auto"` or one of `"narrowPeak"`, `"broadPeak"`,
#'   `"gappedPeak"`, `"macs2_xls"`, `"seacr"`, `"homer"`, `"homer_annotated"`,
#'   `"diffbind"`, `"bed"`.
#' @param minReps With several files: minimum number of files a consensus
#'   peak must overlap. Defaults to 2 (or 1 for a single file).
#' @param blacklist Optional blacklist (path or `GRanges`); overlapping peaks
#'   are removed.
#' @param verbose Print what was detected.
#' @return A `GRanges` with `peak_id`, `peak_format`, and the format's extra
#'   columns (for consensus peaks: `n_reps`).
#' @examples
#' f <- system.file("extdata", "toy_peaks.bed", package = "LocusPick")
#' readPeaks(f)
#' @export
readPeaks <- function(x, format = "auto", minReps = NULL, blacklist = NULL,
                      verbose = TRUE) {
    multi <- (is.character(x) && length(x) > 1) || (is.list(x) && !is.data.frame(x))
    if (multi) {
        files <- if (is.list(x)) x else as.list(x)
        reps <- lapply(files, readPeaks, format = format, verbose = verbose)
        if (is.null(minReps)) minReps <- min(2L, length(reps))
        gr <- .consensus(reps, minReps)
        if (verbose) message(sprintf("Consensus: %d peaks present in >= %d of %d replicates",
                                     length(gr), minReps, length(reps)))
    } else {
        gr <- .readOnePeakSet(x, format, verbose)
    }
    if (!is.null(blacklist)) {
        bl <- if (is.character(blacklist)) rtracklayer::import(blacklist) else blacklist
        bl <- .matchSeqStyle(bl, gr)
        drop <- IRanges::overlapsAny(gr, bl, ignore.strand = TRUE)
        if (verbose) message("Blacklist: removed ", sum(drop), " of ", length(gr), " peaks")
        gr <- gr[!drop]
    }
    gr$peak_id <- .peakId(gr)
    gr
}

.consensus <- function(reps, minReps) {
    fmt <- unique(unlist(lapply(reps, function(r) unique(r$peak_format))))
    plain <- lapply(reps, function(r) GenomicRanges::granges(r))
    plain <- lapply(plain, function(r) .matchSeqStyle(r, plain[[1]]))
    merged <- GenomicRanges::reduce(do.call(c, unname(plain)), ignore.strand = TRUE)
    n <- Reduce(`+`, lapply(plain, function(r)
        IRanges::overlapsAny(merged, r, ignore.strand = TRUE)))
    out <- merged[n >= minReps]
    out$n_reps <- as.integer(n[n >= minReps])
    out$peak_format <- paste0("consensus (", paste(fmt, collapse = "/"), ")")
    out
}

.readOnePeakSet <- function(x, format, verbose) {
    label <- if (is.character(x)) basename(x) else class(x)[1]
    if (is(x, "GRanges") || is.data.frame(x)) {
        cols <- if (is(x, "GRanges")) colnames(S4Vectors::mcols(x)) else names(x)
        if (format == "auto")
            format <- if (all(c("Fold", "FDR") %in% cols)) "diffbind" else "ranges"
        gr <- if (is(x, "GRanges")) x else asPeaks(x)
        if (format == "diffbind") gr <- .tidyDiffBind(gr)
        gr$peak_format <- format
        if (verbose) message(label, ": ", format, ", ", length(gr), " peaks")
        return(gr)
    }
    if (!is.character(x) || length(x) != 1) stop("Unsupported peak input.", call. = FALSE)
    if (!file.exists(x)) stop("Peak file not found: ", x, call. = FALSE)
    if (format == "auto") format <- .detectPeakFormat(x)
    gr <- switch(format,
        narrowPeak = .readBedCols(x, c("chr", "start", "end", "name", "score", "strand",
                                       "signal", "pvalue", "qvalue", "summit")),
        broadPeak  = .readBedCols(x, c("chr", "start", "end", "name", "score", "strand",
                                       "signal", "pvalue", "qvalue")),
        gappedPeak = .readBedCols(x, c("chr", "start", "end", "name", "score", "strand",
                                       "thick_start", "thick_end", "rgb", "block_count",
                                       "block_sizes", "block_starts", "signal", "pvalue",
                                       "qvalue")),
        seacr      = .readBedCols(x, c("chr", "start", "end", "total_signal",
                                       "max_signal", "max_region")),
        macs2_xls  = .readMacsXls(x),
        homer      = .readHomer(x, annotated = FALSE),
        homer_annotated = .readHomer(x, annotated = TRUE),
        diffbind   = .tidyDiffBind(asPeaks(utils::read.delim(x, stringsAsFactors = FALSE))),
        bed        = .readBedCols(x, NULL),
        stop("Unknown peak format '", format, "'.", call. = FALSE))
    gr$peak_format <- format
    if (verbose) message(label, ": ", format, ", ", length(gr), " peaks")
    gr
}

.firstLines <- function(path, n = 50) {
    con <- file(path, "r"); on.exit(close(con))
    readLines(con, n = n, warn = FALSE)
}

.detectPeakFormat <- function(path) {
    b <- tolower(basename(path))
    b <- sub("\\.gz$", "", b)
    if (grepl("narrowpeak$", b)) return("narrowPeak")
    if (grepl("broadpeak$", b)) return("broadPeak")
    if (grepl("gappedpeak$", b)) return("gappedPeak")
    if (grepl("_peaks\\.xls$", b)) return("macs2_xls")
    ln <- .firstLines(path)
    if (any(grepl("^PeakID \\(cmd=", ln))) return("homer_annotated")
    if (any(grepl("^#?PeakID\\t", ln))) return("homer")
    if (any(grepl("^# Command line: macs", ln))) return("macs2_xls")
    hdr <- strsplit(ln[!grepl("^(#|track|browser)", ln)][1], "\t")[[1]]
    if (all(c("Fold", "FDR") %in% hdr)) return("diffbind")
    data <- ln[!grepl("^(#|track|browser)", ln)]
    f <- strsplit(data[1], "\t")[[1]]
    if (length(f) == 6 && grepl("^[^:]+:[0-9]+-[0-9]+$", f[6])) return("seacr")
    if (length(f) == 10 && grepl("^-?[0-9]+$", f[10]) && grepl("^[0-9.eE+-]+$", f[7]))
        return("narrowPeak")
    "bed"
}

.readBedCols <- function(path, names) {
    d <- utils::read.delim(path, header = FALSE, comment.char = "#",
                           stringsAsFactors = FALSE, quote = "")
    d <- d[!grepl("^(track|browser)", d[[1]]), , drop = FALSE]
    if (!is.numeric(d[[2]])) d <- d[-1, , drop = FALSE]   # header row
    if (is.null(names)) names <- c("chr", "start", "end", "name")
    names <- names[seq_len(min(length(names), ncol(d)))]
    d <- d[, seq_along(names), drop = FALSE]
    names(d) <- names
    d$start <- as.integer(d$start); d$end <- as.integer(d$end)
    if ("strand" %in% names(d)) d$strand[!d$strand %in% c("+", "-")] <- "*"
    if ("summit" %in% names(d)) d$summit <- d$start + 1L + as.integer(d$summit)  # absolute, 1-based
    GenomicRanges::makeGRangesFromDataFrame(d, keep.extra.columns = TRUE,
                                            starts.in.df.are.0based = TRUE)
}

.readMacsXls <- function(path) {
    d <- utils::read.delim(path, comment.char = "#", stringsAsFactors = FALSE)
    names(d) <- sub("^X\\.", "", names(d))
    ren <- c(abs_summit = "summit", pileup = "pileup", log10.pvalue. = "pvalue",
             fold_enrichment = "fold_enrichment", log10.qvalue. = "qvalue", name = "name")
    keep <- intersect(names(ren), names(d))
    out <- data.frame(chr = d$chr, start = as.integer(d$start), end = as.integer(d$end))
    out[ren[keep]] <- d[keep]
    GenomicRanges::makeGRangesFromDataFrame(out, keep.extra.columns = TRUE)
}

.readHomer <- function(path, annotated) {
    ln <- readLines(path, warn = FALSE)
    hdrLine <- if (annotated) grep("^PeakID \\(cmd=", ln)[1] else grep("^#?PeakID\\t", ln)[1]
    d <- utils::read.delim(text = ln[hdrLine:length(ln)], header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE, quote = "")
    names(d)[1] <- "PeakID"
    nm <- tolower(names(d))
    col <- function(...) { i <- match(tolower(c(...)), nm); i <- i[!is.na(i)][1]
                           if (is.na(i)) NULL else d[[i]] }
    out <- data.frame(chr = col("chr"), start = as.integer(col("start")),
                      end = as.integer(col("end")), name = d$PeakID,
                      stringsAsFactors = FALSE)
    st <- col("strand"); if (!is.null(st)) out$strand <- ifelse(st %in% c("+", "-"), st, "*")
    if (annotated) {
        out$homer_gene <- col("Gene Name")
        out$homer_annotation <- sub(" \\(.*$", "", col("Annotation"))
        out$homer_distance_to_tss <- as.numeric(col("Distance to TSS"))
    } else {
        sc <- col("findPeaks Score", "Normalized Tag Count")
        if (!is.null(sc)) out$homer_score <- as.numeric(sc)
    }
    GenomicRanges::makeGRangesFromDataFrame(out, keep.extra.columns = TRUE)
}

.tidyDiffBind <- function(gr) {
    mc <- S4Vectors::mcols(gr)
    new <- S4Vectors::DataFrame(diff_fold = as.numeric(mc$Fold),
                                diff_fdr = as.numeric(mc$FDR))
    if ("Conc" %in% names(mc)) new$diff_conc <- as.numeric(mc$Conc)
    conc <- grep("^Conc_", names(mc), value = TRUE)
    for (cc in conc) new[[paste0("diff_", tolower(cc))]] <- as.numeric(mc[[cc]])
    S4Vectors::mcols(gr) <- new
    gr
}
