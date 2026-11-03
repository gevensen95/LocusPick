#' Check that inputs are compatible before building candidates
#'
#' Catches the mismatches that otherwise fail silently: different chromosome
#' naming (`chr1` vs `1`), a different genome build, peaks outside the signal
#' tracks, and gene groups whose symbols are missing from the gene models.
#'
#' @param peaks Peaks (anything accepted by [asPeaks()]).
#' @param genes Gene models from [readGeneModels()].
#' @param signal Optional named list of signal sources.
#' @param groups Optional gene groups from [defineGroups()] (or any
#'   `data.frame` with a `gene` column).
#' @param maxDist Distance used to count peaks near any gene.
#' @param sample Number of peaks sampled when querying signal tracks.
#' @return Invisibly, a `data.frame` of checks with columns `check`, `value`
#'   and `status` (`ok`, `warn` or `fail`). Results are also printed.
#' @examples
#' ex <- locusPickExample()
#' checkInputs(ex$peaks, ex$genes, ex$signal, ex$groups)
#' @export
checkInputs <- function(peaks, genes, signal = NULL, groups = NULL,
                        maxDist = 5000, sample = 500) {
    pk <- asPeaks(peaks)
    gm <- .asGeneModels(genes)
    res <- list()
    add <- function(check, value, status) {
        res[[length(res) + 1]] <<- data.frame(check = check, value = value,
                                              status = status,
                                              stringsAsFactors = FALSE)
    }
    style <- function(x) tryCatch(paste(GenomeInfoDb::seqlevelsStyle(x), collapse = "/"),
                                  error = function(e) "unknown")

    pStyle <- style(pk); gStyle <- style(gm$genes)
    add("chromosome naming: peaks / genes", paste(pStyle, "/", gStyle),
        if (identical(pStyle, gStyle)) "ok" else "warn")

    g <- .matchSeqStyle(gm$genes, pk)
    pchr <- unique(as.character(GenomicRanges::seqnames(pk)))
    shared <- mean(pchr %in% as.character(GenomicRanges::seqnames(g)))
    add("peak chromosomes found in gene models", sprintf("%.0f%%", 100 * shared),
        if (shared == 1) "ok" else if (shared >= 0.8) "warn" else "fail")

    gx <- suppressWarnings(GenomicRanges::trim(
        GenomicRanges::resize(g, GenomicRanges::width(g) + 2 * maxDist, fix = "center")))
    near <- mean(IRanges::overlapsAny(pk, gx, ignore.strand = TRUE))
    add(sprintf("peaks within %s bp of a gene", format(maxDist, big.mark = ",")),
        sprintf("%.0f%%", 100 * near),
        if (near >= 0.5) "ok" else if (near >= 0.2) "warn" else "fail")

    if (length(signal)) {
        set.seed(1)
        sp <- if (length(pk) > sample) pk[sort(sample(length(pk), sample))] else pk
        for (nm in names(signal)) {
            s <- tryCatch(signalAtPeaks(sp, signal[nm]), error = function(e) e)
            if (inherits(s, "error")) {
                add(paste0("signal '", nm, "' readable"), conditionMessage(s), "fail")
                next
            }
            cov <- s[[paste0(nm, "_cov")]]
            any <- mean(cov > 0)
            add(paste0("signal '", nm, "': peaks with data"), sprintf("%.0f%%", 100 * any),
                if (any >= 0.95) "ok" else if (any >= 0.8) "warn" else "fail")
        }
    }

    if (!is.null(groups)) {
        .stopIfMissingCols(groups, "gene", "groups")
        inM <- mean(unique(groups$gene) %in% g$gene_name)
        add("grouped genes found in gene models", sprintf("%.0f%%", 100 * inM),
            if (inM >= 0.98) "ok" else if (inM >= 0.9) "warn" else "fail")
    }

    out <- do.call(rbind, res)
    hint <- c(
        "chromosome naming: peaks / genes" =
            "LocusPick converts styles automatically, but check that both use the same genome build.",
        "peak chromosomes found in gene models" =
            "Peaks on contigs missing from the GTF (e.g. random/Un/MT) cannot be assigned.",
        "grouped genes found in gene models" =
            "Symbols differ between sources; see harmonizeSymbols().")
    w <- max(nchar(out$check))
    for (i in seq_len(nrow(out))) {
        cat(sprintf("[%-4s] %-*s %s\n", out$status[i], w, out$check[i], out$value[i]))
        h <- hint[out$check[i]]
        if (out$status[i] != "ok" && !is.na(h)) cat("       ", h, "\n")
    }
    if (any(out$status == "fail" & grepl("within", out$check)))
        cat("        Few peaks are near genes: are peaks and gene models from the same genome build?\n")
    if (any(out$status != "ok" & grepl("^signal", out$check)))
        cat("        Peaks without signal: check build/naming of the track, or gaps in its coverage.\n")
    invisible(out)
}
