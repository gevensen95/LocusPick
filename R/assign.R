#' Assign peaks to genes
#'
#' Links each peak to every gene whose body lies within `maxDist` bp, and
#' classifies the position of the peak relative to that gene.
#'
#' @param peaks Peaks (anything accepted by [asPeaks()]).
#' @param genes Gene models from [readGeneModels()] (or a GTF path).
#' @param maxDist Maximum gap in bp between the peak and the gene body.
#' @param promoter Length-2 numeric: bp upstream and downstream of the TSS
#'   counted as promoter.
#' @param mode `"all"` keeps every gene within `maxDist`; `"nearest"` keeps
#'   only the closest gene(s) per peak (ties are kept).
#' @return A `data.frame` with one row per peak-gene pair: `peak_id`,
#'   `seqnames`, `start`, `end`, `gene`, `distance` (0 when the peak overlaps
#'   the gene body), `tss_distance` (signed, upstream negative), and
#'   `feature` (`promoter`, `gene body`, `upstream`, `downstream`).
#' @examples
#' gm <- readGeneModels(system.file("extdata", "toy.gtf", package = "LocusPick"))
#' pk <- asPeaks(system.file("extdata", "toy_peaks.bed", package = "LocusPick"))
#' head(assignPeaks(pk, gm))
#' @export
assignPeaks <- function(peaks, genes, maxDist = 5000, promoter = c(1000, 100),
                        mode = c("all", "nearest")) {
    mode <- match.arg(mode)
    pk <- asPeaks(peaks)
    gm <- .asGeneModels(genes)
    g <- .matchSeqStyle(gm$genes, pk)

    gx <- GenomicRanges::resize(g, GenomicRanges::width(g) + 2 * maxDist, fix = "center")
    suppressWarnings(gx <- GenomicRanges::trim(gx))
    hits <- GenomicRanges::findOverlaps(pk, gx, ignore.strand = TRUE)
    qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
    if (!length(qh)) return(.emptyAssignment())

    p <- pk[qh]; gg <- g[sh]
    dist <- GenomicRanges::distance(p, gg, ignore.strand = TRUE)
    dist[is.na(dist)] <- 0L
    mid <- (GenomicRanges::start(p) + GenomicRanges::end(p)) %/% 2
    minus <- as.character(GenomicRanges::strand(gg)) == "-"
    tssd <- ifelse(minus, gg$tss - mid, mid - gg$tss)

    feature <- ifelse(tssd >= -promoter[1] & tssd <= promoter[2], "promoter",
               ifelse(dist == 0, "gene body",
               ifelse(tssd < 0, "upstream", "downstream")))
    out <- data.frame(
        peak_id = p$peak_id,
        seqnames = as.character(GenomicRanges::seqnames(p)),
        start = GenomicRanges::start(p), end = GenomicRanges::end(p),
        gene = gg$gene_name, distance = as.integer(dist),
        tss_distance = as.integer(tssd), feature = feature,
        stringsAsFactors = FALSE)
    out <- out[order(out$peak_id, out$distance, abs(out$tss_distance)), ]
    out <- out[!duplicated(out[, c("peak_id", "gene")]), ]
    if (mode == "nearest") {
        best <- stats::ave(out$distance, out$peak_id, FUN = min)
        out <- out[out$distance == best, ]
    }
    rownames(out) <- NULL
    out
}

.emptyAssignment <- function() {
    data.frame(peak_id = character(), seqnames = character(), start = integer(),
               end = integer(), gene = character(), distance = integer(),
               tss_distance = integer(), feature = character(),
               stringsAsFactors = FALSE)
}

#' Check an existing peak-to-gene annotation against gene models
#'
#' Annotation tools and databases sometimes assign peaks to genes that are far
#' away (for example enhancer-target predictions) or use gene coordinates that
#' disagree with the gene models you plot. This function measures the actual
#' distance from each peak to its assigned gene in `genes`.
#'
#' @param peaks Peaks with a gene column (anything accepted by [asPeaks()]).
#' @param genes Gene models from [readGeneModels()].
#' @param geneCol Column of `peaks` holding the assigned gene symbol.
#' @param maxDist Distance in bp beyond which an assignment is flagged.
#' @return `data.frame` with `peak_id`, `gene`, `gene_in_models`, `distance`
#'   (to the gene body on the same chromosome; `NA` if not found) and
#'   `near_gene`.
#' @examples
#' gm <- readGeneModels(system.file("extdata", "toy.gtf", package = "LocusPick"))
#' pk <- asPeaks(system.file("extdata", "toy_peaks.bed", package = "LocusPick"))
#' pk$gene <- "geneA"
#' checkAssignment(pk, gm, geneCol = "gene")
#' @export
checkAssignment <- function(peaks, genes, geneCol = "gene", maxDist = 5000) {
    pk <- asPeaks(peaks)
    gm <- .asGeneModels(genes)
    if (!geneCol %in% colnames(S4Vectors::mcols(pk)))
        stop("Peaks have no column '", geneCol, "'.", call. = FALSE)
    g <- .matchSeqStyle(gm$genes, pk)
    gname <- as.character(S4Vectors::mcols(pk)[[geneCol]])
    pchr <- as.character(GenomicRanges::seqnames(pk))
    gchr <- as.character(GenomicRanges::seqnames(g))

    dist <- rep(NA_real_, length(pk))
    idx <- split(seq_along(g), paste(g$gene_name, gchr, sep = "\r"))
    key <- paste(gname, pchr, sep = "\r")
    for (i in which(key %in% names(idx))) {
        d <- GenomicRanges::distance(pk[i], g[idx[[key[i]]]], ignore.strand = TRUE)
        dist[i] <- min(d, na.rm = TRUE)
    }
    data.frame(peak_id = pk$peak_id, gene = gname,
               gene_in_models = gname %in% g$gene_name,
               distance = dist,
               near_gene = !is.na(dist) & dist <= maxDist,
               stringsAsFactors = FALSE)
}
