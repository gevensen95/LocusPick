#' Read gene models
#'
#' Builds the gene-level ranges used for peak assignment, distance checks and
#' locus plots.
#'
#' @param x Path to a GTF/GFF file, or a `GRanges` with a feature-type column
#'   (`type`) and a gene name column, e.g. the result of
#'   `rtracklayer::import()` on a GTF.
#' @param geneNameCol Attribute holding the gene symbol. Falls back to
#'   `gene_id` for features where it is missing.
#' @param harmonize Optional function applied to gene symbols, for example
#'   `function(s) harmonizeSymbols(s, org.Mm.eg.db::org.Mm.eg.db)`.
#' @return A list of class `locusGeneModels` with
#'   * `genes`: one range per gene and chromosome (`gene_name`, `tss`);
#'   * `exons`: exons collapsed per gene (`gene_name`).
#' @examples
#' gtf <- system.file("extdata", "toy.gtf", package = "LocusPick")
#' gm <- readGeneModels(gtf)
#' gm$genes
#' @export
readGeneModels <- function(x, geneNameCol = "gene_name", harmonize = NULL) {
    gr <- if (is.character(x)) rtracklayer::import(x) else x
    if (!is(gr, "GRanges")) stop("x must be a GTF path or GRanges.", call. = FALSE)
    mc <- S4Vectors::mcols(gr)
    if (!"type" %in% colnames(mc))
        stop("Gene models need a 'type' column (gene/transcript/exon).", call. = FALSE)

    nm <- if (geneNameCol %in% colnames(mc)) as.character(mc[[geneNameCol]]) else
        rep(NA_character_, length(gr))
    if ("gene_id" %in% colnames(mc)) {
        miss <- is.na(nm) | nm == ""
        nm[miss] <- as.character(mc$gene_id)[miss]
    }
    if (all(is.na(nm))) stop("No gene names found in column '", geneNameCol, "'.",
                             call. = FALSE)
    if (!is.null(harmonize)) nm <- harmonize(nm)
    gr$gene_name <- nm
    gr <- gr[!is.na(gr$gene_name)]

    type <- as.character(gr$type)
    exons <- gr[type == "exon"]
    body  <- gr[type %in% c("gene", "transcript")]
    if (!length(body)) body <- exons
    if (!length(body)) stop("No gene, transcript or exon features found.", call. = FALSE)

    # one span per gene name + chromosome + strand
    key <- paste(body$gene_name, GenomicRanges::seqnames(body),
                 GenomicRanges::strand(body), sep = "\r")
    spans <- unlist(range(GenomicRanges::split(body, key)), use.names = TRUE)
    spans$gene_name <- sub("\r.*", "", names(spans))
    names(spans) <- NULL
    plus <- as.character(GenomicRanges::strand(spans)) != "-"
    spans$tss <- ifelse(plus, GenomicRanges::start(spans), GenomicRanges::end(spans))
    S4Vectors::mcols(spans) <- S4Vectors::mcols(spans)[, c("gene_name", "tss")]

    if (length(exons)) {
        ek <- paste(exons$gene_name, GenomicRanges::seqnames(exons), sep = "\r")
        ex <- unlist(GenomicRanges::reduce(GenomicRanges::split(exons, ek)),
                     use.names = TRUE)
        ex$gene_name <- sub("\r.*", "", names(ex))
        names(ex) <- NULL
    } else {
        ex <- spans[0]
    }
    structure(list(genes = spans, exons = ex), class = "locusGeneModels")
}

#' @export
print.locusGeneModels <- function(x, ...) {
    cat("locusGeneModels:", length(unique(x$genes$gene_name)), "genes,",
        length(x$exons), "collapsed exons\n")
    invisible(x)
}

.asGeneModels <- function(genes) {
    if (inherits(genes, "locusGeneModels")) return(genes)
    readGeneModels(genes)
}
