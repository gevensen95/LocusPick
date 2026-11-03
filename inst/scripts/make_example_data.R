# Generates the synthetic example data in inst/extdata.
# Run from the package root:  Rscript inst/scripts/make_example_data.R
#
# One 60 kb chromosome with five genes. Each gene illustrates one case:
#   geneA  activated, condition down   promoter peak, signal A > B
#   geneB  repressed, condition up     promoter peak (minus strand)
#   geneC  activated, condition up     gene-body peak
#   geneD  repressed, condition down   promoter peak in a low-mappability region,
#                                      plus a downstream peak (44 kb) in a gap of
#                                      the signal tracks
#   geneE  activated, condition down   promoter peak
#   geneF  not differentially expressed, has a peak
suppressPackageStartupMessages({
    library(GenomicRanges); library(rtracklayer)
})
set.seed(42)
out <- "inst/extdata"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
chrLen <- 60000L
si <- Seqinfo("chr1", chrLen, genome = "toy")

## genes ---------------------------------------------------------------------
gdef <- data.frame(
    gene  = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF"),
    start = c(5000, 15000, 25000, 35000, 48000, 54000),
    end   = c(12000, 20000, 30000, 42000, 52000, 57000),
    strand = c("+", "-", "+", "+", "+", "-"))
gtf <- do.call(c, lapply(seq_len(nrow(gdef)), function(i) {
    g <- gdef[i, ]
    len <- g$end - g$start
    ex <- data.frame(s = g$start + c(0, round(len * 0.45), len - 600),
                     e = g$start + c(400, round(len * 0.45) + 250, len))
    ids <- list(gene_id = sprintf("G%03d", i), gene_name = g$gene)
    gene <- GRanges("chr1", IRanges(g$start, g$end), strand = g$strand,
                    type = "gene", source = "toy", gene_id = ids$gene_id,
                    gene_name = ids$gene_name, transcript_id = NA_character_)
    tx <- gene; tx$type <- "transcript"; tx$transcript_id <- paste0(ids$gene_id, ".1")
    exn <- GRanges("chr1", IRanges(ex$s, ex$e), strand = g$strand, type = "exon",
                   source = "toy", gene_id = ids$gene_id, gene_name = ids$gene_name,
                   transcript_id = paste0(ids$gene_id, ".1"))
    c(gene, tx, exn)
}))
export(gtf, file.path(out, "toy.gtf"), format = "gtf")

## peaks ---------------------------------------------------------------------
pk <- data.frame(
    chr = "chr1",
    start = c(4800, 19800, 27000, 34600, 47800, 44000, 56600, 58500),
    end   = c(5300, 20300, 27400, 35100, 48200, 44400, 57100, 58900),
    name  = paste0("peak", 1:8),
    peak_set = c("Shared", "Shared", "Young-only", "Shared", "Shared",
                 "Young-only", "Shared", "Shared"))
write.table(transform(pk, start = start - 1L)[, c("chr", "start", "end", "name")],
            file.path(out, "toy_peaks.bed"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

## signal bigWigs (20 bp bins) -----------------------------------------------
bins <- seq(1L, chrLen - 19L, by = 20L)
mid <- bins + 10
bump <- function(center, height, sd = 120) height * exp(-0.5 * ((mid - center) / sd)^2)
centers <- (pk$start + pk$end) / 2
hA <- c(1.0, 0.8, 0.6, 0.9, 0.7, 0.5, 0.6, 0.4)   # condition A (e.g. young)
hB <- c(0.4, 0.5, 0.5, 0.6, 0.4, 0.3, 0.6, 0.4)   # condition B (e.g. old)
mk <- function(h) {
    s <- 0.02 + abs(rnorm(length(bins), 0, 0.01))
    for (i in seq_along(centers)) s <- s + bump(centers[i], h[i])
    gr <- GRanges("chr1", IRanges(bins, width = 20), seqinfo = si, score = s)
    gr[!(start(gr) >= 43900 & end(gr) <= 44500)]       # coverage hole at peak6
}
export(mk(hA), file.path(out, "toy_A.bw"))
export(mk(hB), file.path(out, "toy_B.bw"))

## mappability (Umap-like: low over 34-35.5 kb) --------------------------------
mp <- GRanges("chr1", IRanges(c(1, 34000, 35501), c(33999, 35500, chrLen)),
              seqinfo = si, score = c(1, 0.3, 1))
export(mp, file.path(out, "toy_map.bw"))

## differential expression ---------------------------------------------------
de <- data.frame(
    gene      = c("geneA", "geneB", "geneC", "geneD", "geneE", "geneF"),
    ko_lfc    = c(-2.0,  1.5, -1.0,  1.0, -1.2,  0.1),
    ko_padj   = c(1e-6, 1e-5, 1e-3, 1e-4, 1e-4,  0.8),
    cond_lfc  = c(-1.0,  1.2,  0.8, -0.9, -0.7,  0.2),
    cond_padj = c(1e-4, 1e-3, 1e-2, 1e-3, 1e-3,  0.6))
write.table(de, file.path(out, "toy_de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## expression (long format, 4 replicates per group) ---------------------------
base <- c(geneA = 200, geneB = 50, geneC = 120, geneD = 80, geneE = 150, geneF = 90)
ex <- do.call(rbind, lapply(names(base), function(g) {
    r <- de[de$gene == g, ]
    rbind(
        data.frame(gene = g, panel = "Knockout", group = rep(c("Ctrl", "KO"), each = 4),
                   value = base[[g]] * 2^(rep(c(0, r$ko_lfc), each = 4) + rnorm(8, 0, 0.15))),
        data.frame(gene = g, panel = "Condition", group = rep(c("Young", "Old"), each = 4),
                   value = base[[g]] * 2^(rep(c(0, r$cond_lfc), each = 4) + rnorm(8, 0, 0.15))))
}))
ex$value <- round(ex$value, 1)
write.table(ex, file.path(out, "toy_expression.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
message("Wrote example data to ", out)
