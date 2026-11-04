# Generates the synthetic single-cell example data in inst/extdata.
# Run from the package root:  Rscript inst/scripts/make_singlecell_example.R
#
# Design: 8 samples = 2 genotypes (WT, KO) x 2 ages (young, old) x 2 replicates,
# two cell types (Hepatocyte, Kupffer), 25 cells per sample and cell type.
# Uses the same toy chromosome, genes and peaks as the bulk example.
#   Hepatocytes: geneA activated by the factor and lost with age;
#                geneB repressed by the factor and up with age;
#                accessibility at the geneA/geneB promoter peaks drops with age.
#   Kupffer:     no genotype or age effects; open chromatin only at geneF.
suppressPackageStartupMessages({ library(Matrix); library(Rsamtools) })
set.seed(11)
out <- "inst/extdata"

samples <- expand.grid(rep = 1:2, age = c("young", "old"), genotype = c("WT", "KO"),
                       stringsAsFactors = FALSE)
samples$sample <- with(samples, paste(genotype, age, rep, sep = "_"))
cellTypes <- c("Hepatocyte", "Kupffer")
nCells <- 25
meta <- merge(samples[, c("sample", "genotype", "age")],
              data.frame(cell_type = cellTypes), by = NULL)
meta <- meta[rep(seq_len(nrow(meta)), each = nCells), ]
meta <- meta[order(meta$sample, meta$cell_type), ]
barcode <- function(n) vapply(seq_len(n), function(i)
    paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = ""), "")
meta$barcode <- paste0(make.unique(barcode(nrow(meta))), "-1")
rownames(meta) <- meta$barcode

## RNA counts ---------------------------------------------------------------
genes <- c(paste0("gene", LETTERS[1:6]), sprintf("filler%02d", 1:60))
base <- c(geneA = 8, geneB = 3, geneC = 4, geneD = 4, geneE = 5, geneF = 2,
          setNames(rgamma(60, 2, 0.5), sprintf("filler%02d", 1:60)))
lfc <- function(g, ct, genotype, age) {
    if (ct != "Hepatocyte") return(if (g == "geneF") 2 else 0)
    ko <- genotype == "KO"; old <- age == "old"
    switch(g,
           geneA = -1.5 * ko - 1.0 * old,
           geneB =  1.2 * ko + 1.0 * old,
           geneC = -1.0 * ko + 0.8 * old,
           geneE = -0.8 * ko - 0.7 * old,
           0)
}
mu <- sapply(seq_len(nrow(meta)), function(j) {
    m <- meta[j, ]
    base * 2^vapply(genes, lfc, 0, ct = m$cell_type, genotype = m$genotype, age = m$age) *
        rgamma(1, 20, 20)
})
counts <- matrix(rnbinom(length(mu), mu = mu, size = 5), nrow = length(genes),
                 dimnames = list(genes, meta$barcode))
saveRDS(Matrix(counts, sparse = TRUE), file.path(out, "toy_sc_rna.rds"), compress = "xz")
write.table(meta[, c("barcode", "sample", "genotype", "age", "cell_type")],
            file.path(out, "toy_sc_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ATAC fragments -------------------------------------------------------------
pk <- read.delim(file.path(out, "toy_peaks.bed"), header = FALSE)
centers <- (pk$V2 + pk$V3) / 2
openness <- function(ct, age) {
    if (ct == "Kupffer") return(c(0, 0, 0, 0, 0, 0, 1.2, 0))
    if (age == "young") c(1.0, 0.8, 0.6, 0.6, 0.6, 0.3, 0.2, 0.2)
    else                c(0.4, 0.4, 0.5, 0.5, 0.5, 0.3, 0.2, 0.2)
}
frags <- do.call(rbind, lapply(seq_len(nrow(meta)), function(j) {
    m <- meta[j, ]
    o <- openness(m$cell_type, m$age)
    nPeak <- rpois(length(o), 6 * o)
    mid <- c(unlist(mapply(function(cn, k) round(rnorm(k, cn, 80)), centers, nPeak)),
             sample(500:59500, 15))                       # background fragments
    len <- sample(80:300, length(mid), replace = TRUE)
    data.frame(chr = "chr1", start = pmax(0L, as.integer(mid - len / 2)),
               end = as.integer(mid + len / 2), barcode = m$barcode, count = 1L)
}))
frags <- frags[order(frags$start, frags$end), ]
tsv <- file.path(tempdir(), "toy_sc_fragments.tsv")
writeLines(c("# id=toy", "# pipeline=synthetic"), tsv)
write.table(frags, tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE,
            append = TRUE)
gz <- bgzip(tsv, file.path(out, "toy_sc_fragments.tsv.gz"), overwrite = TRUE)
indexTabix(gz, seq = 1, start = 2, end = 3, comment = "#", zeroBased = TRUE)

## peak-gene links (as a multiome analysis would produce) --------------------------
links <- data.frame(
    peak  = c("chr1-4800-5300", "chr1-19800-20300", "chr1-27000-27400",
              "chr1-47800-48200", "chr1-27000-27400"),
    gene  = c("geneA", "geneB", "geneC", "geneE", "geneA"),
    score = c(0.62, 0.48, 0.35, 0.30, 0.12),
    pvalue = c(1e-8, 1e-6, 1e-4, 1e-3, 0.03))
write.table(links, file.path(out, "toy_sc_links.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
message("Wrote single-cell example data to ", out)
