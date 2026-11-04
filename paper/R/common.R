## =============================================================================
## Helpers shared by the paper scripts. Sourced after config.R.
## =============================================================================

suppressPackageStartupMessages({
    library(LocusPick)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
})
if (packageVersion("LocusPick") < "0.4.2")
    stop("The paper scripts need LocusPick >= 0.4.2 (scaleSource, replicates). Reinstall from the repo.",
         call. = FALSE)

say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

## Fail early, listing every missing input rather than the first one.
needFiles <- function(paths, what = "input") {
    paths <- path.expand(unlist(paths))
    miss <- paths[!file.exists(paths)]
    if (length(miss))
        stop("Missing ", what, " file(s):\n  ", paste(miss, collapse = "\n  "),
             if (!TEST) "\nRun the download/processing steps in README.md first.",
             call. = FALSE)
    invisible(paths)
}

figDir <- function(fig) {
    d <- file.path(RESULTS, fig)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    d
}

writeTSV <- function(x, path) {
    utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(path)
}

saveFig <- function(p, path, width = 7, height = 5) {
    ggsave(paste0(path, ".pdf"), p, width = width, height = height)
    ggsave(paste0(path, ".png"), p, width = width, height = height, dpi = 200)
    invisible(path)
}

theme_paper <- function(base = 9) {
    theme_classic(base_size = base) +
        theme(strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              legend.key.size = unit(0.35, "cm"))
}

## Okabe-Ito, colour-blind safe
PAL <- c(LocusPick = "#0072B2", `LocusPick +change` = "#56B4E9", `LocusPick +change +reps` = "#009E73",
         top_signal = "#E69F00", top_lfc = "#D55E00",
         top_change = "#CC79A7", random = "#999999", CPM = "#E69F00", spikein = "#0072B2")

## ---- Statistics ---------------------------------------------------------------------
## Wilson score interval for a proportion.
wilson <- function(x, n, conf = 0.95) {
    z <- stats::qnorm(1 - (1 - conf) / 2)
    p <- ifelse(n > 0, x / n, NA_real_)
    den <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / den
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
    data.frame(rate = p, lo = pmax(0, ctr - half), hi = pmin(1, ctr + half))
}

## Cohen's kappa for two raters (unweighted), with a bootstrap CI.
cohenKappa <- function(a, b, B = 2000) {
    ok <- !is.na(a) & !is.na(b)
    a <- as.character(a[ok]); b <- as.character(b[ok])
    k <- function(a, b) {
        lev <- union(a, b)
        t <- table(factor(a, lev), factor(b, lev))
        po <- sum(diag(t)) / sum(t)
        pe <- sum(rowSums(t) * colSums(t)) / sum(t)^2
        if (pe == 1) return(NA_real_)
        (po - pe) / (1 - pe)
    }
    est <- k(a, b)
    bs <- replicate(B, { i <- sample(length(a), replace = TRUE); k(a[i], b[i]) })
    q <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE)
    data.frame(n = length(a), agreement = mean(a == b), kappa = est, lo = q[[1]], hi = q[[2]])
}

jaccard <- function(a, b) if (!length(union(a, b))) NA_real_ else
    length(intersect(a, b)) / length(union(a, b))

## ---- Peak/gene helpers ---------------------------------------------------------------
## The common shortcut: every peak goes to the gene with the nearest TSS, at any
## distance, with no further checks (what most annotation tools report).
nearestTSS <- function(peaks, gm) {
    tss <- resize(gm$genes, 1, fix = "start")
    pk <- peaks; strand(pk) <- "*"
    hit <- distanceToNearest(pk, tss, ignore.strand = TRUE)
    data.frame(peak_id = pk$peak_id[S4Vectors::queryHits(hit)],
               gene = tss$gene_name[S4Vectors::subjectHits(hit)],
               tss_distance = S4Vectors::mcols(hit)$distance, stringsAsFactors = FALSE)
}

## Evaluate any (peak, gene) pairs with the same checks LocusPick applies, plus
## blacklist overlap, so picks from different rules are judged identically.
pairFlags <- function(pairs, peaks, gm, signal, mappability = NULL, blacklist = NULL,
                      maxDist = 5000, minMappability = 0.8, minCoverage = 0.5) {
    pk <- peaks[match(pairs$peak_id, peaks$peak_id)]
    ## distance to the nearest copy of the gene (names repeat: PAR genes etc.),
    ## computed exactly as LocusPick does
    d <- LocusPick:::.distToGene(pk, pairs$gene, gm)
    out <- pairs
    out$gene_distance <- d
    out$near_gene <- !is.na(d) & d <= maxDist
    sig <- signalAtPeaks(unique(pk), signal)
    sig <- sig[match(pairs$peak_id, sig$peak_id), , drop = FALSE]
    covCols <- grep("_cov$", names(sig), value = TRUE)
    out$has_signal <- Reduce(`&`, lapply(covCols, function(cc) !is.na(sig[[cc]]) & sig[[cc]] >= minCoverage))
    if (!is.null(mappability)) {
        m <- mappabilityAtPeaks(unique(pk), mappability)
        out$mappability <- m[match(pairs$peak_id, unique(pk)$peak_id)]
        out$mappable <- !is.na(out$mappability) & out$mappability >= minMappability
    }
    if (!is.null(blacklist)) out$blacklisted <- overlapsAny(pk, blacklist, ignore.strand = TRUE)
    out
}

readBlacklist <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    b <- utils::read.table(path, sep = "\t", header = FALSE, colClasses = c("character", "integer", "integer"),
                           fill = TRUE)[, 1:3]
    GRanges(b[[1]], IRanges(b[[2]] + 1L, b[[3]]))
}

## ---- ENCODE ----------------------------------------------------------------------------
readManifest <- function() {
    needFiles(ENCODE$manifest, "ENCODE manifest (run 00_select_encode.R)")
    m <- utils::read.delim(ENCODE$manifest, stringsAsFactors = FALSE)
    m$file <- ifelse(file.exists(m$file), m$file, file.path(ENCODE$dir, basename(m$file)))
    m
}

encodeFile <- function(m, target, cell, output, replicate = NULL) {
    r <- m[m$target == target & m$cell == cell & m$output == output, ]
    if (!is.null(replicate)) r <- r[r$replicate %in% replicate, ]
    if (!nrow(r)) stop("No ", output, " file for ", target, " in ", cell, " in the manifest.", call. = FALSE)
    needFiles(r$file, paste(target, cell, output))
}

## HepG2 vs K562 differential expression from ENCODE RSEM gene quantifications.
## Gene IDs (GENCODE v29) are mapped to symbols with the GTF used for the peaks,
## so IDs and models always agree. Cached.
encodeRNADE <- function(m, gtf) {
    cache <- file.path(figDir("cache"), "encode_rna_de.tsv")
    if (file.exists(cache)) return(utils::read.delim(cache, stringsAsFactors = FALSE))
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is required.", call. = FALSE)
    r <- m[m$output == "rna_quant", ]
    needFiles(r$file, "ENCODE RNA-seq quantification")
    cnt <- sapply(r$file, function(f) {
        x <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
        stats::setNames(round(x$expected_count), x$gene_id)
    })
    colnames(cnt) <- paste(r$cell, r$replicate, sep = "_")
    cnt <- cnt[!grepl("^(gSpikein|ERCC|tSpikein)", rownames(cnt)), ]
    cd <- data.frame(cell = factor(r$cell, levels = rev(ENCODE$cells)), row.names = colnames(cnt))
    dds <- DESeq2::DESeqDataSetFromMatrix(cnt, cd, ~ cell)
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
    res <- DESeq2::results(dds, contrast = c("cell", ENCODE$cells[["test"]], ENCODE$cells[["ref"]]))
    g <- rtracklayer::import(gtf, feature.type = "gene")
    id2sym <- stats::setNames(g$gene_name, sub("\\.[0-9]+(_PAR_Y)?$", "", g$gene_id))
    out <- data.frame(gene = unname(id2sym[sub("\\.[0-9]+(_PAR_Y)?$", "", rownames(res))]),
                      baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                      padj = res$padj, stringsAsFactors = FALSE)
    out <- out[!is.na(out$gene) & !is.na(out$padj), ]
    out <- out[order(out$padj), ]
    out <- out[!duplicated(out$gene), ]
    writeTSV(out, cache)
    out
}

## ---- Hnf4a ------------------------------------------------------------------------------
## chromDict .rds files become chromDictSource(); anything else is passed through.
hnf4aSource <- function(path) {
    path <- path.expand(path)
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        if (!requireNamespace("chromDicts", quietly = TRUE))
            stop("chromDicts must be installed to read ", basename(path), call. = FALSE)
        chromDictSource(readRDS(path))
    } else path
}

hnf4aDE <- function() {
    needFiles(c(HNF4A$age_de, HNF4A$targets), "Hnf4a DE")
    age <- utils::read.delim(HNF4A$age_de, stringsAsFactors = FALSE)
    tg  <- utils::read.delim(HNF4A$targets, stringsAsFactors = FALSE)
    list(age = readDE(age, gene = "Gene", lfc = "log2FoldChange", padj = "padj",
                      direction = "Old/Young", verbose = FALSE),
         ko  = readDE(tg, gene = "Gene", lfc = "RNA_KO_log2FC", padj = "RNA_KO_FDR",
                      direction = "KO/WT", verbose = FALSE))
}

hnf4aPeaks <- function() {
    needFiles(HNF4A$peaks, "Hnf4a consensus peaks")
    pk <- utils::read.delim(HNF4A$peaks, header = FALSE, stringsAsFactors = FALSE)
    pk <- stats::setNames(pk[, 1:5], c("chr", "start", "end", "feature", "annot_gene"))
    if (!is.numeric(pk$start))
        stop(basename(HNF4A$peaks), " is not chr/start/end/feature/gene in columns 1-5.", call. = FALSE)
    gr <- GRanges(pk$chr, IRanges(pk$start + 1L, pk$end))   # BED is 0-based
    gr$feature <- pk$feature
    gr$annot_gene <- pk$annot_gene
    asPeaks(gr)
}

## ---- Test-mode note on every figure ------------------------------------------------------
testCaption <- function() if (TEST) labs(caption = "TEST MODE: simulated data, numbers are meaningless")
