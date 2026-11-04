## =============================================================================
## Test mode: simulated stand-ins for every real input of the paper scripts.
##
## Sourced by config.R when LOCUSPICK_PAPER_TEST=1. Builds (once) a small
## two-chromosome genome with files laid out like the real downloads and
## pipeline outputs, then points REF, ENCODE, CHIPRX, PBMC and HNF4A at them.
##
## The simulation exists only to check that the scripts run end to end. Its
## numbers mean nothing; do not report them.
##
## Built-in pitfalls, so every check in the scripts has something to find:
##   * a low-mappability "paralog" block with inflated signal (chr2:2.00-2.30 Mb)
##   * a blacklisted region with a very strong artefact peak (chr1:5.50-5.52 Mb)
##   * a coverage gap with no signal at all (chr1:3.00-3.10 Mb)
##   * distal peaks far from every gene (wrong nearest-TSS assignments)
##   * a genome-wide loss of signal in the ChIP-Rx "treated" and Hnf4a "Old" samples
## =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

.SIM_SL <- c(chr1 = 6e6, chr2 = 4e6)
.SIM_BIN <- 100L

.simGenes <- function(n = 400) {
    set.seed(SEED)
    out <- list()
    for (chr in names(.SIM_SL)) {
        k <- round(n * .SIM_SL[[chr]] / sum(.SIM_SL))
        starts <- sort(sample(seq(5e4, .SIM_SL[[chr]] - 6e4, by = 1000), k))
        len <- round(stats::runif(k, 3e3, 4e4))
        ok <- c(TRUE, diff(starts) > 1.2e4)          # leave room between genes
        out[[chr]] <- GRanges(chr, IRanges(starts[ok], width = len[ok]),
                              strand = sample(c("+", "-"), sum(ok), TRUE))
    }
    g <- unlist(GRangesList(unname(out)))
    g$gene_name <- sprintf("GENE%03d", seq_along(g))
    g$gene_name[overlapsAny(g, GRanges("chr2:2000000-2300000"))] <-
        sprintf("PARA%d", seq_len(sum(overlapsAny(g, GRanges("chr2:2000000-2300000")))))
    g$gene_id <- sprintf("ENSG%011d.1", seq_along(g))
    seqlengths(g) <- .SIM_SL[seqlevels(g)]
    g
}

.writeGTF <- function(g, path) {
    tx <- g; tx$type <- "transcript"; tx$transcript_id <- paste0(g$gene_id, "_t1")
    gn <- g; gn$type <- "gene"; gn$transcript_id <- NA
    ex <- do.call(c, lapply(seq_along(g), function(i) {
        n <- sample(3:7, 1); w <- width(g)[i]
        s <- sort(unique(c(0, sample(seq(500, w - 600, by = 100), n - 1))))
        e <- GRanges(seqnames(g)[i], IRanges(start(g)[i] + s, width = pmin(400, w - s)),
                     strand = strand(g)[i])
        e$gene_name <- g$gene_name[i]; e$gene_id <- g$gene_id[i]
        e$transcript_id <- paste0(g$gene_id[i], "_t1"); e$type <- "exon"
        e
    }))
    all <- c(gn, tx, ex)
    all$source <- "sim"
    all <- sort(all)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    plain <- sub("\\.gz$", "", path)
    export(all, plain, format = "gtf")
    if (grepl("\\.gz$", path)) { R.utils_gzip(plain, path) }
}

R.utils_gzip <- function(from, to) {
    con <- gzfile(to, "w"); writeLines(readLines(from), con); close(con); unlink(from)
}

.tss <- function(g) resize(g, 1, fix = "start")

## Peaks: promoter, gene-proximal and distal; plus artefact peaks.
.simPeaks <- function(g, n = 1500, prefix = "pk") {
    k <- round(n * c(0.6, 0.25, 0.15))
    tss <- .tss(g)
    p1 <- shift(tss[sample(length(g), k[1], TRUE)], round(stats::rnorm(k[1], 0, 150)))
    gi <- sample(length(g), k[2], TRUE)
    p2 <- GRanges(seqnames(g)[gi],
                  IRanges(pmax(1, start(g)[gi] + round(stats::runif(k[2], -4000, width(g)[gi]))), width = 1))
    chr <- sample(names(.SIM_SL), k[3], TRUE, prob = .SIM_SL)
    p3 <- GRanges(chr, IRanges(round(stats::runif(k[3], 1e4, .SIM_SL[chr] - 1e4)), width = 1))
    p3 <- p3[distance(p3, g[nearest(p3, g)]) > 2e4 | is.na(nearest(p3, g))]
    art <- GRanges(c("chr1", "chr2", "chr2", "chr2", "chr1", "chr1"),
                   IRanges(c(5510000, 2050000, 2150000, 2250000, 3040000, 3070000), width = 1))
    p <- c(granges(p1), granges(p2), granges(p3), art)
    strand(p) <- "*"
    p <- resize(p, round(stats::runif(length(p), 200, 600)), fix = "center")
    p <- reduce(p)
    seqlengths(p) <- .SIM_SL[seqlevels(p)]
    p <- trim(p)
    p$name <- sprintf("%s%05d", prefix, seq_along(p))
    p
}

## Coverage track (100-bp bins) = background + peak heights; gaps removed.
.simTrack <- function(peaks, height, path, bg = 0.3, gaps = GRanges("chr1:3000000-3100000"),
                      artefacts = TRUE) {
    pk <- peaks; strand(pk) <- "*"
    if (artefacts) {
        add <- GRanges(c("chr1", "chr2"), IRanges(c(5500000, 2000000), c(5520000, 2300000)))
        pk <- c(granges(pk), add)
        height <- c(height, 40, 6)
    }
    cov <- coverage(pk, weight = height)
    cov <- cov[names(.SIM_SL)]
    tiles <- tileGenome(.SIM_SL, tilewidth = .SIM_BIN, cut.last.tile.in.chrom = TRUE)
    for (chr in names(.SIM_SL)) cov[[chr]] <- c(cov[[chr]], Rle(0, max(0, .SIM_SL[[chr]] - length(cov[[chr]]))))
    tiles <- binnedAverage(tiles, cov, "score")
    tiles$score <- tiles$score + stats::rgamma(length(tiles), 2, 2 / bg)
    if (!is.null(gaps)) tiles <- tiles[!overlapsAny(tiles, gaps)]
    seqlengths(tiles) <- .SIM_SL
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    export.bw(tiles, path)
    invisible(path)
}

.writeNarrowPeak <- function(p, path, score = stats::runif(length(p), 5, 50)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    d <- data.frame(as.character(seqnames(p)), start(p) - 1L, end(p), p$name, 500L, ".",
                    score, 10, 8, round(width(p) / 2))
    con <- gzfile(path, "w"); utils::write.table(d, con, sep = "\t", quote = FALSE,
                                                 row.names = FALSE, col.names = FALSE); close(con)
}

.simulateAll <- function(dir) {
    set.seed(SEED)
    message("Simulating test data in ", dir, " (once) ...")
    g <- .simGenes()
    .writeGTF(g, REF$hg38$gtf)
    ## blacklist and mappability
    bl <- GRanges("chr1", IRanges(5500001, 5520000))
    con <- gzfile(REF$hg38$blacklist, "w")
    writeLines(sprintf("%s\t%d\t%d\tHigh Signal Region", seqnames(bl), start(bl) - 1L, end(bl)), con)
    close(con)
    tiles <- tileGenome(.SIM_SL, tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
    tiles$score <- pmin(1, stats::rbeta(length(tiles), 20, 1))
    tiles$score[overlapsAny(tiles, GRanges("chr2:2000000-2300000"))] <- 0.3
    seqlengths(tiles) <- .SIM_SL
    export.bw(tiles, REF$hg38$umap)
    for (f in c(REF$mm10$blacklist)) file.copy(REF$hg38$blacklist, f, overwrite = TRUE)

    ## ---- ENCODE-like: 2 TFs + ATAC in two cell lines, RNA-seq quantifications ---
    nG <- length(g)
    mu <- exp(stats::rnorm(nG, 5, 1.5))
    de <- stats::rbinom(nG, 1, 0.3) * stats::rnorm(nG, 0, 2)
    man <- list()
    for (cell in ENCODE$cells) {
        f <- if (cell == ENCODE$cells[["test"]]) 2^(de / 2) else 2^(-de / 2)
        for (r in 1:2) {
            cnt <- stats::rnbinom(nG, mu = mu * f, size = 10)
            path <- file.path(ENCODE$dir, sprintf("rna_%s_rep%d.tsv", cell, r))
            utils::write.table(data.frame(gene_id = g$gene_id, `transcript_id(s)` = paste0(g$gene_id, "_t1"),
                                          length = width(g), effective_length = width(g),
                                          expected_count = cnt, TPM = cnt / sum(cnt) * 1e6,
                                          FPKM = 0, check.names = FALSE),
                               path, sep = "\t", quote = FALSE, row.names = FALSE)
            man[[length(man) + 1]] <- data.frame(target = "RNA", cell = cell, output = "rna_quant",
                                                 replicate = r, accession = sprintf("ENCFFSIMRNA%s%d", cell, r),
                                                 file = path)
        }
    }
    for (tf in c("TFA", "TFB", "ATAC")) {
        pk <- .simPeaks(g, if (tf == "ATAC") 2500 else 1200, prefix = tolower(tf))
        near <- nearest(pk, .tss(g))
        lfcAt <- ifelse(is.na(near), 0, de[near]) * (distance(pk, .tss(g)[pmax(1, near)]) < 2000)
        h0 <- exp(stats::rnorm(length(pk), 1.5, 0.8))
        for (cell in ENCODE$cells) {
            s <- if (cell == ENCODE$cells[["test"]]) 1 else -1
            h <- h0 * 2^(s * 0.35 * lfcAt + stats::rnorm(length(pk), 0, 0.3))
            keep <- h > stats::quantile(h, 0.15)
            .writeNarrowPeak(pk[keep], f <- file.path(ENCODE$dir, sprintf("%s_%s_peaks.bed.gz", tf, cell)))
            man[[length(man) + 1]] <- data.frame(target = tf, cell = cell, output = "peaks", replicate = NA,
                                                 accession = sprintf("ENCFFSIM%s%sPK", tf, cell), file = f)
            reps <- character()
            for (r in 1:2) {
                fr <- file.path(ENCODE$dir, sprintf("%s_%s_rep%d_fc.bw", tf, cell, r))
                .simTrack(pk, h * exp(stats::rnorm(length(pk), 0, 0.25)), fr)
                man[[length(man) + 1]] <- data.frame(target = tf, cell = cell, output = "fc_rep", replicate = r,
                                                     accession = sprintf("ENCFFSIM%s%sR%d", tf, cell, r), file = fr)
            }
            fp <- file.path(ENCODE$dir, sprintf("%s_%s_pooled_fc.bw", tf, cell))
            .simTrack(pk, h, fp)
            man[[length(man) + 1]] <- data.frame(target = tf, cell = cell, output = "fc_pooled", replicate = NA,
                                                 accession = sprintf("ENCFFSIM%s%sFC", tf, cell), file = fp)
        }
    }
    man <- do.call(rbind, man)
    man$experiment <- paste0("ENCSRSIM", match(paste(man$target, man$cell), unique(paste(man$target, man$cell))))
    man$assembly <- "GRCh38"
    utils::write.table(man, ENCODE$manifest, sep = "\t", quote = FALSE, row.names = FALSE)

    ## ---- ChIP-Rx-like: H3K79me2 over 5' gene bodies, global loss on treatment ----
    expressed <- g[mu > stats::median(mu)]
    dom <- resize(expressed, pmin(width(expressed), 10000), fix = "start"); strand(dom) <- "*"
    dom$name <- sprintf("k79_%04d", seq_along(dom))
    trueCtrl <- exp(stats::rnorm(length(dom), 1, 0.5))
    lossFrac <- stats::runif(length(dom), 0.1, 0.5)        # treated keeps 10-50 %
    nf <- file.path(CHIPRX$nfcore, "bwa/mergedLibrary")
    design <- list()
    for (cond in c(CHIPRX$ref, CHIPRX$test)) for (r in 1:2) {
        tru <- trueCtrl * if (cond == CHIPRX$test) lossFrac else 1
        tru <- tru * exp(stats::rnorm(length(tru), 0, 0.15))
        human <- 2e7 * sum(tru * width(dom)) / sum(trueCtrl * width(dom))
        fly <- 2e6 * exp(stats::rnorm(1, 0, 0.05))
        cpm <- tru / sum(tru * width(dom)) * 1e6 * (human / (human + fly)) * 2e3
        id <- sprintf("%s_REP%d", cond, r)
        .simTrack(dom, cpm, file.path(nf, "bigwig", paste0(id, ".bigWig")), bg = 0.05,
                  gaps = NULL, artefacts = FALSE)
        bp <- file.path(nf, "macs2/broadPeak", paste0(id, "_peaks.broadPeak"))
        dir.create(dirname(bp), recursive = TRUE, showWarnings = FALSE)
        kp <- dom[tru > stats::quantile(tru, if (cond == CHIPRX$test) 0.3 else 0.05)]
        utils::write.table(data.frame(as.character(seqnames(kp)), start(kp) - 1L, end(kp), kp$name,
                                      500, ".", 5, 10, 8), bp, sep = "\t", quote = FALSE,
                           row.names = FALSE, col.names = FALSE)
        ix <- file.path(nf, "samtools_stats", paste0(id, ".mLb.clN.sorted.bam.idxstats"))
        dir.create(dirname(ix), recursive = TRUE, showWarnings = FALSE)
        hs <- round(human * .SIM_SL / sum(.SIM_SL))
        utils::write.table(data.frame(c(names(.SIM_SL), "dm6_chr2L", "dm6_chr3R", "*"),
                                      c(.SIM_SL, 23513712, 32079331, 0),
                                      c(hs, round(fly * 0.45), round(fly * 0.55), 0), 0),
                           ix, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        design[[length(design) + 1]] <- data.frame(gsm = sprintf("GSMSIM%d", length(design) + 1),
                                                   srr = sprintf("SRRSIM%d", length(design) + 1),
                                                   sample = cond, replicate = r, condition = cond,
                                                   antibody = "H3K79me2", is_input = FALSE,
                                                   control = "", control_replicate = NA)
    }
    CHIPRX$design <<- file.path(dir, "gse60104_design.csv")
    utils::write.csv(do.call(rbind, design), CHIPRX$design, row.names = FALSE)

    ## ---- Hnf4a-like: C&R / ATAC / H3K27ac young vs old, KO and age DE ---------
    hg <- g; hg$gene_name <- sub("^GENE", "Gene", hg$gene_name)    # mouse-style symbols
    .writeGTF(hg, HNF4A$gtf)
    file.copy(REF$hg38$umap, HNF4A$umap, overwrite = TRUE)
    pk <- .simPeaks(hg, 1500, prefix = "cnr")
    near <- nearest(pk, hg)
    dist <- distance(pk, hg[pmax(1, near)])
    feat <- ifelse(is.na(dist) | dist > 5000, "Intergenic >5kb",
                   ifelse(distance(pk, .tss(hg)[pmax(1, near)]) < 1000, "Promoter-TSS", "intron"))
    annGene <- hg$gene_name[pmax(1, near)]
    swap <- sample(length(pk), 25)                                # wrong annotations
    annGene[swap] <- sample(hg$gene_name, 25)
    dir.create(dirname(HNF4A$peaks), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(data.frame(as.character(seqnames(pk)), start(pk) - 1L, end(pk), feat, annGene),
                       HNF4A$peaks, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    ko  <- stats::rbinom(nG, 1, 0.25) * stats::rnorm(nG, 0, 1.5)
    age <- stats::rbinom(nG, 1, 0.25) * stats::rnorm(nG, 0, 1.2)
    pv <- function(l) ifelse(l == 0, stats::runif(nG, 0.06, 1), 10^-stats::runif(nG, 2, 8))
    utils::write.table(data.frame(Gene = hg$gene_name, baseMean = mu, log2FoldChange = age + stats::rnorm(nG, 0, 0.05),
                                  padj = pv(age)), HNF4A$age_de, sep = "\t", quote = FALSE, row.names = FALSE)
    tgt <- unique(annGene)
    ix <- match(tgt, hg$gene_name)
    utils::write.table(data.frame(Gene = tgt, RNA_KO_log2FC = ko[ix] + stats::rnorm(length(ix), 0, 0.05),
                                  RNA_KO_FDR = pv(ko)[ix]), HNF4A$targets, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    nearG <- match(annGene, hg$gene_name)
    h <- exp(stats::rnorm(length(pk), 1.5, 0.7)) * 2^(-0.4 * pmin(0, ko[nearG]))
    oldF <- 0.6 * 2^(0.3 * age[nearG] + stats::rnorm(length(pk), 0, 0.2))   # global loss with age
    cdDir <- dirname(HNF4A$chromDicts[[1]])
    for (a in c("Yng", "Old")) {
        hh <- if (a == "Old") h * oldF else h
        .simTrack(pk, hh, file.path(cdDir, sprintf("CNR_%s.bw", a)))
        for (r in 1:3) .simTrack(pk, hh * exp(stats::rnorm(length(pk), 0, 0.3)), HNF4A$cnr_reps[[a]][r])
        .simTrack(pk, hh * 0.5 + 1, file.path(cdDir, sprintf("ATAC_%s.bw", a)))
        .simTrack(pk, hh * 0.8 + 0.5, file.path(cdDir, sprintf("H3K27ac_%s.bw", a)))
    }
    smp <- c(paste0("Y", 1:4), paste0("O", 1:4))
    cnt <- sapply(smp, function(s) stats::rnbinom(nG, mu = mu * 2^(age * (substr(s, 1, 1) == "O")), size = 20))
    rownames(cnt) <- hg$gene_name
    utils::write.table(cnt, HNF4A$age_counts, sep = "\t", quote = FALSE)
    utils::write.table(data.frame(X = smp, Age = ifelse(grepl("^O", smp), "Old", "Young")),
                       HNF4A$age_meta, sep = "\t", quote = FALSE, row.names = FALSE)

    ## ---- PBMC-like and liver sn-like single-cell objects -------------------------
    if (requireNamespace("Signac", quietly = TRUE)) {
        so <- suppressMessages(suppressWarnings(LocusPick::locusPickSingleCellExample(seurat = TRUE)$object))
        so$severity <- ifelse(so$age == "old", "decompensated", "stable")
        so$predicted.celltype.l2 <- so$cell_type
        so$donor <- so$sample
        dir.create(dirname(PBMC$object), recursive = TRUE, showWarnings = FALSE)
        saveRDS(so, PBMC$object)
        so$age_lbl <- so$age
        saveRDS(so, HNF4A$sn_rna)
        saveRDS(so, HNF4A$sn_atac)
    }
    writeLines(format(Sys.time()), file.path(dir, ".simulated"))
}

## ---- Point the configuration at the simulated files ---------------------------------
REF$hg38 <- list(gtf = file.path(DATA, "ref/hg38/genes.gtf.gz"),
                 blacklist = file.path(DATA, "ref/hg38/blacklist.bed.gz"),
                 umap = file.path(DATA, "ref/hg38/umap.bw"))
REF$mm10 <- list(blacklist = file.path(DATA, "ref/mm10/blacklist.bed.gz"))
ENCODE$manifest <- file.path(DATA, "encode/manifest.tsv")
ENCODE$dir <- file.path(DATA, "encode")
CHIPRX$nfcore <- file.path(DATA, "gse60104/chipseq")
CHIPRX$design <- file.path(DATA, "gse60104_design.csv")
PBMC$object <- file.path(DATA, "gse215915/pbmc_multiome.rds")
PBMC$condition <- c("severity", "decompensated", "stable")
HNF4A$age_de     <- file.path(DATA, "hnf4a/JvsJ_allGenes.txt")
HNF4A$age_counts <- file.path(DATA, "hnf4a/NormCounts_JvJ.txt")
HNF4A$age_meta   <- file.path(DATA, "hnf4a/meta.txt")
HNF4A$targets    <- file.path(DATA, "hnf4a/CUTNRUN_Targets_Yng.tsv")
HNF4A$ko_expr    <- file.path(DATA, "hnf4a/NormExp.tsv")
HNF4A$peaks      <- file.path(DATA, "hnf4a/ConsensusPeaks_YngAll_annot.bed")
HNF4A$gtf        <- file.path(DATA, "hnf4a/mm10.gtf")
HNF4A$umap       <- file.path(DATA, "hnf4a/umap.bw")
HNF4A$chromDicts <- file.path(DATA, "hnf4a/tracks", paste0(c("CNR_Yng", "CNR_Old", "ATAC_Yng", "ATAC_Old",
                                                             "H3K27ac_Yng", "H3K27ac_Old"), ".bw"))
names(HNF4A$chromDicts) <- c("CNR_Yng", "CNR_Old", "ATAC_Yng", "ATAC_Old", "H3K27ac_Yng", "H3K27ac_Old")
HNF4A$cnr_reps <- list(Yng = file.path(DATA, sprintf("hnf4a/reps/CNR_Yng_rep%d.bw", 1:3)),
                       Old = file.path(DATA, sprintf("hnf4a/reps/CNR_Old_rep%d.bw", 1:3)))
HNF4A$sn_rna  <- file.path(DATA, "hnf4a/sn_rna.rds")
HNF4A$sn_atac <- file.path(DATA, "hnf4a/sn_atac.rds")
HNF4A$sn_cols <- list(celltype_rna = "cell_type", celltype_atac = "cell_type",
                      sample_rna = "sample", sample_atac = "sample",
                      age_rna = "age", age_atac = "age", old = "old", young = "young",
                      rna_assay = "RNA")

if (!file.exists(file.path(DATA, ".simulated"))) {
    for (d in c("ref/hg38", "ref/mm10", "encode", "hnf4a/tracks", "hnf4a/reps"))
        dir.create(file.path(DATA, d), recursive = TRUE, showWarnings = FALSE)
    .simulateAll(DATA)
}
