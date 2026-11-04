## =============================================================================
## Fig. 7. Case study: Hnf4a in the aging mouse liver.
##
## A. Bulk: Hnf4a CUT&RUN (young, old) with ATAC and H3K27ac, Hnf4a-KO
##    RNA-seq (perturbation, restricted to CUT&RUN targets) and young vs old
##    hepatocyte RNA-seq (condition). Four groups:
##      G1 activated_age_down  G2 repressed_age_up  G3 activated_age_up  G4 repressed_age_down
##    Ranking as in the original selection: C&R young signal, |KO log2FC|,
##    |age log2FC|, plus C&R loss with age (young - old) for G1 and G2.
##    Also reported:
##      - how many of the peak file's own gene annotations (HOMER-style column 5)
##        fail LocusPick's distance check (checkAssignment)
##      - whether the age loss at the picks exceeds the genome-wide loss (signalShift)
##      - rank and flags of the genes in the original hand-picked panels (config: manual_picks)
## B. Single nucleus: unpaired snRNA-seq + snATAC-seq (2 young, 2 old) per cell type.
##
## Outputs (results/fig7/): bulk_candidates.tsv, bulk_top.tsv, bulk_qc.tsv,
## bulk_shift.tsv, annotation_check.tsv, manual_picks.tsv, bulk_top_loci.pdf,
## fig7_panels.pdf, sn_qc.tsv, sn_shift.tsv, sn_top_loci.pdf, bulk_vs_sn.tsv
## =============================================================================

source("config.R"); source("R/common.R")
out <- figDir("fig7")
set.seed(SEED)
orgDb <- if (!TEST && requireNamespace("org.Mm.eg.db", quietly = TRUE)) org.Mm.eg.db::org.Mm.eg.db
harm  <- if (!is.null(orgDb)) function(s) harmonizeSymbols(s, orgDb)

## ---- A. Bulk ------------------------------------------------------------------------------------
de <- hnf4aDE()
if (!is.null(orgDb)) { de$age$gene <- harm(de$age$gene); de$ko$gene <- harm(de$ko$gene) }
groups <- defineGroups(de$ko, de$age, perturbationType = "loss", conditionName = "age")
print(table(groups$group))
gm <- readGeneModels(needFiles(HNF4A$gtf), harmonize = harm)
pk <- hnf4aPeaks()
if (!is.null(orgDb)) pk$annot_gene <- harm(pk$annot_gene)
needFiles(c(HNF4A$chromDicts, HNF4A$umap), "Hnf4a tracks")
signal <- lapply(HNF4A$chromDicts, hnf4aSource)
names(signal) <- names(HNF4A$chromDicts)
tracks <- list(`Hnf4a CUT&RUN` = c("CNR_Yng", "CNR_Old"), ATAC = c("ATAC_Yng", "ATAC_Old"),
               H3K27ac = c("H3K27ac_Yng", "H3K27ac_Old"))

print(checkInputs(pk, gm, signal[c("CNR_Yng", "CNR_Old")], groups))

## The peak file's own gene column vs the gene models
ac <- checkAssignment(pk, gm, geneCol = "annot_gene")
acs <- data.frame(peaks = nrow(ac), gene_not_in_models = sum(!ac$gene_in_models, na.rm = TRUE),
                  farther_than_5kb = sum(ac$gene_in_models & !ac$near_gene, na.rm = TRUE))
writeTSV(acs, file.path(out, "annotation_check.tsv"))
print(acs)

## Candidates: LocusPick's own assignment; eligibility requires C&R signal only
cand <- buildCandidates(pk, gm, groups, signal = signal, mappability = HNF4A$umap,
                        requireSignal = c("CNR_Yng", "CNR_Old"))
cand$signal_change <- cand$CNR_Yng_mean - cand$CNR_Old_mean
base <- c(CNR_Yng_mean = 1, abs_lfc_perturbation = 1, abs_lfc_condition = 1)
w <- list(activated_age_down = c(base, signal_change = 1), repressed_age_up = c(base, signal_change = 1),
          activated_age_up = base, repressed_age_down = base)
w <- w[names(w) %in% unique(as.character(cand$group))]
ranked <- rankLoci(cand, w)
top <- topLoci(ranked, n = N_TOP)
qc <- qcSummary(cand)
shift <- signalShift(cand[cand$eligible, ], "CNR_Yng_mean", "CNR_Old_mean",
                     background = signalAtPeaks(pk, signal[c("CNR_Yng", "CNR_Old")]))
writeTSV(ranked, file.path(out, "bulk_candidates.tsv"))
writeTSV(top, file.path(out, "bulk_top.tsv"))
writeTSV(qc, file.path(out, "bulk_qc.tsv"))
writeTSV(shift, file.path(out, "bulk_shift.tsv"))
print(qc); print(shift)

## Expression insets
expr <- NULL
if (file.exists(HNF4A$age_counts) && file.exists(HNF4A$age_meta)) {
    cnt <- as.matrix(utils::read.delim(HNF4A$age_counts, check.names = FALSE))
    meta <- utils::read.delim(HNF4A$age_meta, stringsAsFactors = FALSE)
    s <- if (all(meta$X %in% colnames(cnt))) meta$X else paste0("S_", meta$X)
    if (all(s %in% colnames(cnt))) {
        if (!is.null(orgDb)) rownames(cnt) <- harm(rownames(cnt))
        age <- if ("Age" %in% names(meta)) meta$Age else meta[[2]]
        expr <- expressionFrom(cnt[, s], group = age, panel = "Age (RNA-seq)")
    } else message("Age RNA-seq sample names do not match the counts; inset skipped.")
}
if (file.exists(HNF4A$ko_expr)) {
    ko <- utils::read.delim(HNF4A$ko_expr, check.names = FALSE)
    num <- names(ko)[vapply(ko, is.numeric, TRUE)]
    if (length(num) >= 10) {                     # first 10 = young Rosa26 (5), young Hnf4a KO (5)
        km <- as.matrix(ko[, num[1:10]]); rownames(km) <- if (!is.null(orgDb)) harm(ko$Gene.name) else ko$Gene.name
        km <- km[!duplicated(rownames(km)) & !is.na(rownames(km)), ]
        expr <- rbind(expr, expressionFrom(km, group = rep(c("Rosa26-KO", "Hnf4a-KO"), each = 5),
                                           panel = "Hnf4a KO (young)"))
    }
}

reportLoci(top, file.path(out, "bulk_top_loci.pdf"), signal = signal, genes = gm,
           trackGroups = tracks, expression = expr)

## One figure-style panel per group (best pick), shared C&R axis across groups
best <- top[!duplicated(top$group), ]
lim <- c(0, max(unlist(lapply(seq_len(nrow(best)), function(i) {
    p <- GRanges(best$seqnames[i], IRanges(best$start[i] - 5000, best$end[i] + 5000))
    unlist(lapply(signal[tracks[[1]]], function(s) fetchSignal(s, p)$score)) })), na.rm = TRUE))
pdf(file.path(out, "fig7_panels.pdf"), width = 7, height = 6)
for (i in seq_len(nrow(best)))
    print(plotLocus(best[i, ], signal, gm, trackGroups = tracks, expression = expr, annotate = FALSE,
                    yLimits = stats::setNames(list(lim), names(tracks)[1])) +
              plot_annotation(title = sprintf("%s: %s", best$group[i], best$gene[i]),
                              caption = if (TEST) "TEST MODE"))
dev.off()

## The original hand-picked genes: where do they rank, and do they pass the checks?
manual <- HNF4A$manual_picks
if (TEST) manual <- stats::setNames(c(top$gene[1], ranked$gene[nrow(ranked)]),
                                    c(as.character(top$group[1]), as.character(ranked$group[nrow(ranked)])))
if (length(manual)) {
    mp <- do.call(rbind, lapply(seq_along(manual), function(i) {
        r <- ranked[ranked$gene == manual[[i]], ]
        if (!is.null(names(manual)) && names(manual)[i] != "") r <- r[as.character(r$group) == names(manual)[i], ]
        if (!nrow(r)) return(data.frame(gene = manual[[i]], group = names(manual)[i], rank = NA,
                                        of_ranked_peaks = NA, percentile = NA, eligible = NA,
                                        near_gene = NA, mappable = NA, has_signal = NA,
                                        note = "no candidate (not in group or no peak within maxDist)"))
        r <- r[order(r$rank), ][1, ]
        n <- sum(!is.na(ranked$score) & ranked$group == r$group)   # ranks are over peaks
        data.frame(gene = r$gene, group = r$group, rank = r$rank, of_ranked_peaks = n,
                   percentile = if (is.na(r$rank)) NA else round(1 - (r$rank - 1) / n, 3),
                   eligible = r$eligible,
                   near_gene = r$near_gene, mappable = r$mappable, has_signal = r$has_signal,
                   note = "")
    }))
    writeTSV(mp, file.path(out, "manual_picks.tsv"))
    print(mp)
}

## ---- B. Single nucleus (unpaired) -------------------------------------------------------------------
snOK <- file.exists(path.expand(HNF4A$sn_rna)) && file.exists(path.expand(HNF4A$sn_atac))
if (!snOK) {
    message("snRNA/snATAC objects not found; part B skipped.")
} else {
    suppressPackageStartupMessages({ library(Seurat); library(Signac) })
    sc <- HNF4A$sn_cols
    rna  <- readRDS(path.expand(HNF4A$sn_rna))
    atac <- readRDS(path.expand(HNF4A$sn_atac))
    chrom <- names(Filter(function(a) inherits(a, "ChromatinAssay"), atac@assays))[1]
    for (o in c("rna", "atac")) {
        obj <- get(o)
        obj$lp_celltype <- as.character(obj@meta.data[[sc[[paste0("celltype_", o)]]]])
        obj$lp_sample   <- as.character(obj@meta.data[[sc[[paste0("sample_", o)]]]])
        obj$lp_age      <- as.character(obj@meta.data[[sc[[paste0("age_", o)]]]])
        assign(o, obj)
    }
    DefaultAssay(rna) <- sc$rna_assay
    if (length(grep("^counts", Layers(rna[[sc$rna_assay]]))) > 1) rna <- JoinLayers(rna, assay = sc$rna_assay)
    shared <- intersect(unique(rna$lp_celltype), unique(atac$lp_celltype))
    snGenes <- if (!is.null(Annotation(atac[[chrom]]))) readGeneModels(Annotation(atac[[chrom]])) else
        if (TEST) locusPickSingleCellExample()$genes else gm
    sn <- scLocusPick(rna, atac = atac, cellType = "lp_celltype", sample = "lp_sample",
                      condition = c("lp_age", sc$old, sc$young), genes = snGenes,
                      rnaAssay = sc$rna_assay, atacAssay = chrom, links = NULL, cellTypes = shared,
                      orgDb = orgDb, covariates = if (TEST) "genotype", n = N_TOP)
    print(sn)
    writeTSV(sn$qc, file.path(out, "sn_qc.tsv"))
    if (!is.null(sn$shift)) writeTSV(sn$shift, file.path(out, "sn_shift.tsv"))
    reportLoci(sn, file.path(out, "sn_top_loci.pdf"))

    ## Do bulk age-down / age-up genes reappear among the hepatocyte sn candidates?
    hep <- grep("hepatocyte", shared, ignore.case = TRUE, value = TRUE)
    bvs <- do.call(rbind, lapply(hep, function(h) {
        snc <- sn$candidates[startsWith(as.character(sn$candidates$group), paste0(h, " |")), ]
        do.call(rbind, lapply(c("down", "up"), function(d) {
            bulkG <- unique(top$gene[grepl(paste0("age_", d, "$"), top$group)])
            snG <- unique(snc$gene[grepl(paste0("_", d, "$"), snc$group) & snc$eligible])
            data.frame(cell_type = h, direction = d, bulk_top_genes = length(bulkG),
                       sn_eligible_genes = length(snG), shared = length(intersect(bulkG, snG)),
                       shared_genes = paste(intersect(bulkG, snG), collapse = ","))
        }))
    }))
    if (!is.null(bvs)) { writeTSV(bvs, file.path(out, "bulk_vs_sn.tsv")); print(bvs) }
}
say("Fig. 7 done: ", out)
