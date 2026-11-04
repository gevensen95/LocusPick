## =============================================================================
## Fig. 6. Single-cell multiome: what changes when loci are picked per cell type?
##
## Data: PBMC 10x multiome, GSE215915 (McClain et al. 2024; COVID-19,
## decompensated vs stable), processed by 00_process_pbmc.R into one Seurat +
## Signac object with Azimuth cell types and LinkPeaks() peak-gene links.
##
##   A  Peak-gene links vs nearest gene: scLocusPick() with links = "auto" vs
##      links = NULL; overlap of the top loci and distance of picks to the gene.
##   B  Pseudobulk vs cell-level DE under a null: donor labels are permuted
##      (no true condition effect) and genes at padj < 0.05 are counted with
##      pseudobulkDE() and with Seurat FindMarkers() (each cell a replicate).
##   C  Paired vs unpaired: cells are split at random into an "RNA half" and an
##      "ATAC half", run as unpaired data, and compared with the paired run.
##   D  Donors: top loci with fewer donors per condition vs all donors.
##
## Outputs (results/fig6/): links_vs_nearest.tsv, null_de_counts.tsv,
## paired_vs_unpaired.tsv, donor_downsampling.tsv, fig6.pdf/png, top_loci.pdf
## =============================================================================

source("config.R"); source("R/common.R")
suppressPackageStartupMessages({ library(Seurat); library(Signac) })
out <- figDir("fig6")
set.seed(SEED)

CT  <- "predicted.celltype.l2"   # Azimuth level-2 labels (00_process_pbmc.R)
DON <- "donor"
COND <- PBMC$condition           # c(column, test, reference)
N_PERM  <- if (TEST) 3 else 20
N_DOWN  <- if (TEST) 2 else 5
DOWN_K  <- if (TEST) c(2, 3) else c(3, 5)
MAX_CT_NULL <- 4                 # largest cell types used for the null DE test

so <- readRDS(needFiles(PBMC$object, "PBMC multiome object (run 00_process_pbmc.R)"))
atacAssay <- names(Filter(function(a) inherits(a, "ChromatinAssay"), so@assays))[1]
genes <- if (!is.null(Annotation(so[[atacAssay]]))) readGeneModels(Annotation(so[[atacAssay]])) else
    if (TEST) locusPickSingleCellExample()$genes else readGeneModels(REF$hg38$gtf)
md <- so@meta.data
donors <- unique(md[, c(DON, COND[1])])
say(nrow(donors), " donors: ", paste(names(table(donors[[COND[1]]])), table(donors[[COND[1]]]),
                                     sep = " = ", collapse = ", "))
cellTypes <- if (!is.null(PBMC$cellTypes)) PBMC$cellTypes else {
    tab <- table(md[[CT]], md[[DON]])
    rownames(tab)[rowSums(tab >= 10) >= 2 * 2]          # >= 10 cells in >= 4 donors
}
covs <- if (TEST) "genotype" else NULL                    # the toy data have a second factor

run <- function(obj, atac = obj, links = "auto", donorsKeep = NULL, ...) {
    if (!is.null(donorsKeep)) {
        obj  <- subset(obj, cells = colnames(obj)[obj@meta.data[[DON]] %in% donorsKeep])
        if (!identical(atac, NULL)) atac <- subset(atac, cells = colnames(atac)[atac@meta.data[[DON]] %in% donorsKeep])
    }
    suppressMessages(scLocusPick(obj, atac = atac, cellType = CT, sample = DON, condition = COND,
                                 genes = genes, links = links, cellTypes = cellTypes,
                                 atacAssay = atacAssay, covariates = covs, n = N_TOP, ...))
}

## ---- A. Links vs nearest gene ------------------------------------------------------------
hasLinks <- length(Links(so[[atacAssay]])) > 0
if (!hasLinks) message("No LinkPeaks() links in the object: panel A skipped.")
ref <- run(so, links = if (hasLinks) "auto" else NULL)
nearest <- run(so, links = NULL)
la <- NULL
if (hasLinks) {
    la <- topOverlap(list(links = ref$candidates, nearest = nearest$candidates), n = N_TOP)
    la <- la[la$ranking == "nearest", ]
    dd <- rbind(data.frame(method = "links", distance = ref$top$distance),
                data.frame(method = "nearest gene", distance = nearest$top$distance))
    writeTSV(la, file.path(out, "links_vs_nearest.tsv"))
    writeTSV(dd, file.path(out, "pick_distances.tsv"))
}

## ---- B. Null DE: permuted donor labels ------------------------------------------------------
DefaultAssay(so) <- "RNA"
if (!"data" %in% Layers(so[["RNA"]]) || all(LayerData(so, "data")@x == LayerData(so, "counts")@x))
    so <- NormalizeData(so, verbose = FALSE)
ctNull <- utils::head(names(sort(table(md[[CT]][md[[CT]] %in% cellTypes]), decreasing = TRUE)), MAX_CT_NULL)
null <- list()
for (b in seq_len(N_PERM)) {
    perm <- stats::setNames(sample(donors[[COND[1]]]), donors[[DON]])
    so$perm <- unname(perm[so@meta.data[[DON]]])
    pb <- suppressMessages(pseudobulkDE(so, contrast = c("perm", COND[2], COND[3]), sample = DON,
                                        cellType = CT, assay = "RNA", cellTypes = ctNull,
                                        covariates = covs))
    for (ct in ctNull) {
        sub <- subset(so, cells = colnames(so)[so@meta.data[[CT]] == ct])
        Idents(sub) <- sub$perm
        fm <- tryCatch(suppressWarnings(FindMarkers(sub, ident.1 = COND[2], ident.2 = COND[3],
                                                    assay = "RNA", logfc.threshold = 0, min.pct = 0.1,
                                                    verbose = FALSE)),
                       error = function(e) NULL)
        null[[length(null) + 1]] <- data.frame(
            permutation = b, cell_type = ct,
            pseudobulk = if (!is.null(pb[[ct]])) sum(pb[[ct]]$padj < 0.05, na.rm = TRUE) else NA,
            FindMarkers = if (!is.null(fm)) sum(fm$p_val_adj < 0.05, na.rm = TRUE) else NA)
    }
    say("permutation ", b, "/", N_PERM)
}
null <- do.call(rbind, null)
writeTSV(null, file.path(out, "null_de_counts.tsv"))

## ---- C. Paired vs unpaired ---------------------------------------------------------------------
half <- sample(c(TRUE, FALSE), ncol(so), replace = TRUE)
rnaHalf  <- subset(so, cells = colnames(so)[half])
atacHalf <- subset(so, cells = colnames(so)[!half])
unp <- run(rnaHalf, atac = atacHalf, links = NULL)
pu <- topOverlap(list(paired_nearest = nearest$candidates, unpaired = unp$candidates), n = N_TOP)
pu <- pu[pu$ranking == "unpaired", ]
writeTSV(pu, file.path(out, "paired_vs_unpaired.tsv"))

## ---- D. Donor downsampling ------------------------------------------------------------------------
down <- list()
byCond <- split(donors[[DON]], donors[[COND[1]]])
for (k in DOWN_K[DOWN_K <= min(lengths(byCond))]) for (r in seq_len(N_DOWN)) {
    keep <- unlist(lapply(byCond, sample, k))
    rr <- tryCatch(run(so, links = if (hasLinks) "auto" else NULL, donorsKeep = keep),
                   error = function(e) { message("k = ", k, ": ", conditionMessage(e)); NULL })
    if (is.null(rr)) next
    o <- topOverlap(list(all = ref$candidates, sub = rr$candidates), n = N_TOP)
    o <- o[o$ranking == "sub", ]
    if (nrow(o)) down[[length(down) + 1]] <- cbind(donors_per_condition = k, draw = r, o)
}
down <- do.call(rbind, down)
if (!is.null(down)) writeTSV(down, file.path(out, "donor_downsampling.tsv"))

## ---- Figure ---------------------------------------------------------------------------------------------
panels <- list()
if (!is.null(la)) panels$A <- ggplot(dd, aes(method, pmax(distance, 1))) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.5) +
    geom_jitter(width = 0.1, size = 0.8) + scale_y_log10() +
    labs(x = NULL, y = "Peak-to-gene distance of top picks (bp)",
         title = sprintf("A  Links vs nearest (median Jaccard %.2f)", stats::median(la$jaccard))) + theme_paper()
nl <- rbind(data.frame(method = "pseudobulk", n = null$pseudobulk),
            data.frame(method = "FindMarkers (cells)", n = null$FindMarkers))
panels$B <- ggplot(nl, aes(method, n + 1)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.5) +
    geom_jitter(width = 0.1, size = 0.8) + scale_y_log10() +
    labs(x = NULL, y = "Genes at padj < 0.05 + 1\n(permuted donor labels)", title = "B  DE under the null") +
    theme_paper()
panels$C <- ggplot(pu, aes(group, jaccard)) + geom_col(fill = PAL[["LocusPick"]], width = 0.6) +
    scale_y_continuous(limits = c(0, 1)) + coord_flip() +
    labs(x = NULL, y = "Jaccard, unpaired vs paired (top loci)", title = "C  Paired vs unpaired") + theme_paper()
if (!is.null(down)) panels$D <- ggplot(down, aes(factor(donors_per_condition), jaccard)) +
    geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.5) + geom_jitter(width = 0.1, size = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Donors per condition", y = "Jaccard vs all donors", title = "D  Fewer donors") + theme_paper()
p <- wrap_plots(panels, ncol = 2) +
    plot_annotation(caption = if (TEST) "TEST MODE: simulated data, numbers are meaningless")
saveFig(p, file.path(out, "fig6"), width = 9, height = 7)
reportLoci(ref, file.path(out, "top_loci.pdf"))
say("Fig. 6 done: ", out)
