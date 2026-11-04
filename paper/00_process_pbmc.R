## =============================================================================
## Fig. 6 data: GSE215915 PBMC 10x multiome -> one Seurat + Signac object.
##
## 1. Download the per-sample Cell Ranger ARC outputs from GEO (GSE215915:
##    filtered_feature_bc_matrix.h5 and atac_fragments.tsv.gz + .tbi) and fill
##    templates/gse215915_samples.csv:
##      donor      - sample / donor id (one row per 10x library)
##      severity   - "decompensated" or "stable" (from the GEO sample records)
##      h5, fragments - paths (relative to paper/ or absolute)
##    Check first that GEO really holds fragment files for every sample; if it
##    only has count matrices, Fig. 6 panel A/C cannot be made from GEO alone.
## 2. Rscript 00_process_pbmc.R
##
## Steps: per-sample objects -> common peak set (union of Cell Ranger peaks,
## 20 bp-10 kb) -> merge -> QC (Signac multiome thresholds) -> Azimuth pbmcref
## labels -> peak-gene links (LinkPeaks) for genes differentially expressed in
## any cell type plus 2,000 variable genes -> data/gse215915/pbmc_multiome.rds
##
## Needs: Seurat >= 5, Signac, Azimuth, EnsDb.Hsapiens.v86,
## BSgenome.Hsapiens.UCSC.hg38, LocusPick. Memory: ~64 GB for ~100k cells.
## =============================================================================

source("config.R")
if (TEST) stop("Test mode uses the simulated object; nothing to process.", call. = FALSE)
suppressPackageStartupMessages({
    library(Seurat); library(Signac); library(GenomicRanges); library(LocusPick)
})
for (p in c("Azimuth", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, " first.", call. = FALSE)
set.seed(SEED)
options(future.globals.maxSize = 50 * 1024^3)

ss <- utils::read.csv(PBMC$samples, stringsAsFactors = FALSE)
if (!nrow(ss)) stop("Fill in ", PBMC$samples, " first.", call. = FALSE)
ss$h5 <- path.expand(ss$h5); ss$fragments <- path.expand(ss$fragments)
miss <- c(ss$h5, ss$fragments, paste0(ss$fragments, ".tbi"))
miss <- miss[!file.exists(miss)]
if (length(miss)) stop("Missing files:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)
stopifnot(all(ss$severity %in% PBMC$condition[2:3]), !anyDuplicated(ss$donor))
print(table(ss$severity))

ann <- suppressWarnings(GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86))
seqlevelsStyle(ann) <- "UCSC"
std <- paste0("chr", c(1:22, "X", "Y"))

## ---- 1. Read each sample; union of peaks -------------------------------------------------
raw <- lapply(seq_len(nrow(ss)), function(i) {
    message("reading ", ss$donor[i])
    x <- Read10X_h5(ss$h5[i])
    pk <- StringToGRanges(rownames(x$Peaks), sep = c(":", "-"))
    list(rna = x$`Gene Expression`, peaks = pk[as.character(seqnames(pk)) %in% std])
})
allPk <- reduce(do.call(c, lapply(raw, `[[`, "peaks")))
allPk <- allPk[width(allPk) >= 20 & width(allPk) <= 10000]
allPk <- keepStandardChromosomes(allPk, pruning.mode = "coarse")
blk <- if (file.exists(REF$hg38$blacklist)) rtracklayer::import(REF$hg38$blacklist) else NULL
if (!is.null(blk)) allPk <- allPk[!overlapsAny(allPk, blk)]
message(length(allPk), " peaks in the common set")

objs <- lapply(seq_len(nrow(ss)), function(i) {
    cells <- colnames(raw[[i]]$rna)
    fr <- CreateFragmentObject(ss$fragments[i], cells = cells)
    cnt <- FeatureMatrix(fragments = fr, features = allPk, cells = cells)
    o <- CreateSeuratObject(raw[[i]]$rna, assay = "RNA", project = ss$donor[i])
    o[["ATAC"]] <- CreateChromatinAssay(cnt, fragments = fr, annotation = ann, sep = c("-", "-"))
    o$donor <- ss$donor[i]; o$severity <- ss$severity[i]
    RenameCells(o, add.cell.id = ss$donor[i])
})
rm(raw); gc()
so <- merge(objs[[1]], objs[-1])
rm(objs); gc()
so <- JoinLayers(so, assay = "RNA")

## ---- 2. QC ------------------------------------------------------------------------------------
DefaultAssay(so) <- "ATAC"
so <- NucleosomeSignal(so)
so <- TSSEnrichment(so)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-", assay = "RNA")
n0 <- ncol(so)
so <- subset(so, nCount_RNA > 1000 & nCount_RNA < 25000 & nCount_ATAC > 1000 & nCount_ATAC < 100000 &
                 percent.mt < 20 & nucleosome_signal < 2 & TSS.enrichment > 1)
message("QC kept ", ncol(so), " of ", n0, " cells")

## ---- 3. Cell types (Azimuth pbmcref) ------------------------------------------------------------
DefaultAssay(so) <- "RNA"
so <- Azimuth::RunAzimuth(so, reference = "pbmcref", assay = "RNA")
so <- subset(so, predicted.celltype.l2.score >= 0.5)
print(table(so$predicted.celltype.l2, so$severity))

## ---- 4. Peak-gene links -------------------------------------------------------------------------
so <- NormalizeData(so, assay = "RNA", verbose = FALSE)
so <- FindVariableFeatures(so, assay = "RNA", nfeatures = 2000, verbose = FALSE)
de <- pseudobulkDE(so, contrast = PBMC$condition, sample = "donor", cellType = "predicted.celltype.l2",
                   assay = "RNA")
deGenes <- unique(unlist(lapply(de, function(d) d$gene[!is.na(d$padj) & d$padj < 0.05])))
genesUse <- intersect(unique(c(deGenes, VariableFeatures(so, assay = "RNA"))), rownames(so[["RNA"]]))
message("LinkPeaks for ", length(genesUse), " genes (", length(deGenes), " DE in some cell type)")
DefaultAssay(so) <- "ATAC"
so <- RegionStats(so, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
so <- LinkPeaks(so, peak.assay = "ATAC", expression.assay = "RNA", genes.use = genesUse)
message(length(Links(so[["ATAC"]])), " links")

dir.create(dirname(PBMC$object), recursive = TRUE, showWarnings = FALSE)
saveRDS(so, PBMC$object)
message("Saved ", PBMC$object)
