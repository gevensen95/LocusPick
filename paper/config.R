## =============================================================================
## Paths and settings for the LocusPick manuscript analyses.
##
## Every figure script sources this file. Run scripts from the paper/ folder:
##     cd paper && Rscript 02_fig3_shortcuts.R
## Test mode (simulated data, minutes, no downloads):
##     LOCUSPICK_PAPER_TEST=1 Rscript 02_fig3_shortcuts.R
## =============================================================================

TEST <- identical(Sys.getenv("LOCUSPICK_PAPER_TEST"), "1")
PAPER <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(PAPER, "config.R")))
    stop("Run the paper scripts from the paper/ folder (cd paper).", call. = FALSE)

DATA    <- file.path(PAPER, if (TEST) "test_data" else "data")
RESULTS <- file.path(PAPER, if (TEST) "results_test" else "results")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

SEED    <- 20261103
N_TOP   <- 10          # loci per group reported as "top"
N_RULE  <- 20          # picks per group and rule in Fig. 3
THREADS <- as.integer(Sys.getenv("LOCUSPICK_THREADS", "8"))

## ---- Reference files (downloaded by 00_download.sh) -----------------------------
REF <- list(
    hg38 = list(gtf       = file.path(DATA, "ref/hg38/gencode.v29.primary_assembly.annotation.gtf.gz"),
                blacklist = file.path(DATA, "ref/hg38/hg38-blacklist.v2.bed.gz"),
                umap      = file.path(DATA, "ref/hg38/k50.Umap.MultiTrackMappability.bw")),
    mm10 = list(blacklist = file.path(DATA, "ref/mm10/mm10-blacklist.v2.bed.gz")))

## ---- Fig. 3 and 5: ENCODE HepG2 vs K562 -------------------------------------------
ENCODE <- list(
    manifest = file.path(DATA, "encode/manifest.tsv"),   # written by 00_select_encode.R
    dir      = file.path(DATA, "encode"),
    cells    = c(test = "HepG2", ref = "K562"))           # fold changes = HepG2 / K562

## ---- Fig. 4: ChIP-Rx global loss (GSE60104, Orlando et al. 2014) -----------------------
CHIPRX <- list(
    design   = file.path(PAPER, "templates/gse60104_design.csv"),  # fill in before running
    nfcore   = file.path(DATA, "gse60104/chipseq"),                # nf-core/chipseq --outdir
    spikePrefix = "dm6_",                                          # fly contig prefix in the combined genome
    test = "EPZ5676", ref = "DMSO")                                # condition labels in the design

## ---- Fig. 6: PBMC multiome (GSE215915, McClain et al. 2024) ---------------------------
PBMC <- list(
    samples  = file.path(PAPER, "templates/gse215915_samples.csv"),  # fill in before running
    object   = file.path(DATA, "gse215915/pbmc_multiome.rds"),        # written by 00_process_pbmc.R
    condition = c("severity", "decompensated", "stable"),             # column, test, reference
    cellTypes = NULL)                                                 # NULL = all with enough cells

## ---- Fig. 7 (and 4, 5): Hnf4a liver case study --------------------------------------
KAT   <- "/vast/igc/analyses/garrett/P01/kat"
HNF4A <- list(
    age_de     = file.path(KAT, "JvJ_RNAseq_Redos/Hepatocyte_YO/JvsJ_allGenes.txt"),
    age_counts = file.path(KAT, "JvJ_RNAseq_Redos/Hepatocyte_YO/NormCounts_JvJ.txt"),
    age_meta   = file.path(KAT, "JvJ_RNAseq_Redos/Hepatocyte_YO/meta.txt"),
    targets    = file.path(KAT, "Hnf4a_Epigenetics/YvO_HNF4a_KO_WL/CUTNRUN_Targets_Yng.tsv"),
    ko_expr    = file.path(KAT, "Hnf4a_Epigenetics/YvO_HNF4a_KO_WL/NormExp.tsv"),
    peaks      = file.path(KAT, "Hnf4a_Epigenetics/Hnf4a_CutnRun_Xue/ConsensusPeaks_YngAll_annot.bed"),
    gtf        = file.path(KAT, "genomes/mm10/mm10.gtf"),
    umap       = "/vast/igc/analyses/garrett/P01/ref/mm10/k50.Umap.MultiTrackMappability.bw",
    chromDicts = file.path(KAT, "2026_Epigenomics_ChromDicts", c(
        CNR_Yng     = "Hnf4a_CUTNRUN/results/perCondition_chromDicts/Yng.rds",
        CNR_Old     = "Hnf4a_CUTNRUN/results/perCondition_chromDicts/Old.rds",
        ATAC_Yng    = "ATAC/results/perCondition_chromDicts/Yng.rds",
        ATAC_Old    = "ATAC/results/perCondition_chromDicts/Old.rds",
        H3K27ac_Yng = "H3K27Ac/results/perCondition_chromDicts/yng-H3K27ac.rds",
        H3K27ac_Old = "H3K27Ac/results/perCondition_chromDicts/old-H3K27ac.rds")),
    ## PLACEHOLDER: per-replicate Hnf4a CUT&RUN bigWigs (needed for leave-one-
    ## replicate-out in Fig. 5). Replace with the real paths; scripts skip this
    ## part with a message while the files do not exist.
    cnr_reps = list(
        Yng = file.path(KAT, sprintf("Hnf4a_Epigenetics/CNR_bigwigs/PLACEHOLDER_Yng_rep%d.bw", 1:3)),
        Old = file.path(KAT, sprintf("Hnf4a_Epigenetics/CNR_bigwigs/PLACEHOLDER_Old_rep%d.bw", 1:3))),
    ## Genes shown in the original hand-picked Fig. 2F-G panels, e.g.
    ## c(activated_age_down = "Cyp2e1", repressed_age_up = "Saa1"); Fig. 7
    ## reports their rank and QC flags. NULL to skip.
    manual_picks = NULL,
    sn_rna  = "~/Desktop/Salk/P01/Single/snRNAseq/age.combined_final.rds",
    sn_atac = "~/Desktop/Salk/P01/Single/snATAC/atac_harmony.rds",
    ## column names in the sn objects (as in test_LocusPick_snRNA_snATAC.R)
    sn_cols = list(celltype_rna = "cell_type", celltype_atac = "cell_type",
                   sample_rna = "sample", sample_atac = "sample",
                   age_rna = "age", age_atac = "age", old = "old", young = "young",
                   rna_assay = "RNA"))

## ---- Test mode: point everything at simulated data ------------------------------------
if (TEST) source(file.path(PAPER, "R/simulate.R"))
