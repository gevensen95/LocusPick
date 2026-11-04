## =============================================================================
## Select ENCODE files for Figs. 3 and 5 and write data/encode/manifest.tsv.
##
## Rule (fixed before looking at any result):
##   * cell lines HepG2 and K562 (config.R: ENCODE$cells), assembly GRCh38,
##     status released, no ERROR audit on the experiment
##   * TF ChIP-seq targets profiled in BOTH lines, taken in the order of
##     TF_PRIORITY below, up to N_TF factors
##   * per experiment: IDR thresholded peaks (bed narrowPeak), the pooled
##     fold-change-over-control bigWig and the two per-replicate fold-change
##     bigWigs; when a file type exists more than once (several pipeline
##     runs), the newest released file is used
##   * ATAC-seq in both lines: same file types (peaks: IDR thresholded peaks,
##     else pseudoreplicated peaks)
##   * total RNA-seq in both lines: RSEM gene quantifications, two biological
##     replicates, GENCODE V29
## The script records every choice (experiment and file accessions) in the
## manifest so the Methods can list them.
##
## Needs internet access to www.encodeproject.org. Run once, then
## 00_download.sh to fetch the files.
## =============================================================================

source("config.R")
if (TEST) stop("00_select_encode.R needs the ENCODE portal; test mode uses simulated files.", call. = FALSE)
suppressPackageStartupMessages(library(jsonlite))

N_TF <- 8
TF_PRIORITY <- c("CTCF", "YY1", "MAX", "REST", "NRF1", "JUND", "CEBPB", "ZBTB33", "GABPA", "ELF1",
                 "SP1", "USF1", "USF2", "RAD21", "EGR1", "ATF3", "SRF", "TCF12", "ZNF143", "MAZ")
PORTAL <- "https://www.encodeproject.org"

## Search wrapper: repeated keys (field=...) go in as a character vector
search <- function(params) {
    qs <- paste(mapply(function(k, v) paste0(k, "=", utils::URLencode(v, reserved = TRUE)),
                       rep(names(params), lengths(params)), unlist(params)), collapse = "&")
    url <- paste0(PORTAL, "/search/?", qs, "&format=json&limit=all")
    r <- tryCatch(fromJSON(url, simplifyVector = TRUE, flatten = TRUE), error = function(e) e)
    if (inherits(r, "error")) {
        if (grepl("404", conditionMessage(r))) return(data.frame())   # ENCODE returns 404 for no hits
        stop("ENCODE search failed: ", url, "\n", conditionMessage(r), call. = FALSE)
    }
    r[["@graph"]]
}

experiments <- function(assay, cell) {
    e <- search(list(type = "Experiment", assay_title = assay, status = "released",
                     biosample_ontology.term_name = cell, assembly = "GRCh38",
                     field = c("accession", "target.label", "audit.ERROR",
                               "replicates.biological_replicate_number")))
    if (!nrow(e)) return(e)
    err <- if ("audit.ERROR" %in% names(e)) lengths(e$audit.ERROR) > 0 else rep(FALSE, nrow(e))
    e$n_bio_reps <- vapply(e$replicates, function(r)
        if (is.data.frame(r)) length(unique(r$biological_replicate_number)) else 0L, 0L)
    e <- e[!err, ]
    e$target <- if ("target.label" %in% names(e)) e$target.label else NA
    e[, c("accession", "target", "n_bio_reps")]
}

files <- function(exp) {
    f <- search(list(type = "File", dataset = paste0("/experiments/", exp, "/"), status = "released",
                     frame = "object"))
    if (!nrow(f)) return(f)
    f$n_bio <- lengths(f$biological_replicates)
    f$rep <- vapply(f$biological_replicates, function(x) if (length(x) == 1) as.integer(x) else NA_integer_, 1L)
    f[is.na(f$assembly) | f$assembly == "GRCh38", ]
}

newest <- function(f) if (nrow(f)) f[order(f$date_created, decreasing = TRUE), ][1, ] else f

pickFiles <- function(exp, target, cell, peakTypes) {
    f <- files(exp)
    pk <- NULL
    for (pt in peakTypes) {
        pk <- newest(f[f$output_type == pt & f$file_format == "bed", ])
        if (nrow(pk)) break
    }
    fc <- f[f$output_type == "fold change over control" & f$file_format == "bigWig", ]
    pooled <- newest(fc[fc$n_bio >= 2, ])
    reps <- do.call(rbind, lapply(1:2, function(r) {
        x <- fc[!is.na(fc$rep), ]
        x <- x[x$rep == sort(unique(x$rep))[r], ]
        newest(x)
    }))
    if (!nrow(pk) || !nrow(pooled) || is.null(reps) || nrow(reps) < 2) return(NULL)
    row <- function(x, output, replicate)
        data.frame(target = target, cell = cell, output = output, replicate = replicate,
                   accession = x$accession, experiment = exp, assembly = "GRCh38",
                   url = paste0(PORTAL, x$href), md5sum = x$md5sum,
                   file = file.path(ENCODE$dir, basename(x$href)))
    rbind(row(pk, "peaks", NA), row(pooled, "fc_pooled", NA),
          row(reps[1, ], "fc_rep", 1), row(reps[2, ], "fc_rep", 2))
}

cells <- unname(ENCODE$cells)
man <- list()

## ---- TF ChIP-seq ----------------------------------------------------------------------
tf <- lapply(cells, function(cl) experiments("TF ChIP-seq", cl))
names(tf) <- cells
both <- intersect(tf[[1]]$target, tf[[2]]$target)
cand <- c(intersect(TF_PRIORITY, both), setdiff(sort(both), TF_PRIORITY))
message(length(both), " TFs profiled in both lines; trying in priority order")
chosen <- character()
for (t in cand) {
    if (length(chosen) >= N_TF) break
    rows <- lapply(cells, function(cl) {
        ex <- tf[[cl]][tf[[cl]]$target == t & tf[[cl]]$n_bio_reps >= 2, ]
        for (e in ex$accession) { r <- pickFiles(e, t, cl, c("IDR thresholded peaks",
                                                             "conservative IDR thresholded peaks"))
                                  if (!is.null(r)) return(r) }
        NULL
    })
    if (any(vapply(rows, is.null, TRUE))) { message("  ", t, ": files incomplete, skipped"); next }
    man <- c(man, rows); chosen <- c(chosen, t)
    message("  ", t, ": ", paste(vapply(rows, function(r) r$experiment[1], ""), collapse = " / "))
}
if (length(chosen) < N_TF) warning("Only ", length(chosen), " factors had complete files.", call. = FALSE)

## ---- ATAC-seq ---------------------------------------------------------------------------
for (cl in cells) {
    ex <- experiments("ATAC-seq", cl)
    ex <- ex[ex$n_bio_reps >= 2, ]
    r <- NULL
    for (e in ex$accession) {
        r <- pickFiles(e, "ATAC", cl, c("IDR thresholded peaks", "conservative IDR thresholded peaks",
                                        "pseudoreplicated peaks", "replicated peaks"))
        if (!is.null(r)) break
    }
    if (is.null(r)) stop("No complete ATAC-seq experiment for ", cl, call. = FALSE)
    man <- c(man, list(r))
}

## ---- RNA-seq ------------------------------------------------------------------------------
for (cl in cells) {
    ex <- experiments("total RNA-seq", cl)
    ex <- ex[ex$n_bio_reps >= 2, ]
    got <- NULL
    for (e in ex$accession) {
        f <- files(e)
        q <- f[f$output_type == "gene quantifications" & f$file_format == "tsv" & !is.na(f$rep), ]
        if ("genome_annotation" %in% names(q)) q <- q[is.na(q$genome_annotation) | q$genome_annotation == "V29", ]
        q <- do.call(rbind, lapply(split(q, q$rep), newest))
        if (!is.null(q) && nrow(q) >= 2) { got <- q[1:2, ]; exp <- e; break }
    }
    if (is.null(got)) stop("No total RNA-seq gene quantifications with 2 replicates for ", cl, call. = FALSE)
    man <- c(man, list(data.frame(target = "RNA", cell = cl, output = "rna_quant", replicate = 1:2,
                                  accession = got$accession, experiment = exp, assembly = "GRCh38",
                                  url = paste0(PORTAL, got$href), md5sum = got$md5sum,
                                  file = file.path(ENCODE$dir, basename(got$href)))))
}

man <- do.call(rbind, man)
dir.create(dirname(ENCODE$manifest), recursive = TRUE, showWarnings = FALSE)
utils::write.table(man, ENCODE$manifest, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(man), " files to ", ENCODE$manifest, "; factors: ", paste(chosen, collapse = ", "))
