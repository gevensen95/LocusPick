#' Find signal and peak files in a pipeline output folder
#'
#' Scans a results folder (for example from nf-core/cutandrun, nf-core/chipseq
#' or nf-core/atacseq, or any folder of MACS2/SEACR output) and lists the
#' bigWig, bedGraph, peak and consensus-peak files it contains, with the sample
#' name parsed from each file name. Detection uses file name patterns rather
#' than fixed paths, so it works across pipeline versions.
#'
#' Control samples (IgG, input) are flagged so they are not used as signal or
#' peaks by mistake.
#'
#' @param dir Results folder.
#' @param recursive Search subfolders.
#' @return A `data.frame` with `file`, `type` (`bigwig`, `bedgraph`, `peaks`,
#'   `consensus`), `caller` (`macs2`, `seacr`, `genrich`, `unknown`),
#'   `sample`, and `control`.
#' @seealso [signalFromFiles()] to turn the bigWigs into averaged condition
#'   tracks, and [readPeaks()] for the peak files.
#' @examples
#' d <- file.path(tempdir(), "results")
#' dir.create(file.path(d, "peaks"), recursive = TRUE, showWarnings = FALSE)
#' file.create(file.path(d, "peaks", c("Y1_peaks.narrowPeak", "IgG_peaks.narrowPeak")))
#' findPipelineOutputs(d)
#' @export
findPipelineOutputs <- function(dir, recursive = TRUE) {
    if (!dir.exists(dir)) stop("Folder not found: ", dir, call. = FALSE)
    f <- list.files(dir, recursive = recursive, full.names = TRUE)
    f <- f[!grepl("/(work|\\.nextflow|pipeline_info)/", f)]
    b <- basename(f)
    type <- ifelse(grepl("\\.(bigwig|bw)$", b, ignore.case = TRUE), "bigwig",
            ifelse(grepl("\\.(bedgraph|bdg)(\\.gz)?$", b, ignore.case = TRUE), "bedgraph",
            ifelse(grepl("consensus[._]peaks?\\.(bed|txt)$", b, ignore.case = TRUE), "consensus",
            ifelse(grepl("(_peaks\\.(narrowPeak|broadPeak|gappedPeak|xls)|\\.(stringent|relaxed)\\.bed|\\.peaks.*\\.bed|\\.narrowPeak|\\.broadPeak)(\\.gz)?$",
                         b, ignore.case = TRUE), "peaks", NA))))
    keep <- !is.na(type)
    f <- f[keep]; b <- b[keep]; type <- type[keep]
    if (!length(f)) {
        warning("No signal or peak files found in ", dir, call. = FALSE)
        return(data.frame(file = character(), type = character(), caller = character(),
                          sample = character(), control = logical()))
    }
    path <- tolower(f)
    caller <- ifelse(grepl("seacr|stringent|relaxed", path), "seacr",
              ifelse(grepl("macs", path) | grepl("_peaks\\.(narrowpeak|broadpeak|xls)", path), "macs2",
              ifelse(grepl("genrich", path), "genrich", "unknown")))
    caller[type %in% c("bigwig", "bedgraph")] <- NA
    sample <- sub("(\\.(bigwig|bw|bedgraph|bdg)(\\.gz)?)$", "", b, ignore.case = TRUE)
    sample <- sub("(_peaks\\.(narrowPeak|broadPeak|gappedPeak|xls)|\\.(stringent|relaxed)\\.bed|\\.peaks.*\\.bed|\\.consensus[._]peaks?\\.(bed|txt)|\\.(narrowPeak|broadPeak))(\\.gz)?$",
                  "", sample, ignore.case = TRUE)
    sample <- sub("(\\.(sorted|dedup|filtered|mLb\\.clN|mRp\\.clN|bowtie2|target|markdup|norm|cpm|bpm|seacr|macs2|macs3|genrich))+$",
                  "", sample, ignore.case = TRUE)
    control <- grepl("igg|input|control|^ctrl", sample, ignore.case = TRUE)
    out <- data.frame(file = f, type = type, caller = caller, sample = sample,
                      control = control, stringsAsFactors = FALSE)
    out[order(out$type, out$sample), ]
}

#' Build condition signal tracks from replicate files
#'
#' Groups signal files (bigWig or BAM) by condition using patterns matched
#' against file or sample names, and returns one averaged source per
#' condition, ready for the `signal` argument of [buildCandidates()] and
#' [plotLocus()].
#'
#' @param files Character vector of files, or the output of
#'   [findPipelineOutputs()] (bigWig rows are used; controls are dropped).
#' @param conditions Named character vector of regular expressions, one per
#'   condition, e.g. `c(young = "^Y", old = "^O")`.
#' @param ... Passed to [bamSource()] for BAM files.
#' @return A named list of signal sources.
#' @examples
#' bw <- system.file("extdata", c("toy_A.bw", "toy_B.bw"), package = "LocusPick")
#' signalFromFiles(bw, c(A = "toy_A", B = "toy_B"))
#' @export
signalFromFiles <- function(files, conditions, ...) {
    if (is.data.frame(files)) {
        files <- files[files$type == "bigwig" & !files$control, , drop = FALSE]
        nm <- files$sample; files <- files$file
    } else {
        nm <- basename(files)
    }
    if (is.null(names(conditions)) || any(names(conditions) == ""))
        stop("conditions must be a named vector of patterns.", call. = FALSE)
    out <- lapply(names(conditions), function(cn) {
        hit <- grepl(conditions[[cn]], nm) | grepl(conditions[[cn]], basename(files))
        if (!any(hit)) stop("No files match condition '", cn, "' (pattern '",
                            conditions[[cn]], "').", call. = FALSE)
        message(cn, ": ", paste(basename(files[hit]), collapse = ", "))
        .fileSource(files[hit], ...)
    })
    stats::setNames(out, names(conditions))
}

.fileSource <- function(files, ...) {
    srcs <- lapply(files, function(f) {
        if (grepl("\\.bam$", f, ignore.case = TRUE)) bamSource(f, ...) else bigWigSource(f)
    })
    if (length(srcs) == 1) srcs[[1]] else meanSource(srcs)
}
