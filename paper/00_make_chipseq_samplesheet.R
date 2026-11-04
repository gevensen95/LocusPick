## Build the nf-core/chipseq samplesheet from the GSE60104 design table and the
## nf-core/fetchngs samplesheet. Called by 00_process_gse60104.sh.
## Usage: Rscript 00_make_chipseq_samplesheet.R design.csv fetchngs_samplesheet.csv out.csv
a <- commandArgs(TRUE)
design <- utils::read.csv(a[1], stringsAsFactors = FALSE)
fq <- utils::read.csv(a[2], stringsAsFactors = FALSE)
runs <- strsplit(design$srr, ";")
rows <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    f <- fq[fq$run_accession %in% trimws(runs[[i]]), ]
    if (!nrow(f)) stop("No FASTQ for ", design$gsm[i], " (", design$srr[i], ") in ", a[2], call. = FALSE)
    inp <- as.logical(design$is_input[i])
    data.frame(sample = design$sample[i], fastq_1 = f$fastq_1,
               fastq_2 = ifelse(is.na(f$fastq_2), "", f$fastq_2),
               replicate = design$replicate[i],
               antibody = if (inp) "" else design$antibody[i],
               control = if (inp) "" else design$control[i],
               control_replicate = if (inp || is.na(design$control_replicate[i])) "" else design$control_replicate[i])
}))
utils::write.csv(rows, a[3], row.names = FALSE, quote = FALSE, na = "")
message("Wrote ", nrow(rows), " rows to ", a[3])
