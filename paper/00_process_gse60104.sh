#!/usr/bin/env bash
## =============================================================================
## Fig. 4 data: GSE60104 (Orlando et al. 2014, ChIP-Rx) through nf-core.
##   1. Fill templates/gse60104_design.csv: one row per sequencing run group
##      (gsm, srr, sample, replicate, condition, antibody, is_input, control,
##      control_replicate). SRR IDs: SRA Run Selector for GSE60104, or
##      `pysradb gsm-to-srr GSM...`. nf-core/fetchngs no longer accepts GEO IDs.
##   2. bash 00_process_gse60104.sh   (runs fetchngs, then chipseq, with Singularity)
## Needs: nextflow, singularity; ~200 GB scratch for the combined genome index.
## =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
D=$PWD/data/gse60104
DESIGN=$PWD/templates/gse60104_design.csv
GENOME=$PWD/data/ref/combined
mkdir -p $D
[ "$(wc -l < $DESIGN)" -gt 1 ] || { echo "Fill in $DESIGN first" >&2; exit 1; }
[ -s $GENOME/hg38_dm6.fa ] || { echo "Run 00_download.sh first (combined hg38 + dm6 genome)" >&2; exit 1; }

# 1. FASTQs: one SRR per line (the srr column may hold several runs separated by ;)
Rscript -e 'd <- read.csv(commandArgs(TRUE)[1]); writeLines(unique(trimws(unlist(strsplit(d$srr, ";")))), commandArgs(TRUE)[2])' \
  $DESIGN $D/ids.csv
nextflow run nf-core/fetchngs -r 1.12.0 -profile singularity \
  --input $D/ids.csv --outdir $D/fetchngs -resume

# 2. chipseq samplesheet from the design + fetchngs output
Rscript 00_make_chipseq_samplesheet.R $DESIGN $D/fetchngs/samplesheet/samplesheet.csv $D/chipseq_samplesheet.csv

# 3. nf-core/chipseq on hg38 + dm6. Broad peaks (H3K79me2): --narrow_peak false (default).
#    Multi-mapping reads between human and fly are removed by the pipeline's MAPQ filter.
nextflow run nf-core/chipseq -r 2.1.0 -profile singularity \
  --input $D/chipseq_samplesheet.csv --outdir $D/chipseq \
  --fasta $GENOME/hg38_dm6.fa --gtf $GENOME/genes.gtf \
  --blacklist $PWD/data/ref/hg38/hg38-blacklist.v2.bed.gz \
  --macs_gsize 2.7e9 -resume

# 4. Spike-in counts: idxstats of the filtered merged-library BAMs (made here if the
#    pipeline version did not write them).
for bam in $D/chipseq/bwa/mergedLibrary/*.mLb.clN.sorted.bam; do
  out=$D/chipseq/bwa/mergedLibrary/samtools_stats/$(basename $bam).idxstats
  [ -s "$out" ] || { mkdir -p "$(dirname $out)"; samtools idxstats "$bam" > "$out"; }
done
echo "Done. Next: Rscript 03_fig4_global_shift.R"
