#!/usr/bin/env bash
# One-time: build the mm10 Umap multi-read mappability bigWig used by cnr_example_tracks.R
# Umap distributes mm10 multi-read mappability only as a bedGraph, so download it
# and convert to an indexed bigWig (fast region lookups from R).
#
# Usage:  bash make_umap_bw.sh [k]      k = 24, 36, 50 (default) or 100
# Needs:  bedGraphToBigWig  ->  conda install -c bioconda ucsc-bedgraphtobigwig
# Disk:   the uncompressed bedGraph is several GB while converting (deleted after)

set -euo pipefail
K=${1:-50}
OUT=/vast/igc/analyses/garrett/P01/ref/mm10
mkdir -p "$OUT"
cd "$OUT"

command -v bedGraphToBigWig >/dev/null 2>&1 || {
  echo "bedGraphToBigWig not found. Install with: conda install -c bioconda ucsc-bedgraphtobigwig" >&2
  exit 1
}

[ -s mm10.chrom.sizes ] || wget -q https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
[ -s "k${K}.umap.bedgraph.gz" ] || wget "https://bismap.hoffmanlab.org/raw/mm10/k${K}.umap.bedgraph.gz"

# keep only data lines (numeric start) on chromosomes in the sizes file, then sort
zcat "k${K}.umap.bedgraph.gz" \
  | awk 'NR==FNR { ok[$1]=1; next } ($1 in ok) && $2 ~ /^[0-9]+$/ && NF==4' mm10.chrom.sizes - \
  | LC_ALL=C sort -k1,1 -k2,2n -S 4G -T "$OUT" \
  > "k${K}.umap.sorted.bedgraph"

[ -s "k${K}.umap.sorted.bedgraph" ] || { echo "Sorted bedGraph is empty: check the download and filter" >&2; exit 1; }

# write to a temp name and only rename on success, so a killed or failed run
# never leaves a broken file at the final path
BW="k${K}.Umap.MultiTrackMappability.bw"
rm -f "$BW" "$BW.tmp"
bedGraphToBigWig "k${K}.umap.sorted.bedgraph" mm10.chrom.sizes "$BW.tmp"
if command -v bigWigInfo >/dev/null 2>&1; then
  bigWigInfo "$BW.tmp" >/dev/null || { echo "bigWig failed validation" >&2; exit 1; }
fi
mv "$BW.tmp" "$BW"
rm "k${K}.umap.sorted.bedgraph"

ls -lh "$OUT/k${K}.Umap.MultiTrackMappability.bw"
echo "Done. Delete k${K}.umap.bedgraph.gz if you want the space back."
