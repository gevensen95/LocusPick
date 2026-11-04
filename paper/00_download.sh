#!/usr/bin/env bash
## =============================================================================
## Download reference files and the ENCODE files listed in the manifest.
## Run from paper/:   bash 00_download.sh
## Needs: curl, md5sum (or md5 on macOS), gzip. Re-running skips finished files.
## =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
D=data
mkdir -p $D/ref/hg38 $D/ref/mm10 $D/ref/combined $D/encode

fetch() {  # url dest
  local url=$1 dest=$2
  [ -s "$dest" ] && return 0
  echo "  $url"
  curl -fL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
}
md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

echo "== hg38 / mm10 references"
fetch https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz \
      $D/ref/hg38/gencode.v29.primary_assembly.annotation.gtf.gz          # ENCODE GRCh38 uses GENCODE V29
fetch https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz $D/ref/hg38/hg38-blacklist.v2.bed.gz
fetch https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz $D/ref/mm10/mm10-blacklist.v2.bed.gz
fetch https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Umap.MultiTrackMappability.bw \
      $D/ref/hg38/k50.Umap.MultiTrackMappability.bw || \
  echo "!! hg38 Umap bigWig not found at UCSC: build it from bismap.hoffmanlab.org (see inst/scripts/make_umap_bw.sh)" >&2

echo "== hg38 + dm6 combined genome for the ChIP-Rx spike-in (Fig. 4)"
fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz $D/ref/combined/hg38.fa.gz
fetch https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz $D/ref/combined/dm6.fa.gz
if [ ! -s $D/ref/combined/hg38_dm6.fa ]; then
  # human: primary chromosomes only (chr1-22, X, Y, M); fly: all contigs, prefixed dm6_
  { zcat $D/ref/combined/hg38.fa.gz | awk '/^>/{keep = ($1 ~ /^>chr([0-9]+|X|Y|M)$/)} keep'
    zcat $D/ref/combined/dm6.fa.gz | sed 's/^>/>dm6_/'; } > $D/ref/combined/hg38_dm6.fa.part
  mv $D/ref/combined/hg38_dm6.fa.part $D/ref/combined/hg38_dm6.fa
fi
[ -s $D/ref/combined/genes.gtf ] || zcat $D/ref/hg38/gencode.v29.primary_assembly.annotation.gtf.gz > $D/ref/combined/genes.gtf

echo "== ENCODE files (manifest from 00_select_encode.R)"
M=$D/encode/manifest.tsv
if [ ! -s $M ]; then echo "!! $M missing: run Rscript 00_select_encode.R first" >&2; exit 1; fi
# columns by name: url, md5sum, file
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {print $c["url"] "\t" $c["md5sum"] "\t" $c["file"]}' $M |
while IFS=$'\t' read -r url sum dest; do
  fetch "$url" "$dest"
  got=$(md5 "$dest")
  if [ "$got" != "$sum" ]; then echo "!! md5 mismatch for $dest (expected $sum, got $got); deleting" >&2; rm -f "$dest"; exit 1; fi
done
echo "Done."
