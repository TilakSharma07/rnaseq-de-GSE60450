#!/usr/bin/env bash
# download_data.sh — fetch every input the pipeline needs. Idempotent.
#
# GSE60450 (Fu et al. 2015, Nat Cell Biol) — mouse mammary gland, basal and
# luminal cells at three developmental stages, 2 replicates each.
#
# The sample sheet is BUILT here rather than shipped, because the mapping from
# counts columns to GEO samples is the one place this dataset can go silently
# wrong: the counts file lists basal libraries first (MCL1-DG..DL) while the
# series matrix lists luminal first. The link is the "Sample name: MCL1-xx"
# field in !Sample_description, and build_sample_sheet.R joins on it.
set -euo pipefail
cd "$(dirname "$0")"

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450
mkdir -p data/raw data/annotation results/figures

fetch () {  # url, dest
  if [ -s "$2" ]; then echo "  cached  $2"; return; fi
  echo "  fetch   $2"
  curl -sSfL -o "$2" "$1"
}

echo "GSE60450 counts and metadata:"
fetch "$GEO/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz" data/raw/GSE60450_counts.txt.gz
fetch "$GEO/matrix/GSE60450_series_matrix.txt.gz"           data/raw/GSE60450_series_matrix.txt.gz

# Annotation is taken from primary sources rather than an org.*.eg.db snapshot,
# so the pipeline pins release identity by file rather than by package version.
echo "annotation:"
fetch "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz" \
      data/annotation/Mus_musculus.gene_info.gz
fetch "https://reactome.org/download/current/NCBI2Reactome.txt" \
      data/annotation/NCBI2Reactome.txt

echo "sample sheet:"
Rscript build_sample_sheet.R

echo
echo "inputs ready:"
ls -lh data/raw data/annotation data/sample_sheet.csv | grep -v '^total' | awk '{print "  "$NF" "$5}'
