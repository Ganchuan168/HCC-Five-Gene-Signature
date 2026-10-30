#!/bin/bash
HCC_ROOT="${HCC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
cd ${HCC_ROOT}/data/raw
dl() { # $1 path, $2 url
  for i in 1 2 3 4 5 6; do
    gzip -t "$1" 2>/dev/null && { echo "OK $1"; return 0; }
    rm -f "$1" "$1.aria2"
    aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=error --max-tries=0 --retry-wait=3 -o "$(basename $1)" -d "$(dirname $1)" --allow-overwrite=true "$2" >/dev/null 2>&1
    sleep 2
  done
  echo "FAIL $1"
}
dl geo/GSE14520-GPL3921_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz"
dl geo/GSE76427_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz"
dl geo/GSE62232_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62232/matrix/GSE62232_series_matrix.txt.gz"
dl tcga/TCGA-LIHC.star_counts.tsv.gz "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz"
dl geo/GSE25097_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25097/matrix/GSE25097_series_matrix.txt.gz"
dl geo/GSE149614_count.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
echo ALLDONE > /tmp/dl_fast.done
