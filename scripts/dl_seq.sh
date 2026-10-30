#!/bin/bash
HCC_ROOT="${HCC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
cd ${HCC_ROOT}/data/raw
get() { # $1=path $2=url
  for i in $(seq 1 60); do
    gzip -t "$1" 2>/dev/null && return 0
    curl -s -C - --max-time 600 -o "$1" "$2"
    sleep 3
  done
}
get tcga/TCGA-LIHC.star_counts.tsv.gz "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz"
get geo/GSE14520-GPL3921_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz"
get geo/GSE76427_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz"
get geo/GSE25097_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25097/matrix/GSE25097_series_matrix.txt.gz"
get geo/GSE149614_count.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
echo ALLDONE > /tmp/dl_seq.done
