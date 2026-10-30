#!/bin/bash
HCC_ROOT="${HCC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
F=${HCC_ROOT}/data/raw/tcga/TCGA-LIHC.star_counts.tsv.gz
U=https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz
for i in $(seq 1 200); do
  if gzip -t "$F" 2>/dev/null; then echo "OK $(stat -c %s $F)" > /tmp/tcga_matrix.done; exit 0; fi
  curl -s -C - --max-time 120 -o "$F" "$U"
  sleep 5
done
