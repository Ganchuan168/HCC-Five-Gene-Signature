#!/bin/bash
HCC_ROOT="${HCC_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
cd ${HCC_ROOT}/data/raw/eqtl
zcat GTEx_v8_Liver_eQTL.tsv.gz | awk -F'\t' 'NR==1 || $3<1e-5' > ${HCC_ROOT}/data/processed/GTEx_liver_eQTL_sig1e5.tsv
echo DONE_LIVER > ${HCC_ROOT}/data/processed/.liver_done
zcat GTEx_v8_WholeBlood_eQTL.tsv.gz | awk -F'\t' 'NR==1 || $3<1e-5' > ${HCC_ROOT}/data/processed/GTEx_blood_eQTL_sig1e5.tsv
echo DONE_BLOOD > ${HCC_ROOT}/data/processed/.blood_done
