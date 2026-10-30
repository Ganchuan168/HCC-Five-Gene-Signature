# A five-gene signature for hepatocellular carcinoma

Analysis code accompanying the manuscript *"A five-gene signature for hepatocellular carcinoma
derived from co-expression network analysis and multi-algorithm consensus ranking: discrimination,
prognosis and association with the immune microenvironment"* (submitted to *Biomarker Research*).

The repository contains every script needed to reproduce the differential expression, co-expression
network, multi-algorithm consensus ranking, validation, survival, immune, Mendelian randomisation
and single-cell analyses, and all eight figures of the paper.

## What the code does

| Step | Scripts | Output |
|---|---|---|
| Data retrieval and preprocessing | `01`–`07c` | expression matrices and clinical tables for GSE14520, TCGA-LIHC, GSE76427 and GSE62232 |
| Differential expression and WGCNA | `03`, `06` | 1,019 DEGs; turquoise module; 418 candidate genes |
| Consensus feature selection | `08`–`08d`, `ml_ranks.R` | 10 algorithms applied singly and in all 45 pairs (55 rankings), patient-grouped 5-fold CV, 1-SE rule, five-gene model |
| Validation, calibration, decision curves | `07b`, `10` | AUCs, calibration slope and intercept, Brier score, net benefit |
| Survival, immune infiltration, subtypes | `09`, `09b`, `12` | Cox models, ssGSEA scores, consensus clustering (k = 2) |
| Mendelian randomisation | `10_mr.R`, `23_mr2.R`, `extract_eqtl.sh` | cis-eQTL Wald ratio / IVW, BH-FDR |
| Enrichment | `11` | GO, KEGG and Hallmark results |
| Single-cell | `20`–`22` | 69,471 cells, clustering, per-gene positivity |
| Figures | `13`–`19`, `22` | Figures 1–8 |

## Requirements

- R 4.6.1 with limma, WGCNA, glmnet, randomForest, e1071, survival, survminer, pROC, rms, cluster,
  ConsensusClusterPlus, GSVA, org.Hs.eg.db, clusterProfiler, TwoSampleMR
- Python 3.12 with pandas, numpy, scikit-learn, scanpy 1.12.4, anndata 0.13.4, leidenalg 0.12.0,
  matplotlib, python-docx

## Data

All input data are public; nothing is redistributed here. See `data/README.md` for the accession
numbers and the directory layout the scripts expect. GTEx v8 cis-eQTL and FinnGen R12 summary
statistics must be downloaded from their portals before `10_mr.R` and `23_mr2.R` can run.

## Running the code

```bash
# 1. create the directory layout described in data/README.md and place the downloads there
# 2. run the scripts in numeric order; later scripts read the objects written by earlier ones
Rscript scripts/06_deg_wgcna.R
Rscript scripts/08b_ml_patient_split.R
Rscript scripts/10_calibration_dca.R
Rscript scripts/10_mr.R
python3 scripts/20_sc_analysis.py
```

Scripts write to `results/` and `results/figures/`; both are created on first run. Paths at the top
of each script point to the project root and can be edited in one place.

## Reproducibility notes

- Cross-validation is grouped by patient, so repeated samples from one tumour never straddle the
  training and test folds.
- Expression values are z-scored within each cohort before the linear predictor is computed; the
  cut-point −3.49 is fixed on that scale.
- `set.seed(20261029)` (or `np.random.seed(20261029)`) is set in every script that resamples or
  clusters, so runs from this release are deterministic. The published numbers were produced by the
  same scripts before the seed was fixed, so bootstrap and clustering quantities may differ in the
  third decimal on a rerun.

## Citation

If you use this code, please cite the manuscript above and the archived release:

> Gan C, Gan L, Xia M. A five-gene signature for hepatocellular carcinoma: analysis code.
> Zenodo; 2026. doi:10.5281/zenodo.23804625

## Licence

MIT — see `LICENSE`.

## Contact

Meng Xia, School of Basic Medicine, Guangxi University of Chinese Medicine, Nanning, China.
xiam@gxtcmu.edu.cn · ORCID 0000-0001-8075-9698
