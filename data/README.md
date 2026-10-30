# Input data

No data are redistributed in this repository; every dataset is public and must be downloaded from
its source. The scripts read from this directory using the layout below. Set the environment
variable `HCC_ROOT` to the repository root before running them (the scripts default to the current
working directory).

```
data/
├── raw/
│   ├── geo/                      GSE14520, GSE76427, GSE62232, GSE149614 (GEO supplementary files)
│   ├── tcga/                     TCGA-LIHC STAR counts and clinical XML from the GDC
│   ├── gtex/                     GTEx v8 liver and whole-blood cis-eQTL pairs
│   └── finngen/                  FinnGen R12 HCC (C3_HEPATOCELLULAR_CARCINOMA) summary statistics
└── processed/                    objects written by the preprocessing scripts (created at run time)
```

| Dataset | Role | Accession / portal |
|---|---|---|
| GSE14520 (GPL3921) | discovery cohort, 225 tumours + 220 adjacent tissues | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520 |
| TCGA-LIHC | validation, prognosis, immune microenvironment, subtypes | https://portal.gdc.cancer.gov/projects/TCGA-LIHC |
| GSE76427 | validation | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76427 |
| GSE62232 | validation | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62232 |
| GSE149614 | single-cell RNA sequencing, 10 patients | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614 |
| GTEx v8 | cis-eQTL instruments (liver, whole blood) | https://gtexportal.org |
| FinnGen R12 | HCC outcome summary statistics (947 cases) | https://r12.finngen.fi |

Two preparatory steps are not scripted because they are interactive downloads: the GDC client
transfer of TCGA-LIHC (`gdc-client download -m gdc_manifest.txt`) and the GTEx and FinnGen file
selection on their portals. Everything after that is covered by scripts `01`–`23`.

Expression values are used as supplied by each source (log2 intensities for the microarray
cohorts, log2(STAR counts + 1) for TCGA-LIHC) and are z-scored gene-wise within each cohort
immediately before modelling, so no cross-platform normalisation is applied to the raw matrices.
