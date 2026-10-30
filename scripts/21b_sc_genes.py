#!/usr/bin/env python3
"""从完整单细胞矩阵中重新提取最终特征基因的表达，并与既有细胞注释合并"""
import os, gc
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.environ.get("HCC_ROOT", os.getcwd())
PROC = os.path.join(ROOT, "data/processed")
RES = os.path.join(ROOT, "results/tables")
GENES = ["CDKN2A", "CCT6A", "CSTB", "ACLY", "EPRS1"]
ALIAS = {"EPRS1": "EPRS"}  # GSE149614 矩阵沿用旧符号 EPRS

a = sc.read_h5ad(os.path.join(PROC, "GSE149614_sc.h5ad"))
print("loaded", a.shape, flush=True)

a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
sc.pp.filter_cells(a, min_genes=200)
a = a[(a.obs["n_genes_by_counts"] < 6000) & (a.obs["pct_counts_mt"] < 20)].copy()
gc.collect()
print("after QC", a.shape, flush=True)

sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

lookup = {g: (g if g in a.var_names else ALIAS.get(g)) for g in GENES}
lookup = {k: v for k, v in lookup.items() if v and v in a.var_names}
missing = [g for g in GENES if g not in lookup]
genes = list(lookup.values())
print("可用基因:", genes, "缺失:", missing, flush=True)

E = pd.DataFrame(a[:, genes].X.toarray(), index=a.obs_names,
                 columns=[k for k in GENES if k in lookup])

meta = pd.read_csv(os.path.join(RES, "29_sc_cell_meta.csv.gz"), index_col=0)
E = E.join(meta[["celltype", "group", "sample"]], how="inner")
print("merged", E.shape, flush=True)
print(E["celltype"].value_counts(), flush=True)
E.to_csv(os.path.join(RES, "26_sc_gene_expression.csv.gz"), compression="gzip")

summ = []
for g in genes:
    t = E.groupby(["celltype", "group"])[g].agg(
        mean="mean", pct=lambda x: (x > 0).mean() * 100, n="size").reset_index()
    t["gene"] = g
    summ.append(t)
summ = pd.concat(summ)
summ.to_csv(os.path.join(RES, "28_sc_gene_by_celltype.csv"), index=False)
print(summ[summ.group == "Tumor"].round(3).to_string(), flush=True)
