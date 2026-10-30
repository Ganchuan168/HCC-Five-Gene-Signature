#!/usr/bin/env python3
"""GSE149614：质控后聚类、注释与特征基因定位（内存友好版）"""
import os, gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

ROOT = os.environ.get("HCC_ROOT", os.getcwd())
np.random.seed(20261029)
PROC = os.path.join(ROOT, "data/processed")
RES = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures")
sc.settings.verbosity = 2
sc.settings.n_jobs = 4

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

genes = ["TP53I3", "CDKN2A", "CSTB", "IDI1", "CCT6A"]
genes = [g for g in genes if g in a.var_names]
E = pd.DataFrame(a[:, genes].X.toarray(), index=a.obs_names, columns=genes)
E["celltype"] = "NA"; E["group"] = a.obs["group"].values; E["sample"] = a.obs["sample"].values
E.to_csv(os.path.join(RES, "26_sc_gene_expression.csv.gz"), compression="gzip")
print("gene expression saved", E.shape, flush=True)

sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="sample")
hv = a.var_names[a.var["highly_variable"]].tolist()
b = a[:, hv].copy()
del a; gc.collect()
print("HVG subset", b.shape, flush=True)

sc.tl.pca(b, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(b, n_neighbors=15, n_pcs=30)
sc.tl.leiden(b, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(b)
print(b.obs.leiden.value_counts(), flush=True)

markers = {
    "Hepatocyte":    ["ALB", "APOA1", "TTR", "TF", "APOB", "HP"],
    "T cell":        ["CD3D", "CD3E", "CD2", "TRAC"],
    "CD8+ T cell":   ["CD8A", "CD8B", "GZMK"],
    "NK cell":       ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "B cell":        ["MS4A1", "CD79A", "CD19"],
    "Plasma cell":   ["MZB1", "JCHAIN", "IGHG1"],
    "Myeloid":       ["CD68", "CD14", "LYZ", "AIF1"],
    "Endothelial":   ["PECAM1", "VWF", "CDH5", "ENG"],
    "Fibroblast":    ["COL1A1", "COL1A2", "ACTA2", "DCN"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM"],
    "Cycling":       ["MKI67", "TOP2A", "CENPF"],
}
avail = {k: [g for g in v if g in b.var_names] for k, v in markers.items()}
for k, v in avail.items():
    if v:
        sc.tl.score_genes(b, v, score_name="s_" + k)
score_cols = ["s_" + k for k, v in avail.items() if v]
S = b.obs.groupby("leiden", observed=True)[score_cols].mean()
S.columns = [c[2:] for c in S.columns]
print(S.round(3), flush=True)
lab = S.idxmax(axis=1)
print(lab, flush=True)
b.obs["celltype"] = b.obs["leiden"].map(lab).astype(str)
E["celltype"] = b.obs["celltype"].reindex(E.index).values
E.to_csv(os.path.join(RES, "26_sc_gene_expression.csv.gz"), compression="gzip")

comp = b.obs.groupby(["sample", "group", "celltype"], observed=True).size().reset_index(name="n")
comp["frac"] = comp["n"] / comp.groupby("sample")["n"].transform("sum")
comp.to_csv(os.path.join(RES, "27_sc_celltype_composition.csv"), index=False)
print(comp.groupby(["group", "celltype"])["frac"].mean().unstack().round(3), flush=True)

summ = []
for g in genes:
    tmp = pd.DataFrame({"v": E[g].values, "ct": E["celltype"].values, "grp": E["group"].values})
    s = tmp.groupby(["ct", "grp"])["v"].agg(mean="mean", pct=lambda x: (x > 0).mean() * 100, n="size")
    s = s.reset_index(); s["gene"] = g
    summ.append(s)
summ = pd.concat(summ)
summ.to_csv(os.path.join(RES, "28_sc_gene_by_celltype.csv"), index=False)
print(summ[summ.grp == "Tumor"].round(3).head(40), flush=True)

b.obs[["leiden", "celltype", "group", "sample"]].to_csv(os.path.join(RES, "29_sc_cell_meta.csv.gz"),
                                                        compression="gzip")
b.write_h5ad(os.path.join(PROC, "GSE149614_sc_hvg.h5ad"))
print("done", flush=True)
