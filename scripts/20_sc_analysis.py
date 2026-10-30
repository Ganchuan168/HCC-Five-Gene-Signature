#!/usr/bin/env python3
"""GSE149614 单细胞转录组：质控、聚类、注释与特征基因定位"""
import os, sys, gc, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 1
ROOT = os.environ.get("HCC_ROOT", os.getcwd())
np.random.seed(20261029)
RAW = os.path.join(ROOT, "data/raw/geo/GSE149614_sc.count.txt.gz")
PROC = os.path.join(ROOT, "data/processed")
FIG = os.path.join(ROOT, "results/figures")
RES = os.path.join(ROOT, "results/tables")
H5 = os.path.join(PROC, "GSE149614_sc.h5ad")

os.makedirs(FIG, exist_ok=True); os.makedirs(RES, exist_ok=True)

if not os.path.exists(H5):
    t0 = time.time()
    cells = pd.read_csv(RAW, sep="\t", nrows=0).columns.tolist()
    print("cells:", len(cells), flush=True)

    blocks, genes = [], []
    for chunk in pd.read_csv(RAW, sep="\t", index_col=0, chunksize=1500,
                             dtype=dict.fromkeys(cells, np.int32), engine="c"):
        genes.extend(chunk.index.tolist())
        blocks.append(sp.csr_matrix(chunk.values.astype(np.float32)))
        del chunk; gc.collect()
    print("blocks:", len(blocks), "genes:", len(genes), "%.1fs" % (time.time() - t0), flush=True)

    M = sp.vstack(blocks, format="csr")         # genes × cells
    del blocks; gc.collect()
    keep = np.asarray((M > 0).sum(axis=1)).ravel() >= 10
    M = M[keep]; genes = [g for g, k in zip(genes, keep) if k]
    print("after gene filter:", M.shape, flush=True)

    a = ad.AnnData(X=sp.csr_matrix(M.T), obs=pd.DataFrame(index=cells),
                   var=pd.DataFrame(index=genes))
    a.obs["sample"] = [c.split("_")[0] for c in cells]
    a.obs["patient"] = [c.split("_")[0][:5] for c in cells]
    a.obs["group"] = ["Tumor" if c.split("_")[0].endswith("T") else "Normal" for c in cells]
    del M; gc.collect()
    a.write_h5ad(H5)
    print("written", H5, "%.1fs" % (time.time() - t0), flush=True)

a = sc.read_h5ad(H5)
print(a)

a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
print(a.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe(), flush=True)

sc.pp.filter_cells(a, min_genes=200)
a = a[a.obs["n_genes_by_counts"] < 6000].copy()
a = a[a.obs["pct_counts_mt"] < 20].copy()
print("after QC:", a.shape, flush=True)
print(a.obs.group.value_counts(), flush=True)

sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key="sample")
a.raw = a
a = a[:, a.var.highly_variable].copy()
sc.pp.scale(a, max_value=10)
sc.tl.pca(a, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(a, n_neighbors=15, n_pcs=30)
sc.tl.leiden(a, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2)
sc.tl.umap(a)
print(a.obs.leiden.value_counts(), flush=True)

## 细胞类型注释
markers = {
    "Hepatocyte":     ["ALB", "APOA1", "TTR", "TF", "APOB", "HP"],
    "T cell":         ["CD3D", "CD3E", "CD2", "TRAC"],
    "CD8+ T cell":    ["CD8A", "CD8B", "GZMK"],
    "NK cell":        ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "B cell":         ["MS4A1", "CD79A", "CD19"],
    "Plasma cell":    ["MZB1", "JCHAIN", "IGHG1"],
    "Myeloid":        ["CD68", "CD14", "LYZ", "AIF1"],
    "Endothelial":    ["PECAM1", "VWF", "CDH5", "ENG"],
    "Fibroblast":     ["COL1A1", "COL1A2", "ACTA2", "DCN"],
    "Cholangiocyte":  ["KRT19", "KRT7", "EPCAM"],
    "Cycling":        ["MKI67", "TOP2A", "CENPF"],
}
avail = {k: [g for g in v if g in a.raw.var_names] for k, v in markers.items()}
sc.tl.score_genes(a, avail["Hepatocyte"], score_name="s_Hep")
for k, v in avail.items():
    sc.tl.score_genes(a, v, score_name="s_" + k)

score_cols = ["s_" + k for k in avail]
S = a.obs.groupby("leiden", observed=True)[score_cols].mean()
S.columns = [c[2:] for c in S.columns]
lab = S.idxmax(axis=1)
print(S.round(3), flush=True)
print(lab, flush=True)
a.obs["celltype"] = a.obs["leiden"].map(lab).astype(str)

## 特征基因表达（原始 log 数据）
genes = ["TP53I3", "CDKN2A", "CSTB", "IDI1", "CCT6A"]
genes = [g for g in genes if g in a.raw.var_names]
print("genes found:", genes, flush=True)
E = pd.DataFrame(a.raw[:, genes].X.toarray(), index=a.obs_names, columns=genes)
E["celltype"] = a.obs["celltype"].values
E["group"] = a.obs["group"].values
E["sample"] = a.obs["sample"].values
E.to_csv(os.path.join(RES, "26_sc_gene_expression.csv.gz"), compression="gzip")

## 汇总：细胞类型比例
comp = a.obs.groupby(["sample", "group", "celltype"], observed=True).size().reset_index(name="n")
comp["frac"] = comp["n"] / comp.groupby("sample")["n"].transform("sum")
comp.to_csv(os.path.join(RES, "27_sc_celltype_composition.csv"), index=False)
print(comp.groupby(["group", "celltype"])["frac"].mean().unstack().round(3), flush=True)

## 特征基因在各细胞类型的平均表达与阳性率
summ = []
for g in genes:
    tmp = pd.DataFrame({"v": E[g].values, "ct": E["celltype"].values, "grp": E["group"].values})
    s = tmp.groupby(["ct", "grp"])["v"].agg(mean="mean", pct=lambda x: (x > 0).mean() * 100, n="size")
    s = s.reset_index(); s["gene"] = g
    summ.append(s)
summ = pd.concat(summ)
summ.to_csv(os.path.join(RES, "28_sc_gene_by_celltype.csv"), index=False)
print(summ[summ.grp == "Tumor"].round(3).head(30), flush=True)

a.write_h5ad(os.path.join(PROC, "GSE149614_sc_processed.h5ad"))
print("done", flush=True)
