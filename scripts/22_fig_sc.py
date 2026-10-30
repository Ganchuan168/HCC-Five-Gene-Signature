#!/usr/bin/env python3
"""图8 单细胞转录组：细胞图谱与特征基因定位"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.environ.get("HCC_ROOT", os.getcwd())
PROC = os.path.join(ROOT, "data/processed")
RES = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures")

b = sc.read_h5ad(os.path.join(PROC, "GSE149614_sc_hvg.h5ad"))
E = pd.read_csv(os.path.join(RES, "26_sc_gene_expression.csv.gz"), index_col=0)
genes = ["CDKN2A", "CCT6A", "CSTB", "ACLY", "EPRS1"]

order = ["Hepatocyte", "Cholangiocyte", "Fibroblast", "Endothelial", "Myeloid",
         "T cell", "NK cell", "B cell", "Plasma cell"]
order = [c for c in order if c in set(b.obs["celltype"])]
colors = dict(zip(order, ["#D62728", "#8C564B", "#BCBD22", "#17BECF", "#7F7F7F",
                          "#1F77B4", "#2CA02C", "#9467BD", "#E377C2"]))

fig = plt.figure(figsize=(13.5, 10.5), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.28, wspace=0.26)

# A UMAP
ax = fig.add_subplot(gs[0, :])
um = b.obsm["X_umap"]
for c in order:
    m = (b.obs["celltype"] == c).values
    ax.scatter(um[m, 0], um[m, 1], s=0.35, c=colors[c], label=f"{c} ({(m).sum():,})",
               linewidths=0, rasterized=True)
ax.set_xlabel("UMAP1", fontsize=9); ax.set_ylabel("UMAP2", fontsize=9)
ax.set_title("A  Single-cell atlas of HCC (GSE149614, 69,471 cells)", fontsize=10, fontweight="bold", loc="left")
ax.set_xticks([]); ax.set_yticks([])
leg = ax.legend(markerscale=12, fontsize=7.5, loc="center left", bbox_to_anchor=(1.005, 0.5),
                frameon=False)
plt.setp(leg.get_texts(), fontsize=7.5)

# B dot plot
ax2 = fig.add_subplot(gs[1, 0])
sub = E[E["celltype"].isin(order)]
X, Y, S, C = [], [], [], []
for i, g in enumerate(genes):
    for j, c in enumerate(order):
        v = sub.loc[sub["celltype"] == c, g].values
        X.append(j); Y.append(i); S.append((v > 0).mean() * 100); C.append(v.mean())
ax2.scatter(X, Y, s=np.array(S) * 1.6, c=C, cmap="Reds", edgecolors="grey", linewidths=0.3,
            vmin=0, vmax=max(C) * 1.05)
ax2.set_xticks(range(len(order))); ax2.set_xticklabels(order, rotation=45, ha="right", fontsize=7.5)
ax2.set_yticks(range(len(genes))); ax2.set_yticklabels(genes, fontsize=8.5, style="italic")
ax2.set_xlim(-0.6, len(order) - 0.4); ax2.set_ylim(len(genes) - 0.5, -0.6)
ax2.set_title("B  Expression of signature genes across cell types", fontsize=10, fontweight="bold", loc="left")
h = [Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markersize=np.sqrt(s * 1.6),
            label=f"{s}%") for s in [10, 40, 80]]
ax2.legend(handles=h, title="% expressing", fontsize=6.5, title_fontsize=7, loc="upper right",
           bbox_to_anchor=(1.02, 1.02), frameon=False, labelspacing=0.9)

# C hepatocytes tumour vs adjacent
ax3 = fig.add_subplot(gs[1, 1])
hep = E[E["celltype"] == "Hepatocyte"]
w = 0.36
xs = np.arange(len(genes))
for k, (grp, lab, col) in enumerate([("Normal", "adjacent", "#4DAF4A"), ("Tumor", "tumour", "#E41A1C")]):
    vals = [(hep.loc[hep["group"] == grp, g] > 0).mean() * 100 for g in genes]
    ax3.bar(xs + (k - 0.5) * w, vals, width=w, color=col, label=f"{lab}-derived hepatocytes",
            edgecolor="none")
    for x, v in zip(xs + (k - 0.5) * w, vals):
        ax3.text(x, v + 1.5, f"{v:.0f}", ha="center", fontsize=6.5)
ax3.set_xticks(xs); ax3.set_xticklabels(genes, fontsize=8.5, style="italic")
ax3.set_ylabel("% of hepatocytes expressing", fontsize=8.5)
ax3.set_ylim(0, 105)
ax3.legend(fontsize=7.5, frameon=False, loc="upper left")
ax3.set_title("C  Signature genes in malignant vs adjacent hepatocytes", fontsize=10, fontweight="bold", loc="left")
for s in ("top", "right"):
    ax3.spines[s].set_visible(False)
ax3.tick_params(labelsize=7.5)

plt.savefig(os.path.join(FIG, "Fig8_singlecell.png"), bbox_inches="tight", facecolor="white")
print("Fig8 saved")
