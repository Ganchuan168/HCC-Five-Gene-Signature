#!/usr/bin/env python3
"""图1 研究流程图：三栏（数据来源 → 分析流程 → 产出），栏间为单一大箭头"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.family"] = "DejaVu Sans"
W, H = 10.2, 5.8
fig, ax = plt.subplots(figsize=(W, H), dpi=300)
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

C_DATA, C_AN, C_OUT = "#DEEBF7", "#FDE0D5", "#E5F5E0"
E_DATA, E_AN, E_OUT = "#2B6CB0", "#C0392B", "#2F855A"
FS_BODY, FS_TITLE = 9.0, 10.5


def band(x, y, w, h, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.10",
                                linewidth=0, facecolor=fc, alpha=0.45, zorder=0))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.10",
                                linewidth=1.0, facecolor="none", edgecolor=ec,
                                linestyle=(0, (4, 2)), zorder=1))


def box(x, y, w, h, text, ec, fs=FS_BODY):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.08",
                                linewidth=1.0, facecolor="white", edgecolor=ec, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            linespacing=1.42, zorder=3)


def barrow(x1, y, x2):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>", mutation_scale=22,
                                 color="#444444", linewidth=2.2, shrinkA=0, shrinkB=0, zorder=4))


def varrow(x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>", mutation_scale=10,
                                 color="#7A7A7A", linewidth=1.1, shrinkA=0, shrinkB=0, zorder=4))


# ---------------- 三个栏 ----------------
band(0.10, 0.25, 2.85, 5.00, C_DATA, E_DATA)
band(3.28, 0.25, 3.67, 5.00, C_AN, E_AN)
band(7.28, 0.25, 2.82, 5.00, C_OUT, E_OUT)

ax.text(1.525, 5.50, "Data sources", ha="center", fontsize=FS_TITLE, fontweight="bold", color=E_DATA)
ax.text(5.115, 5.50, "Analysis pipeline", ha="center", fontsize=FS_TITLE, fontweight="bold", color=E_AN)
ax.text(8.69, 5.50, "Outputs", ha="center", fontsize=FS_TITLE, fontweight="bold", color=E_OUT)

barrow(3.00, 2.75, 3.23)
barrow(7.02, 2.75, 7.23)

# ---------------- 数据来源 ----------------
box(0.25, 4.09, 2.55, 0.95,
    "Discovery cohort\nGSE14520 (GPL3921)\n225 tumour / 220 adjacent", E_DATA)
box(0.25, 2.91, 2.55, 0.95,
    "Validation cohorts\nTCGA-LIHC 371 / 50\nGSE76427, GSE62232", E_DATA)
box(0.25, 1.73, 2.55, 0.95,
    "Causal resources\nGTEx v8 cis-eQTL (liver, blood)\nFinnGen R12 HCC (947 cases)", E_DATA, fs=8.4)
box(0.25, 0.55, 2.55, 0.95,
    "Single-cell atlas\nGSE149614 (10 patients)", E_DATA)

# ---------------- 分析流程 ----------------
steps = [
    "1  Differential expression (limma)\n    1,019 DEGs (466 up / 553 down)",
    "2  WGCNA co-expression network\n    turquoise module → 418 candidates",
    "3  GO / KEGG / Hallmark enrichment\n    cell cycle, DNA replication, MYC/E2F",
    "4  Multi-algorithm consensus ranking\n    10 algorithms, 55 combinations → 5 genes",
    "5  Discrimination and prognosis\n    AUC, calibration, DCA, Cox, subtypes",
    "6  cis-eQTL Mendelian randomisation\n    Wald ratio, F > 10, BH-FDR",
]
h, gap = 0.66, 0.15
ys = [0.45 + i * (h + gap) for i in range(6)]          # 自下而上
for y, s in zip(ys[::-1], steps):                        # 第 1 步在最上方
    box(3.40, y, 3.43, h, s, E_AN, fs=8.6)
for i in range(5):
    varrow(5.115, ys[i + 1] - 0.005, ys[i] + h + 0.005)

# ---------------- 产出 ----------------
box(7.40, 4.09, 2.58, 0.95,
    "Five-gene signature\nCDKN2A, CCT6A, CSTB,\nACLY, EPRS1", E_OUT)
box(7.40, 2.91, 2.58, 0.95,
    "Two molecular subtypes (k = 2)\nC1: higher NK/endothelial\nC2: higher risk, worse OS", E_OUT, fs=8.4)
box(7.40, 1.73, 2.58, 0.95,
    "Independent prognostic value\nHR 1.67 per SD\nC-index 0.66", E_OUT)
box(7.40, 0.55, 2.58, 0.95,
    "Causal screening (cis-eQTL MR)\nno gene survives FDR;\nCSTB nominal in liver", E_OUT, fs=8.4)

ROOT = os.environ.get("HCC_ROOT", os.getcwd())
FIG = os.path.join(ROOT, "results/figures")
os.makedirs(FIG, exist_ok=True)
base = os.path.join(FIG, "Fig1_workflow")
fig.savefig(base + ".png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(base + ".pdf", bbox_inches="tight", facecolor="white")
print("Fig1 done")
