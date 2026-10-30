#!/usr/bin/env python3
"""整合 TCGA-LIHC 与 GEO 系列矩阵，输出表达矩阵与样本分组。"""
import os, sys, gzip, glob, re
import pandas as pd
import numpy as np

ROOT = os.environ.get("HCC_ROOT", os.getcwd())
RAW = os.path.join(ROOT, "data/raw")
OUT = os.path.join(ROOT, "data/processed")
os.makedirs(OUT, exist_ok=True)


def prep_tcga():
    files = glob.glob(f"{RAW}/tcga/counts/*.tsv")
    man = pd.read_csv(f"{RAW}/tcga/file_manifest.tsv", sep="\t")
    fname2meta = {r.file_name: r for r in man.itertuples()}
    mats, groups = {}, {}
    for f in files:
        name = os.path.basename(f)
        if name not in fname2meta:
            continue
        df = pd.read_csv(f, sep="\t", comment="#", index_col=0)
        if "unstranded" not in df.columns:
            continue
        s = df["unstranded"].iloc[4:]          # 跳过 N_* 汇总行
        s = s[~s.index.str.startswith("N_")]
        mats[fname2meta[name].sample] = s
        groups[fname2meta[name].sample] = fname2meta[name].sample_type
    m = pd.DataFrame(mats)
    m.to_csv(f"{OUT}/TCGA_LIHC_counts.tsv.gz", sep="\t", compression="gzip")
    pd.Series(groups, name="sample_type").to_csv(f"{OUT}/TCGA_LIHC_groups.tsv", sep="\t")
    print("TCGA:", m.shape)


def parse_series_matrix(path):
    """返回 (表达矩阵, 样本注释 DataFrame)。"""
    meta, lines = {}, []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                parts = [p.strip('"\n') for p in line.rstrip("\n").split("\t")]
                k, vals = parts[0], parts[1:]
                meta.setdefault(k, []).append(vals)
            if line.startswith(("!Series_title", "!Series_summary")):
                pass
    df = pd.read_csv(path, sep="\t", skiprows=_skip_to_table(path), index_col=0)
    df = df.dropna(how="all")
    ann = pd.DataFrame({k.replace("!Sample_", ""): v for k, v in meta.items()
                        if len(v) == 1 and len(v[0]) == df.shape[1]})
    return df, ann


def _skip_to_table(path):
    n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            n += 1
            if line.startswith("!series_matrix_table_begin"):
                return n
    return n


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "tcga":
        prep_tcga()
    else:
        for g in sys.argv[1:]:
            p = f"{RAW}/geo/{g}"
            df, ann = parse_series_matrix(p)
            print(g, df.shape)
            print(ann.head())
