#!/usr/bin/env Rscript
## 03 差异表达分析 + WGCNA（TCGA-LIHC）
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(limma); library(WGCNA); library(jsonlite)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
PROC <- file.path(ROOT, "data/processed")
RES  <- file.path(ROOT, "results/tables")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
options(stringsAsFactors = FALSE)

## ---------- 1. 读入表达矩阵 ----------
expr <- fread(file.path(ROOT, "data/raw/tcga/TCGA-LIHC.star_counts.tsv.gz"))
ens <- sub("\\..*$", "", expr[[1]])
m <- as.matrix(expr[, -1]); rownames(m) <- ens
storage.mode(m) <- "double"

## Xena star_counts 为 log2(count+1)；若不是则转换
if (max(m, na.rm = TRUE) > 30) m <- log2(m + 1)

## 样本分组：01 = 肿瘤，11 = 癌旁
grp <- ifelse(substr(colnames(m), 14, 15) == "01", "Tumor", "Normal")
grp[!substr(colnames(m), 14, 15) %in% c("01", "11")] <- NA
keep_s <- !is.na(grp); m <- m[, keep_s]; grp <- grp[keep_s]
## 同一患者多处取样时保留第一个
pid <- substr(colnames(m), 1, 12)
m <- m[, !duplicated(pid)]; grp <- grp[!duplicated(pid)]
cat("表达矩阵:", dim(m), " 肿瘤:", sum(grp == "Tumor"), " 癌旁:", sum(grp == "Normal"), "\n")

## 过滤低表达基因
keep_g <- rowMeans(m) > 1 & apply(m, 1, function(x) sum(x > 0) > 0.3 * length(x))
m <- m[keep_g, ]
cat("过滤后基因数:", nrow(m), "\n")

## ---------- 2. 差异表达 ----------
design <- model.matrix(~ 0 + factor(grp, levels = c("Normal", "Tumor")))
colnames(design) <- c("Normal", "Tumor")
fit <- lmFit(m, design)
cont <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg$gene <- rownames(deg)
fwrite(deg, file.path(RES, "TCGA_DEG_all.csv"))

sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]
cat("DEG 数:", nrow(sig), " 上调:", sum(sig$logFC > 0), " 下调:", sum(sig$logFC < 0), "\n")
fwrite(sig, file.path(RES, "TCGA_DEG_sig.csv"))

## ---------- 3. WGCNA ----------
dat <- t(m)
gsg <- goodSamplesGenes(dat, verbose = 0)
dat <- dat[gsg$goodSamples, gsg$goodGenes]
cat("WGCNA 输入:", dim(dat), "\n")

## 取变异最大的前 25% 基因
vars <- apply(dat, 2, var)
dat <- dat[, order(vars, decreasing = TRUE)[1:floor(0.25 * ncol(dat))]]

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed hybrid", verbose = 0)
saveRDS(sft, file.path(PROC, "sft.rds"))
fwrite(sft$powerEstimate_, file.path(RES, "soft_threshold.csv"))
