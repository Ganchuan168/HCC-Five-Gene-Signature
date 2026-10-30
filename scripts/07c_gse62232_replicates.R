#!/usr/bin/env Rscript
## 07c GSE62232 癌旁样本的重复性核查与敏感性分析
## GEO 样本标题显示 10 份非肿瘤样本来自 5 例患者（CHC591 与 CHC934 各 3 份、CHC932 2 份、
## CHC469 与 CHC562 各 1 份），.rep2/.rep3 为同一组织的技术重复。主分析按患者去重（n=5），
## 此处保留全部 10 份作为敏感性分析，并给出重复样本间的相关系数。
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(pROC); library(hgu133plus2.db) })
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")

src <- file.path(ROOT, "scripts/04_prep_geo.R")
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con)); meta <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || startsWith(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_")) {
      p <- strsplit(line, "\t")[[1]]; k <- sub("!Sample_", "", p[1])
      meta[[k]] <- c(meta[[k]], list(gsub('"', "", p[-1])))
    }
  }
  tab <- read.delim(con, header = TRUE, check.names = FALSE)
  rownames(tab) <- tab[[1]]; tab <- tab[, -1, drop = FALSE]
  list(expr = as.matrix(tab), meta = meta)
}
sm <- parse_series_matrix(file.path(ROOT, "data/raw/geo/GSE62232_series_matrix.txt.gz"))
p2s <- suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(sm$expr),
        keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
s <- p2s[rownames(sm$expr)]; ok <- !is.na(s) & s != ""
m <- sm$expr[ok, , drop = FALSE]; s <- s[ok]
ord <- order(rowMeans(m, na.rm = TRUE), decreasing = TRUE)
m <- m[ord, ]; s <- s[ord]; keep <- !duplicated(s)
m <- m[keep, , drop = FALSE]; rownames(m) <- s[keep]
if (max(m, na.rm = TRUE) > 30) m <- log2(m + 1)
m <- normalizeBetweenArrays(m, method = "quantile")
title <- unlist(sm$meta$title)
grp <- ifelse(grepl("T$", title), "Tumor", "Normal")

## 重复样本相关性（全基因组）
rep_i <- grep("\\.rep", title)
rep_t <- title[rep_i]; base_t <- sub("\\.rep\\d+$", "", rep_t)
base_i <- match(base_t, title)
cors <- sapply(seq_along(rep_i), function(k) cor(m[, rep_i[k]], m[, base_i[k]], method = "pearson"))
rp <- data.table(sample = rep_t, patient = base_t, r = cors)
fwrite(rp, file.path(RES, "07c_gse62232_replicate_cor.csv"))
cat("技术重复与方法内样本的相关系数：\n"); print(rp)

## 主分析（去重 n=5）与敏感性分析（全部 n=10）
sf <- readRDS(file.path(ROOT, "data/processed/signature_final.rds"))
g <- sf$genes; b <- coef(sf$model)
dedup <- !duplicated(sub("\\.rep\\d+$", "", title))
out <- list()
for (nm in c("deduplicated (n=5 controls)", "all controls (n=10)")) {
  idx <- if (grepl("dedup", nm)) dedup else rep(TRUE, length(title))
  e <- m[g, idx, drop = FALSE]; y <- as.integer(grp[idx] == "Tumor")
  z <- t(scale(t(e))); z[is.na(z)] <- 0
  lp <- as.numeric(t(z) %*% b[g]) + b["(Intercept)"]
  r <- roc(y, lp, quiet = TRUE); ci <- as.numeric(ci.auc(r))
  out[[nm]] <- data.table(scheme = nm, n = length(y), n_tumor = sum(y == 1), n_normal = sum(y == 0),
                          auc = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3])
}
out <- rbindlist(out); fwrite(out, file.path(RES, "07c_gse62232_sensitivity.csv"))
print(out)
