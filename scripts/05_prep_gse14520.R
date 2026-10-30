#!/usr/bin/env Rscript
## 05 GSE14520（发现集）预处理
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(limma)})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())

parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  meta <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || startsWith(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_")) {
      p <- strsplit(line, "\t")[[1]]
      meta[[sub("!Sample_", "", p[1])]] <- gsub('"', "", p[-1])
    }
  }
  tab <- read.delim(con, header = TRUE, check.names = FALSE)
  rownames(tab) <- tab[[1]]; tab <- tab[, -1, drop = FALSE]
  list(expr = as.matrix(tab), meta = meta)
}

collapse_probes <- function(mat, probe2sym) {
  s <- probe2sym[rownames(mat)]
  ok <- !is.na(s) & s != ""
  mat <- mat[ok, , drop = FALSE]; s <- s[ok]
  ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]; s <- s[ord]
  keep <- !duplicated(s)
  mat <- mat[keep, , drop = FALSE]; rownames(mat) <- s[keep]
  mat
}

sm <- parse_series_matrix(file.path(ROOT, "data/raw/geo/GSE14520-GPL3921_series_matrix.txt.gz"))
cat("原始维度:", dim(sm$expr), "\n")
suppressPackageStartupMessages(library(hgu133a2.db))
p2s <- suppressMessages(AnnotationDbi::mapIds(hgu133a2.db, keys = rownames(sm$expr),
      keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
m <- collapse_probes(sm$expr, p2s)
cat("注释后:", dim(m), "\n")

tis <- sm$meta$characteristics_ch1
# characteristics 有多行，取含 "Tissue:" 的那一行
idx <- grep("^Tissue:", tis)
if (length(idx) >= 1) {
  tis <- tis[[idx[1]]]
} else {
  tis <- sm$meta$title
}
grp <- ifelse(grepl("Non-Tumor|non-tumor|Non tumor", tis), "Normal", "Tumor")
cat("分组:\n"); print(table(grp))
cat("样本名示例:", head(sm$meta$title, 3), "\n")

saveRDS(list(expr = m, group = grp, title = sm$meta$title, gsm = colnames(sm$expr)),
        file.path(ROOT, "data/processed/GSE14520.rds"))
