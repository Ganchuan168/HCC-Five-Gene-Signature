#!/usr/bin/env Rscript
## 04 GEO 验证集预处理
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(limma)})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); PROC <- file.path(ROOT, "data/processed")

parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  meta <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || startsWith(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_")) {
      p <- strsplit(line, "\t")[[1]]
      k <- sub("!Sample_", "", p[1])
      meta[[k]] <- c(meta[[k]], list(gsub('"', "", p[-1])))
    }
  }
  tab <- read.delim(con, header = TRUE, check.names = FALSE)
  rownames(tab) <- tab[[1]]; tab <- tab[, -1, drop = FALSE]
  list(expr = as.matrix(tab), meta = meta)
}

pick_chr <- function(meta, pattern) {
  for (v in meta$characteristics_ch1)
    if (any(grepl(pattern, v))) return(v)
  NULL
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

## ---------- GSE62232（GPL570） ----------
p <- file.path(ROOT, "data/raw/geo/GSE62232_series_matrix.txt.gz")
if (file.exists(p)) {
  sm <- parse_series_matrix(p)
  suppressPackageStartupMessages(library(hgu133plus2.db))
  p2s <- suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db, keys = rownames(sm$expr),
        keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
  m <- collapse_probes(sm$expr, p2s)
  if (max(m, na.rm = TRUE) > 30) m <- log2(m + 1)
  m <- normalizeBetweenArrays(m, method = "quantile")
  title <- unlist(sm$meta$title)
  grp <- ifelse(grepl("T$", title), "Tumor", "Normal")
  keep <- !duplicated(sub("\\.rep.*$", "", title))
  m <- m[, keep]; grp <- grp[keep]; title <- title[keep]
  saveRDS(list(expr = m, group = grp, sample = title), file.path(PROC, "GSE62232.rds"))
  cat("GSE62232:", dim(m), "\n"); print(table(grp))
}

## ---------- GSE76427（GPL10558，Illumina HT-12 v4） ----------
p2 <- file.path(ROOT, "data/raw/geo/GSE76427_series_matrix.txt.gz")
if (file.exists(p2)) {
  sm2 <- parse_series_matrix(p2)
  suppressPackageStartupMessages(library(illuminaHumanv4.db))
  map2 <- suppressMessages(AnnotationDbi::mapIds(illuminaHumanv4.db, keys = rownames(sm2$expr),
          keytype = "PROBEID", column = "SYMBOL", multiVals = "first"))
  m2 <- collapse_probes(sm2$expr, map2)
  m2[!is.finite(m2)] <- NA
  if (max(m2, na.rm = TRUE) > 30) m2 <- log2(m2 + 1)
  m2 <- normalizeBetweenArrays(m2, method = "quantile")
  tisline <- pick_chr(sm2$meta, "^tissue:")
  grp2 <- ifelse(grepl("non-tumor", tisline), "Normal", "Tumor")
  ev <- pick_chr(sm2$meta, "^event_os")
  tm <- pick_chr(sm2$meta, "^time_os|^time:")
  saveRDS(list(expr = m2, group = grp2, title = unlist(sm2$meta$title),
               event = as.numeric(sub("^[a-z_]+: ", "", ev)),
               time = as.numeric(sub("^[a-z_]+: ", "", tm))),
          file.path(PROC, "GSE76427.rds"))
  cat("GSE76427:", dim(m2), "\n"); print(table(grp2))
}
