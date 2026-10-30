#!/usr/bin/env Rscript
## 11 功能富集：GO-BP / KEGG（候选基因）+ Hallmark GSEA（全转录组排序）
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
GMT <- file.path(ROOT, "data/raw/msigdb")

read_gmt <- function(f) {
  ln <- readLines(f); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln, "\t")
  data.frame(term = rep(sapply(sp, `[`, 1), lengths(sp) - 2),
             gene = unlist(lapply(sp, function(x) x[-(1:2)])))
}

cand <- fread(file.path(RES, "05_candidate_genes.csv"))
genes <- unique(cand$gene)
cat("候选基因数:", length(genes), "\n")

eg <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes,
        keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
eg <- unique(na.omit(eg))
cat("可映射 ENTREZ:", length(eg), "\n")

## ---- GO-BP ----
go <- enrichGO(gene = eg, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
               pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
god <- as.data.frame(go)
fwrite(god, file.path(RES, "17_GO_BP.csv"))
cat("GO BP 显著条目:", nrow(god), "\n")

## ---- KEGG（KEGG REST）----
kk <- tryCatch(enrichKEGG(eg, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
               error = function(e) NULL)
if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  kd <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kdf <- as.data.frame(kd)
  fwrite(kdf, file.path(RES, "18_KEGG.csv"))
  cat("KEGG 显著通路:", nrow(kdf), "\n")
}

## ---- Hallmark GSEA ----
deg <- fread(file.path(RES, "01_DEG_all.csv"))
deg <- deg[!is.na(logFC) & !is.na(P.Value)]
rk <- deg$logFC; names(rk) <- deg$gene
rk <- sort(rk, decreasing = TRUE)
hl <- read_gmt(file.path(GMT, "h.all.v2023.2.Hs.symbols.gmt"))
set.seed(20261029)
g <- GSEA(rk, TERM2GENE = hl, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
          eps = 0, seed = TRUE, verbose = FALSE)
gd <- as.data.frame(g)
gd$ID <- sub("^HALLMARK_", "", gd$ID)
fwrite(gd, file.path(RES, "19_hallmark_gsea.csv"))
cat("Hallmark 通路总数:", nrow(gd), " FDR<0.25:", sum(gd$p.adjust < 0.25), "\n")
print(head(gd[order(gd$p.adjust), c("ID", "setSize", "NES", "pvalue", "p.adjust")], 12))
