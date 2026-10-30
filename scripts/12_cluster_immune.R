#!/usr/bin/env Rscript
## 12 基于枢纽基因的一致性聚类分型 + 免疫微环境（ssGSEA）
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ConsensusClusterPlus); library(GSVA)
  library(survival)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

tc <- readRDS(file.path(PROC, "TCGA_LIHC.rds"))
sig <- readRDS(file.path(PROC, "signature_final.rds"))
genes <- sig$genes
mt <- tc$expr[, tc$group == "Tumor"]
cat("肿瘤样本:", ncol(mt), " 特征基因:", paste(genes, collapse = ","), "\n")

## ---------- 一致性聚类（k-means，500 次重抽样）----------
X <- t(scale(t(mt[genes, ])))
cc <- ConsensusClusterPlus(X, maxK = 4, reps = 500, pItem = 0.8, pFeature = 1,
                           clusterAlg = "km", distance = "euclidean",
                           title = file.path(FIG, "consensus"), plot = "png",
                           seed = 20261029, verbose = FALSE)
kt <- rbindlist(lapply(2:4, function(k) {
  cm <- cc[[k]]$consensusMatrix; lt <- cm[lower.tri(cm)]
  data.table(k = k, PAC = mean(lt > 0.1 & lt < 0.9),
             mean_consensus_in = mean(sapply(unique(cc[[k]]$consensusClass), function(cl) {
               idx <- which(cc[[k]]$consensusClass == cl); if (length(idx) < 2) return(NA)
               mean(cm[idx, idx][lower.tri(cm[idx, idx])])
             }), na.rm = TRUE),
             size = paste(sort(table(cc[[k]]$consensusClass), decreasing = TRUE), collapse = "/"))
}))
fwrite(kt, file.path(RES, "20_consensus_k.csv")); print(kt)

K <- 2L
cl <- cc[[K]]$consensusClass
cat("K=2 各簇样本数:\n"); print(table(cl))

## ---------- 免疫微环境 ssGSEA ----------
read_gmt <- function(f) {
  ln <- readLines(f); ln <- ln[nzchar(ln)]
  sp <- strsplit(ln, "\t")
  setNames(lapply(sp, function(x) x[-(1:2)]), sapply(sp, `[`, 1))
}
hm <- read_gmt(file.path(ROOT, "data/raw/msigdb/h.all.v2023.2.Hs.symbols.gmt"))
imm <- hm[grepl("IMMUNE|INFLAMMATORY|INTERFERON|ALLOGRAFT|COMPLEMENT|IL6|TNFA", names(hm))]

panels <- list(
  `CD8+ T cell`       = c("CD8A", "CD8B", "GZMK", "CD3D"),
  `CD4+ T cell`       = c("CD4", "IL7R", "CD40LG"),
  `Regulatory T cell` = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  `NK cell`           = c("KLRD1", "NKG7", "GNLY", "PRF1", "KLRF1"),
  `B cell`            = c("CD19", "MS4A1", "CD79A", "CD79B"),
  `Macrophage`        = c("CD68", "CD14", "CSF1R", "FCGR3A", "MRC1"),
  `Neutrophil`        = c("FCGR3B", "CSF3R", "S100A8", "S100A9"),
  `Dendritic cell`    = c("ITGAX", "CD1C", "LAMP3", "CLEC9A"),
  `Mast cell`         = c("TPSAB1", "CPA3", "MSA2"),
  `Endothelial cell`  = c("PECAM1", "VWF", "CDH5"),
  `Fibroblast`        = c("COL1A1", "COL1A2", "ACTA2", "FAP")
)
panels[["Mast cell"]] <- c("TPSAB1", "CPA3", "MS4A2")
panels <- lapply(panels, function(g) intersect(g, rownames(mt)))
stopifnot(all(lengths(panels) >= 2))

sets <- c(imm, panels)
ss <- gsva(ssgseaParam(as.matrix(mt), sets, normalize = TRUE), verbose = FALSE)
ss <- t(scale(t(ss)))
colnames(ss) <- sub("^HALLMARK_", "", colnames(ss))
fwrite(data.table(sample = colnames(ss), as.data.table(t(ss))), file.path(RES, "22_immune_ssgsea.csv"))

df <- data.table(sample = colnames(mt), cluster = paste0("C", cl))
comp <- rbindlist(lapply(rownames(ss), function(s) {
  v <- ss[s, ]; a <- v[df$cluster == "C1"]; b <- v[df$cluster == "C2"]
  w <- wilcox.test(a, b)
  data.table(score = s, mean_C1 = mean(a), mean_C2 = mean(b), diff = mean(a) - mean(b), p = w$p.value)
}))
comp[, fdr := p.adjust(p, "BH")]; comp <- comp[order(p)]
fwrite(comp, file.path(RES, "23_cluster_immune_diff.csv"))
print(head(comp, 12))

## ---------- 簇与风险评分、分期、生存 ----------
sv <- readRDS(file.path(PROC, "survival.rds"))
d <- merge(df, sv$df[, .(sample = barcode, riskScore, OS.time, OS, stage_simple)], by = "sample")
cat("合并样本:", nrow(d), "\n")
print(tapply(d$riskScore, d$cluster, mean))
cat("风险评分 Wilcoxon P =", signif(wilcox.test(riskScore ~ cluster, d)$p.value, 3), "\n")
print(table(d$cluster, d$stage_simple, useNA = "ifany"))
sd <- survdiff(Surv(OS.time, OS) ~ cluster, data = d)
cat("簇间 log-rank P =", signif(1 - pchisq(sd$chisq, 1), 3), "\n")

fwrite(d, file.path(RES, "21_cluster_assignment.csv"))
saveRDS(list(cc = cc, k = K, cluster = cl, ss = ss, panels = panels), file.path(PROC, "cluster.rds"))
