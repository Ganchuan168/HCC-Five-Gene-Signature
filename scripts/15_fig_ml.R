#!/usr/bin/env Rscript
## 图4 机器学习特征筛选、诊断效能与枢纽基因表达
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pROC); library(gridExtra); library(reshape2)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")

sig <- readRDS(file.path(PROC, "signature.rds"))
genes <- sig$genes; model <- sig$model
cat("特征基因:", paste(genes, collapse = ", "), "\n")

## ---------- A 12 种算法排名热图（前 20 基因）----------
rk <- sig$ranks
ful <- names(rk)[sapply(rk, length) >= 300]   # 仅纳入给出完整基因排序的算法
rkf <- rk[ful]
allg <- unique(unlist(lapply(rkf, function(r) head(r, 20))))
allg <- allg[!is.na(allg)]
sc <- sapply(rkf, function(r) match(allg, r))
rownames(sc) <- allg
sc[is.na(sc)] <- 25
avg <- rowMeans(sc)
top <- names(sort(avg))[1:20]
mdf <- as.data.table(melt(sc[top, ]))
setnames(mdf, c("gene", "algorithm", "rank"))
mdf[, algorithm := factor(algorithm, levels = names(rkf))]
pA <- ggplot(mdf, aes(algorithm, gene, fill = rank)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient(low = "#B2182B", high = "#F7F7F7", trans = "reverse",
                      name = "Rank") +
  labs(x = NULL, y = NULL, title = sprintf("A  Consensus ranking across %d algorithms (top 20)", length(rkf))) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10, face = "bold"))

## ---------- B 四队列 ROC ----------
tc <- readRDS(file.path(PROC, "TCGA_LIHC.rds"))
coh <- list(`GSE14520 (training)` = readRDS(file.path(PROC, "discovery_wgcna.rds")),
            `TCGA-LIHC` = tc,
            `GSE76427` = readRDS(file.path(PROC, "GSE76427.rds")),
            `GSE62232` = readRDS(file.path(PROC, "GSE62232.rds")))
rocdf <- list(); aucs <- list()
for (nm in names(coh)) {
  e <- coh[[nm]]; gs <- intersect(genes, rownames(e$expr))
  Xe <- scale(t(e$expr[gs, , drop = FALSE])); ye <- as.integer(e$group == "Tumor")
  full <- matrix(0, nrow = nrow(Xe), ncol = length(genes), dimnames = list(NULL, genes))
  for (g in genes) full[, g] <- if (g %in% gs) Xe[, g] else 0
  pr <- predict(model, newdata = data.frame(full, check.names = FALSE), type = "response")
  ok <- !is.na(pr); r <- roc(ye[ok], pr[ok], quiet = TRUE)
  rocdf[[nm]] <- data.table(FPR = 1 - r$specificities, TPR = r$sensitivities, cohort = nm)
  aucs[[nm]] <- data.table(cohort = nm, auc = as.numeric(auc(r)), n = sum(ok))
}
rocdf <- rbindlist(rocdf); aucs <- rbindlist(aucs)
aucs[, label := sprintf("%s (AUC=%.3f)", cohort, auc)]
lbl <- setNames(aucs$label, aucs$cohort)
pB <- ggplot(rocdf, aes(FPR, TPR, color = cohort)) +
  geom_line(size = 0.8) + geom_abline(linetype = "dashed", color = "grey60", size = 0.3) +
  scale_color_manual(values = c("#B2182B", "#2166AC", "#4DAF4A", "#984EA3"), labels = lbl) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity", color = NULL,
                       title = "B  Diagnostic ROC in four cohorts") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 7))

## ---------- C 枢纽基因表达（四队列）----------
box <- list()
for (nm in names(coh)) {
  e <- coh[[nm]]; gs <- intersect(genes, rownames(e$expr))
  z <- e$expr[gs, , drop = FALSE]
  z <- t(scale(t(z)))
  box[[nm]] <- data.table(gene = rep(gs, times = ncol(z)),
                          value = as.numeric(z),
                          group = rep(e$group, each = length(gs)), cohort = nm)
}
box <- rbindlist(box)
box[, group := factor(group, levels = c("Normal", "Tumor"))]
pC <- ggplot(box, aes(gene, value, fill = group)) +
  geom_boxplot(outlier.size = 0.25, size = 0.25, width = 0.7) +
  facet_wrap(~cohort, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "#4DAF4A", Tumor = "#E41A1C")) +
  labs(x = NULL, y = "Expression (z-score)", fill = NULL,
       title = "C  Hub gene expression across cohorts") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

png(file.path(FIG, "Fig4_ML_signature.png"), width = 11, height = 9, units = "in", res = 300)
grid.arrange(pA, pB, pC, layout_matrix = rbind(c(1, 2), c(3, 3)),
             heights = c(1.05, 1))
dev.off()
fwrite(aucs, file.path(RES, "24_roc_auc_all.csv"))
cat("Fig4 done\n"); print(aucs)
