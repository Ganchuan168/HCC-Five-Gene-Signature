#!/usr/bin/env Rscript
## 图6 分子分型与免疫微环境
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(pheatmap); library(gridExtra)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")

cl <- readRDS(file.path(PROC, "cluster.rds"))
cc <- cl$cc; K <- cl$k; grp <- paste0("C", cl$cluster); ss <- cl$ss

## A 共识矩阵
cm <- cc[[K]]$consensusMatrix
ord <- order(grp)
ph <- pheatmap(cm[ord, ord], cluster_rows = FALSE, cluster_cols = FALSE,
               show_rownames = FALSE, show_colnames = FALSE,
               color = colorRampPalette(c("white", "#4393C3", "#053061"))(100),
               border_color = NA, silent = TRUE,
               main = sprintf("A  Consensus matrix (k = %d)", K))

## B CDF
cdf <- rbindlist(lapply(2:4, function(k) {
  m <- cc[[k]]$consensusMatrix
  data.table(k = paste0("k=", k), x = seq(0, 1, 0.02),
             y = sapply(seq(0, 1, 0.02), function(q) mean(m[lower.tri(m)] <= q)))
}))
kt <- fread(file.path(RES, "20_consensus_k.csv"))
pB <- ggplot(cdf, aes(x, y, color = k)) + geom_line(size = 0.6) +
  scale_color_manual(values = c("#2166AC", "#4DAF4A", "#B2182B")) +
  labs(x = "Consensus index", y = "CDF", color = NULL, title = "B  Cumulative distribution") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

## C KM
d <- fread(file.path(RES, "21_cluster_assignment.csv"))
fit <- survfit(Surv(OS.time, OS) ~ cluster, data = d)
pC <- ggsurvplot(fit, data = d, pval = TRUE, pval.size = 3, conf.int = FALSE,
                 palette = c("#2166AC", "#B2182B"), legend.labs = c("C1", "C2"),
                 legend.title = "Subtype", xlab = "Time (days)", ylab = "Overall survival",
                 ggtheme = theme_bw(base_size = 9), risk.table = FALSE, size = 0.6)
pC <- pC$plot + labs(title = "C  Survival by subtype") +
  theme(plot.title = element_text(size = 10, face = "bold"))

## D 风险评分
pD <- ggplot(d, aes(cluster, riskScore, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4, size = 0.3, width = 0.6) +
  scale_fill_manual(values = c(C1 = "#2166AC", C2 = "#B2182B")) +
  labs(x = NULL, y = "Risk score", title = "D  Risk score by subtype") +
  theme_bw(base_size = 9) + theme(legend.position = "none",
    plot.title = element_text(size = 10, face = "bold"))

## E 免疫评分
comp <- fread(file.path(RES, "23_cluster_immune_diff.csv"))
sel <- comp[order(p)][1:8]$score
m <- melt(as.data.table(t(ss[sel, , drop = FALSE])), variable.name = "score", value.name = "val")
m[, cluster := rep(grp, times = length(sel))]
pE <- ggplot(m, aes(score, val, fill = cluster)) +
  geom_boxplot(outlier.size = 0.2, size = 0.25, width = 0.7) +
  scale_fill_manual(values = c(C1 = "#2166AC", C2 = "#B2182B")) +
  labs(x = NULL, y = "ssGSEA score (z)", fill = "Subtype",
       title = "E  Immune infiltration between subtypes") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 10, face = "bold"))

png(file.path(FIG, "Fig6_subtype_immune.png"), width = 11, height = 8.6, units = "in", res = 300)
grid.arrange(grid::grobTree(ph$gtable), pB, pC, pD, pE,
             layout_matrix = rbind(c(1, 1, 2), c(3, 3, 4), c(5, 5, 5)),
             heights = c(1.15, 1, 0.9))
dev.off()
cat("Fig6 done\n")
