#!/usr/bin/env Rscript
## 图2 差异表达与 WGCNA
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pheatmap); library(grid); library(gridExtra)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

## ---------- A 火山图 ----------
deg <- fread(file.path(RES, "01_DEG_all.csv"))
deg[, sig := fifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                     fifelse(logFC > 1, "Up", "Down"), "NS")]
deg[, sig := factor(sig, levels = c("Up", "Down", "NS"))]
lab <- deg[sig != "NS"][order(-abs(logFC))][1:12]
pA <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey75")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", size = 0.3, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3, color = "grey40") +
  geom_text(data = lab, aes(label = gene), size = 2.4, vjust = -0.4, color = "black") +
  labs(x = expression(log[2]~"fold change"), y = expression(-log[10]~"adjusted P"),
       title = "A  Differential expression (GSE14520)") +
  theme_bw(base_size = 9) + theme(legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0), legend.title = element_blank(),
    plot.title = element_text(size = 10, face = "bold"))

## ---------- B 热图（前 50 基因）----------
d <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
top <- deg[sig != "NS"][order(adj.P.Val)][1:50]$gene
top <- intersect(top, rownames(d$expr))
hm <- d$expr[top, ]
hm <- t(scale(t(hm)))
ann <- data.frame(Group = d$group); rownames(ann) <- colnames(hm)
ann$Group <- ifelse(ann$Group == "Tumor", "Tumour", "Adjacent")
ph <- pheatmap(hm, annotation_col = ann, show_colnames = FALSE, cluster_cols = TRUE,
               cluster_rows = TRUE, fontsize_row = 6, treeheight_row = 15,
               treeheight_col = 10, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               annotation_colors = list(Group = c(Adjacent = "#4DAF4A", Tumour = "#E41A1C")),
               silent = TRUE, main = "B  Top 50 DEGs")

## ---------- C 软阈值 ----------
sf <- fread(file.path(RES, "03_soft_threshold.csv"))
pC <- ggplot(sf, aes(Power, SFT.R.sq)) +
  geom_point(size = 1.4, color = "#B2182B") + geom_line(color = "grey40", size = 0.3) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "blue", size = 0.3) +
  geom_text(aes(label = Power), size = 2, vjust = -0.7) +
  scale_x_continuous(breaks = sf$Power) +
  labs(x = "Soft threshold (power)", y = expression("Scale-free topology"~R^2),
       title = "C  Soft threshold") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

## ---------- D 模块-性状相关（肿瘤 vs 癌旁）----------
mt <- fread(file.path(RES, "04_module_trait.csv"))
mt[, module := sub("^ME", "", module)]
mt[, stars := ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))]
mt[, module := factor(module, levels = module[order(cor)])]
pD <- ggplot(mt, aes(module, cor, fill = cor)) +
  geom_col(width = 0.7, color = "grey30", size = 0.2) +
  geom_text(aes(label = stars), vjust = ifelse(mt$cor > 0, -0.4, 1.2), size = 3) +
  coord_flip() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
  labs(x = NULL, y = "Module-trait correlation (tumour vs. adjacent)",
       title = "D  Module-trait relationships") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"),
                                  legend.position = "none")

png(file.path(FIG, "Fig2_DEG_WGCNA.png"), width = 11, height = 8, units = "in", res = 300)
grid.arrange(pA, pC, grid::grobTree(ph$gtable), pD, ncol = 2,
             widths = c(1, 1.15), heights = c(1, 1))
dev.off()
cat("Fig2 done\n")
