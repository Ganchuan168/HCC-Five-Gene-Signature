#!/usr/bin/env Rscript
## 图3 功能富集
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(gridExtra)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
FIG <- file.path(ROOT, "results/figures")

go <- fread(file.path(RES, "17_GO_BP.csv"))
go <- go[order(p.adjust)][1:12]
go[, Description := factor(Description, levels = rev(Description))]
pA <- ggplot(go, aes(Count, Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "#92C5DE", high = "#B2182B") +
  labs(x = "Gene count", y = NULL, color = expression(-log[10]~"FDR"),
       title = "A  GO biological process") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

kk <- fread(file.path(RES, "18_KEGG.csv"))
kk <- kk[order(p.adjust)][1:12]
kk[, Description := factor(Description, levels = rev(Description))]
pB <- ggplot(kk, aes(-log10(p.adjust), Description, fill = Count)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(low = "#92C5DE", high = "#B2182B") +
  labs(x = expression(-log[10]~"FDR"), y = NULL, title = "B  KEGG pathway") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

hm <- fread(file.path(RES, "19_hallmark_gsea.csv"))
hm <- hm[order(-NES)]
sel <- rbind(hm[1:8], hm[(.N - 7):.N])
sel[, ID := gsub("_", " ", ID)]
sel[, ID := factor(ID, levels = ID[order(NES)])]
sel[, dir := fifelse(NES > 0, "Activated in tumour", "Suppressed in tumour")]
pC <- ggplot(sel, aes(NES, ID, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_fill_manual(values = c("Activated in tumour" = "#B2182B",
                               "Suppressed in tumour" = "#2166AC")) +
  labs(x = "Normalised enrichment score (NES)", y = NULL, fill = NULL,
       title = "C  Hallmark pathways (GSEA, FDR < 0.25)") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "bottom")

png(file.path(FIG, "Fig3_enrichment.png"), width = 12, height = 7.2, units = "in", res = 300)
grid.arrange(pA, pB, pC, layout_matrix = rbind(c(1, 1, 2), c(3, 3, 3)),
             heights = c(1, 1.15))
dev.off()
cat("Fig3 done\n")
