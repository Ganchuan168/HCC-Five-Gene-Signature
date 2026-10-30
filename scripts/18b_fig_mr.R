#!/usr/bin/env Rscript
## 图7 cis-eQTL 孟德尔随机化（限定于 418 个候选基因，回文 SNP 已剔除，多 SNP 基因采用 IVW）
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(gridExtra)})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
FIG <- file.path(ROOT, "results/figures"); PROC <- file.path(ROOT, "data/processed")

sig <- readRDS(file.path(PROC, "signature_final.rds"))$genes
m <- fread(file.path(RES, "17b_mr_combined.csv"))
m[, logOR := log(OR)]
m <- m[is.finite(logOR) & is.finite(p) & p > 0]
m[, hi := symbol %in% sig]
lv <- m[tissue == "Liver"]; bd <- m[tissue == "Blood"]

mk_vol <- function(d, ttl) {
  d <- as.data.frame(d)
  ggplot(d, aes(x = logOR, y = -log10(p), colour = hi)) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25, color = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey40") +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"),
                        labels = c(`TRUE` = "Signature gene", `FALSE` = "Other candidate")) +
    labs(x = expression(log~OR), y = expression(-log[10]~"P"), colour = NULL, title = ttl) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
}
pA <- mk_vol(lv, sprintf("A  Liver cis-eQTL (%d genes)", length(unique(lv$symbol))))
pB <- mk_vol(bd, sprintf("B  Whole blood cis-eQTL (%d genes)", length(unique(bd$symbol))))

top <- m[order(p)][!duplicated(paste(symbol, tissue))][1:14]
top[, lab := sprintf("%s (%s)", symbol, tissue)]
top[, lab := factor(lab, levels = rev(lab))]
top[, sig := symbol %in% sig]
pC <- ggplot(top, aes(OR, lab, color = tissue)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, color = "grey40") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.18, linewidth = 0.35) +
  geom_point(aes(shape = sig), size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("Candidate", "Signature gene")) +
  scale_x_log10(limits = c(0.1, 20), oob = scales::squish) +
  scale_color_manual(values = c(Liver = "#B2182B", Blood = "#2166AC")) +
  labs(x = "Odds ratio (95% CI, log scale)", y = NULL, color = "eQTL tissue", shape = NULL,
       title = "C  Top nominal causal estimates") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

w <- dcast(m[, .(symbol, tissue, logOR)], symbol ~ tissue, value.var = "logOR")
w <- w[complete.cases(w)]
r <- cor(w$Liver, w$Blood, method = "spearman")
pD <- ggplot(w, aes(Liver, Blood)) +
  geom_point(size = 0.9, alpha = 0.7, color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#B2182B") +
  geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(x = "log OR (liver)", y = "log OR (whole blood)",
       title = sprintf("D  Cross-tissue concordance (rho = %.2f, n = %d)", r, nrow(w))) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

png(file.path(FIG, "Fig7_MR.png"), width = 11, height = 8.4, units = "in", res = 300)
grid.arrange(pA, pB, pC, pD, layout_matrix = rbind(c(1, 2), c(3, 4)))
dev.off()
cat("Fig7 done; rho =", round(r, 3), "n =", nrow(w), "\n")
fwrite(w, file.path(RES, "25_mr_cross_tissue.csv"))
