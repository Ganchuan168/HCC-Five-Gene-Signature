#!/usr/bin/env Rscript
## 图5 风险评分分层与预后
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(pROC); library(gridExtra)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")

sv <- readRDS(file.path(PROC, "survival.rds"))
d <- as.data.table(sv$df)
d <- d[order(riskScore)]
d[, idx := .I]
d[, status := fifelse(OS == 1, "Dead", "Alive")]

pa <- ggplot(d, aes(idx, riskScore, color = riskGroup)) +
  geom_point(size = 0.7) +
  geom_vline(xintercept = sum(d$riskGroup == "Low") + 0.5, linetype = "dashed", size = 0.3) +
  scale_color_manual(values = c(High = "#B2182B", Low = "#2166AC")) +
  labs(x = NULL, y = "Risk score", color = NULL, title = "A  Risk score distribution") +
  theme_bw(base_size = 9) + theme(legend.position = "top",
    plot.title = element_text(size = 10, face = "bold"))

pb <- ggplot(d, aes(idx, OS.time, color = status)) +
  geom_point(size = 0.7) +
  geom_vline(xintercept = sum(d$riskGroup == "Low") + 0.5, linetype = "dashed", size = 0.3) +
  scale_color_manual(values = c(Dead = "#B2182B", Alive = "#2166AC")) +
  labs(x = "Patients (ranked by risk score)", y = "Survival time (days)", color = NULL,
       title = "B  Survival status") +
  theme_bw(base_size = 9) + theme(legend.position = "top",
    plot.title = element_text(size = 10, face = "bold"))

fit <- survfit(Surv(OS.time, OS) ~ riskGroup, data = d)
pcm <- ggsurvplot(fit, data = d, pval = TRUE, pval.size = 3, conf.int = FALSE,
                 palette = c("#B2182B", "#2166AC"), legend.labs = c("High risk", "Low risk"),
                 legend.title = "", xlab = "Time (days)", ylab = "Overall survival",
                 ggtheme = theme_bw(base_size = 9), risk.table = TRUE, risk.table.height = 0.26,
                 break.time.by = 365, tables.theme = theme_cleantable(), size = 0.6)
pc <- pcm$plot + labs(title = "C  Kaplan-Meier (TCGA-LIHC)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

tr <- sv$roc
rocd <- data.table(FPR = 1 - tr$FP[, 1], TPR = tr$TP[, 1], year = "1-year")
rocd <- rbind(rocd, data.table(FPR = 1 - tr$FP[, 2], TPR = tr$TP[, 2], year = "3-year"))
pd <- ggplot(rocd, aes(FPR, TPR, color = year)) + geom_line(size = 0.7) +
  geom_abline(linetype = "dashed", color = "grey60", size = 0.3) + coord_equal() +
  scale_color_manual(values = c("#2166AC", "#B2182B"),
                     labels = sprintf("%s (AUC=%.3f)", c("1-year", "3-year"), tr$AUC)) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL,
       title = "D  Time-dependent ROC") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA))

cox <- fread(file.path(RES, "09_cox_multivariate.csv"))
cox[, variable := c(riskScore = "Risk score", `stage_simpleIII-IV` = "Stage III-IV",
                    age = "Age", gendermale = "Male sex")[variable]]
cox[, variable := factor(variable, levels = rev(variable))]
pe <- ggplot(cox, aes(HR, variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.3, color = "grey40") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.15, size = 0.4) +
  geom_point(size = 1.8, color = "#B2182B") +
  scale_x_log10() +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), P%s", HR, low, high,
        ifelse(p < 0.001, "<0.001", paste0("=", formatC(p, format = "f", digits = 3))))),
        vjust = -0.9, size = 2.4) +
  labs(x = "Hazard ratio (log scale)", y = NULL,
       title = "E  Multivariable Cox regression") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 10, face = "bold"))

png(file.path(FIG, "Fig5_survival.png"), width = 11, height = 8.5, units = "in", res = 300)
if (!is.null(pcm$table)) {
  cplot <- arrangeGrob(grobs = list(ggplotGrob(pc), ggplotGrob(pcm$table)), heights = c(0.76, 0.24))
} else {
  cplot <- ggplotGrob(pc)
}
grid.arrange(pa, pb, cplot, pd, pe,
             layout_matrix = rbind(c(1, 1, 3, 3), c(2, 2, 3, 3), c(4, 4, 5, 5)),
             heights = c(1, 1, 1.1))
dev.off()
cat("Fig5 done\n")
