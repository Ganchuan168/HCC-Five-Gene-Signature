#!/usr/bin/env Rscript
## 图4 多算法共识排序、五基因特征的判别效能、校准与决策曲线
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pROC); library(gridExtra); library(reshape2)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
PROC <- file.path(ROOT, "data/processed"); FIG <- file.path(ROOT, "results/figures")

sig <- readRDS(file.path(PROC, "signature_final.rds"))
genes <- sig$genes; model <- sig$model; co <- coef(model)
cat("特征基因:", paste(genes, collapse = ", "), "\n")

## ---------- A 共识排序热图 ----------
rk <- sig$ranks
ful <- names(rk)[sapply(rk, length) >= 300]
rkf <- rk[ful]
allg <- unique(unlist(lapply(rkf, function(r) head(r, 20))))
allg <- allg[!is.na(allg)]
sc <- sapply(rkf, function(r) match(allg, r)); rownames(sc) <- allg
sc[is.na(sc)] <- 25
top <- names(sort(rowMeans(sc)))[1:20]
mdf <- as.data.table(melt(sc[top, ]))
setnames(mdf, c("gene", "algorithm", "rank"))
mdf[, algorithm := factor(algorithm, levels = names(rkf))]
pA <- ggplot(mdf, aes(algorithm, factor(gene, levels = rev(top)), fill = rank)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#B2182B", high = "#F7F7F7", trans = "reverse", name = "Rank") +
  labs(x = NULL, y = NULL,
       title = sprintf("A  Consensus ranking across %d algorithms (top 20 genes)", length(rkf))) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10, face = "bold"))

## ---------- B ROC ----------
coh <- list(`GSE14520` = readRDS(file.path(PROC, "discovery_wgcna.rds")),
            `TCGA-LIHC` = readRDS(file.path(PROC, "TCGA_LIHC.rds")),
            `GSE76427` = readRDS(file.path(PROC, "GSE76427.rds")),
            `GSE62232` = readRDS(file.path(PROC, "GSE62232.rds")))
rocdf <- list(); aucs <- list()
for (nm in names(coh)) {
  e <- coh[[nm]]
  z <- t(scale(t(e$expr[genes, , drop = FALSE]))); z[is.na(z)] <- 0
  lp <- as.numeric(t(z) %*% co[genes]) + co["(Intercept)"]
  ye <- as.integer(e$group == "Tumor")
  r <- roc(ye, lp, quiet = TRUE)
  rocdf[[nm]] <- data.table(FPR = 1 - r$specificities, TPR = r$sensitivities, cohort = nm)
  aucs[[nm]] <- data.table(cohort = nm, auc = as.numeric(auc(r)), n = length(ye))
}
rocdf <- rbindlist(rocdf); aucs <- rbindlist(aucs)
aucs[, label := sprintf("%s (AUC=%.3f)", cohort, auc)]
lbl <- setNames(aucs$label, aucs$cohort)
pB <- ggplot(rocdf, aes(FPR, TPR, color = cohort)) +
  geom_line(linewidth = 0.8) + geom_abline(linetype = "dashed", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("#B2182B", "#2166AC", "#4DAF4A", "#984EA3"), labels = lbl) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity", color = NULL,
                       title = "B  Discrimination in four cohorts") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 7))

## ---------- C 表达 ----------
box <- list()
for (nm in names(coh)) {
  e <- coh[[nm]]
  z <- t(scale(t(e$expr[genes, , drop = FALSE])))
  box[[nm]] <- data.table(gene = rep(genes, times = ncol(z)), value = as.numeric(z),
                          group = rep(e$group, each = length(genes)), cohort = nm)
}
box <- rbindlist(box)
box[, group := factor(group, levels = c("Normal", "Tumor"), labels = c("Adjacent", "Tumour"))]
box[, cohort := factor(cohort, levels = names(coh))]
pC <- ggplot(box, aes(gene, value, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.7) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_manual(values = c(Adjacent = "#4DAF4A", Tumour = "#E41A1C")) +
  labs(x = NULL, y = "Expression (cohort-wise z-score)", fill = NULL,
       title = "C  Signature gene expression") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

## ---------- D 校准曲线 ----------
cal <- rbindlist(lapply(c("GSE14520", "TCGALIHC", "GSE76427"), function(tag) {
  f <- file.path(RES, paste0("30_calib_", tag, ".csv"))
  if (file.exists(f)) fread(f) else NULL
}))
cal <- cal[cohort != "GSE14520"]
cal[, cohort := factor(cohort, levels = c("TCGA-LIHC", "GSE76427"))]
cal[, `:=`(lo = qbeta(0.025, obs * n + 1, (1 - obs) * n + 1),
           hi = qbeta(0.975, obs * n + 1, (1 - obs) * n + 1))]
pD <- ggplot(cal, aes(pred, obs, color = cohort)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.35) +
  geom_point(size = 0.9) + geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("#2166AC", "#4DAF4A")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability (decile mean)", y = "Observed proportion",
       color = NULL, title = "D  Calibration in external cohorts") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 7))

## ---------- E 决策曲线 ----------
pr <- readRDS(file.path(PROC, "prob.rds"))
ths <- seq(0.02, 0.8, by = 0.02)
nb <- function(y, p, t) {
  tp <- mean(p >= t & y == 1); fp <- mean(p >= t & y == 0)
  tp - fp * t / (1 - t)
}
dca <- rbindlist(lapply(c("GSE14520", "TCGA-LIHC"), function(nm) {
  d <- pr[[nm]]
  lp <- qlogis(pmin(pmax(d$p, 1e-6), 1 - 1e-6))
  rec <- plogis(coef(glm(d$y ~ lp, family = binomial))[1] + coef(glm(d$y ~ lp, family = binomial))[2] * lp)
  data.table(cohort = nm, threshold = ths,
             signature = sapply(ths, function(t) nb(d$y, rec, t)),
             AFP = sapply(ths, function(t) nb(d$y, d$afp, t)),
             all = sapply(ths, function(t) mean(d$y) - (1 - mean(d$y)) * t / (1 - t)))
}))
dcal <- as.data.table(melt(dca, id.vars = c("cohort", "threshold"), variable.name = "model", value.name = "nb"))
dcal[, cohort := factor(cohort, levels = c("GSE14520", "TCGA-LIHC"))]
pE <- ggplot(dcal, aes(threshold, nb, color = model, linetype = model)) +
  geom_line(linewidth = 0.6) + facet_wrap(~cohort) +
  scale_color_manual(values = c(signature = "#B2182B", AFP = "#2166AC", all = "grey50"),
                     labels = c("Five-gene signature", "AFP", "Treat all")) +
  scale_linetype_manual(values = c(signature = "solid", AFP = "solid", all = "dashed"),
                        labels = c("Five-gene signature", "AFP", "Treat all")) +
  coord_cartesian(ylim = c(-0.05, 1)) +
  labs(x = "Threshold probability", y = "Net benefit", color = NULL, linetype = NULL,
       title = "E  Decision curve analysis") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "bottom",
        legend.text = element_text(size = 7))
fwrite(dca, file.path(RES, "32_dca_figure.csv"))

png(file.path(FIG, "Fig4_ML_signature.png"), width = 11, height = 13, units = "in", res = 300)
grid.arrange(pA, pB, pC, pD, pE, layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 5)),
             heights = c(1.15, 1, 1))
dev.off()
fwrite(aucs, file.path(RES, "24_roc_auc_all.csv"))
cat("Fig4 done\n"); print(aucs)
