#!/usr/bin/env Rscript
## 10 判别效能补充：AFP 对照、校准与决策曲线分析
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

sig <- readRDS(file.path(PROC, "signature_final.rds"))
genes <- sig$genes; co <- coef(sig$model)

disc <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
sets <- list(
  "GSE14520" = list(expr = disc$expr[, ], group = disc$group),
  "TCGA-LIHC" = { d <- readRDS(file.path(PROC, "TCGA_LIHC.rds")); list(expr = d$expr, group = d$group) },
  "GSE76427" = { d <- readRDS(file.path(PROC, "GSE76427.rds")); list(expr = d$expr, group = d$group) },
  "GSE62232" = { d <- readRDS(file.path(PROC, "GSE62232.rds")); list(expr = d$expr, group = d$group) }
)

res <- list(); prob <- list()
for (nm in names(sets)) {
  e <- sets[[nm]]; y <- as.integer(e$group == "Tumor")
  ## 特征评分（队列内标准化）
  z <- t(scale(t(e$expr[genes, ]))); z[is.na(z)] <- 0
  lp <- as.numeric(t(z) %*% co[genes]) + co["(Intercept)"]
  p <- plogis(lp)
  ## AFP 单基因
  af <- as.numeric(scale(as.numeric(e$expr["AFP", ])))
  af[is.na(af)] <- 0
  ra <- roc(y, af, quiet = TRUE); rs <- roc(y, p, quiet = TRUE)
  dl <- roc.test(rs, ra, method = "delong")
  ## 校准：logistic 重校准的斜率与截距
  cal <- suppressWarnings(glm(y ~ lp, family = binomial))          # 自由斜率：校准斜率
  citl <- suppressWarnings(coef(glm(y ~ offset(lp), family = binomial))[1])  # 斜率固定为 1：校准截距
  brier <- mean((p - y)^2)
  ## 十分位校准
  q <- cut(rank(p, ties.method = "first"), breaks = 10, labels = FALSE)
  cd <- data.table(cohort = nm, decile = q, p = p, y = y)[, .(pred = mean(p), obs = mean(y), n = .N), by = .(cohort, decile)]
  prob[[nm]] <- data.table(cohort = nm, y = y, p = p, afp = plogis(coef(glm(y ~ af, family = binomial))[1] +
                                                                      coef(glm(y ~ af, family = binomial))[2] * af))
  res[[nm]] <- data.table(cohort = nm, n = length(y), auc_sig = as.numeric(auc(rs)),
                          auc_afp = as.numeric(auc(ra)), p_delong = dl$p.value,
                          cal_slope = coef(cal)[2], cal_citl = citl,
                          brier = brier, mean_pred = mean(p), obs_rate = mean(y))
  fwrite(cd, file.path(RES, paste0("30_calib_", gsub("[^A-Za-z0-9]", "", nm), ".csv")))
  cat(sprintf("%-10s AUC 特征=%.3f AFP=%.3f DeLong P=%.3g 校准斜率=%.2f CITL=%.2f Brier=%.3f 平均预测=%.2f 实际=%.2f\n",
              nm, as.numeric(auc(rs)), as.numeric(auc(ra)), dl$p.value, coef(cal)[2], citl,
              brier, mean(p), mean(y)))
}
res <- rbindlist(res); fwrite(res, file.path(RES, "31_afp_calibration.csv"))
saveRDS(prob, file.path(PROC, "prob.rds"))

## ---------- 决策曲线（净获益）----------
nb <- function(y, p, pt) {
  tp <- mean(p >= pt & y == 1); fp <- mean(p >= pt & y == 0)
  tp - fp * pt / (1 - pt)
}
ths <- seq(0.01, 0.9, by = 0.01)
for (nm in c("GSE14520", "TCGA-LIHC", "GSE76427")) {
  d <- prob[[nm]]
  dd <- data.table(threshold = ths,
                   signature = sapply(ths, function(t) nb(d$y, d$p, t)),
                   AFP = sapply(ths, function(t) nb(d$y, d$afp, t)),
                   all = sapply(ths, function(t) mean(d$y) - (1 - mean(d$y)) * t / (1 - t)))
  fwrite(dd, file.path(RES, paste0("32_dca_", gsub("[^A-Za-z0-9]", "", nm), ".csv")))
}
cat("saved DCA\n")
