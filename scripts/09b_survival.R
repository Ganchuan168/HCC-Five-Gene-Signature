#!/usr/bin/env Rscript
## 09b TCGA-LIHC：五基因风险评分的预后价值（含缺失分期处理与敏感性分析）
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC); library(pROC)
})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

tc <- readRDS(file.path(PROC, "TCGA_LIHC.rds"))
sig <- readRDS(file.path(PROC, "signature_final.rds"))
genes <- sig$genes
coefs <- coef(sig$model)[-1]
cat("特征基因:", paste(genes, collapse = ", "), "\n")

m <- tc$expr; tum <- tc$group == "Tumor"; mt <- m[, tum]
z <- t(scale(t(mt[genes, ])))
risk <- as.numeric(t(z) %*% coefs[genes])
names(risk) <- colnames(mt)
cat("风险评分范围:", round(range(risk), 2), " SD:", round(sd(risk), 2), "\n")

clin <- as.data.table(tc$clin)
barc <- colnames(mt); pid <- substr(barc, 1, 12)
df <- clin[data.table(barcode = barc, patient = pid, riskScore = risk), on = "patient", nomatch = 0]
df <- df[!is.na(OS.time) & !is.na(OS) & OS.time > 0]
df[, stage_simple := fifelse(grepl("^Stage (I|II)[ABC]?$", stage), "I-II",
                     fifelse(grepl("^Stage (III|IV)", stage), "III-IV", NA))]
df[, stage_cat := fifelse(is.na(stage_simple), "Unknown", stage_simple)]
cat("可分析样本:", nrow(df), " 事件:", sum(df$OS == 1), "\n")
print(df[, .N, by = stage_cat])
cat("中位随访(天):", median(df$OS.time), " 范围:", range(df$OS.time), "\n")

## 随访时间：以事件者中位生存时间与 censored 者随访时间共同描述
cat("死亡者中位生存:", median(df[OS == 1, OS.time]), " 存活者中位随访:", median(df[OS == 0, OS.time]), "\n")

df[, riskGroup := fifelse(riskScore > median(riskScore), "High", "Low")]
df[, riskSD := riskScore / sd(riskScore)]
fwrite(df, file.path(RES, "08_tcga_riskscore.csv"))

## KM
sd1 <- survdiff(Surv(OS.time, OS) ~ riskGroup, data = df)
p_lr <- 1 - pchisq(sd1$chisq, df = 1)
cat("Log-rank P =", signif(p_lr, 3), "\n")

## Cox
f0 <- coxph(Surv(OS.time, OS) ~ riskScore, data = df)
f1 <- coxph(Surv(OS.time, OS) ~ riskScore + age + gender, data = df)
f2 <- coxph(Surv(OS.time, OS) ~ riskScore + stage_simple + age + gender, data = df[!is.na(stage_simple)])
f3 <- coxph(Surv(OS.time, OS) ~ riskScore + stage_cat + age + gender, data = df)
f4 <- coxph(Surv(OS.time, OS) ~ riskSD, data = df)
f5 <- coxph(Surv(OS.time, OS) ~ riskGroup, data = df)
for (nm in c("f0", "f1", "f2", "f3", "f4", "f5")) {
  f <- get(nm); s <- summary(f)
  cat("\n==", nm, " n =", f$n, " events =", f$nevent, " C-index =", round(s$concordance[1], 4), "\n")
  print(round(s$conf.int[, c(1, 3, 4)], 4)); print(signif(s$coefficients[, 5], 3))
}
s2 <- summary(f2)
fwrite(data.table(variable = rownames(s2$conf.int), HR = s2$conf.int[, 1],
                  low = s2$conf.int[, 3], high = s2$conf.int[, 4],
                  p = s2$coefficients[, 5], n = f2$n, events = f2$nevent,
                  cindex = s2$concordance[1]),
       file.path(RES, "09_cox_multivariate.csv"))

## 时间依赖 ROC（随访中位数约 1.6 年，5 年估计不稳定，主分析取 1 年与 3 年）
tr <- timeROC(T = df$OS.time, delta = df$OS, marker = df$riskScore, cause = 1,
              times = c(365, 1095), iid = TRUE)
cat("\n时间依赖 AUC:\n"); print(round(tr$AUC, 3)); if (!is.null(tr$AUC.se)) print(round(tr$AUC.se, 3))

## 免疫检查点
ckpt <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD27", "CD40", "IDO1", "LGALS9")
ck <- intersect(ckpt, rownames(mt))
ckdf <- data.table(gene = ck,
                   rho = sapply(ck, function(g) cor(mt[g, ], risk, method = "spearman")),
                   p = sapply(ck, function(g) cor.test(mt[g, ], risk, method = "spearman")$p.value))
ckdf[, fdr := p.adjust(p, "BH")]
setorder(ckdf, -rho)
fwrite(ckdf, file.path(RES, "10_checkpoint_correlation.csv"))
print(ckdf)

saveRDS(list(df = df, risk = risk, p_lr = p_lr, roc = tr, fits = list(f0 = f0, f1 = f1, f2 = f2, f3 = f3, f4 = f4, f5 = f5)),
        file.path(PROC, "survival.rds"))
cat("saved survival.rds\n")
