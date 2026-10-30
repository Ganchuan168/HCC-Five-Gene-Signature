#!/usr/bin/env Rscript
## 09 TCGA：特征评分、预后、免疫微环境
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(survival); library(survminer); library(timeROC)
  library(GSVA); library(msigdbr); library(ggplot2)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

tc <- readRDS(file.path(PROC, "TCGA_LIHC.rds"))
sig <- readRDS(file.path(PROC, "signature.rds"))
genes <- sig$genes
cat("特征基因:", paste(genes, collapse = ", "), "\n")

m <- tc$expr
tum <- tc$group == "Tumor"
mt <- m[, tum]
cat("肿瘤样本:", sum(tum), "\n")

## 风险评分：z 标准化后按逻辑回归系数加权
z <- t(scale(t(mt[genes, ])))
coefs <- coef(sig$model)[-1]
risk <- as.numeric(t(z) %*% coefs[names(coefs)])
names(risk) <- colnames(mt)
cat("风险评分范围:", round(range(risk), 3), "\n")

## 结合临床与生存
clin <- as.data.table(tc$clin)
barc <- colnames(mt); pid <- substr(barc, 1, 12)
df <- clin[data.table(barcode = barc, patient = pid, riskScore = risk), on = "patient", nomatch = 0]
df <- df[!is.na(OS.time) & !is.na(OS) & OS.time > 0]
df[, stage_simple := fifelse(grepl("^Stage (I|II)[AB]?$", stage), "I-II",
                     fifelse(grepl("^Stage (III|IV)", stage), "III-IV", NA))]
cat("可分析样本:", nrow(df), " 事件数:", sum(df$OS == 1), "\n")

df[, riskGroup := fifelse(riskScore > median(riskScore), "High", "Low")]
fwrite(df, file.path(RES, "08_tcga_riskscore.csv"))

## KM
fit <- survfit(Surv(OS.time, OS) ~ riskGroup, data = df)
sd <- survdiff(Surv(OS.time, OS) ~ riskGroup, data = df)
p_lr <- 1 - pchisq(sd$chisq, df = 1)
cat("Log-rank P =", signif(p_lr, 3), "\n")
print(survdiff(Surv(OS.time, OS) ~ riskGroup, data = df))

## Cox
cox1 <- coxph(Surv(OS.time, OS) ~ riskScore, data = df)
cox2 <- coxph(Surv(OS.time, OS) ~ riskScore + stage_simple + age + gender, data = df)
s1 <- summary(cox1); s2 <- summary(cox2)
cat("单因素 HR:", round(s1$conf.int[1, 1], 3), "(", round(s1$conf.int[1, 3], 3), "-",
    round(s1$conf.int[1, 4], 3), ") P =", signif(s1$coefficients[1, 5], 3), "\n")
print(s2$conf.int)
fwrite(data.table(variable = rownames(s2$conf.int), HR = s2$conf.int[, 1],
                  low = s2$conf.int[, 3], high = s2$conf.int[, 4],
                  p = s2$coefficients[, 5]), file.path(RES, "09_cox_multivariate.csv"))

## 时间依赖 ROC
tr <- timeROC(T = df$OS.time, delta = df$OS, marker = df$riskScore,
              cause = 1, times = c(365, 1095, 1825), iid = TRUE)
print(tr$AUC)

## 免疫检查点
ckpt <- c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD27", "CD40", "IDO1", "LGALS9")
ck <- intersect(ckpt, rownames(mt))
ckdf <- data.table(gene = ck, rho = sapply(ck, function(g) cor(mt[g, ], risk, method = "spearman")),
                   p = sapply(ck, function(g) cor.test(mt[g, ], risk, method = "spearman")$p.value))
ckdf <- ckdf[order(-abs(rho))]
fwrite(ckdf, file.path(RES, "10_checkpoint_correlation.csv"))
print(ckdf)

saveRDS(list(df = df, risk = risk, km = fit, p_lr = p_lr, roc = tr), file.path(PROC, "survival.rds"))
