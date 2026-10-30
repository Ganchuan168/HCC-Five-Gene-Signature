#!/usr/bin/env Rscript
## 08d 依 1-SE 准则确定特征数，拟合最终五基因模型并保存
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

sf <- readRDS(file.path(PROC, "signature_full.rds"))
res <- sf$res

## 1-SE 准则：AUC 在最优值 1 个标准误以内的最简约 k
best <- res[1]
se1 <- best$cv_sd / sqrt(5)
ok <- res[cv_auc >= best$cv_auc - se1]
k_sel <- min(ok$k)
cat("最优 AUC:", round(best$cv_auc, 5), " 1SE:", round(se1, 5), " 选定 k =", k_sel, "\n")

## 在 1-SE 范围内按“跨组合入选频率”打破平局，降低对单一算法组合的依赖
cand <- res[k == k_sel & cv_auc >= best$cv_auc - se1]
freq <- table(unlist(strsplit(res$genes, ",")))
cand[, stab := sapply(genes, function(s) sum(freq[strsplit(s, ",")[[1]]]))]
setorder(cand, -stab, -cv_auc)
print(head(cand[, .(combo, k, cv_auc, stab, genes)], 6))
sel <- cand[1]
genes <- strsplit(sel$genes, ",")[[1]]
cat("组合:", sel$combo, " k =", sel$k, "\n特征基因:", paste(genes, collapse = ", "), "\n")
cat("该组合跨组合出现频率之和:", sel$stab, "\n")
cat("该组合 CV AUC:", round(sel$cv_auc, 5), "±", round(sel$cv_sd, 5), "\n")

disc <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
X <- t(disc$expr[genes, ]); y <- as.integer(disc$group == "Tumor")
Xz <- scale(X)
mu <- attr(Xz, "scaled:center"); sdv <- attr(Xz, "scaled:scale")
fit <- glm(y ~ ., data = data.frame(y = y, Xz, check.names = FALSE), family = binomial)
print(round(coef(fit), 4))
pr <- predict(fit, type = "response")
r <- roc(y, pr, quiet = TRUE)
thr <- as.numeric(coords(r, "best", ret = "threshold", transpose = FALSE)[1, 1])
cat("表观 AUC:", round(as.numeric(auc(r)), 4), " 阈值:", round(thr, 4), "\n")

## 各基因在肿瘤/癌旁间的平均差（log2 尺度）
delta <- sapply(genes, function(g) mean(disc$expr[g, y == 1]) - mean(disc$expr[g, y == 0]))
cat("肿瘤-癌旁均值差:\n"); print(round(delta, 3))

saveRDS(list(genes = genes, model = fit, mu = mu, sd = sdv, threshold = thr,
             combo = sel$combo, k = sel$k, cv_auc = sel$cv_auc, cv_sd = sel$cv_sd, stability = sel$stab,
             apparent_auc = as.numeric(auc(r)), delta = delta,
             res = res, cons = sf$cons, ranks = sf$ranks,
             optimism = sf$optimism, auc_corrected = sf$auc_corrected),
        file.path(PROC, "signature_final.rds"))
cat("saved signature_final.rds\n")
