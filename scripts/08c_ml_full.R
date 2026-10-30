#!/usr/bin/env Rscript
## 08c 在全部发现集样本（233 例患者、445 个样本）上进行多算法共识排序，
## 并以患者分组的 5 折交叉验证选择特征数目 k，避免原先随机划分导致的乐观偏倚。
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(randomForest); library(e1071)
  library(xgboost); library(gbm); library(MASS); library(pROC)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")
source(file.path(ROOT, "scripts/ml_ranks.R"))

disc <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
cand <- disc$genes
X <- t(disc$expr[cand, ]); y <- as.integer(disc$group == "Tumor")
d0 <- readRDS(file.path(PROC, "GSE14520.rds"))
pid <- sub("[AB]$", "", sub(".*(LCS-[0-9]+).*", "\\1", d0$title))
stopifnot(length(pid) == nrow(X))
cat("候选基因:", ncol(X), " 样本:", nrow(X), " 患者:", length(unique(pid)), "\n")

## ---------- 1. 全部发现集的共识排序 ----------
Xtr <- X; ytr <- y
ranks <- list()
for (nm in names(algos)) {
  r <- tryCatch(algos[[nm]](Xtr, ytr), error = function(e) NULL)
  if (is.null(r) || length(r) < 5) { cat("  跳过", nm, "\n"); next }
  ranks[[nm]] <- r
  cat("  ", nm, "->", length(r), "基因\n")
}
rk <- sapply(ranks, function(r) match(colnames(X), r))
rownames(rk) <- colnames(X)
rk[is.na(rk)] <- ncol(X) + 1

## 每个基因在 10 种算法中的平均名次（越小越靠前）
cons <- data.table(gene = rownames(rk), mean_rank = rowMeans(rk),
                   n_top10 = rowSums(rk <= 10), n_top20 = rowSums(rk <= 20))
setorder(cons, mean_rank)
fwrite(cons, file.path(RES, "06a_consensus_rank_full.csv"))

## ---------- 2. 患者分组 5 折交叉验证选择 k ----------
set.seed(20261029)
upid <- unique(pid)
fold_of_pid <- integer(length(upid)); names(fold_of_pid) <- upid
fold_of_pid[sample(upid)] <- rep(1:5, length.out = length(upid))
fold <- fold_of_pid[pid]

combos <- c(as.list(names(ranks)), combn(names(ranks), 2, simplify = FALSE))
ks <- c(3, 4, 5, 6, 8, 10, 15, 20)
res <- list()
for (cb in combos) {
  sc <- rowMeans(rk[, cb, drop = FALSE])
  ord <- names(sort(sc))
  for (k in ks) {
    g <- ord[1:k]
    aucs <- numeric(5)
    for (f in 1:5) {
      tr <- fold != f; te <- fold == f
      fit <- suppressWarnings(glm(y ~ ., family = binomial,
                                  data = data.frame(y = y[tr], X[tr, g, drop = FALSE], check.names = FALSE)))
      pr <- suppressWarnings(predict(fit, newdata = data.frame(X[te, g, drop = FALSE], check.names = FALSE),
                                     type = "response"))
      aucs[f] <- as.numeric(pROC::auc(pROC::roc(y[te], pr, quiet = TRUE)))
    }
    res[[length(res) + 1]] <- data.table(combo = paste(cb, collapse = "+"), k = k,
                                         cv_auc = mean(aucs), cv_sd = sd(aucs),
                                         genes = paste(g, collapse = ","))
  }
}
res <- rbindlist(res)[order(-cv_auc)]
fwrite(res, file.path(RES, "06b_ml_combos_cv.csv"))
print(head(res[, .(combo, k, cv_auc, cv_sd)], 12))

best <- res[1]
genes <- strsplit(best$genes, ",")[[1]]
cat("最优组合:", best$combo, " k=", best$k, " 交叉验证AUC=", round(best$cv_auc, 4),
    " (SD", round(best$cv_sd, 4), ")\n", sep = "")
cat("特征基因:", paste(genes, collapse = ", "), "\n")

## ---------- 3. 最终模型 ----------
Xz <- scale(X[, genes, drop = FALSE])
mu <- attr(Xz, "scaled:center"); sdv <- attr(Xz, "scaled:scale")
fit_final <- suppressWarnings(glm(y ~ ., data = data.frame(y = y, Xz, check.names = FALSE),
                                  family = binomial))
cat("系数:\n"); print(round(coef(fit_final), 4))

pr <- predict(fit_final, type = "response")
r <- pROC::roc(y, pr, quiet = TRUE)
thr <- as.numeric(pROC::coords(r, "best", ret = "threshold", transpose = FALSE)[1, 1])
cat("表观AUC:", round(as.numeric(pROC::auc(r)), 4), " 阈值:", round(thr, 4), "\n")

## ---------- 4. Bootstrap 乐观度校正（Harrell） ----------
B <- 300
opt <- numeric(B)
for (b in 1:B) {
  idx <- sample(seq_along(upid), replace = TRUE)
  bp <- upid[idx]
  keep <- pid %in% bp
  ## 同一患者可能被抽中多次，按患者重采样（重复患者整体复制会导致重复样本，故仅取一次）
  fitb <- suppressWarnings(glm(y ~ ., family = binomial,
                               data = data.frame(y = y[keep], Xz[keep, , drop = FALSE], check.names = FALSE)))
  pb <- predict(fitb, type = "response")
  a_app <- as.numeric(pROC::auc(pROC::roc(y[keep], pb, quiet = TRUE)))
  pt <- suppressWarnings(predict(fitb, newdata = data.frame(Xz, check.names = FALSE), type = "response"))
  a_test <- as.numeric(pROC::auc(pROC::roc(y, pt, quiet = TRUE)))
  opt[b] <- a_app - a_test
}
corrected <- as.numeric(pROC::auc(r)) - mean(opt, na.rm = TRUE)
cat("Bootstrap 平均乐观度:", round(mean(opt, na.rm = TRUE), 4),
    " 校正后AUC:", round(corrected, 4), "\n")

saveRDS(list(genes = genes, model = fit_final, mu = mu, sd = sdv, threshold = thr,
             combo = best$combo, k = best$k, cv_auc = best$cv_auc, cv_sd = best$cv_sd,
             apparent_auc = as.numeric(pROC::auc(r)), optimism = mean(opt, na.rm = TRUE),
             auc_corrected = corrected, ranks = ranks, cons = cons, res = res),
        file.path(PROC, "signature_full.rds"))
cat("saved signature_full.rds\n")
