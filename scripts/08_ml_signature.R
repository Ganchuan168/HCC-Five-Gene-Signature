#!/usr/bin/env Rscript
## 08 多算法特征筛选与外部验证
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(randomForest); library(e1071)
  library(xgboost); library(gbm); library(MASS); library(pROC); library(caret)
})
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

disc <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
cand <- disc$genes
X <- t(disc$expr[cand, ]); y <- as.integer(disc$group == "Tumor")
cat("候选基因:", length(cand), " 样本:", nrow(X), "\n")

## 训练/内部验证划分（分层 7:3）
set.seed(20261029)
idx <- createDataPartition(y, p = 0.7, list = FALSE)[, 1]
Xtr <- X[idx, ]; ytr <- y[idx]; Xte <- X[-idx, ]; yte <- y[-idx]

## ---------- 12 种算法的基因排序 ----------
rank_lasso <- function(X, y) {
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)
  cf <- as.matrix(coef(cv, s = "lambda.min"))[-1, ]
  names(sort(-abs(cf)))
}
rank_ridge <- function(X, y) {
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 5)
  cf <- as.matrix(coef(cv, s = "lambda.min"))[-1, ]
  names(sort(-abs(cf)))
}
rank_enet <- function(X, y) {
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0.5, nfolds = 5)
  cf <- as.matrix(coef(cv, s = "lambda.min"))[-1, ]
  names(sort(-abs(cf)))
}
rank_stepglm <- function(X, y) {
  ## 逐步 GLM：先按单因素 p 值预筛前 60 个基因，再做双向逐步回归
  pv <- apply(X, 2, function(x) tryCatch(summary(glm(y ~ x, family = binomial))$coef[2, 4], error = function(e) 1))
  pre <- names(sort(pv))[1:min(60, length(pv))]
  d <- data.frame(y = y, X[, pre, drop = FALSE], check.names = TRUE)
  fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
  st <- suppressWarnings(step(fit, direction = "both", trace = 0, k = 2))
  cf <- coef(st)[-1]
  if (length(cf) == 0) return(names(sort(pv)))
  names(sort(-abs(cf)))
}
rank_svm <- function(X, y) {
  fit <- svm(x = X, y = factor(y), kernel = "linear", scale = TRUE, probability = TRUE)
  w <- t(fit$coefs) %*% fit$SV
  v <- abs(w[1, ]); names(v) <- colnames(X)
  names(sort(-v))
}
rank_glmboost <- function(X, y) {
  suppressPackageStartupMessages(library(mboost))
  d <- data.frame(y = y, X, check.names = TRUE)
  fit <- suppressWarnings(glmboost(y ~ ., data = d, family = Binomial(), control = boost_control(mstop = 100)))
  v <- abs(coef(fit)[-1]); names(v) <- names(coef(fit))[-1]
  names(sort(-v))
}
rank_rf <- function(X, y) {
  set.seed(1)
  fit <- randomForest(x = X, y = factor(y), ntree = 500, importance = TRUE)
  v <- importance(fit)[, "MeanDecreaseGini"]
  names(sort(-v))
}
rank_lda <- function(X, y) {
  fit <- lda(x = X, grouping = factor(y))
  v <- abs(fit$scaling[, 1]); names(v) <- rownames(fit$scaling)
  names(sort(-v))
}
rank_xgb <- function(X, y) {
  dt <- xgb.DMatrix(X, label = y)
  fit <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1,
                                 nthread = 4, eval_metric = "auc"),
                   data = dt, nrounds = 100, verbose = 0)
  v <- xgb.importance(model = fit)
  v$Feature[order(-v$Gain)]
}
rank_gbm <- function(X, y) {
  d <- data.frame(y = y, X, check.names = TRUE)
  fit <- gbm(y ~ ., data = d, distribution = "bernoulli", n.trees = 300, interaction.depth = 2,
             shrinkage = 0.05, verbose = FALSE)
  v <- summary(fit, plotit = FALSE)
  v$var[order(-v$rel.inf)]
}
rank_nb <- function(X, y) {
  ## 朴素贝叶斯无原生排序，采用组间标准化均值差的绝对值作为基因排序
  s <- apply(X, 2, function(x) abs(mean(x[y == 1]) - mean(x[y == 0])) / (sd(x) + 1e-8))
  names(sort(-s))
}
rank_pls <- function(X, y) {
  ## 偏最小二乘：以第一主成分载荷绝对值排序
  pl <- prcomp(scale(X), center = FALSE, scale. = FALSE)
  v <- abs(pl$rotation[, 1]); names(v) <- colnames(X)
  names(sort(-v))
}
rank_plsglm <- function(X, y) {
  suppressPackageStartupMessages(library(plsRglm))
  d <- data.frame(y = y, X, check.names = TRUE)
  fit <- suppressWarnings(plsRglm(y ~ ., data = d, modele = "pls-glm-logistic", nk = 3, verbose = FALSE))
  cf <- fit$Coeffs
  if (is.null(cf)) return(rank_pls(X, y))
  v <- abs(cf); names(v) <- colnames(X)
  names(sort(-v))
}

algos <- list(
  Lasso = rank_lasso, Ridge = rank_ridge, Enet = rank_enet, Stepglm = rank_stepglm,
  SVM = rank_svm, glmBoost = rank_glmboost, RF = rank_rf, LDA = rank_lda,
  XGBoost = rank_xgb, GBM = rank_gbm, NB = rank_nb, plsRglm = rank_plsglm
)

cat("运行 12 种算法排序...\n")
ranks <- list()
for (nm in names(algos)) {
  r <- tryCatch(algos[[nm]](Xtr, ytr), error = function(e) NULL)
  if (is.null(r) || length(r) < 5) { cat("  跳过", nm, "\n"); next }
  ranks[[nm]] <- r
  cat("  ", nm, "->", length(r), "基因\n")
}

rk <- sapply(ranks, function(r) match(colnames(Xtr), r))
rownames(rk) <- colnames(Xtr)
rk[is.na(rk)] <- length(ranks[[1]]) + 1

## ---------- 组合评估 ----------
combos <- c(as.list(names(ranks)), combn(names(ranks), 2, simplify = FALSE))
ks <- c(5, 8, 10, 15, 20)
res <- list()
for (cb in combos) {
  sc <- rowMeans(rk[, cb, drop = FALSE])
  ord <- names(sort(sc))
  for (k in ks) {
    g <- ord[1:k]
    fit <- suppressWarnings(glm(ytr ~ ., data = data.frame(ytr = ytr, Xtr[, g, drop = FALSE]), family = binomial))
    pr <- predict(fit, newdata = data.frame(Xte[, g, drop = FALSE]), type = "response")
    a <- as.numeric(pROC::auc(pROC::roc(yte, pr, quiet = TRUE)))
    res[[length(res) + 1]] <- data.table(combo = paste(cb, collapse = "+"), k = k, auc_internal = a,
                                         genes = paste(g, collapse = ","))
  }
}
res <- rbindlist(res)[order(-auc_internal)]
fwrite(res, file.path(RES, "06_ml_combos.csv"))
print(head(res[, .(combo, k, auc_internal)], 10))

## ---------- 外部验证 ----------
best <- res[1]
genes <- strsplit(best$genes, ",")[[1]]
cat("最优组合:", best$combo, " k=", best$k, " 内部AUC=", round(best$auc_internal, 3), "\n")
cat("特征基因:", paste(genes, collapse = ", "), "\n")

Xtr_z <- scale(Xtr[, genes, drop = FALSE])
fit_final <- suppressWarnings(glm(y ~ ., data = data.frame(y = ytr, Xtr_z, check.names = FALSE), family = binomial))

ext <- list()
tc <- readRDS(file.path(PROC, "TCGA_LIHC.rds"))
ext[["TCGA-LIHC"]] <- list(expr = tc$expr, group = tc$group)
for (g in c("GSE76427", "GSE62232")) {
  d <- readRDS(file.path(PROC, paste0(g, ".rds")))
  ext[[g]] <- list(expr = d$expr, group = d$group)
}
out <- list()
for (nm in names(ext)) {
  e <- ext[[nm]]
  gs <- intersect(genes, rownames(e$expr))
  if (length(gs) < 3) { cat(nm, "可用基因不足\n"); next }
  Xe <- scale(t(e$expr[gs, , drop = FALSE])); ye <- as.integer(e$group == "Tumor")
  ## 缺失基因以队列均值补齐
  full <- matrix(0, nrow = nrow(Xe), ncol = length(genes), dimnames = list(NULL, genes))
  for (g in genes) full[, g] <- if (g %in% gs) Xe[, g] else 0
  pr <- predict(fit_final, newdata = data.frame(full, check.names = FALSE), type = "response")
  ok <- !is.na(pr)
  a <- as.numeric(pROC::auc(pROC::roc(ye[ok], pr[ok], quiet = TRUE)))
  ci <- as.numeric(pROC::ci.auc(pROC::roc(ye[ok], pr[ok], quiet = TRUE)))
  ## 最佳阈值下的敏感度/特异度
  co <- pROC::coords(pROC::roc(ye[ok], pr[ok], quiet = TRUE), "best", ret = c("sensitivity", "specificity"))
  out[[nm]] <- data.table(cohort = nm, n = sum(ok), n_tumor = sum(ye[ok] == 1), n_normal = sum(ye[ok] == 0),
                          auc = a, ci_low = ci[1], ci_high = ci[3],
                          sens = co$sensitivity, spec = co$specificity,
                          genes_used = length(gs))
}
out <- rbindlist(out)
fwrite(out, file.path(RES, "07_signature_validation.csv"))
print(out)
saveRDS(list(genes = genes, model = fit_final, combo = best$combo, k = best$k,
             auc_internal = best$auc_internal, ranks = ranks, res = res), file.path(PROC, "signature.rds"))
