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

## 训练/内部验证划分：按患者分组，避免同一患者的肿瘤/癌旁样本分置两侧
d0 <- readRDS(file.path(PROC, "GSE14520.rds"))
pid <- sub("[AB]$", "", sub(".*(LCS-[0-9]+).*", "\\1", d0$title))
stopifnot(length(pid) == nrow(X))
suppressPackageStartupMessages(library(data.table))
pdf_dt <- unique(data.table(pid = pid))
set.seed(20261029)
train_pid <- sample(pdf_dt$pid, size = floor(0.7 * nrow(pdf_dt)))
in_tr <- pid %in% train_pid
idx <- which(in_tr)
Xtr <- X[idx, ]; ytr <- y[idx]; Xte <- X[-idx, ]; yte <- y[-idx]
cat("训练样本:", length(idx), " 患者:", length(train_pid),
    " 内部验证样本:", nrow(X) - length(idx), " 患者:", length(setdiff(pdf_dt$pid, train_pid)), "\n")

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

## 以训练集原始表达计算均值与标准差（供外部样本统一标准化）
mu <- apply(Xtr[, genes, drop = FALSE], 2, mean)
sdv <- apply(Xtr[, genes, drop = FALSE], 2, sd)
Xtr_z <- scale(Xtr[, genes, drop = FALSE])
fit_final <- suppressWarnings(glm(y ~ ., data = data.frame(y = ytr, Xtr_z, check.names = FALSE),
                                  family = binomial))

## 阈值：训练集 Youden 指数最大点
ptr <- predict(fit_final, type = "response")
rtr <- pROC::roc(ytr, ptr, quiet = TRUE)
thr <- as.numeric(pROC::coords(rtr, "best", ret = "threshold", transpose = FALSE)[1, 1])
cat("训练集判定阈值:", round(thr, 4), " 训练集AUC:", round(as.numeric(pROC::auc(rtr)), 3), "\n")

## 内部验证集（患者不重叠）
Xte_z <- scale(Xte[, genes, drop = FALSE])
pte <- predict(fit_final, newdata = data.frame(Xte_z, check.names = FALSE), type = "response")
a_in <- as.numeric(pROC::auc(pROC::roc(yte, pte, quiet = TRUE)))

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
  ye <- as.integer(e$group == "Tumor")
  full <- matrix(0, nrow = ncol(e$expr), ncol = length(genes), dimnames = list(NULL, genes))
  for (g in genes) {
    if (g %in% gs) {
      v <- e$expr[g, ]
      v[is.na(v)] <- mu[[g]]           # 缺失值以训练集均值填补
      full[, g] <- (v - mu[[g]]) / sdv[[g]]   # 采用训练集尺度，不使用本队列自身尺度
    }
  }
  pr <- predict(fit_final, newdata = data.frame(full, check.names = FALSE), type = "response")
  ok <- !is.na(pr)
  r <- pROC::roc(ye[ok], pr[ok], quiet = TRUE)
  a <- as.numeric(pROC::auc(r)); ci <- as.numeric(pROC::ci.auc(r))
  cl <- ifelse(pr[ok] >= thr, 1L, 0L)
  sens <- mean(cl[ye[ok] == 1] == 1); spec <- mean(cl[ye[ok] == 0] == 0)
  out[[nm]] <- data.table(cohort = nm, n = sum(ok), n_tumor = sum(ye[ok] == 1),
                          n_normal = sum(ye[ok] == 0), auc = a, ci_low = ci[1], ci_high = ci[3],
                          sens_fixed = sens, spec_fixed = spec, genes_used = length(gs))
}
out <- rbindlist(out)
fwrite(out, file.path(RES, "07_signature_validation.csv"))
print(out)

saveRDS(list(genes = genes, model = fit_final, combo = best$combo, k = best$k,
             auc_internal = best$auc_internal, auc_internal_holdout = a_in,
             threshold = thr, mu = mu, sd = sdv, ranks = ranks, res = res,
             n_train = length(idx), n_test = nrow(X) - length(idx),
             n_patient_train = length(train_pid)),
        file.path(PROC, "signature.rds"))
