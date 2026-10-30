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
