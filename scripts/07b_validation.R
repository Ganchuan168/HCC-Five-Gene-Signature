#!/usr/bin/env Rscript
## 07b 五基因特征的多队列判别效能
## 主分析：各队列内对特征基因行 z 分数标准化（跨平台量纲可比的常规做法），
##        阈值固定为发现集 Youden 最优点；
## 敏感性分析：以发现集均值/标准差直接外推（不做队列内标准化）。
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); PROC <- file.path(ROOT, "data/processed")

sf <- readRDS(file.path(PROC, "signature_final.rds"))
genes <- sf$genes; mu <- sf$mu; sdv <- sf$sd; thr <- sf$threshold
## 决策规则定义在“队列内标准化线性预测值”尺度上：切点取发现集 Youden 指数最大点在该尺度上的取值
disc0 <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
zd0 <- t(scale(t(disc0$expr[genes, ]))); yd0 <- as.integer(disc0$group == "Tumor")
lpd0 <- as.numeric(t(zd0) %*% coef(sf$model)[-1])
cut_lp <- as.numeric(pROC::coords(pROC::roc(yd0, lpd0, quiet = TRUE), "best",
                                  ret = "threshold")[1])
cat("特征基因:", paste(genes, collapse = ", "), " 阈值(概率):", round(thr, 4),
    " 标准化评分切点:", round(cut_lp, 3), "\n")

cohorts <- list(
  "TCGA-LIHC" = "TCGA_LIHC",
  "GSE76427"  = "GSE76427",
  "GSE62232"  = "GSE62232"
)

out <- list()
for (nm in names(cohorts)) {
  d <- readRDS(file.path(PROC, paste0(cohorts[[nm]], ".rds")))
  e <- d$expr[genes, , drop = FALSE]
  ye <- as.integer(d$group == "Tumor")
  ## 主分析：队列内标准化
  zi <- t(scale(t(e)))
  zi[is.na(zi)] <- 0
  pr <- as.numeric(t(zi) %*% coef(sf$model)[-1])
  p <- plogis(pr + coef(sf$model)[1] * 0)   # 线性预测值，ROC 不受影响
  r <- roc(ye, pr, quiet = TRUE); ci <- as.numeric(ci.auc(r))
  cl <- as.integer(pr >= cut_lp)
  sens <- mean(cl[ye == 1] == 1); spec <- mean(cl[ye == 0] == 0)
  ## 敏感性分析：发现集尺度外推
  zt <- sweep(sweep(e, 1, mu[genes], "-"), 1, sdv[genes], "/")
  zt[is.na(zt)] <- 0
  pt <- as.numeric(t(zt) %*% coef(sf$model)[-1])
  rt <- roc(ye, pt, quiet = TRUE); cit <- as.numeric(ci.auc(rt))
  out[[nm]] <- data.table(cohort = nm, n = ncol(e), n_tumor = sum(ye == 1), n_normal = sum(ye == 0),
                          auc = as.numeric(auc(r)), ci_low = ci[1], ci_high = ci[3],
                          sens = sens, spec = spec,
                          auc_transfer = as.numeric(auc(rt)), t_low = cit[1], t_high = cit[3])
  cat(sprintf("%-10s n=%3d 内标准化 AUC=%.3f (%.3f-%.3f) 敏感度=%.1f%% 特异度=%.1f%% | 训练尺度外推 AUC=%.3f\n",
              nm, ncol(e), as.numeric(auc(r)), ci[1], ci[3], 100 * sens, 100 * spec, as.numeric(auc(rt))))
}
out <- rbindlist(out)
fwrite(out, file.path(RES, "07_signature_validation.csv"))

## 发现集表现与患者分组交叉验证
disc <- readRDS(file.path(PROC, "discovery_wgcna.rds"))
X <- t(disc$expr[genes, ]); y <- as.integer(disc$group == "Tumor")
d0 <- readRDS(file.path(PROC, "GSE14520.rds"))
pid <- sub("[AB]$", "", sub(".*(LCS-[0-9]+).*", "\\1", d0$title))
set.seed(20261029)
upid <- unique(pid); fo <- integer(length(upid)); names(fo) <- upid
fo[sample(upid)] <- rep(1:5, length.out = length(upid))
fold <- fo[pid]
au <- numeric(5)
for (f in 1:5) {
  tr <- fold != f; te <- fold == f
  mu2 <- apply(X[tr, ], 2, mean); sd2 <- apply(X[tr, ], 2, sd)
  Xtr <- sweep(sweep(X[tr, ], 2, mu2, "-"), 2, sd2, "/")
  Xte <- sweep(sweep(X[te, ], 2, mu2, "-"), 2, sd2, "/")
  fit <- glm(y ~ ., data = data.frame(y = y[tr], Xtr, check.names = FALSE), family = binomial)
  pr <- predict(fit, newdata = data.frame(Xte, check.names = FALSE), type = "link")
  au[f] <- as.numeric(auc(roc(y[te], pr, quiet = TRUE)))
}
cat("患者分组 5 折交叉验证 AUC:", round(mean(au), 4), "±", round(sd(au), 4), "\n")
fwrite(data.table(fold = 1:5, auc = au), file.path(RES, "07b_cv_folds.csv"))
saveRDS(list(cv_mean = mean(au), cv_sd = sd(au), folds = au), file.path(PROC, "cv5.rds"))
