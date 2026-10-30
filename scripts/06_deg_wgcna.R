#!/usr/bin/env Rscript
## 06 发现集：差异表达 + WGCNA + 候选基因
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(limma); library(WGCNA)})
allowWGCNAThreads(nThreads = 4)
set.seed(20261029)
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd())
RES <- file.path(ROOT, "results/tables"); dir.create(RES, recursive = TRUE, showWarnings = FALSE)
PROC <- file.path(ROOT, "data/processed")

d <- readRDS(file.path(PROC, "GSE14520.rds"))
m <- d$expr; grp <- factor(d$group, levels = c("Normal", "Tumor"))
cat("输入:", dim(m), "\n")

## ---------- DEG ----------
design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(m, design), makeContrasts(Tumor - Normal, levels = design)))
tt <- topTable(fit, number = Inf, adjust.method = "BH", sort.by = "P")
tt$gene <- rownames(tt)
fwrite(tt, file.path(RES, "01_DEG_all.csv"))
sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, ]
cat("DEG:", nrow(sig), " 上调", sum(sig$logFC > 0), " 下调", sum(sig$logFC < 0), "\n")
fwrite(sig, file.path(RES, "02_DEG_sig.csv"))

## ---------- WGCNA ----------
dat <- t(m)
gsg <- goodSamplesGenes(dat, verbose = 0)
dat <- dat[gsg$goodSamples, gsg$goodGenes]
vars <- apply(dat, 2, var)
dat <- dat[, order(vars, decreasing = TRUE)[1:floor(0.4 * ncol(dat))]]
cat("WGCNA 输入:", dim(dat), "\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed hybrid", verbose = 0)
st <- sft$fitIndices
st$R2 <- -sign(st$slope) * st$SFT.R.sq
fwrite(st, file.path(RES, "03_soft_threshold.csv"))
pow <- st$Power[which(st$R2 >= 0.85)[1]]
if (is.na(pow)) pow <- st$Power[which.max(st$R2)]
cat("软阈值:", pow, " R2:", round(max(st$R2[st$Power == pow]), 3), "\n")

net <- blockwiseModules(dat, power = pow, networkType = "signed hybrid",
                        TOMType = "signed", minModuleSize = 50, mergeCutHeight = 0.25,
                        numericLabels = FALSE, randomSeed = 20261029, verbose = 0)
mods <- net$colors
cat("模块数:", length(unique(mods)), "\n"); print(table(mods))

MEs <- moduleEigengenes(dat, mods)$eigengenes
trait <- data.frame(Tumor = as.numeric(grp == "Tumor"))
rownames(trait) <- rownames(dat)
mt <- cor(MEs, trait$Tumor, use = "p")
mtp <- corPvalueStudent(mt, nrow(dat))
mtab <- data.frame(module = rownames(mt), cor = mt[, 1], p = mtp[, 1])
mtab <- mtab[order(-abs(mtab$cor)), ]
fwrite(mtab, file.path(RES, "04_module_trait.csv"))
print(head(mtab, 8))

best <- mtab$module[1]
best_genes <- names(mods)[mods == sub("^ME", "", best)]
cat("与表型最相关模块:", best, " 基因数:", length(best_genes), "\n")

inter <- intersect(best_genes, sig$gene)
cat("DEG ∩ 模块:", length(inter), "\n")
fwrite(data.frame(gene = inter), file.path(RES, "05_candidate_genes.csv"))
saveRDS(list(expr = m, group = grp, mods = mods, MEs = MEs, module = best,
             genes = inter, power = pow), file.path(PROC, "discovery_wgcna.rds"))
