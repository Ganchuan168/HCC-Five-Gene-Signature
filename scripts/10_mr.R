#!/usr/bin/env Rscript
## 10 孟德尔随机化：肝脏 cis-eQTL → HCC
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table)})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")
args <- commandArgs(TRUE)
ARG  <- ifelse(length(args) >= 1, args[1], "11_deg_eqtl.csv")
ARG2 <- ifelse(length(args) >= 2, args[2], "12_finngen_hits.tsv")
OUTF <- ifelse(length(args) >= 3, args[3], "13_mr_results.csv")

e <- fread(file.path(RES, ARG))
f <- fread(file.path(RES, ARG2), header = FALSE)
setnames(f, c("chr", "pos2", "ref2", "alt2", "rsids", "nearest", "pval", "mlogp",
              "beta", "sebeta", "af", "af_cases", "af_ctrl"))
f[, chr := as.character(chr)]
e[, chr := sub("^chr", "", chr)]
setkey(e, chr, pos, ref, alt); setkey(f, chr, pos2, ref2, alt2)
m <- merge(e, f, by.x = c("chr", "pos", "ref", "alt"), by.y = c("chr", "pos2", "ref2", "alt2"))
cat("匹配到的 SNP-基因对:", nrow(m), "\n")

## F 统计量（近似）：F = (beta/se)^2
m[, F := (slope / slope_se)^2]
m <- m[F > 10]
cat("F>10 的工具变量:", nrow(m), "\n")

## Wald ratio
m[, b := beta / slope]
m[, se := sebeta / abs(slope)]
m[, z := b / se]
m[, p := 2 * pnorm(-abs(z))]
m[, OR := exp(b)]
m[, OR_low := exp(b - 1.96 * se)]
m[, OR_high := exp(b + 1.96 * se)]
m[, fdr := p.adjust(p, method = "BH")]

out <- m[order(p), .(symbol, variant_id, rsids, slope, slope_se, beta, sebeta, F,
                     OR, OR_low, OR_high, p, fdr)]
fwrite(out, file.path(RES, OUTF))
print(head(out[1:12, .(symbol, rsids, OR, OR_low, OR_high, p, fdr)], 12))
cat("\nFDR<0.05 的基因:", length(unique(out[fdr < 0.05]$symbol)), "\n")
cat("P<0.05 的基因:", length(unique(out[p < 0.05]$symbol)), "\n")

## 与机器学习特征基因的交集
sig <- readRDS(file.path(ROOT, "data/processed/signature.rds"))$genes
cat("特征基因在 MR 结果中:", paste(intersect(sig, out$symbol), collapse = ", "), "\n")
sub <- out[symbol %in% sig]
if (nrow(sub) > 0) print(sub)
