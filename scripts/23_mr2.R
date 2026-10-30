#!/usr/bin/env Rscript
## 23 MR（限定于418个候选基因）：等位基因协调、回文SNP剔除、多SNP基因IVW合并
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); RES <- file.path(ROOT, "results/tables")

cand <- fread(file.path(RES, "05_candidate_genes.csv"))$gene

run <- function(eqtl_file, hits_file, out_file, tag) {
  e <- fread(file.path(RES, eqtl_file))
  e <- e[symbol %in% cand]
  f <- fread(file.path(RES, hits_file), header = FALSE)
  setnames(f, c("chr", "pos2", "ref2", "alt2", "rsids", "nearest", "pval", "mlogp",
                "beta", "sebeta", "af", "af_cases", "af_ctrl"))
  f[, chr := sub("^chr", "", as.character(chr))]
  e[, chr := sub("^chr", "", as.character(chr))]
  m <- merge(e, f, by.x = c("chr", "pos", "ref", "alt"), by.y = c("chr", "pos2", "ref2", "alt2"))
  cat(tag, "匹配 SNP-基因对:", nrow(m), "\n")

  ## 等位基因协调：按参照/替代等位基因精确匹配后，效应方向天然一致；
  ## 剔除回文 SNP（A/T、C/G）中次等位基因频率介于 0.4–0.6 者
  m[, palindromic := (toupper(ref) == "A" & toupper(alt) == "T") |
                      (toupper(ref) == "T" & toupper(alt) == "A") |
                      (toupper(ref) == "C" & toupper(alt) == "G") |
                      (toupper(ref) == "G" & toupper(alt) == "C")]
  m[, af := as.numeric(af)]
  amb <- m[palindromic == TRUE & (af > 0.4 & af < 0.6)]
  cat(tag, "剔除回文且MAF接近0.5的SNP:", nrow(amb), "\n")
  m <- m[!(palindromic == TRUE & (af > 0.4 & af < 0.6))]

  m[, F := (slope / slope_se)^2]
  m <- m[F > 10]
  m[, b := beta / slope]
  m[, se := sebeta / abs(slope)]
  m[, p := 2 * pnorm(-abs(b / se))]
  m[, OR := exp(b)]

  ## 逐基因：单 SNP 用 Wald 比值；≥2 个 SNP 用固定效应 IVW
  per_gene <- m[, {
    if (.N == 1) {
      list(n_snp = 1L, b = b, se = se, method = "Wald")
    } else {
      w <- 1 / se^2
      list(n_snp = .N, b = sum(b * w) / sum(w), se = sqrt(1 / sum(w)), method = "IVW")
    }
  }, by = symbol]
  per_gene[, p := 2 * pnorm(-abs(b / se))]
  per_gene[, OR := exp(b)]
  per_gene[, OR_low := exp(b - 1.96 * se)]
  per_gene[, OR_high := exp(b + 1.96 * se)]
  per_gene[, fdr := p.adjust(p, "BH")]
  setorder(per_gene, p)
  fwrite(per_gene, file.path(RES, out_file))

  cat(tag, "基因数:", nrow(per_gene), " 名义显著(P<0.05):", sum(per_gene$p < 0.05),
      " FDR<0.05:", sum(per_gene$fdr < 0.05), " 最小FDR:", signif(min(per_gene$fdr), 3), "\n")
  print(head(per_gene[, .(symbol, n_snp, method, OR, OR_low, OR_high, p, fdr)], 8))
  per_gene[, tissue := tag]
  per_gene
}

lv <- run("11_deg_eqtl.csv", "12_finngen_hits.tsv", "13b_mr_liver.csv", "Liver")
bd <- run("14_blood_eqtl.csv",      "15_blood_finngen.tsv", "16b_mr_blood.csv", "Blood")
fwrite(rbind(lv, bd), file.path(RES, "17b_mr_combined.csv"))
