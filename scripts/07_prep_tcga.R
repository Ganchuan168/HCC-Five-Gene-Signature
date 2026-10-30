#!/usr/bin/env Rscript
## 07 TCGA-LIHC 预处理：表达矩阵、分组、临床与生存
.libPaths(c(file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- Sys.getenv("HCC_ROOT", unset = getwd()); PROC <- file.path(ROOT, "data/processed")

expr <- fread(file.path(ROOT, "data/raw/tcga/TCGA-LIHC.star_counts.tsv.gz"))
ens <- sub("\\..*$", "", expr[[1]])
m <- as.matrix(expr[, -1]); rownames(m) <- ens; storage.mode(m) <- "double"
if (max(m, na.rm = TRUE) > 30) m <- log2(m + 1)
cat("矩阵:", dim(m), " 值域:", round(range(m, na.rm = TRUE), 2), "\n")

barc <- colnames(m)
type <- substr(barc, 14, 15)
grp <- ifelse(type == "01", "Tumor", ifelse(type == "11", "Normal", NA))
keep <- !is.na(grp)
m <- m[, keep]; grp <- grp[keep]; barc <- barc[keep]
pid <- substr(barc, 1, 12)
dup <- duplicated(paste(pid, grp))
cat("重复样本:", sum(dup), "\n")
m <- m[, !dup]; grp <- grp[!dup]; barc <- barc[!dup]; pid <- pid[!dup]
cat("去重后:", dim(m), " 肿瘤", sum(grp == "Tumor"), " 癌旁", sum(grp == "Normal"), "\n")

## Ensembl → SYMBOL
suppressPackageStartupMessages(library(org.Hs.eg.db))
sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(m),
        keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
ok <- !is.na(sym) & sym != ""
m <- m[ok, ]; sym <- sym[ok]
ord <- order(rowMeans(m), decreasing = TRUE)
m <- m[ord, ]; sym <- sym[ord]
keep <- !duplicated(sym)
m <- m[keep, ]; rownames(m) <- sym[keep]
cat("注释后:", dim(m), "\n")

## 临床
cases <- fromJSON(file.path(ROOT, "data/raw/tcga/cases.json"), simplifyVector = FALSE)
clin <- rbindlist(lapply(cases, function(cs) {
  dg <- if (length(cs$diagnoses) > 0) cs$diagnoses[[1]] else list()
  data.table(patient = cs$submitter_id,
             age = ifelse(is.null(cs$demographic$age_at_index), NA, cs$demographic$age_at_index),
             gender = ifelse(is.null(cs$demographic$sex_at_birth), NA, cs$demographic$sex_at_birth),
             stage = ifelse(is.null(dg$ajcc_pathologic_stage), NA, dg$ajcc_pathologic_stage),
             grade = ifelse(is.null(dg$tumor_grade), NA, dg$tumor_grade),
             t_stage = ifelse(is.null(dg$ajcc_pathologic_t), NA, dg$ajcc_pathologic_t),
             vital = ifelse(is.null(cs$demographic$vital_status), NA, cs$demographic$vital_status))
}))
surv <- fread(file.path(ROOT, "data/raw/tcga/TCGA-LIHC.survival.tsv.gz"))
surv[, patient := substr(sample, 1, 12)]
surv <- unique(surv, by = "patient")
clin <- merge(clin[!duplicated(patient)], surv[, .(patient, OS.time, OS)], by = "patient", all.x = TRUE)
cat("临床:", nrow(clin), "\n"); print(table(clin$stage, useNA = "ifany"))

saveRDS(list(expr = m, group = grp, barcode = barc, patient = pid, clin = clin),
        file.path(PROC, "TCGA_LIHC.rds"))
