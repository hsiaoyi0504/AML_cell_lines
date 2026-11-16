# Compare mTOR-High vs mTOR-Low AML cell lines (groups defined by RNA, DepMap 26Q1)
# across every dataset in data/cell_lines/. Base R + xml2 only.
#
# Usage: Rscript scripts/compare_mTOR_RNA_groups.R   (run from project root)

suppressMessages(library(xml2))

root    <- normalizePath(".")
dd      <- file.path(root, "data/cell_lines")
out_dir <- file.path(root, "results/mTOR_RNA_groups")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- groups ----
cl <- read.delim(file.path(dd, "fig2_MYC_mTOR_DepMap_all_cell_line_classification.tsv"))
cl <- cl[grepl("^RNA", cl$modality), ]
model <- read.csv(file.path(dd, "DepMap_26Q3/Model.csv"))
idx <- match(cl$cell_line, model$StrippedCellLineName)
stopifnot(!anyNA(idx))
cl$ModelID  <- model$ModelID[idx]
cl$CCLEName <- paste0(cl$cell_line, "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
cl <- cbind(cl, model[idx, c("OncotreeSubtype", "Sex", "Age", "PrimaryOrMetastasis")])
write.table(cl, file.path(out_dir, "groups_RNA_mTOR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("mTOR groups:", table(cl$mTOR_group), "\n")

# ------------------------------------------------------------- statistics ---
# X: samples x features (rownames = cell_line). Returns one row per feature.
diff_test <- function(X, label, min_n = 3, annot = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  g  <- cl$mTOR_group[match(rownames(X), cl$cell_line)]
  H  <- X[g == "High", , drop = FALSE]
  L  <- X[g == "Low",  , drop = FALSE]
  nH <- colSums(!is.na(H)); nL <- colSums(!is.na(L))
  mH <- colMeans(H, na.rm = TRUE); mL <- colMeans(L, na.rm = TRUE)
  vH <- colSums(sweep(H, 2, mH)^2, na.rm = TRUE) / (nH - 1)
  vL <- colSums(sweep(L, 2, mL)^2, na.rm = TRUE) / (nL - 1)
  keep <- nH >= min_n & nL >= min_n & (vH + vL) > 0
  se  <- sqrt(vH / nH + vL / nL)
  t   <- (mH - mL) / se
  df  <- (vH / nH + vL / nL)^2 / ((vH / nH)^2 / (nH - 1) + (vL / nL)^2 / (nL - 1))
  p   <- 2 * pt(-abs(t), df)
  sdp <- sqrt(((nH - 1) * vH + (nL - 1) * vL) / (nH + nL - 2))
  res <- data.frame(dataset = label, feature = colnames(X), n_High = nH, n_Low = nL,
                    mean_High = mH, mean_Low = mL, diff_High_minus_Low = mH - mL,
                    cohens_d = (mH - mL) / sdp, t = t, p_welch = p,
                    stringsAsFactors = FALSE)[keep, ]
  res$fdr_welch <- p.adjust(res$p_welch, "BH")
  res$p_wilcox <- vapply(res$feature, function(f)
    suppressWarnings(wilcox.test(H[, f], L[, f], exact = FALSE)$p.value), 0)
  res$fdr_wilcox <- p.adjust(res$p_wilcox, "BH")
  res <- res[order(res$p_welch), ]
  if (!is.null(annot)) {                  # extra per-feature columns, e.g. compound names
    a   <- annot[match(res$feature, annot[[1]]), -1, drop = FALSE]
    res <- cbind(res[, 1:2], a, res[, -(1:2)])
  }
  attr(res, "X") <- X                     # per-line values, used for plotting
  write.table(res, file.path(out_dir, paste0("diff_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-22s lines H/L = %2d/%2d  features = %6d  FDR<0.05 = %5d  FDR<0.25 = %5d\n",
              label, sum(g == "High"), sum(g == "Low"), nrow(res),
              sum(res$fdr_welch < 0.05), sum(res$fdr_welch < 0.25)))
  invisible(res)
}

# rows = DepMap ModelIDs -> rows = cell_line
by_model <- function(df, id_col = 1) {
  ids <- df[[id_col]]
  df  <- df[ids %in% cl$ModelID, , drop = FALSE]
  rn  <- cl$cell_line[match(df[[id_col]], cl$ModelID)]
  X   <- as.matrix(df[, -id_col, drop = FALSE]); rownames(X) <- rn
  X
}

# read only the header + rows whose first field is one of our ModelIDs
read_rows <- function(path, ids) {
  hdr <- readLines(path, n = 1)
  tmp <- tempfile()
  pat <- paste0("(^|,)(", paste(ids, collapse = "|"), "),")
  system2("grep", c("-E", shQuote(pat), shQuote(path)), stdout = tmp)
  read.csv(text = c(hdr, readLines(tmp)), check.names = FALSE)
}

# read a features-in-rows table keeping only the needed sample columns
read_cols <- function(path, keep_cols, sep, skip = 0, id_cols = 1) {
  hdr <- strsplit(readLines(path, n = skip + 1)[skip + 1], sep, fixed = TRUE)[[1]]
  hdr <- gsub('^"|"$', "", hdr)
  cc  <- ifelse(hdr %in% keep_cols, "numeric", "NULL")
  cc[id_cols] <- "character"
  read.table(path, sep = sep, header = TRUE, skip = skip, colClasses = cc, quote = "\"",
             check.names = FALSE, na.strings = c("NA", "NaN", ""), strip.white = TRUE,
             comment.char = "")
}

ccle2line <- function(x) sub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", "", x)

results <- list()

# ------------------------------------------------------- 1. Annotations -----
sink(file.path(out_dir, "annotation_tests.txt"))
cat("== OncotreeSubtype x mTOR group\n"); tb <- table(cl$OncotreeSubtype, cl$mTOR_group); print(tb)
print(fisher.test(tb, simulate.p.value = TRUE, B = 1e5))
cat("== Sex\n"); tb <- table(cl$Sex, cl$mTOR_group); print(tb); print(fisher.test(tb))
cat("== KMT2A rearranged\n"); tb <- table(cl$KMT2A_rearranged, cl$mTOR_group); print(tb); print(fisher.test(tb))
cat("== MYC group (RNA)\n"); tb <- table(cl$MYC_group, cl$mTOR_group); print(tb); print(fisher.test(tb))
cat("== MYC score vs mTOR score (Spearman)\n"); print(cor.test(cl$MYC_score, cl$mTOR_score, method = "spearman"))
cat("== Age\n"); print(wilcox.test(Age ~ mTOR_group, data = cl))
sink()

# ------------------------------------------------- 2. Expression (26Q3) -----
ex <- read_rows(file.path(dd, "DepMap_26Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"), cl$ModelID)
ex <- ex[ex$IsDefaultEntryForModel == "Yes", ]
X  <- by_model(ex[, c(4, 7:ncol(ex))])
X  <- X[, colSums(X > 1, na.rm = TRUE) >= 5]          # expressed in >= 5 lines
results$expression <- diff_test(X, "expression_TPMlog")

# ---------------------------------------------- 3. CRISPR Chronos (26Q3) ----
cr <- read_rows(file.path(dd, "DepMap_26Q3/CRISPRGeneEffect.csv"), cl$ModelID)
results$crispr <- diff_test(by_model(cr), "CRISPR_Chronos")

# ------------------------------------------------------- 4. DEMETER2 RNAi ---
d2 <- read_cols(file.path(dd, "DEMETER2_v6/D2_combined_gene_dep_scores.csv"), cl$CCLEName, ",")
X  <- t(as.matrix(d2[, -1])); colnames(X) <- d2[[1]]; rownames(X) <- ccle2line(rownames(X))
results$rnai <- diff_test(X, "RNAi_DEMETER2", min_n = 3)

# ----------------------------------------------------------- 5. Drugs -------
drug_files <- c(PRISM = "PRISM/REPURPOSINGAUCMatrix.csv", CTRP = "CTD2/CTRPAUCMatrix.csv",
                GDSC2 = "GDSC/GDSC2AUCMatrix.csv")
cmpd <- read.csv(file.path(dd, "DepMap_26Q3/PortalCompounds.csv"))[
  , c("CompoundID", "CompoundName", "GeneSymbolOfTargets", "TargetOrMechanism")]
for (s in names(drug_files)) {
  dr <- read_rows(file.path(dd, "DepMap_Portal_Harmonized_Drug_Screens", drug_files[[s]]), cl$ModelID)
  results[[paste0("drug_", s)]] <- diff_test(by_model(dr), paste0("drug_AUC_", s), annot = cmpd)
}

# --------------------------------------------------- 6. RRBS methylation ----
me <- read_cols(file.path(dd, "methylation_RRBS/CCLE_RRBS_TSS_1kb_20180614.txt"),
                cl$CCLEName, "\t", id_cols = 1:7)
X  <- t(as.matrix(me[, -(1:7)])); colnames(X) <- paste(me$gene, me$TSS_id, sep = "|")
rownames(X) <- ccle2line(rownames(X))
X  <- X[, apply(X, 2, sd, na.rm = TRUE) > 0.05 & colSums(!is.na(X)) >= 10]
results$meth <- diff_test(X, "methylation_TSS1kb")

# same RRBS data summarised per CpG cluster near each TSS (finer resolution)
mc <- read_cols(file.path(dd, "methylation_RRBS/CCLE_RRBS_TSS_CpG_clusters_20180614.txt"),
                cl$CCLEName, "\t", id_cols = 1:5)
X  <- t(as.matrix(mc[, -(1:5)])); colnames(X) <- paste(mc$gene_name, mc$cluster_id, sep = "|")
rownames(X) <- ccle2line(rownames(X))
X  <- X[, apply(X, 2, sd, na.rm = TRUE) > 0.05 & colSums(!is.na(X)) >= 10]
results$meth_cpg <- diff_test(X, "methylation_CpGcluster")

# ----------------------------------------------------- 7. miRNA (CCLE) ------
mi <- read_cols(file.path(dd, "CCLE_2019/CCLE_miRNA_20181103.gct"), cl$CCLEName, "\t", skip = 2, id_cols = 1:2)
X  <- log2(t(as.matrix(mi[, -(1:2)])) + 1); colnames(X) <- mi$Description
rownames(X) <- ccle2line(rownames(X))
results$mirna <- diff_test(X, "miRNA_log2")

# -------------------------------------------- 8. Global chromatin (CCLE) ----
ch <- read.csv(file.path(dd, "CCLE_2019/CCLE_GlobalChromatinProfiling_20181130.csv"), check.names = FALSE)
ch <- ch[ch$CellLineName %in% cl$CCLEName, ]
X  <- as.matrix(ch[, -(1:2)]); rownames(X) <- ccle2line(ch$CellLineName)
results$chrom <- diff_test(X, "histone_marks")

# ------------------------------------------- 9. Metabolomics (CCLE 2019) ----
read_xlsx_sheet <- function(dir, sheet_file) {
  ss <- xml_text(xml_find_all(read_xml(file.path(dir, "xl/sharedStrings.xml")), "//d1:si", ns = c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")))
  x  <- read_xml(file.path(dir, "xl/worksheets", sheet_file))
  ns <- c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  cells <- xml_find_all(x, "//d1:sheetData/d1:row/d1:c", ns)
  ref <- xml_attr(cells, "r"); typ <- xml_attr(cells, "t")
  val <- xml_text(xml_find_first(cells, "d1:v", ns))
  val[!is.na(typ) & typ == "s"] <- ss[as.integer(val[!is.na(typ) & typ == "s"]) + 1]
  col <- gsub("[0-9]", "", ref); row <- as.integer(gsub("[A-Z]", "", ref))
  colnum <- vapply(strsplit(col, ""), function(l) sum((match(l, LETTERS)) * 26^(rev(seq_along(l)) - 1)), 0)
  M <- matrix(NA_character_, max(row), max(colnum)); M[cbind(row, colnum)] <- val
  M
}
xdir <- tempfile(); unzip(file.path(dd, "targeted_metabolomics/Supp_Tables.xlsx"),
                          files = c("xl/sharedStrings.xml", "xl/worksheets/sheet3.xml"), exdir = xdir)
M  <- read_xlsx_sheet(xdir, "sheet3.xml")            # sheet "1-clean data"
hdr <- M[1, ]; M <- M[-1, ]
mt_ids <- M[, 1]
X  <- suppressWarnings(matrix(as.numeric(M[, -1]), nrow(M), dimnames = list(ccle2line(mt_ids), hdr[-1])))
X  <- X[rownames(X) %in% cl$cell_line, !is.na(hdr[-1]), drop = FALSE]
results$metab <- diff_test(X, "metabolomics")

# ------------------------------------------ 10. mTOR pathway focus panel ---
panel <- c("MTOR", "RPTOR", "RICTOR", "MLST8", "RHEB", "TSC1", "TSC2", "AKT1", "AKT2", "PIK3CA",
           "PIK3CD", "PTEN", "RPS6KB1", "EIF4EBP1", "EIF4E", "LAMTOR1", "RRAGA", "RRAGC",
           "SLC7A5", "SLC3A2", "MYC", "MAX", "DEPTOR", "PRKAA1", "STK11", "NPRL2", "DEPDC5", "SZT2")
sym <- function(f) sub(" \\(.*$", "", f)
foc <- do.call(rbind, lapply(results[c("expression", "crispr", "rnai")], function(r) r[sym(r$feature) %in% panel, ]))
foc <- foc[order(foc$dataset, foc$p_welch), ]
write.table(foc, file.path(out_dir, "mTOR_pathway_panel.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(results, file.path(out_dir, "all_results.rds"))
cat("done ->", out_dir, "\n")
