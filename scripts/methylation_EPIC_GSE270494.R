# Differential methylation, mTOR-High vs mTOR-Low (RNA groups), using the
# Illumina EPIC array data of GSE270494 (Noguera-Castells et al., average beta).
#
#   * values with detection p > 0.05 set to NA; probes masked by Zhou et al.
#     (M_general), on chrX/chrY/chrM or without a mapping are removed
#   * beta -> M value (log2(beta / (1 - beta)), beta clipped to [0.001, 0.999]);
#     Welch t-test on M values, delta beta reported as the effect size
#   * promoter probes (-1500..+500 bp of a GENCODE v36 TSS) collapsed per gene
#     (Sidak-corrected best probe), then the same enrichment as for RRBS
#   * comparison with the RRBS CpG-cluster results; methylation vs expression
#
# Inputs:  data/cell_lines/GSE270494/*Average_Beta*.csv.gz
#          data/annotation/EPIC.hg38.manifest.gencode.v36.tsv.gz, EPIC.hg38.mask.tsv.gz
#          (https://github.com/zhou-lab/InfiniumAnnotationV1)
# Usage:   Rscript scripts/methylation_EPIC_GSE270494.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"
fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R"); source("scripts/enrichment_utils.R")
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
results <- readRDS(file.path(res_dir, "all_results.rds"))

beta_file <- "data/cell_lines/GSE270494/GSE270494_Noguera-Castells_Average_Beta_Homo_Sapiens.csv.gz"
det_p_max <- 0.05   # 0.01 removes ~18% of values, mostly high-beta ones (weak det. p in this file)

# ------------------------------------------------------ match cell lines ---
hdr  <- scan(pipe(paste("gzip -dc", shQuote(beta_file), "| head -1")), what = "", sep = ",", quiet = TRUE)
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
samp <- hdr[-1][!grepl("Detection Pval$", hdr[-1])]
samp <- samp[norm(samp) %in% cl$cell_line]           # sublines (KG1A, HL-60/MX2) do not match
col_beta <- match(samp, hdr); col_det <- match(paste(samp, "Detection Pval"), hdr)
stopifnot(!anyNA(col_det))
lines_ <- norm(samp)
grp    <- cl$mTOR_group[match(lines_, cl$cell_line)]
cat("GSE270494 lines matched:", length(samp), "| High", sum(grp == "High"), "| Low", sum(grp == "Low"), "\n")

# --------------------------------------------- read only the needed columns
cols <- sort(c(1, col_beta, col_det))
raw  <- read.csv(pipe(paste("gzip -dc", shQuote(beta_file), "| cut -d, -f", paste(cols, collapse = ","))),
                 check.names = FALSE, colClasses = c("character", rep("numeric", length(cols) - 1)))
B <- as.matrix(raw[, samp]); P <- as.matrix(raw[, paste(samp, "Detection Pval")])
rownames(B) <- raw$ID_REF; colnames(B) <- lines_
rm(raw); invisible(gc())
n_raw <- nrow(B); n_failed <- sum(P > det_p_max, na.rm = TRUE)
B[is.na(P) | P > det_p_max] <- NA
rm(P); invisible(gc())

# ----------------------------------------------------------- probe QC -------
man  <- read.delim(gzfile("data/annotation/EPIC.hg38.manifest.gencode.v36.tsv.gz"), colClasses = "character")
mask <- read.delim(gzfile("data/annotation/EPIC.hg38.mask.tsv.gz"), colClasses = "character")
masked <- mask$Probe_ID[mask$M_general == "TRUE"]
man  <- man[match(rownames(B), man$probeID), ]
keep <- !is.na(man$CpG_chrm) & !man$CpG_chrm %in% c("chrX", "chrY", "chrM") & !rownames(B) %in% masked
qc <- c(raw = n_raw, unmapped = sum(is.na(man$CpG_chrm)),
        sex_chrM = sum(man$CpG_chrm %in% c("chrX", "chrY", "chrM")),
        masked = sum(rownames(B) %in% masked & !is.na(man$CpG_chrm)))
B <- B[keep, ]; man <- man[keep, ]
nH <- rowSums(!is.na(B[, grp == "High"])); nL <- rowSums(!is.na(B[, grp == "Low"]))
ok <- nH >= 3 & nL >= 3 & apply(B, 1, sd, na.rm = TRUE) > 0.02     # drop invariant probes
qc <- c(qc, after_QC = sum(keep), too_few_or_invariant = sum(!ok), tested = sum(ok))
B <- B[ok, ]; man <- man[ok, ]

M <- log2(pmin(pmax(B, 0.001), 0.999) / (1 - pmin(pmax(B, 0.001), 0.999)))

# ------------------------------------------------- Welch t on M values -----
welch <- function(X, g) {
  H <- X[, g == "High", drop = FALSE]; L <- X[, g == "Low", drop = FALSE]
  nH <- rowSums(!is.na(H)); nL <- rowSums(!is.na(L))
  mH <- rowMeans(H, na.rm = TRUE); mL <- rowMeans(L, na.rm = TRUE)
  vH <- rowSums((H - mH)^2, na.rm = TRUE) / (nH - 1); vL <- rowSums((L - mL)^2, na.rm = TRUE) / (nL - 1)
  se <- sqrt(vH / nH + vL / nL); t <- (mH - mL) / se
  df <- (vH / nH + vL / nL)^2 / ((vH / nH)^2 / (nH - 1) + (vL / nL)^2 / (nL - 1))
  list(nH = nH, nL = nL, mH = mH, mL = mL, t = t, p = 2 * pt(-abs(t), df))
}
w <- welch(M, grp)
gene1 <- sub(";.*$", "", man$genesUniq)
res <- data.frame(probe = rownames(M), gene = gene1, chr = man$CpG_chrm, pos = man$CpG_beg,
                  CGI = man$CGIposition, n_High = w$nH, n_Low = w$nL,
                  meanM_High = w$mH, meanM_Low = w$mL, diff_M = w$mH - w$mL,
                  meanBeta_High = rowMeans(B[, grp == "High"], na.rm = TRUE),
                  meanBeta_Low  = rowMeans(B[, grp == "Low"],  na.rm = TRUE),
                  t = w$t, p_welch = w$p, stringsAsFactors = FALSE)
res$delta_beta <- res$meanBeta_High - res$meanBeta_Low
res$fdr_welch  <- p.adjust(res$p_welch, "BH")
res <- res[order(res$p_welch), ]
gz <- gzfile(file.path(res_dir, "diff_methylation_EPIC_GSE270494.tsv.gz"), "w")
write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

# ------------------------------------------- promoter probes -> genes -------
# a probe is a promoter probe of every gene with a transcript TSS within -1500..+500 bp
gn <- strsplit(man$geneNames, ";"); dt <- strsplit(man$distToTSS, ";")
pg <- do.call(rbind, lapply(seq_along(gn), function(i) {
  d <- suppressWarnings(as.numeric(dt[[i]]))
  if (!length(d) || all(is.na(d))) return(NULL)
  g <- unique(gn[[i]][!is.na(d) & d >= -1500 & d <= 500])
  if (length(g)) data.frame(probe = rownames(M)[i], gene = g) else NULL
}))
pg <- merge(pg, res[, c("probe", "p_welch", "diff_M", "delta_beta", "t")], by = "probe")
pg <- pg[order(pg$gene, pg$p_welch), ]
k  <- table(pg$gene)
gl <- pg[!duplicated(pg$gene), ]
names(gl)[names(gl) == "probe"] <- "best_probe"
gl$n_probes  <- as.integer(k[gl$gene])
gl$p_gene    <- 1 - (1 - gl$p_welch)^gl$n_probes                  # Sidak
gl$fdr_gene  <- p.adjust(gl$p_gene, "BH")
gl$score     <- -log10(gl$p_gene) * sign(gl$diff_M)               # + = more methylated in High
gl$direction <- ifelse(gl$diff_M > 0, "hyper_in_High", "hypo_in_High")
gl <- gl[order(gl$p_gene), ]

# methylation (best promoter probe) vs expression of the same gene, across matched lines
E <- attr(results$expression, "X"); colnames(E) <- sub(" \\([0-9]+\\)$", "", colnames(E))
ex_rows <- match(colnames(M), rownames(E))
gl$rho_meth_vs_expr <- NA_real_; gl$p_meth_vs_expr <- NA_real_
top_idx <- which(gl$fdr_gene < 0.25 & gl$gene %in% colnames(E))
for (i in top_idx) {
  m <- M[gl$best_probe[i], ]; e <- E[ex_rows, gl$gene[i]]; o <- !is.na(m) & !is.na(e)
  if (sum(o) >= 8) {
    ct <- suppressWarnings(cor.test(m[o], e[o], method = "spearman", exact = FALSE))
    gl$rho_meth_vs_expr[i] <- ct$estimate; gl$p_meth_vs_expr[i] <- ct$p.value
  }
}
write.table(gl, file.path(res_dir, "methylation_EPIC_gene_level.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

hit   <- gl$fdr_gene < 0.25 & abs(gl$delta_beta) >= 0.2
hypo  <- gl$gene[hit & gl$direction == "hypo_in_High"]
hyper <- gl$gene[hit & gl$direction == "hyper_in_High"]
en <- run_enrichment(hypo, hyper, gl$gene, setNames(gl$score, gl$gene),
                     file.path(res_dir, "methylation_EPIC_enrichment"))

# ---------------------------------------------- comparison with RRBS --------
rr <- read.delim(file.path(res_dir, "methylation_CpG_gene_level.tsv"))
sh <- merge(gl[, c("gene", "score", "fdr_gene", "delta_beta")],
            rr[, c("gene", "score", "fdr_gene", "best_diff")], by = "gene", suffixes = c("_epic", "_rrbs"))
both_sig <- sh$fdr_gene_epic < 0.25 & sh$fdr_gene_rrbs < 0.25
goi <- c("STAB1", "AZU1", "ELANE", "FGR", "SIRPA", "MPO", "TAL1", "GATA1", "GP1BB", "ERMAP",
         "KLF1", "PIK3CD", "PIK3CD-AS1", "NPRL3", "MAP1LC3A", "CEBPA", "IRF8", "SPI1")
cmp <- data.frame(gene = goi,
                  epic_dBeta = round(gl$delta_beta[match(goi, gl$gene)], 2),
                  epic_fdr   = signif(gl$fdr_gene[match(goi, gl$gene)], 2),
                  epic_probes = gl$n_probes[match(goi, gl$gene)],
                  rrbs_dBeta = round(rr$best_diff[match(goi, rr$gene)], 2),
                  rrbs_fdr   = signif(rr$fdr_gene[match(goi, rr$gene)], 2))
write.table(cmp, file.path(res_dir, "methylation_EPIC_vs_RRBS_genes_of_interest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- summary ---
sink(file.path(res_dir, "methylation_EPIC_summary.txt"))
cat("GSE270494 EPIC: lines matched", length(samp), "(High", sum(grp == "High"), "| Low", sum(grp == "Low"), ")\n")
cat("  High:", sort(lines_[grp == "High"]), "\n  Low: ", sort(lines_[grp == "Low"]), "\n")
cat(sprintf("Values failing detection p > %.2f: %d (%.2f%%)\n", det_p_max, n_failed,
            100 * n_failed / (n_raw * length(samp))))
cat("Probe QC:\n"); print(qc)
cat(sprintf("\nProbes tested: %d | FDR<0.05: %d | FDR<0.25: %d\n", nrow(res),
            sum(res$fdr_welch < 0.05), sum(res$fdr_welch < 0.25)))
cat(sprintf("  FDR<0.05: more methylated in High %d | less %d\n",
            sum(res$fdr_welch < 0.05 & res$diff_M > 0), sum(res$fdr_welch < 0.05 & res$diff_M < 0)))
cat(sprintf("  FDR<0.05 and |delta beta| >= 0.2: %d\n", sum(res$fdr_welch < 0.05 & abs(res$delta_beta) >= 0.2)))
cat(sprintf("Mean beta per line: High %.3f | Low %.3f\n",
            median(colMeans(B[, grp == "High"], na.rm = TRUE)), median(colMeans(B[, grp == "Low"], na.rm = TRUE))))
cat(sprintf("\nGenes with promoter probes: %d | gene FDR<0.25 & |delta beta|>=0.2: %d (hypo in High %d | hyper %d)\n",
            nrow(gl), sum(hit), length(hypo), length(hyper)))
s_ok <- !is.na(gl$rho_meth_vs_expr[hit])
cat(sprintf("  of those with expression data: %d; silencing pattern (rho < -0.4, p < 0.01): %d\n",
            sum(s_ok), sum(gl$rho_meth_vs_expr[hit] < -0.4 & gl$p_meth_vs_expr[hit] < 0.01, na.rm = TRUE)))
cat(sprintf("\nEPIC vs RRBS CpG clusters, gene level (%d shared genes): Spearman rho of signed score = %.2f\n",
            nrow(sh), cor(sh$score_epic, sh$score_rrbs, method = "spearman")))
cat(sprintf("  genes at FDR<0.25 in both: %d, same direction: %d\n", sum(both_sig),
            sum(sign(sh$score_epic[both_sig]) == sign(sh$score_rrbs[both_sig]))))
cat("\n== Genes of interest: EPIC (GSE270494) vs RRBS CpG clusters\n"); print(cmp, row.names = FALSE)
cat("\n== Top genes (promoter, gene level)\n")
print(head(data.frame(gene = gl$gene, probes = gl$n_probes, probe = gl$best_probe,
                      dBeta = round(gl$delta_beta, 2), dM = round(gl$diff_M, 2),
                      fdr = signif(gl$fdr_gene, 2), rho_expr = round(gl$rho_meth_vs_expr, 2)), 40), row.names = FALSE)
print_enrichment(en)
sink()

# ---------------------------------------------------------------- figures ---
# volcano on M values; non-significant probes thinned for drawing only
set.seed(1)
draw <- res$p_welch < 0.01 | seq_len(nrow(res)) %in% sample(nrow(res), 60000)
vr <- data.frame(feature = paste(res$gene, res$probe, sep = "|"), diff_High_minus_Low = res$diff_M,
                 p_welch = res$p_welch, fdr_welch = res$fdr_welch, n_High = res$n_High, n_Low = res$n_Low)[draw, ]
vr$feature[is.na(res$gene[draw])] <- res$probe[draw][is.na(res$gene[draw])]
volcano_png(vr, list(title = "Methylation EPIC (GSE270494), probe level", x = "M-value difference",
                     dep = FALSE, n_features = nrow(res), n_sig = sum(res$fdr_welch < 0.05),
                     always = "^(STAB1|AZU1|TAL1|GATA1|PIK3CD|NPRL3|MAP1LC3A)\\|"),
            file.path(fig_dir, "volcano_meth_EPIC.png"), n_lab = 14, xpad = 1.6)
plot_enrichment(en$ORA, length(hypo), length(hyper), file.path(fig_dir, "methylation_EPIC_enrichment.png"),
                title_suffix = "EPIC promoter, over-representation")

# EPIC vs RRBS gene-level agreement
png(file.path(fig_dir, "methylation_EPIC_vs_RRBS.png"), width = 1700, height = 1600, res = 220,
    type = "cairo", bg = "#fcfcfb")
par(mar = c(4.2, 4.2, 3.2, 1.2), fg = ink2)
x <- sh$score_rrbs; y <- sh$score_epic
lim <- max(abs(c(x, y)), na.rm = TRUE) * 1.05 * c(-1, 1)
col <- ifelse(both_sig & sign(x) == sign(y), ifelse(y > 0, col_hi, col_lo), col_ns)
plot(NA, xlim = lim, ylim = lim, axes = FALSE, xlab = "", ylab = "")
abline(h = pretty(lim), v = pretty(lim), col = grid_col, lwd = 0.8); abline(h = 0, v = 0, col = ink2, lwd = 0.8)
axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.8, tcl = -0.25, mgp = c(2, 0.4, 0))
axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.8, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
o <- order(col != col_ns)
points(x[o], y[o], pch = 21, bg = col[o], col = "#fcfcfb", lwd = 0.3, cex = ifelse(col[o] == col_ns, 0.45, 0.8))
lab <- which(sh$gene %in% goi)
text(x[lab], y[lab], sh$gene[lab], pos = 4, offset = 0.3, cex = 0.55, col = ink)
points(x[lab], y[lab], pch = 21, bg = NA, col = ink, lwd = 0.8, cex = 0.9)
mtext("Promoter methylation: EPIC array vs RRBS", side = 3, line = 1.6, adj = 0, font = 2, cex = 0.9, col = ink)
mtext(sprintf("%d shared genes · Spearman rho = %.2f · colored: FDR<0.25 in both, same direction (%d)",
              nrow(sh), cor(x, y, method = "spearman"), sum(both_sig & sign(x) == sign(y))),
      side = 3, line = 0.5, adj = 0, cex = 0.6, col = ink2)
mtext("RRBS CpG clusters: signed -log10 p  (+ = more methylated in mTOR-High)", side = 1, line = 2.2, cex = 0.65, col = ink2)
mtext("EPIC (GSE270494): signed -log10 p", side = 2, line = 2.4, cex = 0.65, col = ink2)
dev.off()

cat(readLines(file.path(res_dir, "methylation_EPIC_summary.txt")), sep = "\n")
