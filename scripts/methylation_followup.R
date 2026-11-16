# Follow-up of differential TSS methylation (mTOR-High vs mTOR-Low, RNA groups):
#   * does promoter methylation track expression of the same gene (silencing)?
#   * does it track the continuous mTOR score?
#   * methylation of mTOR-pathway / GATOR1 genes
# Reads results written by scripts/compare_mTOR_RNA_groups.R. Base R only.
#
# Usage: Rscript scripts/methylation_followup.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"
fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))

me <- results$meth; M <- attr(me, "X")
E  <- attr(results$expression, "X"); colnames(E) <- sub(" \\([0-9]+\\)$", "", colnames(E))
cr <- results$crispr; C <- attr(cr, "X"); colnames(C) <- sub(" \\([0-9]+\\)$", "", colnames(C))
me$gene <- sub("\\|.*$", "", me$feature)
score <- setNames(cl$mTOR_score, cl$cell_line)

scor <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  if (sum(ok) < 8) return(c(rho = NA, p = NA, n = sum(ok)))
  t <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE))
  c(rho = unname(t$estimate), p = t$p.value, n = sum(ok))
}

# --------------------------------------- annotate every tested TSS region ---
ann <- t(vapply(seq_len(nrow(me)), function(i) {
  f <- me$feature[i]; g <- me$gene[i]; m <- M[, f]
  ex <- if (g %in% colnames(E)) scor(m, E[match(names(m), rownames(E)), g]) else c(NA, NA, NA)
  sc <- scor(m, score[names(m)])
  c(ex, sc)
}, numeric(6)))
colnames(ann) <- c("rho_meth_vs_expr", "p_meth_vs_expr", "n_meth_expr",
                   "rho_meth_vs_mTORscore", "p_meth_vs_mTORscore", "n_meth_score")
me <- cbind(me, ann)
me$silencing <- !is.na(me$rho_meth_vs_expr) & me$rho_meth_vs_expr < -0.4 & me$p_meth_vs_expr < 0.01
# expression direction implied by methylation (lower methylation -> higher expression)
me$expr_diff_High_minus_Low <- NA
hasE <- me$gene %in% colnames(E)
g_e  <- cl$mTOR_group[match(rownames(E), cl$cell_line)]
me$expr_diff_High_minus_Low[hasE] <- vapply(me$gene[hasE], function(g)
  mean(E[g_e == "High", g]) - mean(E[g_e == "Low", g]), 0)

keep <- c("feature", "gene", "n_High", "n_Low", "mean_High", "mean_Low", "diff_High_minus_Low",
          "p_welch", "fdr_welch", "rho_meth_vs_expr", "p_meth_vs_expr", "silencing",
          "expr_diff_High_minus_Low", "rho_meth_vs_mTORscore", "p_meth_vs_mTORscore")
out <- me[order(me$p_welch), keep]
write.table(out, file.path(res_dir, "methylation_followup_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
hits <- out[out$fdr_welch < 0.25, ]
write.table(hits, file.path(res_dir, "methylation_followup_FDR0.25.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------- mTOR pathway genes ---
pathway <- c("NPRL2", "NPRL3", "DEPDC5", "TSC1", "TSC2", "TBC1D7", "DEPTOR", "PTEN", "STK11",
             "SESN1", "SESN2", "SESN3", "CASTOR1", "SZT2", "KPTN", "ITFG2", "MIOS", "WDR24",
             "WDR59", "SEH1L", "SEC13", "MTOR", "RPTOR", "RICTOR", "MLST8", "AKT1S1", "RHEB",
             "LAMTOR1", "LAMTOR2", "RRAGA", "RRAGB", "RRAGC", "RRAGD", "FLCN", "PIK3CA",
             "PIK3CD", "PIK3CD-AS1", "PIK3CD-AS2", "PIK3R1", "INPP4B", "DDIT4", "EIF4EBP1",
             "RPS6KB1", "ULK1", "MAP1LC3A", "MAP1LC3B", "TFEB", "SLC7A5", "SLC38A9", "MYC")
pw <- out[out$gene %in% pathway, ]
write.table(pw, file.path(res_dir, "methylation_followup_mTOR_pathway.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- report ----
sink(file.path(res_dir, "methylation_followup_summary.txt"))
cat("TSS regions tested:", nrow(out), "| FDR<0.05:", sum(out$fdr_welch < 0.05),
    "| FDR<0.25:", nrow(hits), "\n")
cat("FDR<0.25 regions: more methylated in mTOR-High:", sum(hits$diff_High_minus_Low > 0),
    "| less methylated in mTOR-High:", sum(hits$diff_High_minus_Low < 0), "\n")
cat("FDR<0.25 regions with expression data:", sum(!is.na(hits$rho_meth_vs_expr)),
    "| with silencing pattern (rho < -0.4, p < 0.01):", sum(hits$silencing), "\n")
ok <- !is.na(hits$expr_diff_High_minus_Low)
cat("Methylation and expression move in opposite directions between groups for",
    sum(sign(hits$diff_High_minus_Low[ok]) != sign(hits$expr_diff_High_minus_Low[ok])), "of", sum(ok), "\n\n")
f <- function(d) print(data.frame(gene = d$gene, region = sub("^.*\\|", "", d$feature),
  methH = round(d$mean_High, 2), methL = round(d$mean_Low, 2), fdr = signif(d$fdr_welch, 2),
  rho_expr = round(d$rho_meth_vs_expr, 2), p_expr = signif(d$p_meth_vs_expr, 2),
  dExpr = round(d$expr_diff_High_minus_Low, 2), rho_score = round(d$rho_meth_vs_mTORscore, 2)),
  row.names = FALSE)
cat("== FDR<0.25 regions with a silencing pattern\n"); f(hits[hits$silencing, ])
cat("\n== mTOR pathway / GATOR1 genes (all TSS regions tested, by p)\n"); f(head(pw, 25))
sink()

# ---------------------------------------------------------------- figure ----
col_hi <- "#2a78d6"; col_lo <- "#e34948"; ink <- "#0b0b0b"; ink2 <- "#52514e"
grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
pick <- c("NPRL3", "TAL1", "GP1BB", "ERMAP", "AZU1", "ELANE", "FGR", "MAP1LC3A")
regs <- vapply(pick, function(g) hits$feature[hits$gene == g][1], "")
regs <- regs[!is.na(regs)]
png(file.path(fig_dir, "methylation_vs_expression.png"), width = 3200, height = 1650, res = 200,
    type = "cairo", bg = bg)
layout(rbind(matrix(1:8, 2, byrow = TRUE), 9), heights = c(1, 1, 0.12))
par(mar = c(3.6, 3.8, 3, 1), fg = ink2)
for (k in seq_along(regs)) {
  f0 <- regs[k]; g <- names(regs)[k]
  m <- M[, f0]; e <- E[match(names(m), rownames(E)), g]
  grp <- cl$mTOR_group[match(names(m), cl$cell_line)]
  ok <- !is.na(m) & !is.na(e)
  plot(NA, xlim = c(0, 1), ylim = range(e[ok]) + c(-0.05, 0.08) * diff(range(e[ok])),
       axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(e[ok]), v = seq(0, 1, 0.25), col = grid_col, lwd = 0.8)
  axis(1, at = seq(0, 1, 0.25), col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
  points(m[ok], e[ok], pch = 21, cex = 1.3, lwd = 0.5, col = bg, bg = ifelse(grp[ok] == "High", col_hi, col_lo))
  r <- hits[hits$feature == f0, ]
  mtext(g, side = 3, line = 1.3, adj = 0, font = 2, cex = 0.8, col = ink)
  mtext(sprintf("High vs Low methylation FDR = %.2g · meth~expr rho = %.2f", r$fdr_welch, r$rho_meth_vs_expr),
        side = 3, line = 0.35, adj = 0, cex = 0.58, col = ink2)
  mtext("TSS methylation (beta)", side = 1, line = 2, cex = 0.6, col = ink2)
  mtext("log2(TPM+1)", side = 2, line = 2.4, cex = 0.6, col = ink2)
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", pch = 21, pt.cex = 1.4, cex = 0.85, text.col = ink, col = bg,
       pt.bg = c(col_hi, col_lo), legend = c("mTOR-High", "mTOR-Low"))
dev.off()
cat(readLines(file.path(res_dir, "methylation_followup_summary.txt")), sep = "\n")
