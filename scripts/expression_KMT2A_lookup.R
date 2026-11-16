# Expression of one gene in KMT2A-rearranged vs other AML lines (the 56 RNA-classified
# lines; KMT2A status from fig2_MYC_mTOR_DepMap_all_cell_line_classification.tsv).
# Also repeated within mTOR-High lines, since KMT2A-r lines are almost all mTOR-High.
#
# Run from the project root: Rscript scripts/expression_KMT2A_lookup.R  (set the gene below)

gene <- "MYPOP"                     # gene to test
res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
out_dir <- file.path(res_dir, "gene_lookup"); dir.create(out_dir, showWarnings = FALSE)
source("scripts/plot_utils.R")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
results <- readRDS(file.path(res_dir, "all_results.rds"))
E <- attr(results$expression, "X"); colnames(E) <- sub(" \\([0-9]+\\)$", "", colnames(E))

if (!gene %in% colnames(E)) {   # the stored matrix keeps genes expressed in >= 5 lines; fall back to the raw file
  hdr <- strsplit(readLines("data/cell_lines/DepMap_26Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", n = 1), ",")[[1]]
  j <- which(sub(" \\([0-9]+\\)$", "", hdr) == gene)
  if (!length(j)) stop(gene, " not found in DepMap expression data")
  ex <- read.csv(pipe(paste0("cut -d, -f4,5,6,", j, " data/cell_lines/DepMap_26Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")))
  ex <- ex[ex$IsDefaultEntryForModel == "Yes" & ex$ModelID %in% cl$ModelID, ]
  v <- setNames(ex[[4]], cl$cell_line[match(ex$ModelID, cl$ModelID)])
} else v <- E[, gene]

d <- data.frame(cell_line = names(v), expr = v)
d <- merge(d, cl[, c("cell_line", "KMT2A_rearranged", "KMT2A_partner", "mTOR_group", "OncotreeSubtype")])
d$KMT2A <- ifelse(d$KMT2A_rearranged, "KMT2A-r", "other")

cmp <- function(x, lab) {
  a <- x$expr[x$KMT2A == "KMT2A-r"]; b <- x$expr[x$KMT2A == "other"]
  tt <- t.test(a, b); wt <- suppressWarnings(wilcox.test(a, b, exact = FALSE))
  sdp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  data.frame(comparison = lab, n_KMT2Ar = length(a), n_other = length(b),
             mean_KMT2Ar = mean(a), mean_other = mean(b), median_KMT2Ar = median(a), median_other = median(b),
             diff_log2 = mean(a) - mean(b), cohens_d = (mean(a) - mean(b)) / sdp,
             p_welch = tt$p.value, p_wilcox = wt$p.value)
}
res <- rbind(cmp(d, "all 56 AML lines"), cmp(d[d$mTOR_group == "High", ], "within mTOR-High lines"))
# rank of the gene among all genes in the same comparison (genome-wide context)
gk <- d$KMT2A[match(rownames(E), d$cell_line)]
t_all <- apply(E, 2, function(e) { a <- e[gk == "KMT2A-r"]; b <- e[gk == "other"]
  if (var(a) + var(b) == 0) return(NA); t.test(a, b)$p.value })
fdr_all <- p.adjust(t_all, "BH")
res$genome_wide_fdr <- c(if (gene %in% names(fdr_all)) unname(fdr_all[gene]) else NA, NA)
res$rank_of_genes <- c(if (gene %in% names(t_all)) sum(t_all <= t_all[gene], na.rm = TRUE) else NA, NA)
res$n_genes <- c(sum(!is.na(t_all)), NA)

write.table(d[order(d$KMT2A, -d$expr), ], file.path(out_dir, paste0(gene, "_expression_KMT2A_lines.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(out_dir, paste0(gene, "_expression_KMT2A_test.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- figure ---
draw <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4.2, 3.4, 1), fg = ink2)
  for (k in 1:2) {
    x <- if (k == 1) d else d[d$mTOR_group == "High", ]
    r <- res[k, ]; f <- factor(x$KMT2A, levels = c("other", "KMT2A-r"))
    set.seed(1); xx <- as.integer(f) + runif(nrow(x), -0.2, 0.2)
    yl <- range(d$expr) + c(-0.05, 0.1) * diff(range(d$expr))
    plot(NA, xlim = c(0.4, 2.6), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(yl), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.8, tcl = -0.25)
    axis(1, at = 1:2, labels = c(sprintf("other\nn=%d", sum(f == "other")), sprintf("KMT2A-r\nn=%d", sum(f == "KMT2A-r"))),
         col = NA, col.axis = ink, cex.axis = 0.8, padj = 0.5)
    # fill = mTOR group (consistent with the other figures); KMT2A-r points ringed in ink
    points(xx, x$expr, pch = 21, cex = 1.3, lwd = ifelse(x$KMT2A == "KMT2A-r", 1.4, 0.5),
           col = ifelse(x$KMT2A == "KMT2A-r", ink, "#fcfcfb"), bg = ifelse(x$mTOR_group == "High", col_hi, col_lo))
    kr <- x$KMT2A == "KMT2A-r"
    text(xx[kr], x$expr[kr], paste0(x$cell_line[kr], " (", x$KMT2A_partner[kr], ")"), pos = 4, cex = 0.5, col = ink2, offset = 0.4)
    for (j in 1:2) segments(j - 0.32, median(x$expr[as.integer(f) == j]), j + 0.32, lwd = 2, col = ink)
    mtext(sprintf("%s · %s", gene, r$comparison), side = 3, line = 1.6, adj = 0, font = 2, cex = 0.85, col = ink)
    mtext(sprintf("Δ = %+.2f log2 · Welch p = %.2g · Wilcoxon p = %.2g%s", r$diff_log2, r$p_welch, r$p_wilcox,
                  if (!is.na(r$genome_wide_fdr)) sprintf(" · genome-wide FDR = %.2g", r$genome_wide_fdr) else ""),
          side = 3, line = 0.5, adj = 0, cex = 0.62, col = ink2)
    mtext("log2(TPM+1)", side = 2, line = 2.6, cex = 0.65, col = ink2)
  }
  legend("bottomright", bty = "n", pch = 21, pt.cex = 1.2, cex = 0.62, text.col = ink,
         pt.bg = c(col_hi, col_lo, "#fcfcfb"), col = c("#fcfcfb", "#fcfcfb", ink), pt.lwd = c(0.5, 0.5, 1.4),
         legend = c("mTOR-High", "mTOR-Low", "KMT2A-r (ring)"))
}
fb <- file.path(fig_dir, paste0("gene_", gene, "_expression_KMT2A"))
png(paste0(fb, ".png"), width = 2600, height = 1300, res = 220, type = "cairo", bg = "#fcfcfb"); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2600 / 220, height = 1300 / 220, bg = "#fcfcfb"); draw(); dev.off()

print(res, row.names = FALSE)
cat("\nPer line (KMT2A-r first):\n")
print(d[order(d$KMT2A != "KMT2A-r", -d$expr), c("cell_line", "KMT2A", "KMT2A_partner", "mTOR_group", "expr")], row.names = FALSE)
