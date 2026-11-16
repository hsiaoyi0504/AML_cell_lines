# Volcano plots + top-hit strip plots for mTOR-High vs mTOR-Low (RNA groups).
# Reads results written by scripts/compare_mTOR_RNA_groups.R. Base R only.
#
# Usage: Rscript scripts/plot_mTOR_RNA_groups.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

results <- readRDS(file.path(res_dir, "all_results.rds"))
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))

source("scripts/plot_utils.R")   # palette, spread(), volcano(), volcano_png()

# RNA expression is left out: the mTOR groups were derived from it (circular)
panels <- list(
  crispr            = list(title = "CRISPR (Chronos)",     x = "Gene effect difference", dep = TRUE),
  rnai               = list(title = "RNAi (DEMETER2)",      x = "Gene effect difference", dep = TRUE),
  drug_CTRP          = list(title = "Drug AUC (CTRP)",      x = "AUC difference",         dep = TRUE),
  drug_PRISM         = list(title = "Drug AUC (PRISM)",     x = "AUC difference",         dep = TRUE),
  drug_GDSC2         = list(title = "Drug AUC (GDSC2)",     x = "AUC difference",         dep = TRUE),
  meth               = list(title = "Methylation (TSS 1kb)", x = "Beta difference",       dep = FALSE,
                            always = "^(NPRL3|TAL1|PIK3CD-AS1|MAP1LC3A)\\|"),
  meth_cpg           = list(title = "Methylation (CpG clusters)", x = "Beta difference",  dep = FALSE,
                            always = "^(NPRL3|TAL1|PIK3CD-AS1|MAP1LC3A|GATA1)\\|"),
  mirna              = list(title = "miRNA",                x = "log2 difference",        dep = FALSE),
  chrom              = list(title = "Histone marks",        x = "Abundance difference",   dep = FALSE),
  metab              = list(title = "Metabolomics",         x = "log10 abundance difference", dep = FALSE, fdr = 0.1)
)

legend_row <- function() {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pt.cex = 1.3, cex = 0.8, text.col = ink,
         pch = 21, pt.bg = c(col_hi, tint_hi, col_lo, tint_lo, col_ns), col = "#fcfcfb",
         legend = c("Favours mTOR-High, FDR<0.05", "FDR<0.25",
                    "Favours mTOR-Low, FDR<0.05", "FDR<0.25", "n.s."))
}

draw_grid <- function(keys, file, ncol, w, h) {
  png(file, width = w, height = h, res = 200, type = "cairo", bg = "#fcfcfb")
  nrow <- ceiling(length(keys) / ncol)
  layout(rbind(matrix(seq_len(nrow * ncol), nrow, byrow = TRUE), nrow * ncol + 1),
         heights = c(rep(1, nrow), 0.12))
  par(mar = c(3.6, 3.2, 2.6, 1), family = "sans", fg = ink2)
  for (k in keys) volcano(results[[k]], panels[[k]])
  for (i in seq_len(nrow * ncol - length(keys))) plot.new()
  legend_row()
  dev.off()
}

# 1. overview: all datasets
draw_grid(names(panels), file.path(fig_dir, "volcano_all_datasets.png"), ncol = 3, w = 2900,
          h = 960 * ceiling(length(panels) / 3))

# 2. one large volcano per dataset
for (k in names(panels)) {
  volcano_png(results[[k]], panels[[k]], file.path(fig_dir, sprintf("volcano_%s.png", k)))
}

# 3. strip plots of selected hits, one dot per cell line
{
  mats <- lapply(results, attr, "X")
  picks <- list(c("meth", "STAB1|STAB1_3_52528355_52529355"), c("mirna", "hsa-miR-126"),
                c("crispr", "BCL2L1 (598)"), c("rnai", "MASTL (84930)"),
                c("rnai", "CEBPA (1050)"), c("drug_CTRP", "DPC-007478"),
                c("drug_CTRP", "DPC-005166"), c("mirna", "hsa-miR-155"),
                c("chrom", "H3K4me1"), c("metab", "glucuronate"))
  png(file.path(fig_dir, "top_hits_strip.png"), width = 3600, height = 1600, res = 200,
      type = "cairo", bg = "#fcfcfb")
  par(mfrow = c(2, 5), mar = c(2.6, 3.4, 3, 1), fg = ink2)
  set.seed(1)
  for (pk in picks) {
    m <- mats[[pk[1]]]; v <- m[, pk[2]]
    g <- cl$mTOR_group[match(rownames(m), cl$cell_line)]
    ok <- !is.na(v); v <- v[ok]; g <- factor(g[ok], levels = c("Low", "High"))
    r <- results[[pk[1]]]; rr <- r[r$feature == pk[2], ]
    plot(NA, xlim = c(0.4, 2.6), ylim = range(v) + c(-0.05, 0.12) * diff(range(v)),
         axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(v), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.8, tcl = -0.25)
    mtext(sub(" difference$", "", panels[[pk[1]]]$x), side = 2, line = 2.4, cex = 0.6, col = ink2)
    axis(1, at = 1:2, labels = c(sprintf("mTOR-Low\nn=%d", sum(g == "Low")),
                                 sprintf("mTOR-High\nn=%d", sum(g == "High"))),
         col = NA, col.axis = ink, cex.axis = 0.8, padj = 0.5)
    xs <- as.integer(g) + runif(length(v), -0.18, 0.18)
    points(xs, v, pch = 21, cex = 1.15, lwd = 0.5, col = "#fcfcfb",
           bg = ifelse(g == "High", col_hi, col_lo))
    for (j in 1:2) {
      md <- median(v[as.integer(g) == j])
      segments(j - 0.3, md, j + 0.3, md, lwd = 2, col = ink)
    }
    mtext(sprintf("%s · %s", label_of(rr, 1), panels[[pk[1]]]$title), side = 3, line = 1.2,
          adj = 0, font = 2, cex = 0.72, col = ink)
    mtext(sprintf("Welch p = %.1e · FDR = %.2g", rr$p_welch, rr$fdr_welch), side = 3,
          line = 0.3, adj = 0, cex = 0.6, col = ink2)
  }
  dev.off()
}
cat("figures ->", fig_dir, "\n"); print(list.files(fig_dir))
