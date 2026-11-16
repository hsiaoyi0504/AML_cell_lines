# Metabolomics and global chromatin profiling volcanoes, MYC-High vs mTOR-High, side by side.
# Reads the tables written by scripts/metabolomics_MYC_vs_mTOR.R and scripts/chromatin_MYC_vs_mTOR.R
# (run those first). Positive difference = higher in MYC-High.
#
# Run from the project root: Rscript scripts/volcano_metab_chromatin.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R")
colA <- "#eb6834"; colB <- "#4a3aa7"

as_volcano <- function(st, feature, nA, nB) data.frame(feature = st[[feature]], diff_High_minus_Low = st$diff_MYC_minus_mTOR,
  p_welch = st$p_welch, fdr_welch = st$fdr_welch, n_High = st[[nA]], n_Low = st[[nB]])
met <- read.delim(file.path(res_dir, "metabolomics_MYC_High_vs_mTOR_High.tsv"))
chr <- read.delim(file.path(res_dir, "chromatin_MYC_High_vs_mTOR_High.tsv"))
chr <- chr[!is.na(chr$p_welch), ]

panels <- list(
  list(r = as_volcano(met, "metabolite", "n_MYC", "n_mTOR"), n_lab = sum(met$fdr_welch < 0.25), xpad = 2.2,
       p = list(title = "Metabolomics (225 metabolites)", x = "log10 abundance difference (MYC-High − mTOR-High)",
                dep = FALSE, fdr = 0.1, always = "^acetylcarnitine$", groups = c("MYC-High", "mTOR-High"), cols = c(colA, colB))),
  list(r = as_volcano(chr, "mark", "n_MYC", "n_mTOR"), n_lab = 10, xpad = 1.9,
       p = list(title = "Global chromatin profiling (42 histone marks)", x = "Relative abundance difference (MYC-High − mTOR-High)",
                dep = FALSE, fdr = 0.05, groups = c("MYC-High", "mTOR-High"), cols = c(colA, colB))))

draw <- function() {
  layout(matrix(1:4, 2), heights = c(1, 0.1))
  for (pn in panels) {
    par(mar = c(3.8, 3.4, 2.8, 1.2), fg = ink2)
    volcano(pn$r, pn$p, n_lab = pn$n_lab, xpad = pn$xpad)
    par(mar = c(0, 0, 0, 0)); plot.new()
    q <- vfdr(pn$p)
    legend("center", horiz = TRUE, bty = "n", pt.cex = 1.1, cex = 0.62, text.col = ink, pch = 21,
           pt.bg = c(vcols(pn$p), col_ns), col = "#fcfcfb",
           legend = c(paste0("MYC-High, FDR<", q), "FDR<.25", paste0("mTOR-High, FDR<", q), "FDR<.25", "n.s."))
  }
}
fb <- file.path(fig_dir, "volcano_metabolomics_and_chromatin_MYC_High_vs_mTOR_High")
png(paste0(fb, ".png"), width = 3000, height = 1400, res = 220, type = "cairo", bg = "#fcfcfb"); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3000 / 220, height = 1400 / 220, bg = "#fcfcfb"); draw(); dev.off()
