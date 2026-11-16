# Scheme panel: AML cell-line classification (RNA MYC / mTOR scores) followed by
# multi-omic characterisation of the groups. Line counts and headline findings are
# the ones produced by the analysis scripts in scripts/.
#
# Run from the project root: Rscript scripts/scheme_panel.R

source("scripts/plot_utils.R")
fig_dir <- "results/mTOR_RNA_groups/figures"
bg <- "#fcfcfb"; colM <- "#eb6834"; colT <- "#4a3aa7"; line_col <- "#c9c8c3"; card <- "#ffffff"
soft <- "#f3f2ee"

box <- function(x0, y0, x1, y1, fill = card, border = line_col, lwd = 1.2) {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = lwd)
}
arrow <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.08, angle = 25, col = ink2, lwd = 1.4)

# data layers: title, source, lines per group, comparison, headline
layers <- list(
  list("Transcriptome", "DepMap 26Q3 RNA-seq", "56 lines", "MYC-High vs mTOR-High",
       "Megakaryocytic/erythroid program (TAL1, GFI1B, ITGA2B) vs myelomonocytic program (AZU1, MNDA, IRF8)"),
  list("Mutations", "DepMap 26Q3 damaging mutations", "54 lines", "MYC-High vs mTOR-High",
       "TP53 mutated in both groups at similar rates (15/25 vs 12/23, p = 0.77)"),
  list("Genetic dependency", "CRISPR (Chronos) + RNAi (DEMETER2)", "16/16 · 13/10", "mTOR-High vs mTOR-Low",
       "Concordant in both screens: CEBPA, IRF8, MEF2D (mTOR-High); BCL2L1, GATA1 (mTOR-Low)"),
  list("Drug response", "CTRP · GDSC2 · PRISM (AUC)", "19/14 · 15/8 · 14/8", "mTOR-High vs mTOR-Low",
       "CTRP/GDSC2 agree: 5-FU, MDM2 inhibitors more active in mTOR-High; rapalog signal CTRP only"),
  list("DNA methylation", "CCLE RRBS · EPIC array (GSE270494)", "18/15 · 15/11", "mTOR-High vs mTOR-Low",
       "Myeloid genes unmethylated, erythroid/stem regulators methylated in mTOR-High; STAB1 replicates"),
  list("Histone marks", "CCLE global chromatin profiling", "13/19", "MYC-High vs mTOR-High",
       "Only H3K4me1 differs (higher in MYC-High, FDR 0.006) of 42 marks"),
  list("miRNA", "CCLE miRNA", "20/15", "mTOR-High vs mTOR-Low",
       "miR-155, miR-455-3p up; miR-126, miR-125b, miR-99a down in mTOR-High"),
  list("Metabolome", "CCLE targeted metabolomics", "20/17", "mTOR-High vs mTOR-Low",
       "Little difference: glucuronate only (FDR 0.04)")
)

draw <- function() {
  par(mar = c(0, 0, 0, 0), bg = bg, family = "sans")
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  text(2, 97.2, "Characterising MYC-High and mTOR-High AML cell lines", adj = c(0, 0.5), font = 2, cex = 1.55, col = ink)
  text(2, 94.4, "Classification from RNA signature scores, then multi-omic comparison of the groups", adj = c(0, 0.5),
       cex = 0.95, col = ink2)

  # ---- 1. classification -------------------------------------------------
  text(2, 89.8, "1  Classify", adj = c(0, 0.5), font = 2, cex = 1.1, col = ink)
  box(2, 71, 30, 87.6, fill = soft)
  text(4, 85.4, "56 AML cell lines", adj = c(0, 0.5), font = 2, cex = 1.05, col = ink)
  text(4, 82.4, "DepMap 26Q1 RNA-seq", adj = c(0, 0.5), cex = 0.85, col = ink2)
  text(4, 78.8, "MYC signature score (98 genes)", adj = c(0, 0.5), cex = 0.85, col = ink)
  text(4, 75.8, "mTOR signature score (267 genes)", adj = c(0, 0.5), cex = 0.85, col = ink)
  text(4, 72.8, "KMT2A rearrangement annotated (n = 10)", adj = c(0, 0.5), cex = 0.85, col = ink2)
  arrow(30.5, 79.3, 35.5, 79.3)

  # score plane cartoon
  px0 <- 37; px1 <- 55; py0 <- 71; py1 <- 87.6
  box(px0, py0, px1, py1)
  segments((px0 + px1) / 2, py0 + 1, (px0 + px1) / 2, py1 - 1, col = line_col, lty = 2)
  segments(px0 + 1, (py0 + py1) / 2, px1 - 1, (py0 + py1) / 2, col = line_col, lty = 2)
  set.seed(4)
  pts <- rbind(cbind(runif(11, 38.5, 44.5), runif(11, 80.5, 85.5), 2), cbind(runif(11, 47.5, 53.5), runif(11, 73, 78), 1))
  points(pts[, 1], pts[, 2], pch = 21, cex = 1.15, col = bg, bg = ifelse(pts[, 3] == 1, colM, colT))
  text(53.8, 79.6, "MYC-High", adj = c(1, 0.5), cex = 0.78, col = colM, font = 2)
  text(38.2, 79.0, "mTOR-High", adj = c(0, 0.5), cex = 0.78, col = colT, font = 2)
  text((px0 + px1) / 2, py0 - 1.5, "MYC score →", cex = 0.72, col = ink2)
  text(px0 - 1.2, (py0 + py1) / 2, "mTOR score →", srt = 90, cex = 0.72, col = ink2)
  arrow(55.5, 79.3, 60.5, 79.3)

  # groups
  box(61, 80.2, 98, 87.6, fill = colM, border = NA)
  text(63, 85.3, "MYC-High / mTOR-Low   n = 25", adj = c(0, 0.5), font = 2, cex = 1, col = "#ffffff")
  text(63, 82.3, "mostly megakaryocytic / erythroid / stem-like lines; 1 KMT2A-r (OCI-AML4)", adj = c(0, 0.5), cex = 0.8, col = "#ffffff")
  box(61, 71, 98, 78.4, fill = colT, border = NA)
  text(63, 76.1, "mTOR-High / MYC-Low   n = 25", adj = c(0, 0.5), font = 2, cex = 1, col = "#ffffff")
  text(63, 73.1, "mostly myelomonocytic lines; 9 of 10 KMT2A-r", adj = c(0, 0.5), cex = 0.8, col = "#ffffff")
  text(98, 69.3, "6 lines high or low for both scores left out of MYC vs mTOR comparisons", adj = c(1, 0.5), cex = 0.72,
       col = ink2, font = 3)

  # ---- 2. multi-omic layers ------------------------------------------------
  text(2, 64.6, "2  Characterise with multi-omic data", adj = c(0, 0.5), font = 2, cex = 1.1, col = ink)
  arrow(50, 68.5, 50, 66.2)
  cw <- 23; gap <- 1.33; ch <- 22.5
  for (k in seq_along(layers)) {
    L <- layers[[k]]; col_i <- (k - 1) %% 4; row_i <- (k - 1) %/% 4
    x0 <- 2 + col_i * (cw + gap); y1 <- 62 - row_i * (ch + 1.6); y0 <- y1 - ch
    box(x0, y0, x0 + cw, y1)
    right <- if (grepl("^MYC", L[[4]])) colT else "#c9c8c3"; left <- if (grepl("^MYC", L[[4]])) colM else colT
    rect(x0, y1 - 0.9, x0 + cw / 2, y1, col = left, border = NA); rect(x0 + cw / 2, y1 - 0.9, x0 + cw, y1, col = right, border = NA)
    text(x0 + 1.2, y1 - 3.4, L[[1]], adj = c(0, 0.5), font = 2, cex = 1, col = ink)
    text(x0 + 1.2, y1 - 6.2, L[[2]], adj = c(0, 0.5), cex = 0.74, col = ink2)
    text(x0 + 1.2, y1 - 8.8, sprintf("%s  ·  %s", L[[4]], L[[3]]), adj = c(0, 0.5), cex = 0.68, col = ink2, font = 3)
    wr <- strwrap(L[[5]], width = 40)
    for (i in seq_along(wr)) text(x0 + 1.2, y1 - 12.4 - (i - 1) * 2.7, wr[i], adj = c(0, 0.5), cex = 0.8, col = ink)
  }

  # ---- 3. integration ------------------------------------------------------
  yb1 <- 9.6
  text(2, yb1 + 2.2, "3  Integrate", adj = c(0, 0.5), font = 2, cex = 1.1, col = ink)
  box(2, 1.2, 98, yb1, fill = soft)
  text(4, 7.2, "Cross-screen concordance  ·  methylation vs expression  ·  case studies", adj = c(0, 0.5),
       font = 2, cex = 0.9, col = ink)
  text(4, 3.7, paste("Overall: the groups separate mainly by lineage.  OCI-AML4 (KMT2A-MLLT1) is the outlier:",
                     "EVI1/MECOM-high, stem/megakaryocytic, weak KMT2A dependency.  ZMIZ1: intronic methylation tracks group, not RNA."),
       adj = c(0, 0.5), cex = 0.78, col = ink2)

  # key
  rect(66, 63.6, 67, 65.2, col = colM, border = NA); rect(67, 63.6, 68, 65.2, col = colT, border = NA)
  text(68.6, 64.4, "compares MYC-High vs mTOR-High", adj = c(0, 0.5), cex = 0.72, col = ink2)
  rect(83.5, 63.6, 84.5, 65.2, col = colT, border = NA); rect(84.5, 63.6, 85.5, 65.2, col = "#c9c8c3", border = NA)
  text(86.1, 64.4, "compares mTOR-High vs Low", adj = c(0, 0.5), cex = 0.72, col = ink2)
}
fb <- file.path(fig_dir, "scheme_panel")
png(paste0(fb, ".png"), width = 3600, height = 2700, res = 220, type = "cairo", bg = bg); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3600 / 220, height = 2700 / 220, bg = bg); draw(); dev.off()
