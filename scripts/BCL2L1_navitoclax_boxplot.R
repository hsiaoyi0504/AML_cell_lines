# BCL2L1: dependency (CRISPR Chronos, RNAi DEMETER2) and navitoclax response (CTRP, GDSC2 AUC),
# MYC-High vs mTOR-High. More negative gene effect / lower AUC = more dependent / sensitive.
# p = Welch t-test (BCL2L1 is a pre-specified target, so no multiple-testing or sensitivity adjustment).
#
# Run from the project root: Rscript scripts/BCL2L1_navitoclax_boxplot.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp <- ifelse(cl$MYC_group == "High" & cl$mTOR_group == "Low", "MYC-High",
          ifelse(cl$mTOR_group == "High" & cl$MYC_group == "Low", "mTOR-High", NA))

colA <- "#eb6834"; colB <- "#4a3aa7"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
tint <- function(col, f = 0.75) colorRampPalette(c(col, "#ffffff"))(100)[round(f * 100)]

col_of <- function(key, name) { X <- attr(results[[key]], "X")
  j <- if (grepl("^drug", key)) results[[key]]$feature[which(results[[key]]$CompoundName == name)[1]]
       else grep(paste0("^", name, " \\("), colnames(X), value = TRUE)[1]
  X[, j] }
panels <- list(
  list(v = col_of("crispr", "BCL2L1"), title = "BCL2L1 · CRISPR", y = "Chronos gene effect", ref = -1),
  list(v = col_of("rnai", "BCL2L1"), title = "BCL2L1 · RNAi", y = "DEMETER2 score", ref = NA),
  list(v = col_of("drug_CTRP", "NAVITOCLAX"), title = "Navitoclax · CTRP", y = "AUC", ref = NA),
  list(v = col_of("drug_GDSC2", "NAVITOCLAX"), title = "Navitoclax · GDSC2", y = "AUC", ref = NA))

draw <- function() {
  par(mfrow = c(1, 4), mar = c(3.2, 4.4, 3.6, 0.8), oma = c(0.5, 0.5, 3.6, 0.5), fg = ink2, bg = bg)
  set.seed(1)
  for (p in panels) {
    v <- p$v; g <- factor(cl$grp[match(names(v), cl$cell_line)], levels = c("MYC-High", "mTOR-High"))
    ok <- !is.na(v) & !is.na(g); v <- v[ok]; g <- g[ok]; x <- as.integer(g)
    tt <- t.test(v[x == 1], v[x == 2])
    yr <- range(c(v, p$ref), na.rm = TRUE); yr <- yr + c(-0.06, 0.08) * diff(yr)
    plot(NA, xlim = c(0.4, 2.6), ylim = yr, axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(yr), col = grid_col, lwd = 0.8)
    if (!is.na(p$ref)) abline(h = p$ref, lty = 3, col = ink2, lwd = 0.8)
    boxplot(v ~ g, at = 1:2, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.55, lwd = 1, medlwd = 2,
            col = c(tint(colA), tint(colB)), border = c(colA, colB), whisklty = 1, staplewex = 0.4)
    points(x + runif(length(v), -0.13, 0.13), v, pch = 21, cex = 1.25, lwd = 0.5, col = bg, bg = c(colA, colB)[x])
    axis(2, col = ink, col.ticks = ink, col.axis = ink, las = 1, cex.axis = 0.9, tcl = -0.3, lwd = 0.8)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", levels(g), tabulate(x, 2)), col = NA, col.axis = ink,
         cex.axis = 0.9, padj = 0.6, mgp = c(2, 0.4, 0))
    mtext(p$y, side = 2, line = 2.9, cex = 0.72, col = ink)
    mtext(p$title, side = 3, line = 1.5, adj = 0, font = 2, cex = 0.9, col = ink)
    mtext(sprintf("mean %.2f vs %.2f · p = %.2g", mean(v[x == 1]), mean(v[x == 2]), tt$p.value), side = 3, line = 0.4,
          adj = 0, cex = 0.7, col = if (tt$p.value < 0.05) ink else ink2, font = if (tt$p.value < 0.05) 2 else 1)
  }
  mtext("BCL2L1 (BCL-XL): genetic dependency and navitoclax response, MYC-High vs mTOR-High", outer = TRUE, side = 3,
        line = 2, adj = 0.005, font = 2, cex = 1.05, col = ink)
  mtext("Lower = more dependent / more sensitive · dotted = Chronos −1 (median common-essential gene) · p = Welch t-test",
        outer = TRUE, side = 3, line = 0.6, adj = 0.005, cex = 0.75, col = ink2)
}
fb <- file.path(fig_dir, "BCL2L1_navitoclax_MYC_vs_mTOR_boxplot")
png(paste0(fb, ".png"), width = 3200, height = 1250, res = 220, type = "cairo", bg = bg); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3200 / 220, height = 1250 / 220, bg = bg); draw(); dev.off()
