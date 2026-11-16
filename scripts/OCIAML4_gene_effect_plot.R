# CRISPR gene effect: OCI-AML4 vs the other KMT2A-rearranged lines (n=8 with CRISPR).
# Outlier = outside the others' range, |diff| >= 0.5 and |z| >= 3 (sd floor 0.15), as in
# scripts/OCIAML4_vs_KMT2Ar.R. Calibration: the same test run on each other KMT2A-r line
# against the remaining 7, to see how many outliers one line gets by chance.
#
# Run from the project root: Rscript scripts/OCIAML4_gene_effect_plot.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures"); out_dir <- file.path(res_dir, "OCIAML4")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
C <- attr(results$crispr, "X"); colnames(C) <- sub(" \\([0-9]+\\)$", "", colnames(C))
kr <- intersect(cl$cell_line[cl$KMT2A_rearranged %in% TRUE], rownames(C))
target <- "OCIAML4"; others <- setdiff(kr, target)
partner <- setNames(cl$KMT2A_partner, cl$cell_line)

colA <- "#eb6834"; colB <- "#4a3aa7"; col_ns <- "#c9c8c3"
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"

outlier_call <- function(line, ref) {
  R <- C[ref, , drop = FALSE]; x <- C[line, ]
  m <- colMeans(R, na.rm = TRUE); s <- pmax(apply(R, 2, sd, na.rm = TRUE), 0.15)
  lo <- apply(R, 2, min, na.rm = TRUE); hi <- apply(R, 2, max, na.rm = TRUE)
  d <- x - m; z <- d / s
  out <- (x < lo | x > hi) & abs(d) >= 0.5 & abs(z) >= 3; out[is.na(out)] <- FALSE
  # Crawford-Howell single-case t-test: is x outside what n reference lines predict? (same sd floor)
  n <- colSums(!is.na(R)); t <- d / (s * sqrt(1 + 1 / n)); p <- 2 * pt(-abs(t), n - 1)
  data.frame(gene = colnames(C), value = x, mean_ref = m, diff = d, z = z, outlier = out,
             t_CH = t, p_CH = p, fdr_CH = p.adjust(p, "BH"))
}
oc <- outlier_call(target, others)
write.table(oc[order(oc$p_CH), ], file.path(out_dir, "OCIAML4_vs_KMT2Ar_CRISPR_CrawfordHowell.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# calibration: each other KMT2A-r line vs the remaining 7 (OCI-AML4 left out of the reference)
cnt <- function(o) c(more_dep = sum(o$outlier & o$diff < 0), less_dep = sum(o$outlier & o$diff > 0),
                     p001 = sum(o$p_CH < 0.001, na.rm = TRUE), fdr10 = sum(o$fdr_CH < 0.1, na.rm = TRUE))
cal <- do.call(rbind, lapply(others, function(l)
  data.frame(line = l, partner = partner[l], t(cnt(outlier_call(l, setdiff(others, l)))))))
o7 <- outlier_call(target, others[-1])       # OCI-AML4 against 7 lines, for a like-for-like count
cal <- rbind(cal, data.frame(line = "OCIAML4 (vs 7)", partner = "MLLT1", t(cnt(o7))))
write.table(cal, file.path(out_dir, "OCIAML4_CRISPR_outlier_calibration.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(cal, row.names = FALSE)

key <- c("KMT2A", "MEN1", "DOT1L", "MLLT1", "MLLT3", "AFF4", "HOXA9", "MEIS1", "ZEB2", "IRF8", "SPI1", "MYB",
         "BCL2L1", "BCL2", "MCL1", "MTOR", "RPTOR", "NPRL2")

# ------------------------------------------------------ genome-wide scatter ---
draw_scatter <- function() {
  layout(matrix(1:2), heights = c(1, 0.07))
  par(mar = c(4.4, 4.6, 3.8, 1.2), fg = ink2, bg = bg)
  ok <- !is.na(oc$value) & !is.na(oc$mean_ref); d <- oc[ok, ]
  lim <- range(c(d$value, d$mean_ref)); lim <- lim + c(-0.05, 0.05) * diff(lim)
  plot(NA, xlim = lim, ylim = lim, axes = FALSE, xlab = "", ylab = "", asp = 1)
  abline(h = pretty(lim), v = pretty(lim), col = grid_col, lwd = 0.8)
  abline(0, 1, col = ink2, lwd = 0.8); abline(h = 0, v = 0, col = ink2, lwd = 0.5, lty = 3)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, mgp = c(2, 0.45, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
  col <- ifelse(d$outlier, ifelse(d$diff < 0, colA, colB), col_ns)
  o <- order(d$outlier)
  points(d$mean_ref[o], d$value[o], pch = 21, col = bg, lwd = 0.3, bg = col[o], cex = ifelse(d$outlier[o], 0.8, 0.45))
  # labels: top outliers each way by |z| that are dependencies somewhere, plus key genes
  dep <- pmin(d$value, d$mean_ref) < -0.5
  lab <- unique(c(head(order(d$z * (d$outlier & d$diff < 0 & dep)), 6),
                  head(order(-d$z * (d$outlier & d$diff > 0 & dep)), 6),
                  which(d$gene %in% c("KMT2A", "BCL2L1", "HOXA9", "NPRL2", "RPTOR", "MTOR", "IRF8"))))
  lab <- lab[(d$outlier[lab] | d$gene[lab] %in% key)]
  text(d$mean_ref[lab], d$value[lab], d$gene[lab], pos = ifelse(d$diff[lab] < 0, 2, 4), offset = 0.3,
       cex = 0.62, col = ink, xpd = NA, font = ifelse(d$gene[lab] %in% key, 2, 1))
  mtext("CRISPR gene effect: OCI-AML4 vs other KMT2A-rearranged lines", side = 3, line = 2.2, adj = 0, font = 2, cex = 1, col = ink)
  mtext(sprintf("Each dot = gene · x = mean of %d other KMT2A-r lines · outlier = outside their range, |diff| ≥ 0.5, |z| ≥ 3 · bold = KMT2A / lineage / survival genes",
                length(others)), side = 3, line = 1, adj = 0, cex = 0.66, col = ink2)
  mtext(sprintf("%d / %d outliers (more / less dependent); a typical other KMT2A-r line tested the same way: median %g / %g",
                sum(d$outlier & d$diff < 0), sum(d$outlier & d$diff > 0), median(cal$more_dep[cal$line != "OCIAML4 (vs 7)"]),
                median(cal$less_dep[cal$line != "OCIAML4 (vs 7)"])), side = 3, line = 0, adj = 0, cex = 0.66, col = ink2)
  mtext("Mean gene effect, other KMT2A-r lines", side = 1, line = 2.4, cex = 0.8, col = ink2)
  mtext("OCI-AML4 gene effect", side = 2, line = 2.8, cex = 0.8, col = ink2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = 21, pt.cex = 1.2, cex = 0.78, col = bg, text.col = ink,
         pt.bg = c(colA, colB, col_ns), legend = c("More dependent in OCI-AML4", "Less dependent in OCI-AML4", "Not an outlier"))
}
fb <- file.path(fig_dir, "OCIAML4_CRISPR_vs_KMT2Ar_scatter")
png(paste0(fb, ".png"), width = 2000, height = 2100, res = 220, type = "cairo", bg = bg); draw_scatter(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2000 / 220, height = 2100 / 220, bg = bg); draw_scatter(); dev.off()

# ------------------------------------------------------ key genes strip -------
draw_key <- function() {
  par(mar = c(5.2, 4.6, 4, 1), fg = ink2, bg = bg)
  K <- C[c(others, target), key]
  yr <- range(K, na.rm = TRUE); yr <- yr + c(-0.05, 0.08) * diff(yr)
  plot(NA, xlim = c(0.5, length(key) + 0.5), ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), col = grid_col, lwd = 0.8); abline(h = 0, col = ink2, lwd = 0.6); abline(h = -1, col = ink2, lwd = 0.6, lty = 3)
  abline(v = c(8.5, 12.5, 15.5), col = grid_col, lwd = 1.5)
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.85, tcl = -0.25)
  text(seq_along(key), yr[1] - 0.03 * diff(yr), key, srt = 45, adj = c(1, 1), xpd = NA, cex = 0.82,
       col = ifelse(oc$outlier[match(key, oc$gene)], ink, ink2), font = ifelse(oc$outlier[match(key, oc$gene)], 2, 1))
  set.seed(2)
  for (j in seq_along(key)) {
    v <- K[others, j]
    points(j + runif(length(v), -0.18, 0.18), v, pch = 21, cex = 1.05, lwd = 0.5, col = bg, bg = colB)
    segments(j - 0.25, median(v, na.rm = TRUE), j + 0.25, lwd = 1.6, col = ink)
    points(j, K[target, j], pch = 23, cex = 1.7, lwd = 0.8, col = ink, bg = colA)
  }
  for (b in list(c(1, 8, "KMT2A complex / HOX"), c(9, 12, "Myeloid lineage TFs"), c(13, 15, "BCL2 family"), c(16, 18, "mTOR pathway")))
    text(mean(as.numeric(b[1:2])), yr[2], b[3], cex = 0.75, col = ink2, font = 3)
  mtext("OCI-AML4 vs other KMT2A-rearranged lines: selected dependencies (CRISPR)", side = 3, line = 2.2, adj = 0, font = 2, cex = 1, col = ink)
  mtext(sprintf("Diamond = OCI-AML4 (KMT2A-MLLT1) · dots = %s · bar = median of others · bold gene = outlier",
                paste(sprintf("%s (%s)", others, partner[others]), collapse = ", ")), side = 3, line = 0.9, adj = 0, cex = 0.6, col = ink2)
  mtext("Chronos gene effect", side = 2, line = 2.8, cex = 0.8, col = ink2)
}
fb <- file.path(fig_dir, "OCIAML4_CRISPR_key_dependencies")
png(paste0(fb, ".png"), width = 2900, height = 1500, res = 220, type = "cairo", bg = bg); draw_key(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2900 / 220, height = 1500 / 220, bg = bg); draw_key(); dev.off()

# ------------------------------------------------ volcano (Crawford-Howell) ---
draw_volcano <- function() {
  layout(matrix(1:2), heights = c(1, 0.07))
  par(mar = c(4.4, 4.6, 4.2, 1.2), fg = ink2, bg = bg)
  d <- oc[!is.na(oc$p_CH), ]; y <- -log10(d$p_CH)
  hit <- d$p_CH < 0.001 & abs(d$diff) >= 0.5
  xr <- max(abs(d$diff)) * 1.08 * c(-1, 1); yr <- c(0, max(y) * 1.06)
  plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), v = pretty(xr), col = grid_col, lwd = 0.8); abline(v = 0, col = ink2, lwd = 0.8)
  abline(h = -log10(0.001), v = c(-0.5, 0.5), lty = 2, col = ink2, lwd = 0.6)
  fdr_line <- if (any(d$fdr_CH < 0.1)) max(d$p_CH[d$fdr_CH < 0.1]) else NA
  if (!is.na(fdr_line)) abline(h = -log10(fdr_line), lty = 3, col = ink2, lwd = 0.8)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, mgp = c(2, 0.45, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
  col <- ifelse(hit, ifelse(d$diff < 0, colA, colB), col_ns)
  o <- order(hit)
  points(d$diff[o], y[o], pch = 21, col = bg, lwd = 0.3, bg = col[o], cex = ifelse(hit[o], 0.85, 0.45))
  # labels: every coloured hit, plus key genes (bold) with p < 0.05
  lab <- unique(c(which(hit), which(d$gene %in% key & d$p_CH < 0.05)))
  side <- d$diff[lab] < 0
  for (sd_ in c(TRUE, FALSE)) {
    i <- lab[side == sd_]; if (!length(i)) next
    ly <- y[i]; oo <- order(ly); ys <- ly[oo]; gap <- diff(yr) * 0.032
    for (k in seq_along(ys)[-1]) if (ys[k] - ys[k - 1] < gap) ys[k] <- ys[k - 1] + gap
    ly[oo] <- ys
    lx <- if (sd_) xr[1] + diff(xr) * 0.13 else xr[2] - diff(xr) * 0.13
    segments(d$diff[i], y[i], lx, ly, col = ink2, lwd = 0.4)
    text(lx, ly, d$gene[i], pos = if (sd_) 2 else 4, offset = 0.15, cex = 0.62, col = ink, xpd = NA,
         font = ifelse(d$gene[i] %in% key, 2, 1))
  }
  cal_o <- cal[cal$line != "OCIAML4 (vs 7)", ]
  mtext("CRISPR: OCI-AML4 vs other KMT2A-rearranged lines", side = 3, line = 2.6, adj = 0, font = 2, cex = 1, col = ink)
  mtext(sprintf("Crawford–Howell single-case t-test vs %d lines (df = %d) · coloured = p < 0.001 and |diff| ≥ 0.5 · bold = key genes with p < 0.05 · FDR < 0.1%s",
                length(others), length(others) - 1, if (is.na(fdr_line)) ": none" else " = dotted line"),
        side = 3, line = 1.45, adj = 0, cex = 0.66, col = ink2)
  mtext(sprintf("p < 0.001: OCI-AML4 %d genes · other KMT2A-r lines tested the same way (vs remaining 7): median %g, range %d–%d",
                sum(d$p_CH < 0.001), median(cal_o$p001), min(cal_o$p001), max(cal_o$p001)),
        side = 3, line = 0.4, adj = 0, cex = 0.66, col = ink2)
  mtext("Gene effect difference (OCI-AML4 − mean of other KMT2A-r)", side = 1, line = 2.4, cex = 0.8, col = ink2)
  mtext("−log10 p", side = 2, line = 2.8, cex = 0.8, col = ink2)
  text(xr[1], yr[2], "more dependent in OCI-AML4", adj = c(0, 1), cex = 0.7, font = 2, col = colA)
  text(xr[2], yr[2], "less dependent in OCI-AML4", adj = c(1, 1), cex = 0.7, font = 2, col = colB)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = 21, pt.cex = 1.2, cex = 0.78, col = bg, text.col = ink,
         pt.bg = c(colA, colB, col_ns), legend = c("More dependent in OCI-AML4", "Less dependent in OCI-AML4", "Other"))
}
fb <- file.path(fig_dir, "OCIAML4_CRISPR_vs_KMT2Ar_volcano")
png(paste0(fb, ".png"), width = 2000, height = 1900, res = 220, type = "cairo", bg = bg); draw_volcano(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2000 / 220, height = 1900 / 220, bg = bg); draw_volcano(); dev.off()
