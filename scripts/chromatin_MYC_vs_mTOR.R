# Global chromatin profiling (CCLE 2019, 42 histone-modification combinations by mass spec)
# in MYC-High/mTOR-Low vs mTOR-High/MYC-Low AML lines (RNA-defined groups; lines high or low
# for both scores are left out). Values in the file are already centred across all CCLE lines.
#
# Run from the project root: Rscript scripts/chromatin_MYC_vs_mTOR.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
X  <- attr(readRDS(file.path(res_dir, "all_results.rds"))$chrom, "X")
g  <- cl[match(rownames(X), cl$cell_line), ]
grp <- ifelse(g$MYC_group == "High" & g$mTOR_group == "Low", "MYC-High",
       ifelse(g$mTOR_group == "High" & g$MYC_group == "Low", "mTOR-High", NA))
X <- X[!is.na(grp), ]; g <- g[!is.na(grp), ]; grp <- grp[!is.na(grp)]
a <- grp == "MYC-High"

# --------------------------------------------------------------- stats ------
st <- do.call(rbind, lapply(colnames(X), function(m) {
  x <- X[a, m]; y <- X[!a, m]; x <- x[!is.na(x)]; y <- y[!is.na(y)]
  data.frame(mark = m, n_MYC = length(x), n_mTOR = length(y), mean_MYC = mean(x), mean_mTOR = mean(y),
             diff_MYC_minus_mTOR = mean(x) - mean(y),
             p_welch = if (length(x) >= 3 && length(y) >= 3) t.test(x, y)$p.value else NA,
             p_wilcox = if (length(x) >= 3 && length(y) >= 3) suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value) else NA)
}))
st$fdr_welch <- p.adjust(st$p_welch, "BH"); st$fdr_wilcox <- p.adjust(st$p_wilcox, "BH")
st <- st[order(st$p_welch), ]
write.table(st, file.path(res_dir, "chromatin_MYC_High_vs_mTOR_High.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# residue families, used to group the heatmap rows
family <- function(m) ifelse(grepl("^H3K4", m), "H3K4", ifelse(grepl("^H3K9", m), "H3K9 / K14",
          ifelse(grepl("^H3K18", m), "H3K18 / K23", ifelse(grepl("^H3\\.3", m), "H3.3K27 / K36",
          ifelse(grepl("^H3K27", m), "H3K27 / K36", ifelse(grepl("^H3K56", m), "H3K56", "H3K79"))))))

# -------------------------------------------------------------- heatmap -----
bg <- "#fcfcfb"; colA <- "#eb6834"; colB <- "#4a3aa7"; col_k <- "#1b1b1a"; col_other <- "#e8e7e3"
pal <- colorRampPalette(c("#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a8a6", "#e34948", "#9c2423"))(101)
zmax <- 2
Z <- sweep(X, 2, apply(X, 2, median, na.rm = TRUE))           # each mark minus its median across these lines
Z[Z > zmax] <- zmax; Z[Z < -zmax] <- -zmax
fam_order <- c("H3K4", "H3K9 / K14", "H3K18 / K23", "H3K27 / K36", "H3.3K27 / K36", "H3K56", "H3K79")
marks <- colnames(X)[order(match(family(colnames(X)), fam_order), seq_along(colnames(X)))]
sc  <- rowMeans(Z[, st$mark[st$p_welch < 0.05], drop = FALSE] * sign(st$diff_MYC_minus_mTOR[st$p_welch < 0.05]), na.rm = TRUE)
cols <- c(which(a)[order(-sc[a])], which(!a)[order(-sc[!a])])  # lines ordered by the nominally different marks
Zp <- Z[cols, rev(marks)]; gp <- g[cols, ]
pv <- st$p_welch[match(rev(marks), st$mark)]; fd <- st$fdr_welch[match(rev(marks), st$mark)]
dd <- st$diff_MYC_minus_mTOR[match(rev(marks), st$mark)]
star <- ifelse(fd < 0.05, "**", ifelse(pv < 0.05, "*", ""))

draw_hm <- function() {
  nr <- ncol(Zp); nc <- nrow(Zp)
  par(oma = c(0, 1, 5.4, 1))
  layout(matrix(1:4, 2, byrow = TRUE), widths = c(1, 0.14), heights = c(0.07, 1))
  par(mar = c(0.4, 13, 0.6, 1), fg = ink2)
  plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, 2.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rect(seq_len(nc) - 0.5, 1.55, seq_len(nc) + 0.5, 2.45, border = bg, lwd = 0.6, col = ifelse(a[cols], colA, colB))
  rect(seq_len(nc) - 0.5, 0.55, seq_len(nc) + 0.5, 1.45, border = bg, lwd = 0.6, col = ifelse(gp$KMT2A_rearranged, col_k, col_other))
  axis(2, at = c(2, 1), labels = c("Group", "KMT2A-r"), las = 1, tick = FALSE, col.axis = ink, cex.axis = 1)
  mtext(sprintf("Global chromatin profiling: MYC-High (n=%d) vs mTOR-High (n=%d)", sum(a), sum(!a)),
        side = 3, line = 3.4, adj = 0, outer = TRUE, font = 2, cex = 1.3, col = ink)
  mtext("CCLE 2019 histone-mark mass spec · each mark minus its median across these lines (capped ±2) · number = mean MYC-High − mTOR-High · * Welch p<0.05, ** FDR<0.05",
        side = 3, line = 1.9, adj = 0, outer = TRUE, cex = 0.82, col = ink2)
  mtext(sprintf("Lines ordered within group by the nominally different marks · not profiled: %s",
                paste(setdiff(cl$cell_line[(cl$MYC_group == "High") != (cl$mTOR_group == "High")], rownames(X)), collapse = ", ")),
        side = 3, line = 0.7, adj = 0, outer = TRUE, cex = 0.72, col = ink2, font = 3)
  par(mar = c(0.4, 0.5, 0.6, 0.5)); plot.new()
  legend("left", bty = "n", cex = 0.85, text.col = ink, pch = 15, pt.cex = 1.4, col = c(colA, colB, col_k),
         legend = c("MYC-High", "mTOR-High", "KMT2A-r"))
  par(mar = c(7.5, 13, 0.4, 1))
  image(seq_len(nc), seq_len(nr), Zp, useRaster = TRUE, col = pal, zlim = c(-zmax, zmax), axes = FALSE, xlab = "", ylab = "")
  na <- which(is.na(Zp), arr.ind = TRUE)                                      # missing measurements
  if (nrow(na)) rect(na[, 1] - 0.5, na[, 2] - 0.5, na[, 1] + 0.5, na[, 2] + 0.5, col = "#ffffff", border = "#d6d5d0", lwd = 0.4)
  abline(v = sum(a) + 0.5, col = bg, lwd = 3)
  fm <- family(rev(marks)); br <- which(fm[-1] != fm[-nr]); abline(h = br + 0.5, col = bg, lwd = 3)
  axis(2, at = seq_len(nr), labels = sprintf("%s  %+.2f %s", rev(marks), dd, star), las = 1, tick = FALSE,
       cex.axis = 0.78, col.axis = ink, line = -0.6, font = 1)
  text(seq_len(nc), par("usr")[3] - 0.4, rownames(Zp), srt = 90, adj = 1, xpd = NA, cex = 0.8, col = ifelse(a[cols], colA, colB))
  par(mar = c(7.5, 1.2, 8, 3.8))
  zs <- seq(-zmax, zmax, length.out = 101)
  image(1, zs, t(zs), useRaster = TRUE, col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = -2:2, las = 1, cex.axis = 0.9, col = NA, col.ticks = ink2, col.axis = ink2, tcl = -0.2)
  mtext("vs mark\nmedian", side = 3, line = 0.5, cex = 0.9, col = ink2)
}
fb <- file.path(fig_dir, "chromatin_MYC_High_vs_mTOR_High_heatmap")
png(paste0(fb, ".png"), width = 3000, height = 2800, res = 220, type = "cairo", bg = bg); draw_hm(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3000 / 220, height = 2800 / 220, bg = bg); draw_hm(); dev.off()

# ------------------------------------------ per-line strips, top marks ------
top <- head(st$mark, 6)
draw_strip <- function() {
  par(mfrow = c(2, 3), oma = c(0, 1, 3, 1), mar = c(3.2, 4, 3, 1), fg = ink2)
  set.seed(1)
  for (m in top) {
    v <- X[, m]; ok <- !is.na(v); f <- factor(ifelse(a, "MYC-High", "mTOR-High")[ok], levels = c("MYC-High", "mTOR-High"))
    v <- v[ok]; r <- st[st$mark == m, ]
    plot(NA, xlim = c(0.4, 2.6), ylim = range(v) + c(-0.05, 0.1) * diff(range(v)), axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(v), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.85, tcl = -0.25)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", levels(f), as.vector(table(f))), col = NA, col.axis = ink, cex.axis = 0.85, padj = 0.5)
    points(as.integer(f) + runif(length(v), -0.18, 0.18), v, pch = 21, cex = 1.3, lwd = 0.5, col = bg,
           bg = ifelse(f == "MYC-High", colA, colB))
    for (j in 1:2) segments(j - 0.3, median(v[as.integer(f) == j]), j + 0.3, lwd = 2, col = ink)
    mtext(m, side = 3, line = 1.3, adj = 0, font = 2, cex = 0.85, col = ink)
    mtext(sprintf("Δ = %+.2f · Welch p = %.1e · FDR = %.2g", r$diff_MYC_minus_mTOR, r$p_welch, r$fdr_welch),
          side = 3, line = 0.3, adj = 0, cex = 0.65, col = ink2)
    mtext("relative abundance (CCLE-centred)", side = 2, line = 2.6, cex = 0.6, col = ink2)
  }
  mtext("Top histone marks by p-value: MYC-High vs mTOR-High", side = 3, line = 1, adj = 0, outer = TRUE, font = 2, cex = 1.1, col = ink)
}
fb <- file.path(fig_dir, "chromatin_MYC_High_vs_mTOR_High_top_marks")
png(paste0(fb, ".png"), width = 3000, height = 1900, res = 220, type = "cairo", bg = bg); draw_strip(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3000 / 220, height = 1900 / 220, bg = bg); draw_strip(); dev.off()

cat(sprintf("MYC-High n=%d | mTOR-High n=%d | marks %d | Welch p<0.05: %d | FDR<0.05: %d | FDR<0.25: %d\n",
            sum(a), sum(!a), nrow(st), sum(st$p_welch < 0.05, na.rm = TRUE), sum(st$fdr_welch < 0.05, na.rm = TRUE),
            sum(st$fdr_welch < 0.25, na.rm = TRUE)))
print(transform(head(st, 12), mean_MYC = round(mean_MYC, 2), mean_mTOR = round(mean_mTOR, 2),
                diff_MYC_minus_mTOR = round(diff_MYC_minus_mTOR, 2), p_welch = signif(p_welch, 2), p_wilcox = signif(p_wilcox, 2),
                fdr_welch = signif(fdr_welch, 2), fdr_wilcox = signif(fdr_wilcox, 2)), row.names = FALSE)

# ------------------------------------------------------------- volcano ------
vr <- data.frame(feature = st$mark, diff_High_minus_Low = st$diff_MYC_minus_mTOR, p_welch = st$p_welch,
                 fdr_welch = st$fdr_welch, n_High = st$n_MYC, n_Low = st$n_mTOR)
volcano_png(vr, list(title = "Global chromatin profiling: MYC-High vs mTOR-High", x = "Relative abundance difference (MYC-High − mTOR-High)",
                     dep = FALSE, groups = c("MYC-High", "mTOR-High"), cols = c(colA, colB)),
            file.path(fig_dir, "volcano_chromatin_MYC_High_vs_mTOR_High.png"), n_lab = 10, xpad = 1.9)

# ------------------------------------------------- box plot of one mark -----
boxplot_mark <- function(m) {
  v <- X[, m]; ok <- !is.na(v)
  f <- factor(ifelse(a, "MYC-High", "mTOR-High")[ok], levels = c("MYC-High", "mTOR-High")); v <- v[ok]
  r <- st[st$mark == m, ]
  write.table(data.frame(cell_line = rownames(X)[ok], group = f, value = v)[order(f, -v), ],
              file.path(res_dir, paste0("chromatin_", m, "_per_line.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  draw <- function() {
    par(oma = c(0, 1, 4.6, 1), mar = c(3.6, 4.6, 0.5, 1), fg = ink2)
    yl <- range(v) + c(-0.08, 0.08) * diff(range(v))
    plot(NA, xlim = c(0.4, 2.6), ylim = yl, axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(yl), col = grid_col, lwd = 0.8)
    boxplot(v ~ f, add = TRUE, at = 1:2, boxwex = 0.5, axes = FALSE, outline = FALSE, col = c(tint(colA, 0.7), tint(colB, 0.7)),
            border = ink2, lwd = 1, medlwd = 2.2, whisklty = 1, staplewex = 0.4)
    set.seed(1)
    points(as.integer(f) + runif(length(v), -0.15, 0.15), v, pch = 21, cex = 1.35, lwd = 0.6, col = bg,
           bg = ifelse(f == "MYC-High", colA, colB))
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.9, tcl = -0.25)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", levels(f), as.vector(table(f))), col = NA, col.axis = ink, cex.axis = 0.95, padj = 0.6)
    mtext(sprintf("%s relative abundance (CCLE-centred)", m), side = 2, line = 3.2, cex = 0.85, col = ink2)
    mtext(sprintf("%s: MYC-High vs mTOR-High", m), side = 3, line = 3, adj = 0, outer = TRUE, font = 2, cex = 1.2, col = ink)
    mtext(sprintf("CCLE 2019 global chromatin profiling · Δ mean = %+.2f", r$diff_MYC_minus_mTOR),
          side = 3, line = 1.6, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
    mtext(sprintf("Welch p = %.1e (FDR %.2g over 42 marks) · Wilcoxon p = %.1e", r$p_welch, r$fdr_welch, r$p_wilcox),
          side = 3, line = 0.4, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
  }
  fb <- file.path(fig_dir, paste0("chromatin_", m, "_boxplot"))
  png(paste0(fb, ".png"), width = 1500, height = 1500, res = 220, type = "cairo", bg = bg); draw(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 1500 / 220, height = 1500 / 220, bg = bg); draw(); dev.off()
}
boxplot_mark("H3K4me1")
