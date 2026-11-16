# Bar plot of KMT2A dependency (DepMap 26Q3 CRISPR Chronos gene effect) in the KMT2A-rearranged
# AML lines, with the non-KMT2A AML lines (same 56-line panel) as a reference band.
# More negative = stronger dependency; -1 is the median effect of common-essential genes.
#
# Run from the project root: Rscript scripts/plot_KMT2A_dependency.R  (set the gene below)

gene <- "KMT2A"                     # dependency to plot
res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
C  <- attr(readRDS(file.path(res_dir, "all_results.rds"))$crispr, "X")
colnames(C) <- sub(" \\([0-9]+\\)$", "", colnames(C))
v  <- C[, gene]

kr <- cl$cell_line[cl$KMT2A_rearranged]
d  <- data.frame(cell_line = intersect(kr, names(v)))
d$effect  <- v[d$cell_line]
d$partner <- cl$KMT2A_partner[match(d$cell_line, cl$cell_line)]
d$group   <- cl$mTOR_group[match(d$cell_line, cl$cell_line)]
d <- d[order(d$effect), ]
ref <- v[setdiff(names(v), kr)]                     # non-KMT2A AML lines with CRISPR data
missing <- setdiff(kr, names(v))
write.table(d, file.path(res_dir, paste0("CRISPR_", gene, "_KMT2Ar_lines.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

draw <- function() {
  bg <- "#fcfcfb"; col_t <- "#eb6834"; col_o <- "#4a3aa7"
  n <- nrow(d); y <- rev(seq_len(n))
  xl <- c(min(-1.2, d$effect, quantile(ref, 0.05)) - 0.1, max(0.3, d$effect, quantile(ref, 0.95)) + 0.1)
  par(oma = c(0, 1, 3.8, 1), mar = c(4.4, 12, 0.5, 1.5), fg = ink2)
  plot(NA, xlim = xl, ylim = c(0.4, n + 0.6), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  q <- quantile(ref, c(0.25, 0.5, 0.75))
  rect(q[1], 0.4, q[3], n + 0.6, col = "#ecebe6", border = NA)                         # non-KMT2A IQR
  abline(v = pretty(xl), col = grid_col, lwd = 0.8)
  abline(v = q[2], col = ink2, lwd = 1.2)                                               # non-KMT2A median
  abline(v = c(-0.5, -1), col = ink2, lty = c(3, 2), lwd = 0.8)
  abline(v = 0, col = ink, lwd = 0.8)
  is_t <- d$cell_line == "OCIAML4"
  rect(0, y - 0.32, d$effect, y + 0.32, col = ifelse(is_t, col_t, col_o), border = bg)
  text(d$effect, y, sprintf("%.2f", d$effect), pos = ifelse(d$effect < 0, 2, 4), cex = 0.8, col = ink, offset = 0.3)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.9, tcl = -0.25, mgp = c(2, 0.5, 0))
  # per-label colour: axis() takes one colour, so draw the labels with text()
  text(par("usr")[1], y, sprintf("%s (%s)", d$cell_line, d$partner), pos = 2, xpd = NA, offset = 0.6,
       col = ifelse(is_t, col_t, ink), cex = 0.9, font = ifelse(is_t, 2, 1))
  mtext(sprintf("CRISPR gene effect (Chronos) · ← more dependent"), side = 1, line = 2.6, cex = 0.85, col = ink2)
  mtext(sprintf("%s dependency in KMT2A-rearranged AML lines (n=%d)", gene, n), side = 3, line = 2.2, adj = 0,
        outer = TRUE, font = 2, cex = 1.2, col = ink)
  mtext(sprintf("DepMap 26Q3 CRISPR · shaded band = IQR and solid line = median of non-KMT2A AML lines (n=%d, median %.2f) · dashed = -1 (common-essential median) · dotted = -0.5%s",
                length(ref), q[2], if (length(missing)) sprintf(" · no CRISPR data: %s", paste(missing, collapse = ", ")) else ""),
        side = 3, line = 0.8, adj = 0, outer = TRUE, cex = 0.72, col = ink2)
}
fb <- file.path(fig_dir, paste0("CRISPR_", gene, "_KMT2Ar_barplot"))
png(paste0(fb, ".png"), width = 2600, height = 1500, res = 220, type = "cairo", bg = "#fcfcfb"); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2600 / 220, height = 1500 / 220, bg = "#fcfcfb"); draw(); dev.off()

print(d, row.names = FALSE)
cat(sprintf("\nNon-KMT2A AML lines (n=%d): median %.2f, IQR %.2f to %.2f | KMT2A-r median %.2f | Wilcoxon p = %.2g\n",
            length(ref), median(ref), quantile(ref, 0.25), quantile(ref, 0.75), median(d$effect),
            suppressWarnings(wilcox.test(d$effect, ref, exact = FALSE)$p.value)))
cat("KMT2A-r without CRISPR data:", if (length(missing)) missing else "none", "\n")
