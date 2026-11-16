# MYC vs mTOR RNA signature scores for the 56 AML lines (DepMap 26Q1 RNA classification).
# High/Low = above/below each score's median (dashed lines). Triangle = KMT2A-rearranged.
#
# Run from the project root: Rscript scripts/plot_MYC_mTOR_scores.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp <- with(cl, ifelse(MYC_group == "High" & mTOR_group == "Low", "MYC-High",
                   ifelse(mTOR_group == "High" & MYC_group == "Low", "mTOR-High",
                   ifelse(MYC_group == "High", "Both high", "Both low"))))

colA <- "#eb6834"; colB <- "#4a3aa7"; col_mix <- "#a3a29d"
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
cols <- c("MYC-High" = colA, "mTOR-High" = colB, "Both high" = col_mix, "Both low" = col_mix)
# x = mTOR score, y = MYC score
xv <- cl$mTOR_score; yv <- cl$MYC_score
mx <- median(xv); my <- median(yv)

draw <- function() {
  layout(matrix(1:2), heights = c(1, 0.07))
  par(mar = c(4.4, 4.6, 3.6, 1.2), fg = ink2, bg = bg)
  xr <- range(xv) + c(-0.1, 0.14); yr <- range(yv) + c(-0.1, 0.1)
  plot(NA, xlim = xr, ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), v = pretty(xr), col = grid_col, lwd = 0.8)
  abline(v = mx, h = my, lty = 2, col = ink2, lwd = 0.9)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.9, tcl = -0.25, mgp = c(2, 0.45, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.9, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
  points(xv, yv, pch = ifelse(cl$KMT2A_rearranged %in% TRUE, 24, 21),
         cex = 1.5, lwd = 0.6, col = bg, bg = cols[cl$grp])
  # labels: greedy placement - try spots around each point, keep the first that hits no
  # earlier label, no point and no quadrant title; draw a leader line if it moved far
  lab <- cl$cell_line
  lcex <- 0.62; px <- xv; py <- yv
  w <- strwidth(lab, cex = lcex); h <- strheight("M", cex = lcex) * 1.35
  rx <- diff(xr) * 0.012; ry <- diff(yr) * 0.012                 # point radius, roughly
  q_box <- rbind(c(xr[1], yr[2] - h * 1.8, xr[1] + diff(xr) * 0.2, yr[2]),
                 c(xr[2] - diff(xr) * 0.2, yr[2] - h * 1.8, xr[2], yr[2]),
                 c(xr[1], yr[1], xr[1] + diff(xr) * 0.2, yr[1] + h * 1.8),
                 c(xr[2] - diff(xr) * 0.2, yr[1], xr[2], yr[1] + h * 1.8))
  boxes <- q_box
  hits <- function(b) any(b[1] < boxes[, 3] & b[3] > boxes[, 1] & b[2] < boxes[, 4] & b[4] > boxes[, 2]) ||
    any(px > b[1] - rx & px < b[3] + rx & py > b[2] - ry & py < b[4] + ry)
  lx <- ly <- numeric(nrow(cl)); adjx <- numeric(nrow(cl))
  for (k in order(-abs(py - my) - abs(px - mx))) {
    d <- c(1, 1.8, 2.8, 4, 5.5)
    cand <- do.call(rbind, lapply(d, function(s) rbind(
      c(px[k] + rx * s, py[k], 0), c(px[k] - rx * s, py[k], 1),
      c(px[k] + rx * s, py[k] + h * s * 0.8, 0), c(px[k] + rx * s, py[k] - h * s * 0.8, 0),
      c(px[k] - rx * s, py[k] + h * s * 0.8, 1), c(px[k] - rx * s, py[k] - h * s * 0.8, 1))))
    done <- FALSE
    for (j in seq_len(nrow(cand))) {
      x0 <- if (cand[j, 3] == 0) cand[j, 1] else cand[j, 1] - w[k]
      b <- c(x0, cand[j, 2] - h / 2, x0 + w[k], cand[j, 2] + h / 2)
      if (b[1] < xr[1] || b[3] > xr[2] || b[2] < yr[1] || b[4] > yr[2] || hits(b)) next
      lx[k] <- cand[j, 1]; ly[k] <- cand[j, 2]; adjx[k] <- cand[j, 3]; boxes <- rbind(boxes, b); done <- TRUE; break
    }
    if (!done) { lx[k] <- px[k] + rx; ly[k] <- py[k]; adjx[k] <- 0 }
  }
  far <- sqrt(((lx - px) / rx)^2 + ((ly - py) / ry)^2) > 2.2
  segments(px[far], py[far], lx[far], ly[far], col = ink2, lwd = 0.4)
  for (k in seq_len(nrow(cl))) text(lx[k], ly[k], lab[k], adj = c(adjx[k], 0.5), cex = lcex,
       col = if (cl$grp[k] == "MYC-High") colA else if (cl$grp[k] == "mTOR-High") colB else ink2)
  q <- function(x, y, s, a) text(x, y, s, adj = a, cex = 0.8, font = 2, col = ink2)
  q(xr[2], yr[2], sprintf("Both high (n=%d)", sum(cl$grp == "Both high")), c(1, 1))
  q(xr[1], yr[1], sprintf("Both low (n=%d)", sum(cl$grp == "Both low")), c(0, 0))
  q(xr[2], yr[1], sprintf("mTOR-High (n=%d)", sum(cl$grp == "mTOR-High")), c(1, 0))
  q(xr[1], yr[2], sprintf("MYC-High (n=%d)", sum(cl$grp == "MYC-High")), c(0, 1))
  mtext("MYC vs mTOR RNA signature scores, 56 AML cell lines", side = 3, line = 1.9, adj = 0, font = 2, cex = 1.05, col = ink)
  mtext(sprintf("DepMap 26Q1 RNA · dashed = median of each score (mTOR %.3f, MYC %.3f) · triangle = KMT2A-rearranged · Spearman rho = %.2f",
                mx, my, cor(cl$MYC_score, cl$mTOR_score, method = "spearman")),
        side = 3, line = 0.65, adj = 0, cex = 0.68, col = ink2)
  mtext("mTOR signature score", side = 1, line = 2.5, cex = 0.85, col = ink2)
  mtext("MYC signature score", side = 2, line = 2.9, cex = 0.85, col = ink2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = c(21, 21, 21, 21, 24), pt.cex = 1.3, cex = 0.78, col = c(bg, bg, bg, ink2, ink2),
         text.col = ink, pt.bg = c(colA, colB, col_mix, bg, bg), pt.lwd = c(0.6, 0.6, 0.6, 1, 1),
         legend = c("MYC-High", "mTOR-High", "Both high / low (excluded)", "no KMT2A-r", "KMT2A-rearranged"))
}
fb <- file.path(fig_dir, "MYC_vs_mTOR_scores_scatter")
png(paste0(fb, ".png"), width = 2000, height = 1900, res = 220, type = "cairo", bg = bg); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2000 / 220, height = 1900 / 220, bg = bg); draw(); dev.off()
