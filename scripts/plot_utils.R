# Shared plotting helpers (palette, label placement, volcano) for the
# mTOR-High vs mTOR-Low analyses. Base R only. source() from project root.

# palette: diverging blue (High) / red (Low) poles, validated; gray = n.s.
col_hi  <- "#2a78d6"; col_lo <- "#e34948"
tint_hi <- "#9fc2ee"; tint_lo <- "#f2a8a6"
col_ns  <- "#c9c8c3"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"

rapalog <- "SIROLIMUS|RAPAMYCIN|EVEROLIMUS|TEMSIROLIMUS|RIDAFOROLIMUS"
clean <- function(f) sub("\\|.*$", "", sub(" \\([0-9]+\\)$", "", f))
label_of <- function(r, i) {             # compound name for drugs, symbol otherwise
  if (!is.null(r$CompoundName)) {
    nm <- sub(" \\(.*\\)$", "", r$CompoundName[i])
    ifelse(is.na(nm) | nm == "", r$feature[i], substr(nm, 1, 28))
  } else clean(r$feature[i])
}

# greedy vertical de-overlap for labels on one side
spread <- function(y, gap, top) {
  o <- order(y); ys <- y[o]
  n <- length(ys)
  ys[1] <- max(ys[1], gap / 2)                           # upward pass from just above the axis
  for (i in seq_len(n)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  ys[n] <- min(ys[n], top)                               # downward pass from the top
  for (i in rev(seq_len(n))[-1]) ys[i] <- min(ys[i], ys[i + 1] - gap)
  y[o] <- ys; y
}

# r: one row per feature with diff_High_minus_Low, p_welch, fdr_welch, n_High, n_Low
# p: list(title, x, dep, always = regex of features to always label, n_features = override,
#           groups = c(name of the group on the positive side, name on the negative side) [default High, Low],
#           cols = c(colour for the positive-side group, negative-side group) [default blue, red]
#           fdr = FDR cutoff for the solid colour and dashed line [default 0.05])
tint <- function(col, f = 0.55) colorRampPalette(c(col, "#ffffff"))(100)[round(f * 100)]
vgroups <- function(p) if (is.null(p$groups)) c("High", "Low") else p$groups
vfdr <- function(p) if (is.null(p$fdr)) 0.05 else p$fdr
vcols <- function(p) { if (is.null(p$cols)) return(c(col_hi, tint_hi, col_lo, tint_lo))
                       c(p$cols[1], tint(p$cols[1]), p$cols[2], tint(p$cols[2])) }
volcano <- function(r, p, n_lab = 8, xpad = 1.5) {
  r  <- r[is.finite(r$p_welch) & is.finite(r$diff_High_minus_Low), ]
  x  <- r$diff_High_minus_Low; y <- -log10(r$p_welch)
  right <- x > 0
  # colour = the group the feature favours; for dependency/AUC a lower value
  # means more dependent/sensitive, so the negative side belongs to High
  hi  <- if (p$dep) !right else right
  vc <- vcols(p); gn <- vgroups(p); q <- vfdr(p)
  col <- ifelse(r$fdr_welch < q, ifelse(hi, vc[1], vc[3]),
         ifelse(r$fdr_welch < 0.25, ifelse(hi, vc[2], vc[4]), col_ns))
  xl <- max(abs(x)) * c(-xpad, xpad); yl <- c(0, max(y) * 1.12)
  plot(NA, xlim = xl, ylim = yl, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yl), v = pretty(xl), col = grid_col, lwd = 0.8)
  abline(v = 0, col = ink2, lwd = 0.8)
  if (any(r$fdr_welch < q)) {   # p-value at the FDR boundary
    abline(h = -log10(max(r$p_welch[r$fdr_welch < q])), lty = 2, col = ink2, lwd = 0.8)
  }
  # visible L-shaped axes: lines along the plot edge, ticks and labels in dark ink
  axis(1, at = pretty(xl), col = ink, col.ticks = ink, col.axis = ink, lwd = 0.9, cex.axis = 0.85, tcl = -0.35,
       mgp = c(2, 0.45, 0), pos = yl[1] - 0.02 * diff(yl))
  axis(2, at = pretty(yl), col = ink, col.ticks = ink, col.axis = ink, lwd = 0.9, cex.axis = 0.85, tcl = -0.35,
       las = 1, mgp = c(2, 0.55, 0), pos = xl[1] - 0.02 * diff(xl))
  o <- order(col != col_ns, y)                         # draw significant on top
  points(x[o], y[o], pch = 21, bg = col[o], col = "#fcfcfb", lwd = 0.4,
         cex = ifelse(col[o] == col_ns, 0.55, 0.8))
  # labels: top hits by p, split by side, de-overlapped
  lab <- head(order(r$p_welch), n_lab * 4)
  lab <- head(lab[!duplicated(label_of(r, lab))], n_lab)      # one label per gene/compound
  lab <- lab[r$fdr_welch[lab] < 0.25 | seq_along(lab) <= 4]
  if (!is.null(r$CompoundName)) lab <- unique(c(lab, grep(rapalog, r$CompoundName)))  # always show rapalogs
  if (!is.null(p$always)) {                       # features of interest: best region per gene
    a <- grep(p$always, r$feature); a <- a[!duplicated(clean(r$feature[a]))]   # r is sorted by p
    lab <- unique(c(lab, a))
  }
  for (side in c(TRUE, FALSE)) {
    i <- lab[right[lab] == side]; if (!length(i)) next
    gap <- min(diff(yl) * 0.055, diff(yl) * 0.95 / length(i))
    ly <- spread(y[i], gap, yl[2])
    w  <- max(strwidth(label_of(r, i), cex = 0.62)) * 1.1   # keep the label column inside the plot
    lx <- if (side) rep(xl[2] - w, length(i))
          else      rep(xl[1] + w, length(i))
    segments(x[i], y[i], lx, ly, col = ink2, lwd = 0.5)
    text(lx, ly, label_of(r, i), pos = if (side) 4 else 2, offset = 0.15,
         cex = 0.62, col = ink, xpd = NA)
  }
  n_feat <- if (is.null(p$n_features)) nrow(r) else p$n_features
  n_sig  <- if (is.null(p$n_sig)) sum(r$fdr_welch < q) else p$n_sig
  mtext(p$title, side = 3, line = 1.3, adj = 0, font = 2, cex = 0.85, col = ink)
  mtext(sprintf("%s n=%d · %s n=%d · %d features · %d at FDR<%g", gn[1],
                max(r$n_High), gn[2], max(r$n_Low), n_feat, n_sig, q),
        side = 3, line = 0.35, adj = 0, cex = 0.62, col = ink2)
  mtext(p$x, side = 1, line = 1.9, cex = 0.8, col = ink)
  mtext(expression(-log[10]~italic(p)), side = 2, line = 2.1, cex = 0.8, col = ink)
  hint <- if (p$dep) c(sprintf("← more dependent/sensitive in %s", gn[1]), sprintf("more dependent/sensitive in %s →", gn[2]))
          else       c(sprintf("← higher in %s", gn[2]), sprintf("higher in %s →", gn[1]))
  hcol <- if (p$dep) c(vc[1], vc[3]) else c(vc[3], vc[1])        # colour each side by its group
  mtext(hint[1], side = 1, line = 2.9, adj = 0, cex = 0.68, font = 2, col = hcol[1])
  mtext(hint[2], side = 1, line = 2.9, adj = 1, cex = 0.68, font = 2, col = hcol[2])
}

# single-panel volcano with its own legend row; writes <file>.png and the same plot as .pdf
volcano_png <- function(r, p, file, n_lab = 14, xpad = 2.2, pdf = TRUE) {
  draw <- function() {
    layout(matrix(1:2), heights = c(1, 0.1))
    par(mar = c(3.8, 3.4, 2.8, 1.2), fg = ink2)
    volcano(r, p, n_lab = n_lab, xpad = xpad)
    par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", horiz = TRUE, bty = "n", pt.cex = 1.1, cex = 0.62, text.col = ink, pch = 21,
           pt.bg = c(vcols(p), col_ns), col = "#fcfcfb",
           legend = c(paste0(vgroups(p)[1], ", FDR<", vfdr(p)), "FDR<.25", paste0(vgroups(p)[2], ", FDR<", vfdr(p)), "FDR<.25", "n.s."))
  }
  png(file, width = 1500, height = 1400, res = 220, type = "cairo", bg = "#fcfcfb"); draw(); dev.off()
  if (pdf) {   # same 1500 x 1400 px at 220 dpi -> 6.8 x 6.4 in, so text sizes match the PNG
    cairo_pdf(sub("\\.png$", ".pdf", file), width = 1500 / 220, height = 1400 / 220, bg = "#fcfcfb")
    draw(); dev.off()
  }
}
