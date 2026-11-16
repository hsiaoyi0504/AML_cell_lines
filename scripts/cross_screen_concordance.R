# Cross-screen concordance of mTOR-High vs mTOR-Low differences (RNA-defined groups):
#   * CRISPR (Chronos) vs RNAi (DEMETER2)         -> cross_CRISPR_RNAi_*
#   * Drug AUC CTRP vs GDSC2 (largest overlap)     -> cross_CTRP_GDSC2_*
#   * Drug AUC CTRP vs PRISM, GDSC2 vs PRISM       -> cross_CTRP_PRISM_*, cross_GDSC2_PRISM_*
# Reads results written by scripts/compare_mTOR_RNA_groups.R. Base R only.
#
# Usage: Rscript scripts/cross_screen_concordance.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"
fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))

col_hi <- "#2a78d6"; col_lo <- "#e34948"; col_ns <- "#c9c8c3"
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"

# a = first screen, b = second; `active` = per-screen threshold on the more
# sensitive group's mean (feature is a real dependency / active drug somewhere)
configs <- list(
  list(prefix = "cross_CRISPR_RNAi", a = "crispr", b = "rnai", la = "CRISPR", lb = "RNAi",
       unit_a = "Chronos gene effect", unit_b = "DEMETER2 score",
       active_a = -0.5, active_b = -0.3, active_name = "is_dependency",
       title = "CRISPR vs RNAi: differential dependency", word = "dependent", n_lab = 16),
  list(prefix = "cross_CTRP_GDSC2", a = "drug_CTRP", b = "drug_GDSC2", la = "CTRP", lb = "GDSC2",
       unit_a = "CTRP AUC", unit_b = "GDSC2 AUC",
       active_a = 0.8, active_b = 0.8, active_name = "is_active",
       title = "CTRP vs GDSC2: differential drug sensitivity", word = "sensitive", n_lab = 25),
  list(prefix = "cross_CTRP_PRISM", a = "drug_CTRP", b = "drug_PRISM", la = "CTRP", lb = "PRISM",
       unit_a = "CTRP AUC", unit_b = "PRISM AUC",
       active_a = 0.8, active_b = 0.8, active_name = "is_active",
       title = "CTRP vs PRISM: differential drug sensitivity", word = "sensitive", n_lab = 25),
  list(prefix = "cross_GDSC2_PRISM", a = "drug_GDSC2", b = "drug_PRISM", la = "GDSC2", lb = "PRISM",
       unit_a = "GDSC2 AUC", unit_b = "PRISM AUC",
       active_a = 0.8, active_b = 0.8, active_name = "is_active",
       title = "GDSC2 vs PRISM: differential drug sensitivity", word = "sensitive", n_lab = 25)
)

rapalog <- "SIROLIMUS|RAPAMYCIN|EVEROLIMUS|TEMSIROLIMUS|RIDAFOROLIMUS"

spread <- function(y, gap, lo, hi) {
  o <- order(y); ys <- y[o]
  for (i in seq_along(ys)[-1]) if (ys[i] - ys[i - 1] < gap) ys[i] <- ys[i - 1] + gap
  if (max(ys) > hi) ys <- ys - (max(ys) - hi)
  ys <- pmax(ys, lo); y[o] <- ys; y
}

run_cross <- function(cf) {
  A <- results[[cf$a]]; B <- results[[cf$b]]
  cols <- c("feature", "n_High", "n_Low", "mean_High", "mean_Low", "diff_High_minus_Low",
            "cohens_d", "t", "p_welch", "fdr_welch")
  m <- merge(A[, cols], B[, cols], by = "feature", suffixes = c("_a", "_b"))
  names(m) <- sub("_a$", paste0("_", tolower(cf$la)), sub("_b$", paste0("_", tolower(cf$lb)), names(m)))
  ka <- tolower(cf$la); kb <- tolower(cf$lb)
  g <- function(stem, k) m[[paste0(stem, "_", k)]]

  if (!is.null(A$CompoundName)) {          # drugs: name + target, from PortalCompounds
    ann <- A[match(m$feature, A$feature), c("CompoundName", "GeneSymbolOfTargets", "TargetOrMechanism")]
    m <- cbind(ann, m)
    m$label <- ifelse(is.na(m$CompoundName) | m$CompoundName == "", m$feature,
                      substr(sub(" \\(.*\\)$", "", m$CompoundName), 1, 26))
  } else m$label <- sub(" \\([0-9]+\\)$", "", m$feature)

  # signed z: positive = more dependent/sensitive in mTOR-High (lower score in High)
  sz <- function(p, d) -sign(d) * qnorm(p / 2, lower.tail = FALSE)
  m$z_a <- sz(g("p_welch", ka), g("diff_High_minus_Low", ka))
  m$z_b <- sz(g("p_welch", kb), g("diff_High_minus_Low", kb))
  m$z_stouffer   <- (m$z_a + m$z_b) / sqrt(2)
  m$p_stouffer   <- 2 * pnorm(-abs(m$z_stouffer))
  m$fdr_stouffer <- p.adjust(m$p_stouffer, "BH")
  m$same_direction <- sign(m$z_a) == sign(m$z_b)
  m$favours <- ifelse(m$z_stouffer > 0, "mTOR-High", "mTOR-Low")
  m[[cf$active_name]] <- pmin(g("mean_High", ka), g("mean_Low", ka)) < cf$active_a |
                         pmin(g("mean_High", kb), g("mean_Low", kb)) < cf$active_b
  m$concordant_hit <- m$same_direction & g("p_welch", ka) < 0.05 & g("p_welch", kb) < 0.05
  m <- m[order(m$p_stouffer), ]
  names(m)[names(m) == "z_a"] <- paste0("z_", ka); names(m)[names(m) == "z_b"] <- paste0("z_", kb)

  write.table(m, file.path(res_dir, paste0(cf$prefix, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(m[m$concordant_hit, ], file.path(res_dir, paste0(cf$prefix, "_concordant_hits.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ----------------------------------------------------------- summary ----
  za <- m[[paste0("z_", ka)]]; zb <- m[[paste0("z_", kb)]]
  act <- m[[cf$active_name]]
  sink(file.path(res_dir, paste0(cf$prefix, "_summary.txt")))
  la <- rownames(attr(A, "X")); lb <- rownames(attr(B, "X"))
  cat(sprintf("Features in both screens: %d\nCell lines: %s %d | %s %d | shared %d\n\n",
              nrow(m), cf$la, length(la), cf$lb, length(lb), length(intersect(la, lb))))
  cat("== Agreement of signed z, all shared features (Spearman)\n")
  print(cor.test(za, zb, method = "spearman", exact = FALSE))
  cat(sprintf("== Restricted to %s (n = %d)\n", cf$active_name, sum(act)))
  print(cor.test(za[act], zb[act], method = "spearman", exact = FALSE))
  for (thr in c(0.05, 0.01)) {
    s <- g("p_welch", ka) < thr & g("p_welch", kb) < thr
    cat(sprintf("== p < %.2f in both: %d; same direction: %d\n", thr, sum(s), sum(m$same_direction[s])))
    if (sum(s) > 0) print(binom.test(sum(m$same_direction[s]), sum(s), 0.5))
    cat(sprintf("   expected overlap if independent: %.1f\n",
                nrow(m) * mean(g("p_welch", ka) < thr) * mean(g("p_welch", kb) < thr)))
  }
  sig_a <- g("p_welch", ka) < 0.05
  cat(sprintf("== Direction replication: of %d features with %s p < 0.05, %d (%.0f%%) have the same sign in %s\n",
              sum(sig_a), cf$la, sum(m$same_direction[sig_a]), 100 * mean(m$same_direction[sig_a]), cf$lb))
  print(binom.test(sum(m$same_direction[sig_a]), sum(sig_a), 0.5))
  cat("== Concordant hits (same direction, p < 0.05 in both):", sum(m$concordant_hit),
      "| favouring mTOR-High:", sum(m$concordant_hit & m$favours == "mTOR-High"),
      "| mTOR-Low:", sum(m$concordant_hit & m$favours == "mTOR-Low"), "\n")
  cat("   Stouffer p < 0.05:", sum(m$p_stouffer < 0.05), " p < 0.01:", sum(m$p_stouffer < 0.01), "\n\n")
  cat("== Top features by combined (Stouffer) evidence\n")
  tab <- data.frame(label = m$label, favours = m$favours, concordant = m$concordant_hit,
                    dA = round(g("diff_High_minus_Low", ka), 3), pA = signif(g("p_welch", ka), 2),
                    dB = round(g("diff_High_minus_Low", kb), 3), pB = signif(g("p_welch", kb), 2),
                    p_stouffer = signif(m$p_stouffer, 2))
  if (!is.null(m$TargetOrMechanism)) tab$target <- substr(m$TargetOrMechanism, 1, 40)
  names(tab)[4:7] <- c(paste0("d", cf$la), paste0("p", cf$la), paste0("d", cf$lb), paste0("p", cf$lb))
  print(head(tab, 30), row.names = FALSE)
  sink()

  # ------------------------------------------------------- scatter plot ---
  draw_scatter <- function() {
  layout(matrix(1:2), heights = c(1, 0.09))
  par(mar = c(4.2, 4.2, 3.2, 1.2), fg = ink2)
  lim <- max(abs(c(za, zb))) * 1.08 * c(-1, 1); xlim <- lim * 1.85
  col <- ifelse(m$concordant_hit, ifelse(m$favours == "mTOR-High", col_hi, col_lo), col_ns)
  plot(NA, xlim = xlim, ylim = lim, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(lim), v = pretty(xlim), col = grid_col, lwd = 0.8)
  abline(h = 0, v = 0, col = ink2, lwd = 0.8)
  z05 <- qnorm(0.975)
  abline(h = c(-z05, z05), v = c(-z05, z05), lty = 2, col = ink2, lwd = 0.6)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.8, tcl = -0.25, mgp = c(2, 0.4, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.8, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
  o <- order(m$concordant_hit)
  points(za[o], zb[o], pch = 21, bg = col[o], col = bg, lwd = 0.4,
         cex = ifelse(m$concordant_hit[o], 0.95, if (nrow(m) < 500) 0.75 else 0.5))
  # label concordant hits first (active ones first), then combined (Stouffer) p < 0.05
  lab <- unique(c(which(m$concordant_hit & act), which(m$concordant_hit),
                  if (nrow(m) < 500) which(m$p_stouffer < 0.05)))
  lab <- unique(c(head(lab, cf$n_lab), grep(rapalog, m$label)))   # always show rapalogs
  for (side in c(TRUE, FALSE)) {
    i <- lab[(za[lab] > 0) == side]; if (!length(i)) next
    ly <- spread(zb[i], diff(lim) * 0.045, lim[1] + diff(lim) * 0.06, lim[2] - diff(lim) * 0.06)
    w  <- max(strwidth(m$label[i], cex = 0.6)) * 1.1
    lx <- if (side) rep(xlim[2] - w, length(i)) else rep(xlim[1] + w, length(i))
    segments(za[i], zb[i], lx, ly, col = ink2, lwd = 0.5)
    text(lx, ly, m$label[i], pos = if (side) 4 else 2, offset = 0.15, cex = 0.6, col = ink, xpd = NA)
  }
  s05 <- g("p_welch", ka) < 0.05 & g("p_welch", kb) < 0.05
  mtext(cf$title, side = 3, line = 1.6, adj = 0, font = 2, cex = 0.9, col = ink)
  mtext(sprintf("%d shared · rho = %.2f · p<0.05 in both: %d, %d same direction (%d concordant) · labelled: combined p < 0.05",
                nrow(m), cor(za, zb, method = "spearman"), sum(s05), sum(m$same_direction[s05]),
                sum(m$concordant_hit)), side = 3, line = 0.5, adj = 0, cex = 0.6, col = ink2)
  mtext(sprintf("%s signed z  (+ = more %s in mTOR-High)", cf$la, cf$word), side = 1, line = 2.2, cex = 0.68, col = ink2)
  mtext(sprintf("%s signed z  (+ = more %s in mTOR-High)", cf$lb, cf$word), side = 2, line = 2.4, cex = 0.68, col = ink2)
  text(xlim[2], lim[2], sprintf("both: High more %s", cf$word), adj = c(1, 1), cex = 0.58, col = col_hi, font = 2)
  text(xlim[1], lim[1], sprintf("both: Low more %s", cf$word), adj = c(0, 0), cex = 0.58, col = col_lo, font = 2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = 21, pt.cex = 1.1, cex = 0.62, text.col = ink, col = bg,
         pt.bg = c(col_hi, col_lo, col_ns),
         legend = c(sprintf("Concordant, more %s in mTOR-High", cf$word),
                    sprintf("Concordant, more %s in mTOR-Low", cf$word), "Other"))
  }
  fb <- file.path(fig_dir, paste0(cf$prefix, "_scatter"))
  png(paste0(fb, ".png"), width = 1800, height = 1650, res = 220, type = "cairo", bg = bg); draw_scatter(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 1800 / 220, height = 1650 / 220, bg = bg); draw_scatter(); dev.off()

  # ------------------------------------- per-line strips, top features ---
  top <- head(m[m$concordant_hit, ], 6)
  if (nrow(top) < 6) top <- head(unique(rbind(top, m)), 6)
  XA <- attr(A, "X"); XB <- attr(B, "X")
  draw_strips <- function() {
  par(mfcol = c(2, nrow(top)), mar = c(2.8, 3.6, 3, 0.8), fg = ink2)
  set.seed(1)
  for (k in seq_len(nrow(top))) for (scr in c("a", "b")) {
    X <- if (scr == "a") XA else XB
    v <- X[, top$feature[k]]; grp <- cl$mTOR_group[match(rownames(X), cl$cell_line)]
    ok <- !is.na(v); v <- v[ok]; grp <- factor(grp[ok], levels = c("Low", "High"))
    plot(NA, xlim = c(0.4, 2.6), ylim = range(v) + c(-0.05, 0.08) * diff(range(v)), axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(v), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25)
    axis(1, at = 1:2, labels = c(sprintf("Low\nn=%d", sum(grp == "Low")), sprintf("High\nn=%d", sum(grp == "High"))),
         col = NA, col.axis = ink, cex.axis = 0.75, padj = 0.5)
    points(as.integer(grp) + runif(length(v), -0.18, 0.18), v, pch = 21, cex = 1.1, lwd = 0.5, col = bg,
           bg = ifelse(grp == "High", col_hi, col_lo))
    for (j in 1:2) segments(j - 0.3, median(v[as.integer(grp) == j]), j + 0.3, lwd = 2, col = ink)
    k_scr <- if (scr == "a") ka else kb
    mtext(if (scr == "a") cf$unit_a else cf$unit_b, side = 2, line = 2.5, cex = 0.55, col = ink2)
    mtext(sprintf("%s · %s", top$label[k], if (scr == "a") cf$la else cf$lb), side = 3, line = 1.2,
          adj = 0, font = 2, cex = 0.68, col = ink)
    mtext(sprintf("p = %.1e", top[[paste0("p_welch_", k_scr)]][k]), side = 3, line = 0.3, adj = 0, cex = 0.6, col = ink2)
  }
  }
  fb <- file.path(fig_dir, paste0(cf$prefix, "_top_strip"))
  png(paste0(fb, ".png"), width = 600 * nrow(top), height = 1500, res = 200, type = "cairo", bg = bg); draw_strips(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 600 * nrow(top) / 200, height = 1500 / 200, bg = bg); draw_strips(); dev.off()
  cat(readLines(file.path(res_dir, paste0(cf$prefix, "_summary.txt"))), sep = "\n")
}

for (cf in configs) { cat("\n########", cf$prefix, "########\n"); run_cross(cf) }
