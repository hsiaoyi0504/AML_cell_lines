# CRISPR (Chronos) vs RNAi (DEMETER2) concordance of MYC-High/mTOR-Low vs mTOR-High/MYC-Low
# gene-effect differences (RNA-defined groups; lines high or low for both scores are left out).
# Per-line matrices come from all_results.rds (scripts/compare_mTOR_RNA_groups.R).
#
# Run from the project root: Rscript scripts/cross_CRISPR_RNAi_MYC_vs_mTOR.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp  <- ifelse(cl$MYC_group == "High" & cl$mTOR_group == "Low", "MYC-High",
           ifelse(cl$mTOR_group == "High" & cl$MYC_group == "Low", "mTOR-High", NA))

colA <- "#eb6834"; colB <- "#4a3aa7"; col_ns <- "#c9c8c3"     # MYC-High, mTOR-High, other
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
n_lab <- 18            # concordant hits labelled on the scatter
# genes supported by both lines of evidence: concordant dependency AND a drug against it that is selectively
# active in the same group (scripts/drugs_MYC_vs_mTOR_screen.R). A target nominated here is one pre-specified
# hypothesis, so the drug test uses the raw Welch p (no FDR, no general-sensitivity adjustment); p < 0.05.
drug_file <- file.path(res_dir, "drugs_MYC_High_vs_mTOR_High_CTRP_GDSC2.tsv")
dual <- if (file.exists(drug_file)) {
  dd <- read.delim(drug_file); dd <- dd[dd$dep_targets != "" & dd$consistent %in% TRUE & dd$p_welch < 0.05, ]
  dd <- dd[order(dd$p_welch), ]; dd[!duplicated(dd$dep_targets), c("dep_targets", "drug", "screen", "p_welch")]
} else data.frame(dep_targets = character(0), drug = character(0), screen = character(0), p_welch = numeric(0))
dep_crispr <- -0.5; dep_rnai <- -0.3   # "real dependency" in at least one group

# Welch t-test MYC-High vs mTOR-High for every column of X
welch <- function(X, min_n = 3) {
  grp <- cl$grp[match(rownames(X), cl$cell_line)]
  A <- X[grp %in% "MYC-High", , drop = FALSE]; B <- X[grp %in% "mTOR-High", , drop = FALSE]
  nA <- colSums(!is.na(A)); nB <- colSums(!is.na(B))
  mA <- colMeans(A, na.rm = TRUE); mB <- colMeans(B, na.rm = TRUE)
  vA <- colSums(sweep(A, 2, mA)^2, na.rm = TRUE) / (nA - 1)
  vB <- colSums(sweep(B, 2, mB)^2, na.rm = TRUE) / (nB - 1)
  se <- sqrt(vA / nA + vB / nB); t <- (mA - mB) / se
  df <- se^4 / ((vA / nA)^2 / (nA - 1) + (vB / nB)^2 / (nB - 1))
  r <- data.frame(feature = colnames(X), n_MYC = nA, n_mTOR = nB, mean_MYC = mA, mean_mTOR = mB,
                  diff_MYC_minus_mTOR = mA - mB, p_welch = 2 * pt(-abs(t), df))
  r <- r[nA >= min_n & nB >= min_n & vA + vB > 0, ]
  r$fdr_welch <- p.adjust(r$p_welch, "BH")
  attr(r, "grp") <- grp; r
}
XC <- attr(results$crispr, "X"); XR <- attr(results$rnai, "X")
C <- welch(XC); R <- welch(XR)
gC <- attr(C, "grp"); gR <- attr(R, "grp")

m <- merge(C, R, by = "feature", suffixes = c("_crispr", "_rnai"))
m$gene <- sub(" \\([0-9]+\\)$", "", m$feature)
# signed z: positive = more dependent in MYC-High (lower gene effect in MYC-High)
sz <- function(p, d) -sign(d) * qnorm(p / 2, lower.tail = FALSE)
m$z_crispr <- sz(m$p_welch_crispr, m$diff_MYC_minus_mTOR_crispr)
m$z_rnai   <- sz(m$p_welch_rnai,   m$diff_MYC_minus_mTOR_rnai)
m$z_stouffer   <- (m$z_crispr + m$z_rnai) / sqrt(2)
m$p_stouffer   <- 2 * pnorm(-abs(m$z_stouffer))
m$fdr_stouffer <- p.adjust(m$p_stouffer, "BH")
m$same_direction <- sign(m$z_crispr) == sign(m$z_rnai)
m$favours <- ifelse(m$z_stouffer > 0, "MYC-High", "mTOR-High")   # group that is more dependent
m$is_dependency <- pmin(m$mean_MYC_crispr, m$mean_mTOR_crispr) < dep_crispr |
                   pmin(m$mean_MYC_rnai, m$mean_mTOR_rnai) < dep_rnai
m$concordant_hit <- m$same_direction & m$p_welch_crispr < 0.05 & m$p_welch_rnai < 0.05
m <- m[order(m$p_stouffer), ]
pre <- "cross_CRISPR_RNAi_MYC_vs_mTOR"
write.table(m, file.path(res_dir, paste0(pre, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(m[m$concordant_hit, ], file.path(res_dir, paste0(pre, "_concordant_hits.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------- summary ------
zc <- m$z_crispr; zr <- m$z_rnai
sink(file.path(res_dir, paste0(pre, "_summary.txt")))
both <- intersect(rownames(XC)[!is.na(gC)], rownames(XR)[!is.na(gR)])
cat(sprintf("Genes in both screens: %d\nCRISPR lines MYC/mTOR = %d/%d | RNAi lines MYC/mTOR = %d/%d | in both screens: %d\n\n",
            nrow(m), sum(gC == "MYC-High", na.rm = TRUE), sum(gC == "mTOR-High", na.rm = TRUE),
            sum(gR == "MYC-High", na.rm = TRUE), sum(gR == "mTOR-High", na.rm = TRUE), length(both)))
cat("Per-screen FDR < 0.05: CRISPR", sum(C$fdr_welch < 0.05), "| RNAi", sum(R$fdr_welch < 0.05), "\n\n")
cat("== Spearman of signed z, all shared genes\n"); print(cor.test(zc, zr, method = "spearman", exact = FALSE))
cat(sprintf("== Restricted to dependencies (n = %d)\n", sum(m$is_dependency)))
print(cor.test(zc[m$is_dependency], zr[m$is_dependency], method = "spearman", exact = FALSE))
for (thr in c(0.05, 0.01)) {
  s <- m$p_welch_crispr < thr & m$p_welch_rnai < thr
  cat(sprintf("== p < %.2f in both: %d; same direction: %d; expected overlap if independent: %.1f\n", thr, sum(s),
              sum(m$same_direction[s]), nrow(m) * mean(m$p_welch_crispr < thr) * mean(m$p_welch_rnai < thr)))
  if (sum(s) > 0) print(binom.test(sum(m$same_direction[s]), sum(s), 0.5))
}
sig <- m$p_welch_crispr < 0.05
cat(sprintf("== Of %d genes with CRISPR p < 0.05, %d (%.0f%%) have the same sign in RNAi\n",
            sum(sig), sum(m$same_direction[sig]), 100 * mean(m$same_direction[sig])))
print(binom.test(sum(m$same_direction[sig]), sum(sig), 0.5))
cat("== Concordant hits (same direction, p < 0.05 in both):", sum(m$concordant_hit),
    "| more dependent in MYC-High:", sum(m$concordant_hit & m$favours == "MYC-High"),
    "| in mTOR-High:", sum(m$concordant_hit & m$favours == "mTOR-High"), "\n")
cat("   of which dependencies:", sum(m$concordant_hit & m$is_dependency), "\n")
cat("   Stouffer FDR < 0.05:", sum(m$fdr_stouffer < 0.05), " FDR < 0.25:", sum(m$fdr_stouffer < 0.25), "\n\n")
cat("== Concordant hits (dCRISPR/dRNAi = MYC-High minus mTOR-High; negative = more dependent in MYC-High)\n")
h <- m[m$concordant_hit, ]
print(data.frame(gene = h$gene, more_dep_in = h$favours, dependency = h$is_dependency,
                 dCRISPR = round(h$diff_MYC_minus_mTOR_crispr, 2), pCRISPR = signif(h$p_welch_crispr, 2),
                 dRNAi = round(h$diff_MYC_minus_mTOR_rnai, 2), pRNAi = signif(h$p_welch_rnai, 2),
                 fdr_stouffer = signif(h$fdr_stouffer, 2)), row.names = FALSE)
sink()

# -------------------------------------------------------- scatter plot ------
draw_scatter <- function() {
  layout(matrix(1:2), heights = c(1, 0.11))
  par(mar = c(4.4, 4.4, 3.4, 9.5), fg = ink2, bg = bg)
  lim <- max(abs(c(zc, zr))) * 1.08 * c(-1, 1); xlim <- lim * 1.5
  col <- ifelse(m$concordant_hit, ifelse(m$favours == "MYC-High", colA, colB), col_ns)
  plot(NA, xlim = xlim, ylim = lim, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(lim), v = pretty(xlim), col = grid_col, lwd = 0.8)
  abline(h = 0, v = 0, col = ink2, lwd = 0.8)
  z05 <- qnorm(0.975); abline(h = c(-z05, z05), v = c(-z05, z05), lty = 2, col = ink2, lwd = 0.6)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.9, tcl = -0.25, mgp = c(2, 0.45, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.9, tcl = -0.25, las = 1, mgp = c(2, 0.5, 0))
  o <- order(m$concordant_hit)
  points(zc[o], zr[o], pch = 21, bg = col[o], col = bg, lwd = 0.4, cex = ifelse(m$concordant_hit[o], 1, 0.5))
  # label concordant dependencies first, then other concordant hits, by combined evidence
  lab <- head(unique(c(which(m$concordant_hit & m$is_dependency), which(m$concordant_hit))), n_lab)
  for (side in c(TRUE, FALSE)) {
    i <- lab[(zc[lab] > 0) == side]; if (!length(i)) next
    gap <- diff(lim) * 0.042; ly <- zr[i]; oo <- order(ly); ys <- ly[oo]
    for (k in seq_along(ys)[-1]) if (ys[k] - ys[k - 1] < gap) ys[k] <- ys[k - 1] + gap
    if (max(ys) > lim[2] - gap) ys <- ys - (max(ys) - lim[2] + gap)
    ly[oo] <- ys
    w  <- max(strwidth(m$gene[i], cex = 0.75, font = 2)) * 1.1
    lx <- if (side) rep(xlim[2] - w, length(i)) else rep(xlim[1] + w, length(i))
    segments(zc[i], zr[i], lx, ly, col = ink2, lwd = 0.5)
    text(lx, ly, m$gene[i], pos = if (side) 4 else 2, offset = 0.15, cex = 0.75, col = ink, xpd = NA,
         font = ifelse(m$is_dependency[i], 2, 1))
    k <- match(m$gene[i], dual$dep_targets); d <- which(!is.na(k))      # drug note right of the gene name
    if (length(d)) text(lx[d] + strwidth(m$gene[i][d], cex = 0.75, font = 2) + strwidth("  ", cex = 0.75), ly[d],
                        sprintf("\u25c2 %s (%s), p = %.2g", tolower(dual$drug[k[d]]), dual$screen[k[d]], dual$p_welch[k[d]]), adj = c(0, 0.5),
                        cex = 0.66, font = 3, col = ink, xpd = NA)
  }
  # dual-evidence genes: black ring around the point
  h <- match(dual$dep_targets, m$gene)
  points(zc[h], zr[h], pch = 21, cex = 2.1, lwd = 1.8, col = ink, bg = NA)
  s05 <- m$p_welch_crispr < 0.05 & m$p_welch_rnai < 0.05
  mtext("CRISPR vs RNAi: MYC-High vs mTOR-High gene effect", side = 3, line = 1.7, adj = 0, font = 2, cex = 1.05, col = ink)
  mtext(sprintf("%d shared genes · rho = %.2f · p<0.05 in both: %d, %d same direction · bold = dependency in either group",
                nrow(m), cor(zc, zr, method = "spearman"), sum(s05), sum(m$same_direction[s05])),
        side = 3, line = 0.55, adj = 0, cex = 0.7, col = ink2)
  mtext("CRISPR signed z  (+ = more dependent in MYC-High)", side = 1, line = 2.4, cex = 0.8, col = ink2)
  mtext("RNAi signed z  (+ = more dependent in MYC-High)", side = 2, line = 2.6, cex = 0.8, col = ink2)
  text(xlim[2], lim[2], "both: MYC-High more dependent", adj = c(1, 1), cex = 0.7, col = colA, font = 2)
  text(xlim[1], lim[1], "both: mTOR-High more dependent", adj = c(0, 0), cex = 0.7, col = colB, font = 2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", ncol = 2, bty = "n", pt.cex = 1.2, cex = 0.72, text.col = ink, col = c(bg, bg, bg, ink),
         pch = c(21, 21, 21, 1), pt.bg = c(colA, colB, col_ns, NA), pt.lwd = c(0.4, 0.4, 0.4, 1.8),
         legend = c("Concordant, MYC-High more dependent", "Concordant, mTOR-High more dependent", "Other",
                    "Also supported by drug response (ring; drug, screen, Welch p)"))
}
fb <- file.path(fig_dir, paste0(pre, "_scatter"))
png(paste0(fb, ".png"), width = 2250, height = 1750, res = 220, type = "cairo", bg = bg); draw_scatter(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2250 / 220, height = 1750 / 220, bg = bg); draw_scatter(); dev.off()

# ------------------------------- per-line strips, top concordant dependencies --
top <- m[m$concordant_hit & m$is_dependency, ]
top <- rbind(head(top[top$favours == "MYC-High", ], 4), head(top[top$favours == "mTOR-High", ], 4))
draw_strips <- function() {
  par(mfcol = c(2, nrow(top)), mar = c(3, 3.8, 3.2, 0.8), fg = ink2, bg = bg)
  set.seed(1)
  for (k in seq_len(nrow(top))) for (scr in c("CRISPR", "RNAi")) {
    X <- if (scr == "CRISPR") XC else XR; gg <- if (scr == "CRISPR") gC else gR
    v <- X[, top$feature[k]]; ok <- !is.na(v) & !is.na(gg)
    v <- v[ok]; g <- factor(gg[ok], levels = c("MYC-High", "mTOR-High"))
    plot(NA, xlim = c(0.4, 2.6), ylim = range(v) + c(-0.05, 0.1) * diff(range(v)), axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(v), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.85, tcl = -0.25)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", c("MYC", "mTOR"), as.vector(table(g))),
         col = NA, col.axis = ink, cex.axis = 0.85, padj = 0.5)
    points(as.integer(g) + runif(length(v), -0.18, 0.18), v, pch = 21, cex = 1.2, lwd = 0.5, col = bg,
           bg = ifelse(g == "MYC-High", colA, colB))
    for (j in 1:2) segments(j - 0.3, median(v[as.integer(g) == j]), j + 0.3, lwd = 2, col = ink)
    mtext(if (scr == "CRISPR") "Chronos gene effect" else "DEMETER2 score", side = 2, line = 2.6, cex = 0.62, col = ink2)
    mtext(sprintf("%s · %s", top$gene[k], scr), side = 3, line = 1.3, adj = 0, font = 2, cex = 0.8, col = ink)
    mtext(sprintf("p = %.1e", top[[paste0("p_welch_", tolower(scr))]][k]), side = 3, line = 0.35, adj = 0,
          cex = 0.68, col = ink2)
  }
}
if (nrow(top)) {
  fb <- file.path(fig_dir, paste0(pre, "_top_strip"))
  png(paste0(fb, ".png"), width = 520 * nrow(top), height = 1500, res = 200, type = "cairo", bg = bg); draw_strips(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 520 * nrow(top) / 200, height = 1500 / 200, bg = bg); draw_strips(); dev.off()
}
cat(readLines(file.path(res_dir, paste0(pre, "_summary.txt"))), sep = "\n")
