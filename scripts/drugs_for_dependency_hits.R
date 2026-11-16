# Drug response (AUC) MYC-High vs mTOR-High for drugs hitting the concordant CRISPR/RNAi
# dependencies from scripts/cross_CRISPR_RNAi_MYC_vs_mTOR.R. Single agents only; combinations
# are skipped. Lower AUC = more sensitive.
#
# Run from the project root: Rscript scripts/drugs_for_dependency_hits.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp  <- ifelse(cl$MYC_group == "High" & cl$mTOR_group == "Low", "MYC-High",
           ifelse(cl$mTOR_group == "High" & cl$MYC_group == "Low", "mTOR-High", NA))
dep     <- read.delim(file.path(res_dir, "cross_CRISPR_RNAi_MYC_vs_mTOR_all.tsv"))

colA <- "#eb6834"; colB <- "#4a3aa7"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"

# dependency gene -> drugs; `link` says how directly the drug hits the gene. `dep` sets the more-dependent
# group for genes that are not concordant CRISPR/RNAi hits (JAK2, KIT: CRISPR dependency + hotspot mutations
# in MYC-High only; see heatmap annotation rows)
panel <- list(
  list(gene = "BCL2L1", link = "direct",    drugs = c("NAVITOCLAX", "ABT-737", "WEHI-539")),
  list(gene = "BCL2 (control)", link = "control", drugs = c("VENETOCLAX")),
  list(gene = "MCL1 (control)", link = "control", drugs = c("AZD5991")),
  list(gene = "PSMB5",  link = "direct",    drugs = c("BORTEZOMIB", "CARFILZOMIB")),
  list(gene = "SOS1 / PTPN11 / GRB2", link = "pathway proxy: MEK, not a direct test", drugs = c("TRAMETINIB", "SELUMETINIB", "COBIMETINIB")),
  list(gene = "CCND3",  link = "partner (CDK4/6)", drugs = c("PALBOCICLIB", "RIBOCICLIB", "ABEMACICLIB")),
  list(gene = "DNM1L",  link = "direct",    drugs = c("MDIVI-1")),
  list(gene = "JAK2", link = "direct; CRISPR + V617F", dep = "MYC-High",
       drugs = c("RUXOLITINIB", "NVP-BSK805", "FEDRATINIB", "MOMELOTINIB", "AZ 960")),
  list(gene = "KIT", link = "direct; CRISPR + N822K/D816V", dep = "MYC-High",
       drugs = c("IMATINIB", "NILOTINIB", "AB1010"))
)
screens <- c(CTRP = "drug_CTRP", GDSC2 = "drug_GDSC2")   # PRISM left out

rows <- list(); vals <- list()
for (pn in panel) for (dn in pn$drugs) for (s in names(screens)) {
  r <- results[[screens[[s]]]]; i <- which(r$CompoundName == dn); if (!length(i)) next
  X <- attr(r, "X"); v <- X[, r$feature[i[1]]]; g <- cl$grp[match(rownames(X), cl$cell_line)]
  a <- v[g %in% "MYC-High" & !is.na(v)]; b <- v[g %in% "mTOR-High" & !is.na(v)]
  if (length(a) < 3 || length(b) < 3) next
  tt <- t.test(a, b)
  rows[[length(rows) + 1]] <- data.frame(dep_gene = pn$gene, link = pn$link, dep_set = if (is.null(pn$dep)) NA else pn$dep, drug = dn, screen = s,
    n_MYC = length(a), n_mTOR = length(b), AUC_MYC = mean(a), AUC_mTOR = mean(b),
    diff_MYC_minus_mTOR = mean(a) - mean(b), ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
    p_welch = tt$p.value, p_wilcox = suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value))
  vals[[paste(dn, s)]] <- list(a = a, b = b)
}
tab <- do.call(rbind, rows)
gd <- sub(" .*", "", tab$dep_gene)
k <- match(gd, dep$gene)
tab$dep_more_dependent_in <- ifelse(tab$link == "control", NA, dep$favours[k])
tab$dep_more_dependent_in[!is.na(tab$dep_set)] <- tab$dep_set[!is.na(tab$dep_set)]
tab$dep_p_CRISPR <- ifelse(tab$link == "control", NA, signif(dep$p_welch_crispr[k], 2))
tab$dep_p_RNAi   <- ifelse(tab$link == "control", NA, signif(dep$p_welch_rnai[k], 2))
write.table(tab, file.path(res_dir, "drugs_for_dependency_hits_MYC_vs_mTOR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(data.frame(tab[, c("dep_gene", "drug", "screen", "n_MYC", "n_mTOR")],
                 AUC_MYC = round(tab$AUC_MYC, 3), AUC_mTOR = round(tab$AUC_mTOR, 3),
                 diff = round(tab$diff_MYC_minus_mTOR, 3), p = signif(tab$p_welch, 2),
                 p_wilcox = signif(tab$p_wilcox, 2)), row.names = FALSE)

# --------------------------------------------------- forest plot ---------
draw_forest <- function() {
  n <- nrow(tab)
  blk <- cumsum(c(TRUE, tab$dep_gene[-1] != tab$dep_gene[-n]))     # one header line per gene block
  y <- -(seq_len(n) + 1.2 * blk); y <- y - min(y) + 1
  layout(matrix(1:2), heights = c(1, 0.05))
  par(mar = c(4.2, 13, 4.2, 6), fg = ink2, bg = bg)
  xr <- range(c(tab$ci_lo, tab$ci_hi, 0)); xr <- xr + c(-0.05, 0.05) * diff(xr)
  plot(NA, xlim = xr, ylim = c(0.4, max(y) + 1.3), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  abline(v = pretty(xr), col = grid_col, lwd = 0.8); abline(v = 0, col = ink2, lwd = 1)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, mgp = c(2, 0.45, 0))
  col <- ifelse(tab$diff_MYC_minus_mTOR < 0, colA, colB)
  sig <- tab$p_welch < 0.05
  segments(tab$ci_lo, y, tab$ci_hi, y, col = col, lwd = 2)
  points(tab$diff_MYC_minus_mTOR, y, pch = 21, cex = 1.3, lwd = 1.2, col = col, bg = ifelse(sig, col, bg))
  text(xr[1] - 0.01 * diff(xr), y, sprintf("%s  ·  %s", tab$drug, tab$screen), adj = c(1, 0.5), xpd = NA,
       cex = 0.78, col = ink)
  text(xr[2] + 0.02 * diff(xr), y, sprintf("p = %s", formatC(tab$p_welch, format = "g", digits = 2)),
       adj = c(0, 0.5), xpd = NA, cex = 0.72, col = ifelse(sig, ink, ink2), font = ifelse(sig, 2, 1))
  # gene block headers
  for (b in unique(blk)) {
    i <- which(blk == b); g <- tab$dep_gene[i[1]]; ln <- tab$link[i[1]]; dm <- tab$dep_more_dependent_in[i[1]]
    lab <- if (ln == "control") g else sprintf("%s  (%s; dependency: %s)", g, ln, dm)
    x0 <- grconvertX(0.03, "ndc", "user")
    text(x0, max(y[i]) + 1.05, lab, adj = c(0, 0.5), xpd = NA,
         cex = 0.78, font = 2, col = if (ln == "control") ink2 else if (dm == "MYC-High") colA else colB)
  }
  x0 <- grconvertX(0.03, "ndc", "user")
  text(x0, grconvertY(0.985, "ndc", "user"), "Drugs hitting MYC-High / mTOR-High dependencies",
       adj = c(0, 1), xpd = NA, font = 2, cex = 1.05, col = ink)
  text(x0, grconvertY(0.965, "ndc", "user"),
       "Mean AUC difference MYC-High minus mTOR-High, 95% CI · filled = p < 0.05 (Welch) · colour = more sensitive group",
       adj = c(0, 1), xpd = NA, cex = 0.7, col = ink2)
  mtext("AUC difference (MYC-High − mTOR-High)", side = 1, line = 2.3, cex = 0.8, col = ink2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = 21, pt.cex = 1.2, cex = 0.78, text.col = ink, col = c(colA, colB),
         pt.bg = c(colA, colB), legend = c("MYC-High more sensitive (lower AUC)", "mTOR-High more sensitive"))
}
fb <- file.path(fig_dir, "drugs_for_dependency_hits_forest")
h <- 700 + 80 * nrow(tab) + 90 * length(unique(tab$dep_gene))
png(paste0(fb, ".png"), width = 2000, height = h, res = 220, type = "cairo", bg = bg); draw_forest(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2000 / 220, height = h / 220, bg = bg); draw_forest(); dev.off()

# ------------------------- boxplots: MYC-High and mTOR-High AUC per row --
tint <- function(col, f = 0.75) colorRampPalette(c(col, "#ffffff"))(100)[round(f * 100)]
draw_box <- function() {
  n <- nrow(tab)
  blk <- cumsum(c(TRUE, tab$dep_gene[-1] != tab$dep_gene[-n]))
  y <- -(seq_len(n) + 1.0 * blk); y <- y - min(y) + 1
  layout(matrix(1:2), heights = c(1, 0.05))
  par(mar = c(4.2, 13, 4.2, 6), fg = ink2, bg = bg)
  xr <- c(0, 1.05)
  plot(NA, xlim = xr, ylim = c(0.45, max(y) + 1.1), axes = FALSE, xlab = "", ylab = "", yaxs = "i", xaxs = "i")
  abline(v = pretty(xr), col = grid_col, lwd = 0.8)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25, mgp = c(2, 0.45, 0))
  set.seed(1)
  for (k in seq_len(n)) {
    vv <- vals[[paste(tab$drug[k], tab$screen[k])]]
    at <- y[k] + c(0.2, -0.2)
    boxplot(list(vv$a, vv$b), at = at, horizontal = TRUE, add = TRUE, axes = FALSE, outline = FALSE,
            boxwex = 0.34, col = c(tint(colA), tint(colB)), border = c(colA, colB), lwd = 1,
            medlwd = 2, staplewex = 0.4, whisklty = 1)
    points(vv$a, at[1] + runif(length(vv$a), -0.08, 0.08), pch = 21, cex = 0.75, lwd = 0.4, col = bg, bg = colA)
    points(vv$b, at[2] + runif(length(vv$b), -0.08, 0.08), pch = 21, cex = 0.75, lwd = 0.4, col = bg, bg = colB)
  }
  sig <- tab$p_welch < 0.05
  text(xr[1] - 0.01 * diff(xr), y, sprintf("%s  ·  %s", tab$drug, tab$screen), adj = c(1, 0.5), xpd = NA,
       cex = 0.78, col = ink)
  text(xr[2] + 0.02 * diff(xr), y, sprintf("p = %s", formatC(tab$p_welch, format = "g", digits = 2)),
       adj = c(0, 0.5), xpd = NA, cex = 0.72, col = ifelse(sig, ink, ink2), font = ifelse(sig, 2, 1))
  x0 <- grconvertX(0.03, "ndc", "user")
  for (b in unique(blk)) {
    i <- which(blk == b); g <- tab$dep_gene[i[1]]; ln <- tab$link[i[1]]; dm <- tab$dep_more_dependent_in[i[1]]
    lab <- if (ln == "control") g else sprintf("%s  (%s; dependency: %s)", g, ln, dm)
    text(x0, max(y[i]) + 0.75, lab, adj = c(0, 0.5), xpd = NA,
         cex = 0.78, font = 2, col = if (ln == "control") ink2 else if (dm == "MYC-High") colA else colB)
  }
  text(x0, grconvertY(0.985, "ndc", "user"), "Drugs hitting MYC-High / mTOR-High dependencies",
       adj = c(0, 1), xpd = NA, font = 2, cex = 1.05, col = ink)
  text(x0, grconvertY(0.965, "ndc", "user"),
       "AUC per cell line (lower = more sensitive) · p = Welch t-test, MYC-High vs mTOR-High",
       adj = c(0, 1), xpd = NA, cex = 0.7, col = ink2)
  mtext("AUC", side = 1, line = 2.3, cex = 0.8, col = ink2)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", pch = 22, pt.cex = 1.6, cex = 0.78, text.col = ink, col = c(colA, colB),
         pt.bg = c(tint(colA), tint(colB)), legend = c("MYC-High", "mTOR-High"))
}
fb <- file.path(fig_dir, "drugs_for_dependency_hits_boxplot")
h <- 700 + 115 * nrow(tab) + 90 * length(unique(tab$dep_gene))
png(paste0(fb, ".png"), width = 2000, height = h, res = 220, type = "cairo", bg = bg); draw_box(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2000 / 220, height = h / 220, bg = bg); draw_box(); dev.off()

# ------------------------------------------- per-line strips, BCL-XL block --
key <- tab[tab$dep_gene %in% c("BCL2L1", "BCL2 (control)", "MCL1 (control)"), ]
draw_strips <- function() {
  par(mfrow = c(1, nrow(key)), mar = c(3, 3.6, 3.4, 0.6), fg = ink2, bg = bg)
  set.seed(1)
  for (k in seq_len(nrow(key))) {
    vv <- vals[[paste(key$drug[k], key$screen[k])]]; v <- c(vv$a, vv$b); g <- rep(1:2, c(length(vv$a), length(vv$b)))
    plot(NA, xlim = c(0.4, 2.6), ylim = range(v) + c(-0.05, 0.1) * diff(range(v)), axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(v), col = grid_col, lwd = 0.8)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.8, tcl = -0.25)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", c("MYC", "mTOR"), c(length(vv$a), length(vv$b))),
         col = NA, col.axis = ink, cex.axis = 0.8, padj = 0.5)
    points(g + runif(length(v), -0.18, 0.18), v, pch = 21, cex = 1.2, lwd = 0.5, col = bg, bg = c(colA, colB)[g])
    for (j in 1:2) segments(j - 0.3, median(v[g == j]), j + 0.3, lwd = 2, col = ink)
    if (k == 1) mtext("AUC", side = 2, line = 2.5, cex = 0.65, col = ink2)
    mtext(sprintf("%s · %s", key$drug[k], key$screen[k]), side = 3, line = 1.4, adj = 0, font = 2, cex = 0.66, col = ink)
    mtext(sprintf("%s · p = %.2g", key$dep_gene[k], key$p_welch[k]), side = 3, line = 0.4, adj = 0, cex = 0.6, col = ink2)
  }
}
fb <- file.path(fig_dir, "drugs_for_dependency_hits_BCL2family_strip")
png(paste0(fb, ".png"), width = 430 * nrow(key), height = 1050, res = 200, type = "cairo", bg = bg); draw_strips(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 430 * nrow(key) / 200, height = 1050 / 200, bg = bg); draw_strips(); dev.off()
