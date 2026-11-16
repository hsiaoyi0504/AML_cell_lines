# MTOR and mTOR-complex gene effect (CRISPR Chronos, RNAi DEMETER2) across the four RNA-defined
# categories: MYC-High, mTOR-High, both high, both low. Test = Welch MYC-High vs mTOR-High.
# More negative = stronger dependency (Chronos -1 = median common-essential gene).
#
# Run from the project root: Rscript scripts/MTOR_gene_effect.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$cat <- with(cl, ifelse(MYC_group == "High" & mTOR_group == "Low", "MYC-High",
                   ifelse(mTOR_group == "High" & MYC_group == "Low", "mTOR-High",
                   ifelse(MYC_group == "High", "Both high", "Both low"))))
cats <- c("MYC-High", "mTOR-High", "Both high", "Both low")

genes <- c("MTOR", "RPTOR", "RICTOR", "MLST8", "RHEB")
screens <- list(CRISPR = list(X = attr(results$crispr, "X"), unit = "Chronos gene effect"),
                RNAi   = list(X = attr(results$rnai, "X"),   unit = "DEMETER2 score"))

colA <- "#eb6834"; colB <- "#4a3aa7"; col_mix <- "#a3a29d"
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
cat_col <- c("MYC-High" = colA, "mTOR-High" = colB, "Both high" = col_mix, "Both low" = col_mix)

get <- function(X, g) { j <- grep(paste0("^", g, " \\("), colnames(X)); if (length(j)) X[, j[1]] else NULL }

rows <- list()
for (s in names(screens)) for (g in genes) {
  v <- get(screens[[s]]$X, g); if (is.null(v)) next
  k <- cl$cat[match(names(v), cl$cell_line)]; ok <- !is.na(v) & !is.na(k)
  a <- v[ok & k == "MYC-High"]; b <- v[ok & k == "mTOR-High"]
  rows[[length(rows) + 1]] <- data.frame(screen = s, gene = g, n_MYC = length(a), n_mTOR = length(b),
    mean_MYC = mean(a), mean_mTOR = mean(b), diff_MYC_minus_mTOR = mean(a) - mean(b),
    p_welch = if (length(a) > 2 && length(b) > 2) t.test(a, b)$p.value else NA,
    mean_both_high = mean(v[ok & k == "Both high"]), mean_both_low = mean(v[ok & k == "Both low"]))
}
tab <- do.call(rbind, rows)
write.table(tab, file.path(res_dir, "MTOR_gene_effect_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(data.frame(tab[, 1:4], round(tab[, 5:7], 2), p = signif(tab$p_welch, 2),
                 both_high = round(tab$mean_both_high, 2), both_low = round(tab$mean_both_low, 2)), row.names = FALSE)

# Figure: MTOR CRISPR only, MYC-High vs mTOR-High (table above keeps all genes / screens / categories)
tint <- function(col, f = 0.75) colorRampPalette(c(col, "#ffffff"))(100)[round(f * 100)]
draw <- function() {
  par(mar = c(3.6, 4.6, 4.6, 1), fg = ink2, bg = bg)
  v <- get(screens$CRISPR$X, "MTOR")
  k <- factor(cl$cat[match(names(v), cl$cell_line)], levels = c("MYC-High", "mTOR-High"))
  ok <- !is.na(v) & !is.na(k); v <- v[ok]; k <- k[ok]; x <- as.integer(k)
  yr <- range(c(v, -0.8)); yr <- yr + c(-0.06, 0.08) * diff(yr)
  plot(NA, xlim = c(0.4, 2.6), ylim = yr, axes = FALSE, xlab = "", ylab = "")
  abline(h = pretty(yr), col = grid_col, lwd = 0.8); abline(h = -1, col = ink2, lwd = 0.6, lty = 3)
  boxplot(v ~ k, at = 1:2, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.5, lwd = 1, medlwd = 2,
          col = c(tint(colA), tint(colB)), border = c(colA, colB), whisklty = 1, staplewex = 0.4)
  set.seed(1)
  points(x + runif(length(v), -0.12, 0.12), v, pch = 21, cex = 1.2, lwd = 0.5, col = bg, bg = c(colA, colB)[x])
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.9, tcl = -0.25)
  axis(1, at = 1:2, labels = sprintf("%s\nn = %d", levels(k), tabulate(x, 2)), col = NA, col.axis = ink,
       cex.axis = 0.9, padj = 0.6)
  r <- tab[tab$screen == "CRISPR" & tab$gene == "MTOR", ]
  mtext("MTOR dependency (CRISPR)", side = 3, line = 2.6, adj = 0, font = 2, cex = 1.1, col = ink)
  mtext(sprintf("Mean %.2f vs %.2f · Welch p = %.2g", r$mean_MYC, r$mean_mTOR, r$p_welch),
        side = 3, line = 1.4, adj = 0, cex = 0.8, col = ink2)
  mtext("More negative = stronger dependency · dotted = −1 (median common-essential gene)", side = 3, line = 0.4,
        adj = 0, cex = 0.68, col = ink2)
  mtext("Chronos gene effect", side = 2, line = 2.9, cex = 0.85, col = ink2)
}
fb <- file.path(fig_dir, "MTOR_gene_effect_by_group")
png(paste0(fb, ".png"), width = 1300, height = 1500, res = 220, type = "cairo", bg = bg); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 1300 / 220, height = 1500 / 220, bg = bg); draw(); dev.off()
