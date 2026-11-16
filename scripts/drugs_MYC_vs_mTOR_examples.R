# Boxplots of selected drugs, MYC-High vs mTOR-High (AUC; lower = more sensitive). p values from
# scripts/drugs_MYC_vs_mTOR_screen.R: raw Welch, and adjusted for each line's general drug sensitivity.
#
# Run from the project root: Rscript scripts/drugs_MYC_vs_mTOR_examples.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp <- ifelse(cl$MYC_group == "High" & cl$mTOR_group == "Low", "MYC-High",
          ifelse(cl$mTOR_group == "High" & cl$MYC_group == "Low", "mTOR-High", NA))
st <- read.delim(file.path(res_dir, "drugs_MYC_High_vs_mTOR_High_CTRP_GDSC2.tsv"))

colA <- "#eb6834"; colB <- "#4a3aa7"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e8e7e3"; bg <- "#fcfcfb"
tint <- function(col, f = 0.75) colorRampPalette(c(col, "#ffffff"))(100)[round(f * 100)]

# drug, screen, short label for what it hits
ex <- list(
  c("NAVITOCLAX", "CTRP", "BCL-XL/BCL2"), c("BORTEZOMIB", "CTRP", "proteasome (PSMB5)"),
  c("AZD1480", "CTRP", "JAK1/2"), c("RUXOLITINIB", "CTRP", "JAK1/2"),
  c("NVP-BSK805", "CTRP", "JAK2"), c("VORINOSTAT", "CTRP", "HDAC"),
  c("VORINOSTAT", "GDSC2", "HDAC"), c("TACEDINALINE", "CTRP", "HDAC1"),
  c("SNS-032", "CTRP", "CDK2/7/9"), c("PHA 793887", "CTRP", "CDK"),
  c("DINACICLIB", "CTRP", "CDK1/2/5/9"), c("DOXORUBICIN", "CTRP", "TOP2"),
  c("TANESPIMYCIN", "CTRP", "HSP90"), c("TRAMETINIB", "CTRP", "MEK1/2"),
  c("TRAMETINIB", "GDSC2", "MEK1/2"), c("FLUOROURACIL", "GDSC2", "TYMS"))

draw <- function() {
  par(mfrow = c(4, 4), mar = c(2.6, 3.8, 4, 0.8), oma = c(0.5, 0.5, 4.2, 0.5), fg = ink2, bg = bg)
  set.seed(1)
  for (e in ex) {
    r <- results[[paste0("drug_", e[2])]]; i <- which(r$CompoundName == e[1])[1]
    v <- attr(r, "X")[, r$feature[i]]; g <- factor(cl$grp[match(names(v), cl$cell_line)], levels = c("MYC-High", "mTOR-High"))
    ok <- !is.na(v) & !is.na(g); v <- v[ok]; g <- g[ok]; x <- as.integer(g)
    s <- st[st$drug == e[1] & st$screen == e[2], ]
    yr <- range(v); yr <- yr + c(-0.06, 0.1) * diff(yr)
    plot(NA, xlim = c(0.4, 2.6), ylim = yr, axes = FALSE, xlab = "", ylab = "")
    abline(h = pretty(yr), col = grid_col, lwd = 0.8)
    boxplot(v ~ g, at = 1:2, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.55, lwd = 1, medlwd = 2,
            col = c(tint(colA), tint(colB)), border = c(colA, colB), whisklty = 1, staplewex = 0.4)
    points(x + runif(length(v), -0.13, 0.13), v, pch = 21, cex = 1.05, lwd = 0.5, col = bg, bg = c(colA, colB)[x])
    axis(2, col = ink, col.ticks = ink, col.axis = ink, las = 1, cex.axis = 0.8, tcl = -0.3, lwd = 0.8)
    axis(1, at = 1:2, labels = sprintf("%s\nn=%d", c("MYC-High", "mTOR-High"), tabulate(x, 2)), col = NA, col.axis = ink,
         cex.axis = 0.78, padj = 0.6, mgp = c(2, 0.4, 0))
    mtext("AUC", side = 2, line = 2.5, cex = 0.62, col = ink2)
    sel <- if (s$adj_diff_MYC_minus_mTOR < 0) colA else colB
    mtext(sprintf("%s · %s", if (grepl("[0-9]", e[1])) e[1] else tools::toTitleCase(tolower(e[1])), e[2]), side = 3, line = 2.2, adj = 0, font = 2, cex = 0.78, col = ink)
    mtext(e[3], side = 3, line = 1.25, adj = 0, cex = 0.66, col = sel, font = 2)
    mtext(sprintf("p = %.2g · adjusted p = %.2g", s$p_welch, s$p_adj), side = 3, line = 0.3, adj = 0, cex = 0.62,
          col = if (s$p_adj < 0.05) ink else ink2)
  }
  mtext("Drugs with group-selective response: MYC-High vs mTOR-High", outer = TRUE, side = 3, line = 2.6, adj = 0.005,
        font = 2, cex = 1.15, col = ink)
  mtext("Rows 1–2: more active in MYC-High (orange subtitle) · rows 3–4: more active in mTOR-High (purple) · lower AUC = more sensitive · adjusted p = ANCOVA, AUC ~ group + each line's general drug sensitivity",
        outer = TRUE, side = 3, line = 1.2, adj = 0.005, cex = 0.74, col = ink2)
  mtext("mTOR-High lines are more sensitive to most drugs (73%), so the adjusted p shows drug-specific selectivity",
        outer = TRUE, side = 3, line = 0.1, adj = 0.005, cex = 0.74, col = ink2, font = 3)
}
fb <- file.path(fig_dir, "drugs_MYC_vs_mTOR_selective_examples")
png(paste0(fb, ".png"), width = 3000, height = 3300, res = 220, type = "cairo", bg = bg); draw(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 3000 / 220, height = 3300 / 220, bg = bg); draw(); dev.off()
