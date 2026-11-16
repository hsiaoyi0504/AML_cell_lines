# OCI-AML4 vs the other KMT2A-rearranged AML lines (n=9; KMT2A status from the
# classification file). OCI-AML4 is the only KMT2A-r line in the MYC-High/mTOR-Low group.
# One line vs nine: descriptive outlier analysis (z vs the other nine, outside their range),
# no formal test. Data with OCI-AML4: DepMap RNA and CRISPR only.
#
# Usage: Rscript scripts/OCIAML4_vs_KMT2Ar.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
out_dir <- file.path(res_dir, "OCIAML4"); dir.create(out_dir, showWarnings = FALSE)
source("scripts/plot_utils.R"); source("scripts/enrichment_utils.R")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
results <- readRDS(file.path(res_dir, "all_results.rds"))
sym <- function(x) sub(" \\([0-9]+\\)$", "", x)
target <- "OCIAML4"
kr <- cl$cell_line[cl$KMT2A_rearranged]; others <- setdiff(kr, target)
partner <- setNames(cl$KMT2A_partner, cl$cell_line)

# outlier stats of one line vs a reference set, per feature (columns of X)
outlier <- function(X, ref, sd_floor) {
  R <- X[intersect(ref, rownames(X)), , drop = FALSE]; x <- X[target, ]
  m <- colMeans(R, na.rm = TRUE); s <- pmax(apply(R, 2, sd, na.rm = TRUE), sd_floor)
  data.frame(feature = sym(colnames(X)), value_OCIAML4 = x, mean_otherKMT2Ar = m,
             min_other = apply(R, 2, min, na.rm = TRUE), max_other = apply(R, 2, max, na.rm = TRUE),
             diff = x - m, z = (x - m) / s,
             outside_range = x > apply(R, 2, max, na.rm = TRUE) | x < apply(R, 2, min, na.rm = TRUE),
             n_ref = nrow(R), row.names = NULL)
}

# ---------------------------------------------------------------- RNA -------
E <- attr(results$expression, "X"); colnames(E) <- sym(colnames(E))
rna <- outlier(E, others, sd_floor = 0.25)          # floor stops near-constant genes inflating z
mh <- cl$cell_line[cl$MYC_group == "High" & cl$mTOR_group == "Low" & !cl$KMT2A_rearranged]
th <- cl$cell_line[cl$mTOR_group == "High" & cl$MYC_group == "Low" & !cl$KMT2A_rearranged]
rna$mean_MYCHigh_nonKMT2A  <- colMeans(E[mh, ])[rna$feature]
rna$mean_mTORHigh_nonKMT2A <- colMeans(E[th, ])[rna$feature]
rna$outlier <- rna$outside_range & abs(rna$diff) >= 2 & abs(rna$z) >= 3
rna <- rna[order(-abs(rna$z) * rna$outlier, -abs(rna$diff)), ]
write.table(rna, file.path(out_dir, "OCIAML4_vs_KMT2Ar_RNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
up <- rna$feature[rna$outlier & rna$diff > 0]; dn <- rna$feature[rna$outlier & rna$diff < 0]

# global similarity on the 2000 most variable genes
vg <- names(sort(apply(E, 2, var), decreasing = TRUE))[1:2000]
cc <- cor(t(E[, vg]), method = "spearman")[target, ]; cc <- sort(cc[names(cc) != target], decreasing = TRUE)
grp_of <- function(l) ifelse(l %in% others, "KMT2A-r", ifelse(l %in% mh, "MYC-High (non-KMT2A)",
                      ifelse(l %in% th, "mTOR-High (non-KMT2A)", "mixed")))
sim <- data.frame(cell_line = names(cc), spearman = round(cc, 3), group = grp_of(names(cc)))
write.table(sim, file.path(out_dir, "OCIAML4_similarity_to_lines.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# marker panel: KMT2A program, monocytic, granulocytic, progenitor, erythroid/megakaryocytic
markers <- list(
  "KMT2A program" = c("HOXA9", "HOXA10", "MEIS1", "PBX3", "MEF2C", "MECOM", "MYB", "ZEB2", "JMJD1C"),
  "Monocytic"     = c("CD14", "ITGAM", "CSF1R", "LYZ", "CD68", "MAFB", "IRF8", "KLF4"),
  "Granulocytic"  = c("MPO", "ELANE", "AZU1", "PRTN3", "CEBPE"),
  "Progenitor"    = c("CD34", "KIT", "PROM1", "MECOM", "HLF"),
  "Erythroid/MK"  = c("GATA1", "TAL1", "KLF1", "GFI1B", "ITGA2B", "GP1BB", "EPOR"),
  "MYC / mTOR"    = c("MYC", "MTOR", "RPTOR", "SLC7A5"))
mk <- do.call(rbind, lapply(names(markers), function(n) {
  g <- intersect(markers[[n]], rna$feature); cbind(panel = n, rna[match(g, rna$feature), ])
}))
mk <- mk[!duplicated(paste(mk$panel, mk$feature)), ]
write.table(mk, file.path(out_dir, "OCIAML4_marker_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# enrichment of OCI-AML4-specific genes (universe = expressed genes)
en <- list()
for (lb in names(enrich_libs)) {
  sets <- read_gmt(file.path("data/gene_sets", paste0(enrich_libs[[lb]], ".gmt")))
  en[[lb]] <- rbind(ora(sets, up, rna$feature, lb, "higher_in_OCIAML4"), ora(sets, dn, rna$feature, lb, "lower_in_OCIAML4"))
}
en <- do.call(rbind, en)
write.table(en[order(en$library, en$direction, en$p), ], file.path(out_dir, "OCIAML4_specific_genes_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- CRISPR ----
C <- attr(results$crispr, "X")
cr <- outlier(C, others, sd_floor = 0.15)
cr$outlier <- cr$outside_range & abs(cr$diff) >= 0.5 & abs(cr$z) >= 3; cr$outlier[is.na(cr$outlier)] <- FALSE
cr <- cr[order(-abs(cr$z) * cr$outlier, -abs(cr$diff)), ]
write.table(cr, file.path(out_dir, "OCIAML4_vs_KMT2Ar_CRISPR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
kmt2a_deps <- c("KMT2A", "MEN1", "PSIP1", "DOT1L", "MLLT1", "MLLT3", "MLLT10", "AFF4", "ELL", "KAT6A", "KAT7",
                "BRD4", "HOXA9", "HOXA10", "MEIS1", "PBX3", "MYB", "ZEB2", "MEF2C", "MEF2D", "IRF8", "CEBPA",
                "SPI1", "LMO2", "GATA1", "TAL1", "BCL2", "BCL2L1", "MCL1", "MTOR", "RPTOR", "MYC")
crk <- cr[cr$feature %in% kmt2a_deps, ]; crk <- crk[match(intersect(kmt2a_deps, crk$feature), crk$feature), ]
write.table(crk, file.path(out_dir, "OCIAML4_CRISPR_key_dependencies.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- summary ---
r2 <- function(d, cols) { d <- d[, cols]; for (c in cols) if (is.numeric(d[[c]])) d[[c]] <- round(d[[c]], 2); d }
sink(file.path(out_dir, "OCIAML4_summary.txt"))
cat("OCI-AML4 (KMT2A-", partner[target], ") vs other KMT2A-r lines: ", paste0(others, " (", partner[others], ")", collapse = ", "), "\n", sep = "")
cat("Note: MOLM13/MOLM14 and MONOMAC1/MONOMAC6 are pairs from the same patients (~7 independent donors).\n")
cat("Datasets with OCI-AML4: RNA (9 reference lines), CRISPR (", sum(others %in% rownames(C)), " reference lines)\n\n", sep = "")
cat("== Scores\n"); print(cl[cl$KMT2A_rearranged, c("cell_line", "KMT2A_partner", "MYC_score", "MYC_group", "mTOR_score", "mTOR_group", "OncotreeSubtype")], row.names = FALSE)
cat("\n== Global similarity (Spearman, 2000 most variable genes)\n")
cat(sprintf("  mean rho to other KMT2A-r: %.3f | MYC-High non-KMT2A: %.3f | mTOR-High non-KMT2A: %.3f\n",
            mean(cc[others]), mean(cc[intersect(mh, names(cc))]), mean(cc[intersect(th, names(cc))])))
cat("  10 most similar lines:\n"); print(head(sim, 10), row.names = FALSE)
cat(sprintf("\n== RNA outlier genes (outside the range of the other 9, |log2 diff| >= 2, |z| >= 3): higher %d | lower %d (of %d genes)\n",
            length(up), length(dn), nrow(rna)))
cols <- c("feature", "value_OCIAML4", "mean_otherKMT2Ar", "min_other", "max_other", "diff", "z", "mean_MYCHigh_nonKMT2A", "mean_mTORHigh_nonKMT2A")
cat("  Top higher in OCI-AML4:\n"); print(head(r2(rna[rna$outlier & rna$diff > 0, ], cols), 25), row.names = FALSE)
cat("  Top lower in OCI-AML4:\n");  print(head(r2(rna[rna$outlier & rna$diff < 0, ], cols), 25), row.names = FALSE)
cat("\n== Marker genes\n"); print(r2(mk, c("panel", cols, "outside_range")), row.names = FALSE)
cat("\n== Enrichment of OCI-AML4-specific genes (top 6 per library and direction, FDR shown)\n")
for (d in c("higher_in_OCIAML4", "lower_in_OCIAML4")) for (lb in names(enrich_libs)) {
  x <- head(en[en$library == lb & en$direction == d, ], 6)
  cat(sprintf("  [%s | %s] ", d, lb), paste(sprintf("%s (%d/%d, FDR %.2g)", substr(x$term, 1, 45), x$overlap, x$set_size, x$fdr), collapse = "; "), "\n")
}
cat(sprintf("\n== CRISPR outliers (outside range of other KMT2A-r, |diff| >= 0.5, |z| >= 3): more dependent in OCI-AML4 %d | less %d\n",
            sum(cr$outlier & cr$diff < 0, na.rm = TRUE), sum(cr$outlier & cr$diff > 0, na.rm = TRUE)))
ccols <- c("feature", "value_OCIAML4", "mean_otherKMT2Ar", "min_other", "max_other", "diff", "z", "n_ref")
cat("  More dependent in OCI-AML4:\n"); print(head(r2(cr[cr$outlier & cr$diff < 0, ], ccols), 20), row.names = FALSE)
cat("  Less dependent in OCI-AML4:\n"); print(head(r2(cr[cr$outlier & cr$diff > 0, ], ccols), 20), row.names = FALSE)
cat("  KMT2A-complex / lineage / survival dependencies:\n"); print(r2(crk, c(ccols, "outside_range")), row.names = FALSE)
sink()

# ---------------------------------------------------------------- figures ---
bg <- "#fcfcfb"; col_t <- "#eb6834"; col_o <- "#4a3aa7"
pal <- colorRampPalette(c("#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a8a6", "#e34948", "#9c2423"))(101)
# (1) heatmap: top OCI-AML4-specific genes across the 10 KMT2A-r lines + group means for context
sel <- c(head(up, 25), head(dn, 25))
lines10 <- c(target, others[order(-cc[others])])            # others ordered by similarity to OCI-AML4
ctx <- rbind(`MYC-High mean (non-KMT2A)` = colMeans(E[mh, sel]), `mTOR-High mean (non-KMT2A)` = colMeans(E[th, sel]))
M <- rbind(E[lines10, sel], ctx)
Z <- sweep(M, 2, colMeans(E[c(target, others), sel])); Z <- sweep(Z, 2, apply(E[c(target, others), sel], 2, sd), "/")
Z[Z > 2.5] <- 2.5; Z[Z < -2.5] <- -2.5
draw_hm <- function() {
  par(oma = c(0, 1, 5.6, 1)); layout(matrix(1:2, 1), widths = c(1, 0.12))
  par(mar = c(13.5, 9, 0.5, 1), fg = ink2)
  nc <- nrow(Z); nr <- ncol(Z)
  image(seq_len(nc), seq_len(nr), Z[, rev(seq_len(nr))], useRaster = TRUE, col = pal, zlim = c(-2.5, 2.5),
        axes = FALSE, xlab = "", ylab = "")
  abline(v = c(1.5, length(lines10) + 0.5), col = bg, lwd = 4); abline(h = length(dn) + 0.5, col = bg, lwd = 4)
  gl <- rev(sel); d <- rna$diff[match(gl, rna$feature)]
  axis(2, at = seq_len(nr), labels = sprintf("%s  %+.1f", gl, d), las = 1, tick = FALSE, cex.axis = 0.85, col.axis = ink, line = -0.6)
  xl <- c(paste0(lines10, " (", partner[lines10], ")"), rownames(ctx))
  text(seq_len(nc), par("usr")[3] - 0.6, xl, srt = 90, adj = 1, xpd = NA, cex = 0.85,
       col = c(col_t, rep(col_o, length(others)), ink2, ink2), font = c(2, rep(1, length(others)), 3, 3))
  mtext(sprintf("higher in OCI-AML4 (%d)", length(head(up, 25))), side = 2, line = 7.6, at = length(dn) + length(head(up, 25)) / 2, cex = 0.9, col = ink2)
  mtext(sprintf("lower in OCI-AML4 (%d)", length(head(dn, 25))), side = 2, line = 7.6, at = length(dn) / 2, cex = 0.9, col = ink2)
  mtext("OCI-AML4-specific genes vs other KMT2A-rearranged lines (RNA)", side = 3, line = 3.4, adj = 0, outer = TRUE, font = 2, cex = 1.3, col = ink)
  mtext("Genes outside the range of the other 9 KMT2A-r lines, |log2 diff| \u2265 2, |z| \u2265 3 \u00b7 top 25 each way (number = OCI-AML4 \u2212 mean of other KMT2A-r)",
        side = 3, line = 1.8, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
  mtext("Row z-score across the 10 KMT2A-r lines (capped \u00b12.5); the two group-mean columns use the same scale for context", side = 3, line = 0.6, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
  par(mar = c(13.5, 1.2, 8, 3.8))
  image(1, seq(-2.5, 2.5, length.out = 101), t(seq(-2.5, 2.5, length.out = 101)), useRaster = TRUE, col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = -2:2, las = 1, cex.axis = 0.9, col = NA, col.ticks = ink2, col.axis = ink2, tcl = -0.2)
  mtext("z-score", side = 3, line = 0.5, cex = 0.9, col = ink2)
}
fb <- file.path(fig_dir, "OCIAML4_vs_KMT2Ar_heatmap")
png(paste0(fb, ".png"), width = 2600, height = 3000, res = 220, type = "cairo", bg = bg); draw_hm(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2600 / 220, height = 3000 / 220, bg = bg); draw_hm(); dev.off()

# (2) markers: OCI-AML4 vs other KMT2A-r, with group means for reference
draw_mk <- function() {
  par(oma = c(0, 1, 3.8, 1), mar = c(4.2, 9, 0.5, 1), fg = ink2)
  d <- mk; n <- nrow(d); y <- rev(seq_len(n))
  xr <- range(c(d$min_other, d$max_other, d$value_OCIAML4, d$mean_MYCHigh_nonKMT2A, d$mean_mTORHigh_nonKMT2A), na.rm = TRUE)
  plot(NA, xlim = xr + c(-0.3, 0.3), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "", yaxs = "i")
  abline(v = pretty(xr), col = grid_col, lwd = 0.8)
  pb <- which(c(TRUE, d$panel[-1] != d$panel[-n])); abline(h = y[pb] + 0.5, col = ink2, lwd = 0.6)
  segments(d$min_other, y, d$max_other, y, col = "#c9c3ee", lwd = 6, lend = 1)       # range of other KMT2A-r
  points(d$mean_MYCHigh_nonKMT2A, y, pch = 124, cex = 1.1, col = "#e34948")
  points(d$mean_mTORHigh_nonKMT2A, y, pch = 124, cex = 1.1, col = "#2a78d6")
  points(d$mean_otherKMT2Ar, y, pch = 21, bg = col_o, col = bg, cex = 1.1)
  points(d$value_OCIAML4, y, pch = 23, bg = col_t, col = ink, cex = 1.4, lwd = 0.8)
  axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.85, tcl = -0.25)
  axis(2, at = y, labels = d$feature, las = 1, tick = FALSE, cex.axis = 0.82, col.axis = ink,
       font = ifelse(d$outside_range, 2, 1))
  text(xr[1] - 0.3, y[pb], d$panel[pb], adj = c(0, -0.6), cex = 0.75, col = ink2, font = 3, xpd = NA)
  mtext("log2(TPM+1)", side = 1, line = 2.4, cex = 0.85, col = ink2)
  mtext("OCI-AML4 marker profile vs other KMT2A-rearranged lines", side = 3, line = 2, adj = 0, outer = TRUE, font = 2, cex = 1.2, col = ink)
  mtext("bar = range of the other 9 KMT2A-r · purple dot = their mean · orange diamond = OCI-AML4 (bold gene = outside range) · ticks: red = MYC-High mean, blue = mTOR-High mean (non-KMT2A)",
        side = 3, line = 0.6, adj = 0, outer = TRUE, cex = 0.72, col = ink2)
}
fb <- file.path(fig_dir, "OCIAML4_marker_genes")
png(paste0(fb, ".png"), width = 2600, height = 2900, res = 220, type = "cairo", bg = bg); draw_mk(); dev.off()
cairo_pdf(paste0(fb, ".pdf"), width = 2600 / 220, height = 2900 / 220, bg = bg); draw_mk(); dev.off()

# (3) one gene across all 56 AML lines, ranked: OCI-AML4 vs other KMT2A-r vs other AML
plot_gene_rank <- function(gene, label = gene) {
  v <- sort(E[, gene], decreasing = TRUE)
  kind <- ifelse(names(v) == target, "OCI-AML4", ifelse(names(v) %in% others, "other KMT2A-r", "other AML"))
  fill <- c("OCI-AML4" = col_t, "other KMT2A-r" = col_o, "other AML" = "#c9c8c3")[kind]
  write.table(data.frame(cell_line = names(v), log2TPM1 = round(v, 3), group = kind,
                         KMT2A_partner = partner[names(v)], mTOR_group = cl$mTOR_group[match(names(v), cl$cell_line)]),
              file.path(out_dir, paste0(gene, "_expression_all_AML.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  rk <- which(names(v) == target)
  draw <- function() {
    n <- length(v)
    par(oma = c(0, 1, 4.6, 1), mar = c(7.2, 4.6, 0.5, 1), fg = ink2)
    yl <- c(0, max(v) * 1.08)
    plot(NA, xlim = c(0.4, n + 0.6), ylim = yl, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    abline(h = pretty(yl), col = grid_col, lwd = 0.8)
    rect(seq_len(n) - 0.4, 0, seq_len(n) + 0.4, v, col = fill, border = NA)
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.9, tcl = -0.25)
    text(seq_len(n), -yl[2] * 0.015, ifelse(kind == "other AML", names(v), paste0(names(v), " (", partner[names(v)], ")")),
         srt = 90, adj = 1, xpd = NA, cex = 0.62, col = ifelse(kind == "other AML", ink2, fill), font = ifelse(kind == "other AML", 1, 2))
    text(rk, v[rk], sprintf("%.2f", v[rk]), pos = 3, cex = 0.8, col = ink, font = 2, offset = 0.3)
    mtext("log2(TPM+1)", side = 2, line = 3.2, cex = 0.85, col = ink2)
    legend("topright", bty = "n", pch = 15, pt.cex = 1.6, cex = 0.85, text.col = ink,
           col = c(col_t, col_o, "#c9c8c3"), legend = c("OCI-AML4 (KMT2A-r)", "other KMT2A-r", "other AML"))
    mtext(sprintf("%s expression: OCI-AML4 vs other AML lines", label), side = 3, line = 3, adj = 0, outer = TRUE,
          font = 2, cex = 1.2, col = ink)
    ok <- names(v) %in% others; oa <- kind == "other AML"
    mtext(sprintf("DepMap 26Q3 RNA (gene level) · OCI-AML4 = %.2f, rank %d of %d · other KMT2A-r median %.2f (max %.2f) · other AML median %.2f",
                  v[rk], rk, n, median(v[ok]), max(v[ok]), median(v[oa])),
          side = 3, line = 1.6, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
    mtext(sprintf("KMT2A fusion partner in brackets · %s", if (gene == "MECOM") "MECOM gene-level TPM includes the EVI1 and MDS1-EVI1 isoforms" else ""),
          side = 3, line = 0.4, adj = 0, outer = TRUE, cex = 0.8, col = ink2)
  }
  fb <- file.path(fig_dir, paste0("OCIAML4_", gene, "_expression_ranked"))
  png(paste0(fb, ".png"), width = 3200, height = 1500, res = 220, type = "cairo", bg = bg); draw(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 3200 / 220, height = 1500 / 220, bg = bg); draw(); dev.off()
}
plot_gene_rank("MECOM", "EVI1 (MECOM)")

cat(readLines(file.path(out_dir, "OCIAML4_summary.txt")), sep = "\n")
