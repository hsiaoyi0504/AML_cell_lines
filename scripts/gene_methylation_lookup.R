# All methylation evidence for one gene, mTOR-High vs mTOR-Low (RNA groups):
#   EPIC GSE270494 (every probe on the gene), RRBS TSS-1kb and CpG clusters,
#   plus expression, methylation-expression correlation and CRISPR/RNAi effect.
# Needs the outputs of compare_mTOR_RNA_groups.R and methylation_EPIC_GSE270494.R.
#
# Run from the project root: Rscript scripts/gene_methylation_lookup.R  (set the gene below)

gene <- "ZMIZ1"                     # gene to look up
res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
out_dir <- file.path(res_dir, "gene_lookup"); dir.create(out_dir, showWarnings = FALSE)
source("scripts/plot_utils.R")
cl      <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
results <- readRDS(file.path(res_dir, "all_results.rds"))
sym <- function(x) sub(" \\([0-9]+\\)$", "", x)

# ------------------------------------------------------- EPIC probes -------
man <- read.delim(gzfile("data/annotation/EPIC.hg38.manifest.gencode.v36.tsv.gz"), colClasses = "character")
man <- man[vapply(strsplit(man$genesUniq, ";"), function(g) gene %in% g, TRUE), ]
# distance to this gene's own TSSs only (geneNames and distToTSS are parallel lists)
d_own <- mapply(function(g, d) { d <- suppressWarnings(as.numeric(strsplit(d, ";")[[1]]))
                                 d[strsplit(g, ";")[[1]] == gene] },
                man$geneNames, man$distToTSS, SIMPLIFY = FALSE)
man$dist_to_TSS <- vapply(d_own, function(d) if (length(d)) d[which.min(abs(d))] else NA_real_, 0)
man$region <- ifelse(!is.na(man$dist_to_TSS) & man$dist_to_TSS >= -1500 & man$dist_to_TSS <= 500,
                     "promoter", ifelse(!is.na(man$dist_to_TSS) & man$dist_to_TSS < -1500, "upstream", "gene body"))

ep_all <- read.delim(gzfile(file.path(res_dir, "diff_methylation_EPIC_GSE270494.tsv.gz")))
ep <- merge(man[, c("probeID", "CpG_chrm", "CpG_beg", "dist_to_TSS", "region", "CGIposition")],
            ep_all[, c("probe", "n_High", "n_Low", "meanBeta_High", "meanBeta_Low", "delta_beta",
                       "diff_M", "p_welch", "fdr_welch")],
            by.x = "probeID", by.y = "probe", all.x = TRUE)
ep$tested <- !is.na(ep$p_welch)
ep <- ep[order(as.numeric(ep$CpG_beg)), ]
# within-gene correction: BH over this gene's tested probes
ep$fdr_within_gene <- NA; ep$fdr_within_gene[ep$tested] <- p.adjust(ep$p_welch[ep$tested], "BH")

# per-line betas for this gene's probes, straight from the array file
beta_file <- "data/cell_lines/GSE270494/GSE270494_Noguera-Castells_Average_Beta_Homo_Sapiens.csv.gz"
hdr  <- scan(pipe(paste("gzip -dc", shQuote(beta_file), "| head -1")), what = "", sep = ",", quiet = TRUE)
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
samp <- hdr[-1][!grepl("Detection Pval$", hdr[-1])]; samp <- samp[norm(samp) %in% cl$cell_line]
pf <- tempfile(); writeLines(paste0("\"", ep$probeID, "\""), pf)
tmp <- tempfile()     # awk hash lookup on the ID column (BSD grep -F -f is very slow on long lines)
system(paste("gzip -dc", shQuote(beta_file), "| awk -F, 'NR==FNR{a[$1]; next} ($1 in a)'", pf, "- >", tmp))
raw <- read.csv(tmp, header = FALSE, col.names = hdr, check.names = FALSE)
Bm  <- as.matrix(raw[, samp]); P <- as.matrix(raw[, paste(samp, "Detection Pval")])
Bm[P > 0.05] <- NA; rownames(Bm) <- raw[[1]]; colnames(Bm) <- norm(samp)

# ------------------------------------------------------- RRBS --------------
rr_tss <- results$meth[sub("\\|.*$", "", results$meth$feature) == gene, ]
rr_cpg <- results$meth_cpg[sub("\\|.*$", "", results$meth_cpg$feature) == gene, ]
keep_cols <- c("feature", "n_High", "n_Low", "mean_High", "mean_Low", "diff_High_minus_Low", "p_welch", "fdr_welch")

# ------------------------------------------- expression / dependencies ----
E <- attr(results$expression, "X"); colnames(E) <- sym(colnames(E))
ex <- results$expression[sym(results$expression$feature) == gene, ]
cr <- results$crispr[sym(results$crispr$feature) == gene, ]
rn <- results$rnai[sym(results$rnai$feature) == gene, ]
scor <- function(a, b, min_n = 8) {
  o <- !is.na(a) & !is.na(b); if (sum(o) < min_n) return(c(NA, NA))
  t <- suppressWarnings(cor.test(a[o], b[o], method = "spearman", exact = FALSE)); c(t$estimate, t$p.value)
}
# methylation vs RNA cross-check, per probe:
#   all lines; within each group (lineage-independent); expression ~ M + group
#   (does methylation explain expression beyond the group difference?); vs RNA mTOR score
cross <- function(b) {
  m <- log2(pmin(pmax(b, 0.001), 0.999) / (1 - pmin(pmax(b, 0.001), 0.999)))
  ln <- names(b); g <- cl$mTOR_group[match(ln, cl$cell_line)]
  e  <- E[match(ln, rownames(E)), gene]; sc <- cl$mTOR_score[match(ln, cl$cell_line)]
  a  <- scor(b, e); h <- scor(b[g == "High"], e[g == "High"], 6); l <- scor(b[g == "Low"], e[g == "Low"], 6)
  s  <- scor(b, sc)
  o  <- !is.na(m) & !is.na(e)
  adj <- if (sum(o) >= 10 && length(unique(g[o])) == 2) {
    f <- summary(lm(e[o] ~ m[o] + factor(g[o])))$coefficients; c(f[2, 1], f[2, 4])
  } else c(NA, NA)
  c(rho_all = a[1], p_all = a[2], rho_within_High = h[1], p_within_High = h[2],
    rho_within_Low = l[1], p_within_Low = l[2], slope_adj_group = adj[1], p_adj_group = adj[2],
    rho_vs_mTOR_score = s[1], p_vs_mTOR_score = s[2], n = sum(o))
}
has_expr <- gene %in% colnames(E)
if (has_expr) {
  cx <- t(vapply(ep$probeID, function(p) if (p %in% rownames(Bm) && isTRUE(ep$tested[ep$probeID == p])) cross(Bm[p, ]) else rep(NA_real_, 11), numeric(11)))
  colnames(cx) <- c("rho_meth_vs_expr", "p_meth_vs_expr", "rho_expr_within_High", "p_expr_within_High",
                    "rho_expr_within_Low", "p_expr_within_Low", "slope_expr_on_M_adj_group", "p_expr_adj_group",
                    "rho_meth_vs_mTOR_score", "p_meth_vs_mTOR_score", "n_lines_expr")
  ep <- cbind(ep, cx)
  tt <- !is.na(ep$p_meth_vs_expr)
  ep$fdr_meth_vs_expr <- NA;  ep$fdr_meth_vs_expr[tt]  <- p.adjust(ep$p_meth_vs_expr[tt], "BH")
  ta <- !is.na(ep$p_expr_adj_group)
  ep$fdr_expr_adj_group <- NA; ep$fdr_expr_adj_group[ta] <- p.adjust(ep$p_expr_adj_group[ta], "BH")
}
write.table(ep, file.path(out_dir, paste0(gene, "_EPIC_probes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------- report ---
sink(file.path(out_dir, paste0(gene, "_methylation_summary.txt")))
cat("=====", gene, "=====\n\n")
cat(sprintf("EPIC (GSE270494): %d probes on %s, %d tested after QC\n", nrow(ep), gene, sum(ep$tested)))
print(table(region = ep$region, tested = ep$tested))
cat(sprintf("  tested probes: genome-wide FDR<0.05: %d | FDR<0.25: %d | nominal p<0.05: %d\n",
            sum(ep$fdr_welch < 0.05, na.rm = TRUE), sum(ep$fdr_welch < 0.25, na.rm = TRUE), sum(ep$p_welch < 0.05, na.rm = TRUE)))
s <- ep$p_welch < 0.05 & ep$tested
cat(sprintf("  of nominal p<0.05 probes: more methylated in High %d | less %d\n",
            sum(s & ep$delta_beta > 0, na.rm = TRUE), sum(s & ep$delta_beta < 0, na.rm = TRUE)))
cat("  by region (tested probes): mean delta beta, n p<0.05\n")
for (rg in unique(ep$region)) { x <- ep[ep$tested & ep$region == rg, ]
  cat(sprintf("    %-9s n=%3d  mean dBeta=%+.3f  p<0.05: %d (up %d / down %d)\n", rg, nrow(x), mean(x$delta_beta),
              sum(x$p_welch < 0.05), sum(x$p_welch < 0.05 & x$delta_beta > 0), sum(x$p_welch < 0.05 & x$delta_beta < 0))) }
cat("\n  Top EPIC probes:\n")
top <- head(ep[ep$tested, ][order(ep$p_welch[ep$tested]), ], 15)
print(data.frame(probe = top$probeID, pos = top$CpG_beg, region = top$region, dTSS = top$dist_to_TSS,
                 CGI = top$CGIposition, betaH = round(top$meanBeta_High, 2), betaL = round(top$meanBeta_Low, 2),
                 dBeta = round(top$delta_beta, 2), p = signif(top$p_welch, 2), fdr = signif(top$fdr_welch, 2),
                 rho_expr = if (!is.null(top$rho_meth_vs_expr)) round(top$rho_meth_vs_expr, 2) else NA), row.names = FALSE)
cat("\nRRBS TSS-1kb regions:\n"); print(rr_tss[, keep_cols], row.names = FALSE)
cat("\nRRBS CpG clusters:\n"); print(rr_cpg[, keep_cols], row.names = FALSE)
cat("\nExpression (log2 TPM+1):\n"); print(ex[, c("mean_High", "mean_Low", "diff_High_minus_Low", "p_welch", "fdr_welch")], row.names = FALSE)
cat("CRISPR:\n"); print(cr[, c("n_High", "n_Low", "mean_High", "mean_Low", "p_welch", "fdr_welch")], row.names = FALSE)
cat("RNAi:\n");   print(rn[, c("n_High", "n_Low", "mean_High", "mean_Low", "p_welch", "fdr_welch")], row.names = FALSE)
if (has_expr) {
  cat("\n===== Methylation vs RNA cross-check (EPIC probes, matched lines) =====\n")
  x <- ep[!is.na(ep$p_meth_vs_expr), ]
  cat(sprintf("Probes with expression link computed: %d (lines with both: %d)\n", nrow(x), max(x$n_lines_expr)))
  cat(sprintf("  all lines, FDR<0.1: %d (negative %d / positive %d) | nominal p<0.05: %d\n",
              sum(x$fdr_meth_vs_expr < 0.1), sum(x$fdr_meth_vs_expr < 0.1 & x$rho_meth_vs_expr < 0),
              sum(x$fdr_meth_vs_expr < 0.1 & x$rho_meth_vs_expr > 0), sum(x$p_meth_vs_expr < 0.05)))
  cat(sprintf("  adjusted for mTOR group, FDR<0.1: %d | nominal p<0.05: %d\n",
              sum(x$fdr_expr_adj_group < 0.1, na.rm = TRUE), sum(x$p_expr_adj_group < 0.05, na.rm = TRUE)))
  for (rg in unique(x$region)) { y <- x[x$region == rg, ]
    cat(sprintf("  %-9s n=%3d  median rho=%+.2f  p<0.05: neg %d / pos %d\n", rg, nrow(y), median(y$rho_meth_vs_expr),
                sum(y$p_meth_vs_expr < 0.05 & y$rho_meth_vs_expr < 0), sum(y$p_meth_vs_expr < 0.05 & y$rho_meth_vs_expr > 0))) }
  show <- function(d) print(data.frame(probe = d$probeID, pos = d$CpG_beg, region = d$region, dBeta = round(d$delta_beta, 2),
    p_group = signif(d$p_welch, 2), rho_all = round(d$rho_meth_vs_expr, 2), fdr_all = signif(d$fdr_meth_vs_expr, 2),
    rho_inHigh = round(d$rho_expr_within_High, 2), rho_inLow = round(d$rho_expr_within_Low, 2),
    p_adj_group = signif(d$p_expr_adj_group, 2), rho_mTORscore = round(d$rho_meth_vs_mTOR_score, 2)), row.names = FALSE)
  cat("\n  Strongest methylation-expression links (all lines):\n"); show(head(x[order(x$p_meth_vs_expr), ], 12))
  cat("\n  Strongest links after adjusting for mTOR group:\n");      show(head(x[order(x$p_expr_adj_group), ], 8))
  cat("\n  Probes that differ between groups (p<0.01):\n");          d1 <- x[which(x$p_welch < 0.01), ]; show(d1[order(d1$p_welch), ])
}
sink()

# ---------------------------------------------------------------- figure ---
bg <- "#fcfcfb"
png(file.path(fig_dir, paste0("gene_", gene, "_methylation.png")), width = 3200, height = 1700, res = 200,
    type = "cairo", bg = bg)
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(1, 1))
# (a) delta beta along the gene, sized by significance
par(mar = c(3.6, 4, 3, 1), fg = ink2)
x <- as.numeric(ep$CpG_beg) / 1e6; y <- ep$delta_beta; tst <- ep$tested
plot(NA, xlim = range(x), ylim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
abline(h = seq(-1, 1, 0.25), col = grid_col, lwd = 0.8); abline(h = 0, col = ink2, lwd = 0.8)
pr <- ep$region == "promoter"
if (any(pr)) rect(x[pr] - 0.002, -1, x[pr] + 0.002, 1, col = "#ecebe6", border = NA)  # +-2 kb per promoter probe (genes can have several TSSs)
axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
cc <- ifelse(!tst, NA, ifelse(ep$fdr_welch < 0.25, ifelse(y > 0, col_hi, col_lo),
             ifelse(ep$p_welch < 0.05, ifelse(y > 0, tint_hi, tint_lo), col_ns)))
points(x[tst], y[tst], pch = 21, bg = cc[tst], col = bg, lwd = 0.4,
       cex = 0.6 + pmin(-log10(ep$p_welch[tst]), 8) / 4)
mtext(sprintf("%s · EPIC probes along the gene (%s)", gene, ep$CpG_chrm[1]), side = 3, line = 1.3, adj = 0,
      font = 2, cex = 0.85, col = ink)
mtext(sprintf("%d tested probes · dark = genome-wide FDR<0.25, light = p<0.05 · shaded = promoter probes (-1500..+500 bp of a TSS) · size ~ -log10 p",
              sum(tst)), side = 3, line = 0.35, adj = 0, cex = 0.6, col = ink2)
mtext(sprintf("%s position (Mb, hg38)", ep$CpG_chrm[1]), side = 1, line = 2, cex = 0.62, col = ink2)
mtext("Δβ (High − Low)", side = 2, line = 2.6, cex = 0.62, col = ink2)
# (b) per-line beta at the top probe, (c) methylation vs expression
tp <- top$probeID[1]; v <- Bm[tp, ]; g <- cl$mTOR_group[match(names(v), cl$cell_line)]
par(mar = c(3.6, 4, 3, 1))
set.seed(1); o <- !is.na(v); gf <- factor(g[o], levels = c("Low", "High"))
plot(NA, xlim = c(0.4, 2.6), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
abline(h = seq(0, 1, 0.25), col = grid_col, lwd = 0.8)
axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25)
axis(1, at = 1:2, labels = c(sprintf("mTOR-Low\nn=%d", sum(gf == "Low")), sprintf("mTOR-High\nn=%d", sum(gf == "High"))),
     col = NA, col.axis = ink, cex.axis = 0.75, padj = 0.5)
points(as.integer(gf) + runif(sum(o), -0.18, 0.18), v[o], pch = 21, cex = 1.2, lwd = 0.5, col = bg,
       bg = ifelse(gf == "High", col_hi, col_lo))
for (j in 1:2) segments(j - 0.3, median(v[o][as.integer(gf) == j]), j + 0.3, lwd = 2, col = ink)
mtext(sprintf("Top probe %s (%s, %+d bp from TSS)", tp, top$region[1], as.integer(top$dist_to_TSS[1])),
      side = 3, line = 1.3, adj = 0, font = 2, cex = 0.75, col = ink)
mtext(sprintf("Δβ = %+.2f · p = %.1e · genome-wide FDR = %.2g", top$delta_beta[1], top$p_welch[1], top$fdr_welch[1]),
      side = 3, line = 0.35, adj = 0, cex = 0.6, col = ink2)
mtext("beta", side = 2, line = 2.6, cex = 0.62, col = ink2)
if (gene %in% colnames(E)) {
  e <- E[match(names(v), rownames(E)), gene]; o2 <- !is.na(v) & !is.na(e)
  plot(NA, xlim = c(0, 1), ylim = range(e[o2]) + c(-0.1, 0.1), axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, 1, 0.25), h = pretty(e[o2]), col = grid_col, lwd = 0.8)
  axis(1, at = seq(0, 1, 0.25), col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
  axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
  points(v[o2], e[o2], pch = 21, cex = 1.2, lwd = 0.5, col = bg, bg = ifelse(g[o2] == "High", col_hi, col_lo))
  text(v[o2], e[o2], names(v)[o2], pos = 4, cex = 0.45, col = ink2, offset = 0.3)
  rho <- cor(v[o2], e[o2], method = "spearman")
  mtext(sprintf("%s expression vs methylation at %s", gene, tp), side = 3, line = 1.3, adj = 0, font = 2, cex = 0.75, col = ink)
  mtext(sprintf("Spearman rho = %.2f · expression High vs Low: Δ = %+.2f, p = %.1e",
                rho, ex$diff_High_minus_Low, ex$p_welch), side = 3, line = 0.35, adj = 0, cex = 0.6, col = ink2)
  mtext("beta", side = 1, line = 2, cex = 0.62, col = ink2); mtext("log2(TPM+1)", side = 2, line = 2.6, cex = 0.62, col = ink2)
}
dev.off()

# ------------------------------------- figure 2: methylation vs RNA ---------
if (has_expr) {
  draw_cross <- function() {
    x <- ep[!is.na(ep$p_meth_vs_expr), ]
    layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE))
    # (a) correlation with expression along the gene
    par(mar = c(3.6, 4, 3, 1), fg = ink2)
    pos <- as.numeric(x$CpG_beg) / 1e6
    plot(NA, xlim = range(as.numeric(ep$CpG_beg) / 1e6), ylim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    pr <- ep$region == "promoter"
    if (any(pr)) rect(as.numeric(ep$CpG_beg[pr]) / 1e6 - 0.002, -1, as.numeric(ep$CpG_beg[pr]) / 1e6 + 0.002, 1,
                      col = "#ecebe6", border = NA)
    abline(h = seq(-1, 1, 0.25), col = grid_col, lwd = 0.8); abline(h = 0, col = ink2, lwd = 0.8)
    axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
    sig <- x$fdr_meth_vs_expr < 0.1
    points(pos, x$rho_meth_vs_expr, pch = 21, col = "#fcfcfb", lwd = 0.4, bg = ifelse(sig, ink, col_ns),
           cex = 0.6 + pmin(-log10(x$p_meth_vs_expr), 8) / 4)
    gd <- x$p_welch < 0.01                    # also differs between mTOR groups
    points(pos[gd], x$rho_meth_vs_expr[gd], pch = 21, bg = NA, col = col_hi, lwd = 1.6,
           cex = 1.1 + pmin(-log10(x$p_meth_vs_expr[gd]), 8) / 4)
    mtext(sprintf("%s · methylation vs %s RNA, per EPIC probe (%d lines)", gene, gene, max(x$n_lines_expr)),
          side = 3, line = 1.3, adj = 0, font = 2, cex = 0.85, col = ink)
    mtext("dark = correlation FDR<0.1 across probes · blue ring = also differs High vs Low (p<0.01) · shaded = promoter probes · size ~ -log10 p",
          side = 3, line = 0.35, adj = 0, cex = 0.6, col = ink2)
    mtext(sprintf("%s position (Mb, hg38)", ep$CpG_chrm[1]), side = 1, line = 2, cex = 0.62, col = ink2)
    mtext("Spearman rho (methylation vs RNA)", side = 2, line = 2.6, cex = 0.62, col = ink2)
    # (b) within-group correlations: lineage-independent links sit off the axes in the same quadrant
    par(mar = c(4, 4, 3, 1))
    plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    abline(h = seq(-1, 1, 0.5), v = seq(-1, 1, 0.5), col = grid_col, lwd = 0.8); abline(h = 0, v = 0, col = ink2, lwd = 0.8)
    axis(1, col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
    points(x$rho_expr_within_Low, x$rho_expr_within_High, pch = 21, col = "#fcfcfb", lwd = 0.4,
           bg = ifelse(sig, ink, col_ns), cex = 0.9)
    lb <- head(order(x$p_meth_vs_expr), 6)
    text(x$rho_expr_within_Low[lb], x$rho_expr_within_High[lb], x$probeID[lb], pos = 4, cex = 0.5, col = ink, offset = 0.3)
    mtext("Within-group correlation with RNA", side = 3, line = 1.3, adj = 0, font = 2, cex = 0.75, col = ink)
    mtext(sprintf("High n=%d, Low n=%d · same sign in both = link not explained by lineage",
                  sum(cl$mTOR_group[match(colnames(Bm), cl$cell_line)] == "High"),
                  sum(cl$mTOR_group[match(colnames(Bm), cl$cell_line)] == "Low")),
          side = 3, line = 0.35, adj = 0, cex = 0.58, col = ink2)
    mtext("rho within mTOR-Low", side = 1, line = 2.2, cex = 0.62, col = ink2)
    mtext("rho within mTOR-High", side = 2, line = 2.6, cex = 0.62, col = ink2)
    # (c) best group-adjusted probe: expression vs beta, one fit line per group
    bp <- x$probeID[which.min(x$p_expr_adj_group)]
    v <- Bm[bp, ]; e <- E[match(names(v), rownames(E)), gene]; g <- cl$mTOR_group[match(names(v), cl$cell_line)]
    o <- !is.na(v) & !is.na(e)
    plot(NA, xlim = c(0, 1), ylim = range(e[o]) + c(-0.1, 0.1), axes = FALSE, xlab = "", ylab = "")
    abline(v = seq(0, 1, 0.25), h = pretty(e[o]), col = grid_col, lwd = 0.8)
    axis(1, at = seq(0, 1, 0.25), col = NA, col.ticks = ink2, col.axis = ink2, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
    axis(2, col = NA, col.ticks = ink2, col.axis = ink2, las = 1, cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.5, 0))
    for (k in c("Low", "High")) { s <- o & g == k
      if (sum(s) >= 4) { f <- lm(e[s] ~ v[s]); xs <- range(v[s])
        lines(xs, coef(f)[1] + coef(f)[2] * xs, col = if (k == "High") col_hi else col_lo, lwd = 2) } }
    points(v[o], e[o], pch = 21, cex = 1.2, lwd = 0.5, col = "#fcfcfb", bg = ifelse(g[o] == "High", col_hi, col_lo))
    text(v[o], e[o], names(v)[o], pos = 4, cex = 0.45, col = ink2, offset = 0.3)
    r <- x[x$probeID == bp, ]
    mtext(sprintf("Best group-adjusted link: %s (%s, %+d bp from TSS)", bp, r$region, as.integer(r$dist_to_TSS)),
          side = 3, line = 1.3, adj = 0, font = 2, cex = 0.75, col = ink)
    mtext(sprintf("RNA ~ M + group: p = %.1e · rho all = %.2f, in High = %.2f, in Low = %.2f",
                  r$p_expr_adj_group, r$rho_meth_vs_expr, r$rho_expr_within_High, r$rho_expr_within_Low),
          side = 3, line = 0.35, adj = 0, cex = 0.58, col = ink2)
    mtext("beta", side = 1, line = 2, cex = 0.62, col = ink2); mtext("log2(TPM+1)", side = 2, line = 2.6, cex = 0.62, col = ink2)
  }
  f2 <- file.path(fig_dir, paste0("gene_", gene, "_methylation_vs_RNA"))
  png(paste0(f2, ".png"), width = 3200, height = 1800, res = 200, type = "cairo", bg = "#fcfcfb"); draw_cross(); dev.off()
  cairo_pdf(paste0(f2, ".pdf"), width = 16, height = 9, bg = "#fcfcfb"); draw_cross(); dev.off()
}
cat(readLines(file.path(out_dir, paste0(gene, "_methylation_summary.txt"))), sep = "\n")
