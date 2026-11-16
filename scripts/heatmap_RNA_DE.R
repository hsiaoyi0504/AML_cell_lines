# Differential RNA expression between two groups of AML lines (DepMap 26Q3
# log2(TPM+1), genes expressed in >= 5 lines) and heatmaps of
#   * the top DE genes in each direction (optionally excluding signature genes), and
#   * the genes of the MYC / mTOR RNA signatures (data/*_signature_RNA.tsv).
#
# Comparisons:
#   MYC_vs_mTOR : MYC-High/mTOR-Low (n=25) vs mTOR-High/MYC-Low (n=25); lines high or
#                 low for both scores are excluded. Both groups are defined from RNA
#                 scores, so this shows the programs behind the groups (circular).
#   KMT2A       : KMT2A-rearranged (n=10) vs all other lines (n=46)
#
# Cell values: "z" = row z-score across the plotted lines (capped +-2.5);
#              "median" = log2(TPM+1) minus the gene's median across the plotted lines
#              (log2 fold vs the gene's median, capped +-4).
#
# Run from the project root: Rscript scripts/heatmap_RNA_DE.R
# The heatmaps produced are listed at the bottom of this file.

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R")
cl  <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
E0  <- attr(readRDS(file.path(res_dir, "all_results.rds"))$expression, "X")
colnames(E0) <- sub(" \\([0-9]+\\)$", "", colnames(E0))
read_genes <- function(f) { x <- readLines(f)[-1]; unique(x[x != ""]) }   # one symbol per line, header skipped

# TP53 damaging mutations (DepMap 26Q3 OmicsSomaticMutationsMatrixDamaging, default entry per model);
# NA = line not in the mutation file. Only TP53's column is read.
mut_file <- "data/cell_lines/DepMap_26Q3/OmicsSomaticMutationsMatrixDamaging.csv"
j <- which(sub(" \\(.*$", "", strsplit(readLines(mut_file, n = 1), ",")[[1]]) == "TP53")
mut <- read.csv(pipe(paste0("cut -d, -f3,5,", j, " ", mut_file)), check.names = FALSE)
mut <- mut[mut$IsDefaultEntryForModel == "Yes", ]
cl$TP53_mut <- mut[[3]][match(cl$ModelID, mut$ModelID)] > 0
col_tp53 <- "#eda100"

# Activating hotspot mutations (DepMap 26Q3 OmicsSomaticMutationsMatrixHotspot, default entry);
# RAS = NRAS/KRAS/PTPN11, JAK2/KIT = JAK2/KIT. NA = line not in the hotspot file.
hot <- read.csv("data/cell_lines/DepMap_26Q3/OmicsSomaticMutationsMatrixHotspot.csv", check.names = FALSE)
hot <- hot[hot$IsDefaultEntryForModel == "Yes", ]
hot_any <- function(genes) {
  j <- which(sub(" \\(.*$", "", names(hot)) %in% genes)
  rowSums(hot[, j, drop = FALSE] > 0, na.rm = TRUE)[match(cl$ModelID, hot$ModelID)] > 0
}
cl$RAS_mut    <- hot_any(c("NRAS", "KRAS", "PTPN11"))
cl$JAKKIT_mut <- hot_any(c("JAK2", "KIT"))
col_ras <- "#1f6fb8"; col_jakkit <- "#a8326e"   # validated with TP53/Mito colours (CVD dE >= 10.7)
col_neg <- "#ffffff"; col_neg_border <- "#c9c8c3"; col_na <- "#a3a29d"   # negative = white, NA = gray

# Mito-high / Mito-low calls from the published heatmap (data/mito_classification.tsv); NA = not in that study
mito <- read.delim("data/mito_classification.tsv")
cl$Mito <- mito$Mito_group[match(cl$cell_line, toupper(gsub("[^A-Za-z0-9]", "", mito$cell_line_as_published)))]
col_mito_hi <- "#008300"; col_mito_lo <- "#e87ba4"   # green vs pink: validated pair (CVD dE 17.6)
pal <- colorRampPalette(c("#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f2a8a6", "#e34948", "#9c2423"))(101)
bg  <- "#fcfcfb"; col_k <- "#1b1b1a"; col_other <- "#e8e7e3"

# ------------------------------------------------- groups + DE table ---------
de_groups <- function(cmp) {
  g <- cl[match(rownames(E0), cl$cell_line), ]
  if (cmp == "MYC_vs_mTOR") {
    grpA <- g$MYC_group == "High" & g$mTOR_group == "Low"; grpB <- g$mTOR_group == "High" & g$MYC_group == "Low"
    info <- list(labA = "MYC-High", labB = "mTOR-High", colA = "#eb6834", colB = "#4a3aa7",
                 note = "groups defined from RNA scores, so differences are partly by construction")
  } else if (cmp == "KMT2A") {
    grpA <- g$KMT2A_rearranged; grpB <- !g$KMT2A_rearranged
    info <- list(labA = "KMT2A-r", labB = "other", colA = "#4a3aa7", colB = "#d6d5d0",
                 note = "KMT2A status from the classification file")
  } else stop("comparison must be MYC_vs_mTOR or KMT2A")
  keep <- grpA | grpB
  E <- E0[keep, ]; a <- grpA[keep]
  A <- E[a, , drop = FALSE]; B <- E[!a, , drop = FALSE]
  mA <- colMeans(A); mB <- colMeans(B); vA <- apply(A, 2, var); vB <- apply(B, 2, var)
  nA <- nrow(A); nB <- nrow(B)
  t  <- (mA - mB) / sqrt(vA / nA + vB / nB)
  df <- (vA / nA + vB / nB)^2 / ((vA / nA)^2 / (nA - 1) + (vB / nB)^2 / (nB - 1))
  de <- data.frame(gene = colnames(E), mean_A = mA, mean_B = mB, log2FC_A_vs_B = mA - mB, t = t,
                   p_welch = 2 * pt(-abs(t), df))
  names(de)[2:3] <- paste0("mean_", c(info$labA, info$labB))
  de <- de[is.finite(de$p_welch), ]
  de$fdr <- p.adjust(de$p_welch, "BH")
  de$p_wilcox <- vapply(de$gene, function(x) suppressWarnings(wilcox.test(A[, x], B[, x], exact = FALSE)$p.value), 0)
  de <- de[order(de$p_welch), ]
  info$tag <- gsub("[^A-Za-z0-9]+", "_", paste(info$labA, "vs", info$labB))
  write.table(de, file.path(res_dir, paste0("DE_RNA_", info$tag, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- de[de$fdr < 0.05 & abs(de$log2FC_A_vs_B) >= 1, ]
  cat(sprintf("\n=== %s (n=%d) vs %s (n=%d); excluded lines: %s\n", info$labA, nA, info$labB, nB,
              paste(setdiff(cl$cell_line, rownames(E)), collapse = ", ")))
  cat(sprintf("Genes tested: %d | FDR<0.05: %d | FDR<0.05 & |log2FC|>=1: %d (higher in %s %d / in %s %d)\n",
              nrow(de), sum(de$fdr < 0.05), nrow(sig), info$labA, sum(sig$log2FC_A_vs_B > 0), info$labB, sum(sig$log2FC_A_vs_B < 0)))
  c(info, list(E = E, a = a, de = de, sig = sig, nA = nA, nB = nB, g = cl[match(rownames(E), cl$cell_line), ]))
}

# ------------------------------------------------------------ heatmap --------
# grp       : output of de_groups()
# genes     : NULL -> top n_top DE genes per direction; or a gene-list file to plot instead
# exclude   : gene-list files whose genes are dropped before picking the top DE genes
# values    : "z" or "median" (see header)
heatmap_rna <- function(grp, n_top = 25, genes = NULL, exclude = NULL, gene_labels = TRUE, values = "z") {
  E <- grp$E; a <- grp$a; de <- grp$de; sig <- grp$sig; g <- grp$g
  labA <- grp$labA; labB <- grp$labB; colA <- grp$colA; colB <- grp$colB
  set_lab <- NA; ex_lab <- NA
  block_lab <- c(sprintf("higher in %s", labA), sprintf("higher in %s", labB))
  multi <- length(genes) > 1                  # several gene lists -> one block of rows per list
  if (multi) {
    nm <- sub("_RNA$", "", sub("\\.[^.]*$", "", basename(genes)))
    sets <- lapply(genes, function(f) intersect(read_genes(f), colnames(E)))
    blocks <- lapply(seq_along(sets), function(i) {        # within each list: genes by log2FC
      s <- setdiff(sets[[i]], unlist(sets[seq_len(i - 1)]))
      s[order(-de$log2FC_A_vs_B[match(s, de$gene)])] })
    for (i in seq_along(blocks)) {
      fc <- de$log2FC_A_vs_B[match(blocks[[i]], de$gene)]; fd <- de$fdr[match(blocks[[i]], de$gene)]
      cat(sprintf("%s: %d genes plotted | higher in %s %d (FDR<0.05: %d) | higher in %s %d (FDR<0.05: %d)\n", nm[i],
                  length(blocks[[i]]), labA, sum(fc > 0), sum(fc > 0 & fd < 0.05), labB, sum(fc <= 0), sum(fc <= 0 & fd < 0.05)))
    }
    stopifnot(length(blocks) == 2)            # drawn as top block (first list) and bottom block (second list)
    upA <- blocks[[1]]; upB <- blocks[[2]]
    set_lab <- paste(nm, collapse = "_and_"); block_lab <- nm
  } else if (is.null(genes)) {
    ex_genes <- character(0)
    if (length(exclude)) {
      ex_genes <- unique(unlist(lapply(exclude, read_genes)))
      ex_lab <- paste(sub("\\.[^.]*$", "", basename(exclude)), collapse = " + ")
      cat(sprintf("Excluding %d genes (%s); %d of them among the %d DE genes passing FDR<0.05 & |log2FC|>=1\n",
                  length(ex_genes), ex_lab, sum(sig$gene %in% ex_genes), nrow(sig)))
    }
    s2 <- sig[!sig$gene %in% ex_genes, ]
    upA <- head(s2$gene[s2$log2FC_A_vs_B > 0], n_top); upB <- head(s2$gene[s2$log2FC_A_vs_B < 0], n_top)
  } else {
    gl <- read_genes(genes); set_lab <- sub("\\.[^.]*$", "", basename(genes))
    absent <- setdiff(gl, colnames(E)); gl <- intersect(gl, colnames(E))
    fc <- de$log2FC_A_vs_B[match(gl, de$gene)]; fd <- de$fdr[match(gl, de$gene)]
    upA <- gl[fc > 0][order(-fc[fc > 0])]; upB <- gl[fc <= 0][order(-fc[fc <= 0])]
    cat(sprintf("%s: %d genes, %d plotted (not plotted: %s)\n", set_lab, length(gl) + length(absent), length(gl),
                if (length(absent)) paste(absent, collapse = ", ") else "none"))
    cat(sprintf("  higher in %s: %d (FDR<0.05: %d) | higher in %s: %d (FDR<0.05: %d)\n", labA, sum(fc > 0),
                sum(fc > 0 & fd < 0.05), labB, sum(fc <= 0), sum(fc <= 0 & fd < 0.05)))
  }

  X <- E[, c(upA, upB)]
  if (values == "z") {
    zmax <- 2.5; Z <- scale(X); val_lab <- "z-score"; val_txt <- "row z-score (capped ±2.5)"
  } else {
    zmax <- 4; Z <- sweep(X, 2, apply(X, 2, median)); val_lab <- "log2 vs\ngene median"
    val_txt <- "log2(TPM+1) minus each gene's median across these lines (capped ±4)"
  }
  Z[Z > zmax] <- zmax; Z[Z < -zmax] <- -zmax

  if (is.null(genes)) {                       # top-DE mode: cluster lines and genes within blocks
    ord_block <- function(r) r[hclust(dist(Z[r, , drop = FALSE]), "average")$order]
    cols <- c(ord_block(which(a)), ord_block(which(!a)))
    row_block <- function(x) x[hclust(dist(t(Z[, x, drop = FALSE])), "average")$order]
    rows <- c(row_block(upB), row_block(upA))                    # bottom: higher in B, top: higher in A
  } else {                                    # gene-list mode: lines by mean value over the list, genes by log2FC
    sc <- if (multi) rowMeans(Z[, upA, drop = FALSE]) - rowMeans(Z[, upB, drop = FALSE]) else rowMeans(Z)
    cols <- c(which(a)[order(-sc[a])], which(!a)[order(-sc[!a])])
    rows <- rev(c(upA, upB))
  }
  Zp <- Z[cols, rows]; gp <- g[cols, ]
  lab_col <- ifelse(gp$KMT2A_rearranged, paste0(rownames(Zp), "*"), rownames(Zp))
  title <- if (is.null(genes)) sprintf("Top differentially expressed genes%s: %s (n=%d) vs %s (n=%d)",
                                       if (is.na(ex_lab)) "" else " (signature genes removed)", labA, grp$nA, labB, grp$nB)
           else if (multi) sprintf("%s (%d) and %s (%d) genes: %s (n=%d) vs %s (n=%d)", block_lab[1], length(upA),
                                   block_lab[2], length(upB), labA, grp$nA, labB, grp$nB)
           else sprintf("%s genes (%d): %s (n=%d) vs %s (n=%d)", set_lab, length(rows), labA, grp$nA, labB, grp$nB)
  subtitle <- if (multi)
    sprintf("DepMap 26Q3 RNA · genes ordered by log2FC within each signature · lines ordered within group by mean %s − mean %s value · %s",
            block_lab[1], block_lab[2], val_txt)
  else if (is.null(genes))
    sprintf("DepMap 26Q3 RNA · Welch t, FDR<0.05 and |log2FC|≥1 · top %d per direction by p%s (number = log2FC) · %s",
            n_top, if (is.na(ex_lab)) "" else sprintf(", excluding genes in %s", ex_lab), val_txt)
  else sprintf("DepMap 26Q3 RNA · genes ordered by log2FC%s · lines ordered by mean value within group · %s",
               if (gene_labels) sprintf(" (number = %s − %s; * = FDR<0.05)", labA, labB) else "", val_txt)

  draw <- function() {
    nr <- ncol(Zp); nc <- nrow(Zp)
    par(oma = c(0, 1, 5.4, 1))
    layout(matrix(1:4, 2, byrow = TRUE), widths = c(1, 0.14), heights = c(0.17, 1))
    lm <- if (gene_labels) 11.5 else 7.5
    par(mar = c(0.4, lm, 0.6, 1), fg = ink2)
    plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(-2.5, 3.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    rect(seq_len(nc) - 0.5, 2.55, seq_len(nc) + 0.5, 3.45, border = bg, lwd = 0.6, col = ifelse(a[cols], colA, colB))
    # annotation cell: positive = colour, negative = white (outlined), NA = gray
    ann_row <- function(y, fill, is_na) {
      fill[is_na] <- col_na
      rect(seq_len(nc) - 0.5, y - 0.45, seq_len(nc) + 0.5, y + 0.45, col = fill, lwd = 0.6,
           border = ifelse(fill == col_neg, col_neg_border, bg))
    }
    ann_row(0, ifelse(gp$Mito %in% "High", col_mito_hi, col_mito_lo), is.na(gp$Mito))
    ann_row(1, ifelse(gp$TP53_mut %in% TRUE, col_tp53, col_neg), is.na(gp$TP53_mut))
    ann_row(2, ifelse(gp$KMT2A_rearranged %in% TRUE, col_k, col_neg), is.na(gp$KMT2A_rearranged))
    ann_row(-1, ifelse(gp$RAS_mut %in% TRUE, col_ras, col_neg), is.na(gp$RAS_mut))
    ann_row(-2, ifelse(gp$JAKKIT_mut %in% TRUE, col_jakkit, col_neg), is.na(gp$JAKKIT_mut))
    axis(2, at = c(3, 2, 1, 0, -1, -2), labels = c("Group", "KMT2A-r", "TP53 mut", "Mito", "RAS mut", "JAK2/KIT mut"),
         las = 1, tick = FALSE, col.axis = ink, cex.axis = 1)
    mtext(title, side = 3, line = 3.4, adj = 0, outer = TRUE, font = 2, cex = 1.3, col = ink)
    mtext(subtitle, side = 3, line = 1.9, adj = 0, outer = TRUE, cex = 0.82, col = ink2)
    mtext(paste0("Note: ", grp$note), side = 3, line = 0.7, adj = 0, outer = TRUE, cex = 0.82, col = ink2, font = 3)
    par(mar = c(0.4, 0.5, 0.6, 0.5)); plot.new()
    legend("left", bty = "n", cex = 0.85, text.col = ink, pt.cex = 1.4, pch = 22, pt.lwd = 0.8,
           col = c(colA, colB, col_k, col_tp53, col_mito_hi, col_mito_lo, col_ras, col_jakkit, col_neg_border, col_na),
           pt.bg = c(colA, colB, col_k, col_tp53, col_mito_hi, col_mito_lo, col_ras, col_jakkit, col_neg, col_na),
           legend = c(labA, labB, "KMT2A-r (*)", "TP53 mutant", "Mito-high", "Mito-low", "NRAS/KRAS/PTPN11", "JAK2/KIT", "negative", "NA / not profiled"))
    par(mar = c(8.2, lm, 0.4, 1))
    image(seq_len(nc), seq_len(nr), Zp, useRaster = TRUE, col = pal, zlim = c(-zmax, zmax), axes = FALSE, xlab = "", ylab = "")
    abline(v = sum(a) + 0.5, col = bg, lwd = 3); abline(h = length(upB) + 0.5, col = bg, lwd = 3)
    if (gene_labels) {
      fc <- de$log2FC_A_vs_B[match(rows, de$gene)]
      star <- if (is.null(genes)) "" else ifelse(de$fdr[match(rows, de$gene)] < 0.05, "*", " ")
      axis(2, at = seq_len(nr), labels = sprintf("%s  %+.1f%s", rows, fc, star), las = 1, tick = FALSE,
           cex.axis = min(0.85, 72 / nr), col.axis = ink, line = -0.6)
    }
    # column labels coloured per line (axis() takes one colour, so draw them with text())
    text(seq_len(nc), par("usr")[3] - 0.6, lab_col, srt = 90, adj = 1, xpd = NA, cex = 0.78,
         col = ifelse(a[cols], colA, colB))
    lab_line <- if (gene_labels) 10 else 1
    if (length(upA) >= max(4, 0.2 * nr)) mtext(sprintf("%s (%d)", block_lab[1], length(upA)), side = 2, line = lab_line,
                                at = length(upB) + length(upA) / 2, cex = 0.9, col = ink2)
    if (length(upB) >= max(4, 0.2 * nr)) mtext(sprintf("%s (%d)", block_lab[2], length(upB)), side = 2, line = lab_line,
                                at = length(upB) / 2, cex = 0.9, col = ink2)
    par(mar = c(8.2, 1.2, 8, 3.8))
    zs <- seq(-zmax, zmax, length.out = 101)
    image(1, zs, t(zs), useRaster = TRUE, col = pal, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = seq(-floor(zmax), floor(zmax)), las = 1, cex.axis = 0.9, col = NA, col.ticks = ink2, col.axis = ink2, tcl = -0.2)
    mtext(val_lab, side = 3, line = 0.5, cex = 0.9, col = ink2)
  }
  fb <- file.path(fig_dir, paste0("heatmap_RNA_", grp$tag, if (is.na(set_lab)) "" else paste0("_", set_lab),
                                  if (is.na(ex_lab)) "" else "_excl_signatures"))
  h_px <- if (gene_labels) 2900 + max(0, length(c(upA, upB)) - 50) * 42 else 2900
  png(paste0(fb, ".png"), width = 3600, height = h_px, res = 220, type = "cairo", bg = bg); draw(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 3600 / 220, height = h_px / 220, bg = bg); draw(); dev.off()
  cat("  ->", paste0(fb, ".png/.pdf"), "\n")
  invisible(list(upA = upA, upB = upB))
}

# ------------------------------------------------ two heatmaps side by side --
# Left: MYC + mTOR signature genes (two blocks). Right: top n_top DE genes per direction with
# signature genes removed. Both median-normalised on one colour scale (±zmax); cell lines in the
# same order in both panels (within group by mean MYC-signature − mean mTOR-signature value).
heatmap_pair <- function(grp, signatures, n_top = 25, zmax = 4) {
  E <- grp$E; a <- grp$a; de <- grp$de; sig <- grp$sig; g <- grp$g
  labA <- grp$labA; labB <- grp$labB; colA <- grp$colA; colB <- grp$colB
  med <- function(X) { Z <- sweep(X, 2, apply(X, 2, median)); Z[Z > zmax] <- zmax; Z[Z < -zmax] <- -zmax; Z }
  fc_of <- function(x) de$log2FC_A_vs_B[match(x, de$gene)]

  # left: signature blocks, genes by log2FC within each list
  nm <- sub("_RNA$", "", sub("\\.[^.]*$", "", basename(signatures)))
  s1 <- intersect(read_genes(signatures[1]), colnames(E))
  s2 <- setdiff(intersect(read_genes(signatures[2]), colnames(E)), s1)
  s1 <- s1[order(-fc_of(s1))]; s2 <- s2[order(-fc_of(s2))]
  ZL <- med(E[, c(s1, s2)])
  # right: top DE genes, signature genes removed, clustered within direction
  ex <- unique(unlist(lapply(signatures, read_genes)))
  sg <- sig[!sig$gene %in% ex, ]
  dA <- head(sg$gene[sg$log2FC_A_vs_B > 0], n_top); dB <- head(sg$gene[sg$log2FC_A_vs_B < 0], n_top)
  ZR <- med(E[, c(dA, dB)])
  row_block <- function(Z, x) x[hclust(dist(t(Z[, x, drop = FALSE])), "average")$order]

  sc <- rowMeans(ZL[, s1]) - rowMeans(ZL[, s2])                 # shared line order
  cols <- c(which(a)[order(-sc[a])], which(!a)[order(-sc[!a])])
  panels <- list(
    list(Z = ZL[cols, rev(c(s1, s2))], nB = length(s2), blocks = c(sprintf("%s (%d)", nm[1], length(s1)), sprintf("%s (%d)", nm[2], length(s2))),
         title = "MYC and mTOR signature genes", labels = FALSE),
    list(Z = ZR[cols, c(row_block(ZR, dB), row_block(ZR, dA))], nB = length(dB),
         blocks = c(sprintf("higher in %s (%d)", labA, length(dA)), sprintf("higher in %s (%d)", labB, length(dB))),
         title = sprintf("Top differentially expressed genes (signature genes removed)"), labels = TRUE))
  gp <- g[cols, ]; ac <- a[cols]; nc <- length(cols)
  lab_col <- ifelse(gp$KMT2A_rearranged, paste0(rownames(E)[cols], "*"), rownames(E)[cols])

  ann_panel <- function(row_labels) {
    plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(-2.5, 3.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    rect(seq_len(nc) - 0.5, 2.55, seq_len(nc) + 0.5, 3.45, border = bg, lwd = 0.6, col = ifelse(ac, colA, colB))
    ann_row <- function(y, fill, is_na) {
      fill[is_na] <- col_na
      rect(seq_len(nc) - 0.5, y - 0.45, seq_len(nc) + 0.5, y + 0.45, col = fill, lwd = 0.6,
           border = ifelse(fill == col_neg, col_neg_border, bg))
    }
    ann_row(2, ifelse(gp$KMT2A_rearranged %in% TRUE, col_k, col_neg), is.na(gp$KMT2A_rearranged))
    ann_row(1, ifelse(gp$TP53_mut %in% TRUE, col_tp53, col_neg), is.na(gp$TP53_mut))
    ann_row(0, ifelse(gp$Mito %in% "High", col_mito_hi, col_mito_lo), is.na(gp$Mito))
    ann_row(-1, ifelse(gp$RAS_mut %in% TRUE, col_ras, col_neg), is.na(gp$RAS_mut))
    ann_row(-2, ifelse(gp$JAKKIT_mut %in% TRUE, col_jakkit, col_neg), is.na(gp$JAKKIT_mut))
    if (row_labels) axis(2, at = c(3, 2, 1, 0, -1, -2), labels = c("Group", "KMT2A-r", "TP53 mut", "Mito", "RAS mut", "JAK2/KIT mut"),
                         las = 1, tick = FALSE, col.axis = ink, cex.axis = 1)
  }
  heat_panel <- function(p, left) {
    nr <- ncol(p$Z)
    image(seq_len(nc), seq_len(nr), p$Z, useRaster = TRUE, col = pal, zlim = c(-zmax, zmax), axes = FALSE, xlab = "", ylab = "")
    abline(v = sum(ac) + 0.5, col = bg, lwd = 3); abline(h = p$nB + 0.5, col = bg, lwd = 3)
    text(seq_len(nc), par("usr")[3] - 0.6 * nr / 100, lab_col, srt = 90, adj = 1, xpd = NA, cex = 0.72, col = ifelse(ac, colA, colB))
    if (p$labels) axis(4, at = seq_len(nr), labels = colnames(p$Z), las = 1, tick = FALSE, cex.axis = min(0.75, 70 / nr),
                       col.axis = ink, line = -0.6)
    side <- if (left) 2 else 4; ln <- if (left) 1 else 5.2
    mtext(p$blocks[1], side = side, line = ln, at = p$nB + (nr - p$nB) / 2, cex = 0.85, col = ink2)
    mtext(p$blocks[2], side = side, line = ln, at = p$nB / 2, cex = 0.85, col = ink2)
  }

  draw <- function() {
    par(oma = c(0, 1, 4.6, 1))
    layout(rbind(c(1, 2, 3), c(4, 5, 6)), widths = c(1, 1, 0.22), heights = c(0.2, 1))
    par(mar = c(0.4, 7.5, 2.2, 0.6), fg = ink2); ann_panel(TRUE)
    mtext(panels[[1]]$title, side = 3, line = 0.5, adj = 0, font = 2, cex = 0.95, col = ink)
    par(mar = c(0.4, 1.5, 2.2, 6.6)); ann_panel(FALSE)
    mtext(panels[[2]]$title, side = 3, line = 0.5, adj = 0, font = 2, cex = 0.95, col = ink)
    par(mar = c(0.4, 0.5, 2.2, 0.5)); plot.new()
    legend("left", bty = "n", cex = 1.05, text.col = ink, pt.cex = 1.8, pch = 22, pt.lwd = 0.8, y.intersp = 1.05,
           col = c(colA, colB, col_k, col_tp53, col_mito_hi, col_mito_lo, col_ras, col_jakkit, col_neg_border, col_na),
           pt.bg = c(colA, colB, col_k, col_tp53, col_mito_hi, col_mito_lo, col_ras, col_jakkit, col_neg, col_na),
           legend = c(labA, labB, "KMT2A-r (*)", "TP53 mutant", "Mito-high", "Mito-low", "NRAS/KRAS/PTPN11", "JAK2/KIT", "negative", "NA / not profiled"))
    par(mar = c(7.6, 7.5, 0.4, 0.6)); heat_panel(panels[[1]], TRUE)
    par(mar = c(7.6, 1.5, 0.4, 6.6)); heat_panel(panels[[2]], FALSE)
    par(mar = c(7.6, 1.5, 12, 4.2))
    zs <- seq(-zmax, zmax, length.out = 101)
    image(1, zs, t(zs), useRaster = TRUE, col = pal, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = seq(-zmax, zmax), las = 1, cex.axis = 0.9, col = NA, col.ticks = ink2, col.axis = ink2, tcl = -0.2)
    mtext("log2 vs\ngene median", side = 3, line = 0.5, cex = 0.85, col = ink2)
    mtext(sprintf("RNA expression: %s (n=%d) vs %s (n=%d)", labA, grp$nA, labB, grp$nB),
          side = 3, line = 3, adj = 0, outer = TRUE, font = 2, cex = 1.3, col = ink)
    mtext(sprintf("DepMap 26Q3 · both panels: log2(TPM+1) minus each gene's median across these lines (capped ±%d), same colour scale · lines in the same order (within group by mean %s − mean %s) · right: Welch FDR<0.05, |log2FC|≥1, top %d per direction",
                  zmax, nm[1], nm[2], n_top), side = 3, line = 1.6, adj = 0, outer = TRUE, cex = 0.75, col = ink2)
    mtext(paste0("Note: ", grp$note), side = 3, line = 0.4, adj = 0, outer = TRUE, cex = 0.75, col = ink2, font = 3)
  }
  fb <- file.path(fig_dir, paste0("heatmap_RNA_", grp$tag, "_signatures_and_topDE"))
  h_px <- 2900 + max(0, 2 * n_top - 50) * 22               # taller when more genes are labelled
  png(paste0(fb, ".png"), width = 7000, height = h_px, res = 220, type = "cairo", bg = bg); draw(); dev.off()
  cairo_pdf(paste0(fb, ".pdf"), width = 7000 / 220, height = h_px / 220, bg = bg); draw(); dev.off()
  cat("  ->", paste0(fb, ".png/.pdf"), "\n")
}

# ================================================================ run =========
myc_mtor <- de_groups("MYC_vs_mTOR")
signatures <- c("data/MYC_signature_RNA.tsv", "data/mTOR_signature_RNA.tsv")

heatmap_rna(myc_mtor)                                         # top 25 DE genes each way
heatmap_rna(myc_mtor, exclude = signatures)                   # same, signature genes removed
heatmap_rna(myc_mtor, genes = signatures[1], gene_labels = FALSE, values = "median")   # MYC signature
heatmap_rna(myc_mtor, genes = signatures[2], gene_labels = FALSE, values = "median")   # mTOR signature
heatmap_rna(myc_mtor, genes = signatures, gene_labels = FALSE, values = "median")       # MYC + mTOR signatures stacked
heatmap_pair(myc_mtor, signatures, n_top = 50)                          # signatures | top DE (excl. signatures), side by side
