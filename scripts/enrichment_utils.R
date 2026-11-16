# Shared gene-set enrichment helpers for the methylation analyses.
# Gene sets: data/gene_sets/*.gmt (Enrichr libraries). Base R only.
# source() from project root; uses the palette from scripts/plot_utils.R for plots.

enrich_libs <- c(Hallmark = "MSigDB_Hallmark_2020", GO_BP = "GO_Biological_Process_2023",
                 KEGG = "KEGG_2021_Human", CellType = "PanglaoDB_Augmented_2021")

read_gmt <- function(f) {
  l <- strsplit(readLines(f, warn = FALSE), "\t")
  setNames(lapply(l, function(x) unique(x[-(1:2)][x[-(1:2)] != ""])), vapply(l, `[`, "", 1))
}

# hypergeometric over-representation of `hits` among `universe` (sets of 10-500 genes)
ora <- function(sets, hits, universe, lib, dir) {
  hits <- intersect(hits, universe)
  r <- do.call(rbind, lapply(names(sets), function(s) {
    g <- intersect(sets[[s]], universe); if (length(g) < 10 || length(g) > 500) return(NULL)
    ov <- intersect(g, hits)
    p  <- phyper(length(ov) - 1, length(g), length(universe) - length(g), length(hits), lower.tail = FALSE)
    data.frame(library = lib, direction = dir, term = s, set_size = length(g), overlap = length(ov),
               expected = length(g) * length(hits) / length(universe), p = p,
               genes = paste(head(ov, 25), collapse = ","))
  }))
  r$fdr <- p.adjust(r$p, "BH"); r[order(r$p), ]
}

# competitive rank test on the signed gene score (ignores inter-gene correlation,
# so p-values are optimistic; use for ranking / direction, not as proof)
rank_test <- function(sets, score, lib) {
  r <- do.call(rbind, lapply(names(sets), function(s) {
    inset <- names(score) %in% sets[[s]]; n <- sum(inset)
    if (n < 10 || n > 500) return(NULL)
    w <- suppressWarnings(wilcox.test(score[inset], score[!inset], exact = FALSE))
    auc <- unname(w$statistic) / (n * sum(!inset))
    data.frame(library = lib, term = s, set_size = n, mean_score = mean(score[inset]),
               auc = auc, direction = ifelse(auc > 0.5, "hyper_in_High", "hypo_in_High"), p = w$p.value)
  }))
  r$fdr <- p.adjust(r$p, "BH"); r[order(r$p), ]
}

# score: named signed gene score (+ = more methylated in mTOR-High)
run_enrichment <- function(hypo, hyper, universe, score, out_prefix) {
  ORA <- list(); RNK <- list()
  for (lb in names(enrich_libs)) {
    sets <- read_gmt(file.path("data/gene_sets", paste0(enrich_libs[[lb]], ".gmt")))
    ORA[[lb]] <- rbind(ora(sets, hypo, universe, lb, "hypo_in_High"), ora(sets, hyper, universe, lb, "hyper_in_High"))
    RNK[[lb]] <- rank_test(sets, score, lb)
  }
  ORA <- do.call(rbind, ORA); RNK <- do.call(rbind, RNK)
  write.table(ORA[order(ORA$library, ORA$direction, ORA$p), ], paste0(out_prefix, "_ORA.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(RNK[order(RNK$library, RNK$p), ], paste0(out_prefix, "_rank.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  list(ORA = ORA, RNK = RNK)
}

print_enrichment <- function(en) {
  for (dir in c("hypo_in_High", "hyper_in_High")) {
    cat(sprintf("\n== ORA, %s (top 8 per library)\n", dir))
    for (lb in names(enrich_libs)) {
      x <- head(en$ORA[en$ORA$library == lb & en$ORA$direction == dir, ], 8)
      print(data.frame(lib = lb, term = substr(x$term, 1, 55), n = x$set_size, ov = x$overlap,
                       exp = round(x$expected, 1), p = signif(x$p, 2), fdr = signif(x$fdr, 2)), row.names = FALSE)
    }
  }
  cat("\n== Rank test (top 10 per library)\n")
  for (lb in names(enrich_libs)) {
    x <- head(en$RNK[en$RNK$library == lb, ], 10)
    print(data.frame(lib = lb, term = substr(x$term, 1, 55), n = x$set_size, auc = round(x$auc, 2),
                     dir = x$direction, fdr = signif(x$fdr, 2)), row.names = FALSE)
  }
}

# diverging bar chart of the top ORA terms per library (blue = more methylated in High)
plot_enrichment <- function(ORA, n_hypo, n_hyper, file, title_suffix = "over-representation") {
  bg <- "#fcfcfb"
  top_terms <- function(lb, dir, n = 6) head(ORA[ORA$library == lb & ORA$direction == dir & ORA$fdr < 0.5, ], n)
  png(file, width = 3200, height = 2300, res = 200, type = "cairo", bg = bg)
  par(mfrow = c(2, 2), mar = c(3.6, 19, 3, 1.5), fg = ink2)
  for (lb in names(enrich_libs)) {
    d <- rbind(top_terms(lb, "hyper_in_High"), top_terms(lb, "hypo_in_High"))
    head_txt <- sprintf("%s  ·  %s", lb, title_suffix)
    if (!nrow(d)) {
      plot.new(); mtext(head_txt, side = 3, line = 1.2, adj = 0, font = 2, cex = 0.85, col = ink)
      text(0.5, 0.5, "no terms at FDR < 0.5", col = ink2, cex = 0.9); next
    }
    v <- ifelse(d$direction == "hyper_in_High", -log10(d$fdr), log10(d$fdr))
    d <- d[order(v), ]; v <- sort(v)
    xl <- max(abs(v), -log10(0.05)) * c(-1.15, 1.15)
    y <- seq_along(v)
    plot(NA, xlim = xl, ylim = c(0.5, length(v) + 0.5), axes = FALSE, xlab = "", ylab = "")
    abline(v = pretty(xl), col = grid_col, lwd = 0.8); abline(v = 0, col = ink2, lwd = 0.8)
    abline(v = c(-1, 1) * -log10(0.05), lty = 2, col = ink2, lwd = 0.6)
    rect(0, y - 0.32, v, y + 0.32, col = ifelse(v > 0, col_hi, col_lo), border = bg, lwd = 1)
    lab <- sprintf("%s  (%d/%d)", substr(d$term, 1, 48), d$overlap, d$set_size)
    axis(2, at = y, labels = lab, las = 1, col = NA, col.axis = ink, cex.axis = 0.72, tick = FALSE)
    axis(1, at = pretty(xl), labels = abs(pretty(xl)), col = NA, col.ticks = ink2, col.axis = ink2,
         cex.axis = 0.75, tcl = -0.25, mgp = c(2, 0.4, 0))
    mtext(head_txt, side = 3, line = 1.2, adj = 0, font = 2, cex = 0.85, col = ink)
    mtext(sprintf("hypo in High: %d genes · hyper in High: %d genes · dashed = FDR 0.05",
                  n_hypo, n_hyper), side = 3, line = 0.3, adj = 0, cex = 0.6, col = ink2)
    mtext("← less methylated in mTOR-High      -log10 FDR      more methylated in mTOR-High →",
          side = 1, line = 2.2, cex = 0.6, col = ink2)
  }
  dev.off()
}
