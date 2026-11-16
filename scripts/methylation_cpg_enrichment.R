# CpG-cluster methylation (CCLE RRBS) for mTOR-High vs mTOR-Low (RNA groups):
#   1. agreement with the TSS-1kb summary for genes of interest
#   2. gene-level aggregation of CpG clusters (Sidak-corrected best cluster)
#   3. gene-set enrichment: over-representation (hypo / hyper in High) and a
#      rank-based test, on MSigDB Hallmark, GO BP, KEGG and PanglaoDB cell types
# Gene sets: data/gene_sets/*.gmt (Enrichr libraries). Base R only.
#
# Usage: Rscript scripts/methylation_cpg_enrichment.R   (run from project root)

res_dir <- "results/mTOR_RNA_groups"
fig_dir <- file.path(res_dir, "figures")
results <- readRDS(file.path(res_dir, "all_results.rds"))
source("scripts/plot_utils.R"); source("scripts/enrichment_utils.R")
cg <- results$meth_cpg; ts <- results$meth
cg$gene <- sub("\\|.*$", "", cg$feature); ts$gene <- sub("\\|.*$", "", ts$feature)

# ------------------------------------------- 1. genes of interest, both files
goi <- c("STAB1", "AZU1", "ELANE", "MPO", "FGR", "SIRPA", "TAL1", "GP1BB", "ERMAP", "GATA1",
         "KLF1", "PIK3CD", "PIK3CD-AS1", "NPRL3", "MAP1LC3A", "CEBPA", "IRF8", "SPI1")
best <- function(d, g) { x <- d[d$gene == g, ]; if (!nrow(x)) return(NULL); x[which.min(x$p_welch), ] }
cmp <- do.call(rbind, lapply(goi, function(g) {
  a <- best(ts, g); b <- best(cg, g)
  data.frame(gene = g,
             tss_diff = if (is.null(a)) NA else a$diff_High_minus_Low,
             tss_fdr  = if (is.null(a)) NA else a$fdr_welch,
             cpg_n_clusters = sum(cg$gene == g),
             cpg_best_cluster = if (is.null(b)) NA else sub("^.*\\|", "", b$feature),
             cpg_diff = if (is.null(b)) NA else b$diff_High_minus_Low,
             cpg_fdr  = if (is.null(b)) NA else b$fdr_welch,
             cpg_n_FDR0.25 = sum(cg$gene == g & cg$fdr_welch < 0.25))
}))
write.table(cmp, file.path(res_dir, "methylation_TSS_vs_CpG_genes_of_interest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# gene-level agreement between the two summaries (best region per gene)
bt <- ts[!duplicated(ts$gene), ]; bc <- cg[!duplicated(cg$gene), ]   # both sorted by p already
sh <- merge(bt[, c("gene", "t")], bc[, c("gene", "t")], by = "gene", suffixes = c("_tss", "_cpg"))

# ------------------------------------------------ 2. gene-level CpG scores --
k   <- table(cg$gene)
gl  <- bc[, c("gene", "feature", "diff_High_minus_Low", "p_welch", "fdr_welch")]
names(gl)[2:5] <- c("best_cluster", "best_diff", "best_p", "best_cluster_fdr")
gl$n_clusters <- as.integer(k[gl$gene])
gl$p_gene   <- 1 - (1 - gl$best_p)^gl$n_clusters        # Sidak: corrects for #clusters
gl$fdr_gene <- p.adjust(gl$p_gene, "BH")
gl$score    <- -log10(gl$p_gene) * sign(gl$best_diff)  # + = more methylated in mTOR-High
gl$direction <- ifelse(gl$best_diff > 0, "hyper_in_High", "hypo_in_High")
gl <- gl[order(gl$p_gene), ]
write.table(gl, file.path(res_dir, "methylation_CpG_gene_level.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

hit <- gl$fdr_gene < 0.25 & abs(gl$best_diff) >= 0.2
hypo  <- gl$gene[hit & gl$direction == "hypo_in_High"]
hyper <- gl$gene[hit & gl$direction == "hyper_in_High"]
universe <- gl$gene

# ------------------------------------------------------- 3. enrichment -----
score <- setNames(gl$score, gl$gene)
en <- run_enrichment(hypo, hyper, universe, score, file.path(res_dir, "methylation_CpG_enrichment"))

# ------------------------------------------------------------- summary ------
sink(file.path(res_dir, "methylation_CpG_summary.txt"))
cat("CpG clusters tested:", nrow(cg), "| FDR<0.05:", sum(cg$fdr_welch < 0.05),
    "| FDR<0.25:", sum(cg$fdr_welch < 0.25), "\n")
cat("  of FDR<0.25 clusters, more methylated in High:", sum(cg$fdr_welch < 0.25 & cg$diff_High_minus_Low > 0),
    "| less:", sum(cg$fdr_welch < 0.25 & cg$diff_High_minus_Low < 0), "\n")
cat("Genes:", nrow(gl), "| gene-level FDR<0.25 & |delta beta|>=0.2:", sum(hit),
    "(hypo in High", length(hypo), "| hyper in High", length(hyper), ")\n")
cat(sprintf("TSS-1kb vs CpG gene-level agreement (best region t, %d shared genes): Spearman rho = %.2f\n\n",
            nrow(sh), cor(sh$t_tss, sh$t_cpg, method = "spearman")))
cat("== Genes of interest: TSS-1kb vs CpG clusters\n")
print(transform(cmp, tss_diff = round(tss_diff, 2), tss_fdr = signif(tss_fdr, 2),
                cpg_diff = round(cpg_diff, 2), cpg_fdr = signif(cpg_fdr, 2)), row.names = FALSE)
cat("\n== Top genes (gene-level)\n")
print(head(transform(gl[, c("gene", "n_clusters", "best_diff", "p_gene", "fdr_gene")],
                     best_diff = round(best_diff, 2), p_gene = signif(p_gene, 2), fdr_gene = signif(fdr_gene, 2)), 40),
      row.names = FALSE)
print_enrichment(en)
sink()

# ---------------------------------------------------------------- figure ----
plot_enrichment(en$ORA, length(hypo), length(hyper), file.path(fig_dir, "methylation_CpG_enrichment.png"))
cat(readLines(file.path(res_dir, "methylation_CpG_summary.txt")), sep = "\n")
