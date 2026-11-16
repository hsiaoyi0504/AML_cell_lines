# CCLE targeted metabolomics (225 metabolites, log10 abundance) in MYC-High/mTOR-Low vs
# mTOR-High/MYC-Low AML lines (RNA-defined groups; lines high or low for both scores left out).
# Positive difference = higher in MYC-High. Volcano significance at FDR < 0.1; labels = all FDR < 0.25.
#
# Run from the project root: Rscript scripts/metabolomics_MYC_vs_mTOR.R

res_dir <- "results/mTOR_RNA_groups"; fig_dir <- file.path(res_dir, "figures")
source("scripts/plot_utils.R")
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
X  <- attr(readRDS(file.path(res_dir, "all_results.rds"))$metab, "X")
g  <- cl[match(rownames(X), cl$cell_line), ]
grp <- ifelse(g$MYC_group == "High" & g$mTOR_group == "Low", "MYC-High",
       ifelse(g$mTOR_group == "High" & g$MYC_group == "Low", "mTOR-High", NA))
X <- X[!is.na(grp), ]; grp <- grp[!is.na(grp)]; a <- grp == "MYC-High"
colA <- "#eb6834"; colB <- "#4a3aa7"

st <- do.call(rbind, lapply(colnames(X), function(m) {
  x <- X[a, m]; y <- X[!a, m]; x <- x[!is.na(x)]; y <- y[!is.na(y)]
  ok <- length(x) >= 3 && length(y) >= 3 && var(c(x, y)) > 0
  data.frame(metabolite = m, n_MYC = length(x), n_mTOR = length(y), mean_MYC = mean(x), mean_mTOR = mean(y),
             diff_MYC_minus_mTOR = mean(x) - mean(y),
             p_welch = if (ok) t.test(x, y)$p.value else NA,
             p_wilcox = if (ok) suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value) else NA)
}))
st <- st[!is.na(st$p_welch), ]
st$fdr_welch <- p.adjust(st$p_welch, "BH"); st$fdr_wilcox <- p.adjust(st$p_wilcox, "BH")
st <- st[order(st$p_welch), ]
write.table(st, file.path(res_dir, "metabolomics_MYC_High_vs_mTOR_High.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("MYC-High n=%d · mTOR-High n=%d · %d metabolites · FDR<0.1: %d · FDR<0.25: %d\n",
            sum(a), sum(!a), nrow(st), sum(st$fdr_welch < 0.1), sum(st$fdr_welch < 0.25)))
print(data.frame(st$metabolite, diff = round(st$diff_MYC_minus_mTOR, 3), p = signif(st$p_welch, 2),
                 fdr = signif(st$fdr_welch, 2))[1:12, ], row.names = FALSE)

vr <- data.frame(feature = st$metabolite, diff_High_minus_Low = st$diff_MYC_minus_mTOR, p_welch = st$p_welch,
                 fdr_welch = st$fdr_welch, n_High = st$n_MYC, n_Low = st$n_mTOR)
volcano_png(vr, list(title = "Metabolomics: MYC-High vs mTOR-High",
                     x = "log10 abundance difference (MYC-High − mTOR-High)", dep = FALSE, fdr = 0.1, always = "^acetylcarnitine$",
                     groups = c("MYC-High", "mTOR-High"), cols = c(colA, colB)),
            file.path(fig_dir, "volcano_metabolomics_MYC_High_vs_mTOR_High.png"), n_lab = sum(st$fdr_welch < 0.25), xpad = 2.2)
