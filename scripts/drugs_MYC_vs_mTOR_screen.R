# Drug response (AUC) MYC-High vs mTOR-High across all single-agent compounds in CTRP and GDSC2,
# and the subset whose annotated targets (PortalCompounds GeneSymbolOfTargets) include a concordant
# CRISPR/RNAi dependency gene (scripts/cross_CRISPR_RNAi_MYC_vs_mTOR.R). Lower AUC = more sensitive.
# "consistent" = the group more dependent on the target gene is also more sensitive to the drug.
#
# Run from the project root: Rscript scripts/drugs_MYC_vs_mTOR_screen.R

res_dir <- "results/mTOR_RNA_groups"
results <- readRDS(file.path(res_dir, "all_results.rds"))
cl <- read.delim(file.path(res_dir, "groups_RNA_mTOR.tsv"))
cl$grp <- ifelse(cl$MYC_group == "High" & cl$mTOR_group == "Low", "MYC-High",
          ifelse(cl$mTOR_group == "High" & cl$MYC_group == "Low", "mTOR-High", NA))
dep <- read.delim(file.path(res_dir, "cross_CRISPR_RNAi_MYC_vs_mTOR_concordant_hits.tsv"))

screen_test <- function(key, s) {
  r <- results[[key]]; X <- attr(r, "X")
  g <- cl$grp[match(rownames(X), cl$cell_line)]
  single <- !grepl(":", r$CompoundName) & !is.na(r$CompoundName) & r$CompoundName != ""
  do.call(rbind, lapply(which(single), function(i) {
    v <- X[, r$feature[i]]; a <- v[g %in% "MYC-High" & !is.na(v)]; b <- v[g %in% "mTOR-High" & !is.na(v)]
    if (length(a) < 4 || length(b) < 4 || var(c(a, b)) == 0) return(NULL)
    data.frame(screen = s, drug = r$CompoundName[i], targets = r$GeneSymbolOfTargets[i], mechanism = r$TargetOrMechanism[i],
               n_MYC = length(a), n_mTOR = length(b), AUC_MYC = mean(a), AUC_mTOR = mean(b),
               diff_MYC_minus_mTOR = mean(a) - mean(b), p_welch = t.test(a, b)$p.value)
  }))
}
all <- rbind(screen_test("drug_CTRP", "CTRP"), screen_test("drug_GDSC2", "GDSC2"))
all$fdr_within_screen <- ave(all$p_welch, all$screen, FUN = function(p) p.adjust(p, "BH"))
all$more_sensitive <- ifelse(all$diff_MYC_minus_mTOR < 0, "MYC-High", "mTOR-High")

# link to dependency genes
tg <- strsplit(ifelse(is.na(all$targets), "", all$targets), ";")
all$dep_targets <- vapply(tg, function(x) paste(intersect(x, dep$gene), collapse = ";"), "")
all$dep_more_dependent <- vapply(tg, function(x) paste(unique(dep$favours[dep$gene %in% x]), collapse = ";"), "")
all$consistent <- ifelse(all$dep_targets == "", NA, all$dep_more_dependent == all$more_sensitive)
all <- all[order(all$p_welch), ]
write.table(all, file.path(res_dir, "drugs_MYC_High_vs_mTOR_High_CTRP_GDSC2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

show <- function(d) print(data.frame(d$screen, drug = substr(d$drug, 1, 22), target = substr(ifelse(is.na(d$targets), "", d$targets), 1, 28),
  sens = d$more_sensitive, n = paste0(d$n_MYC, "/", d$n_mTOR), dAUC = round(d$diff_MYC_minus_mTOR, 3), p = signif(d$p_welch, 2),
  fdr = signif(d$fdr_within_screen, 2), dep = d$dep_targets, cons = d$consistent), row.names = FALSE)

cat(sprintf("Single-agent drugs tested: CTRP %d, GDSC2 %d | p<0.05: %d | FDR<0.25 within screen: %d\n",
            sum(all$screen == "CTRP"), sum(all$screen == "GDSC2"), sum(all$p_welch < 0.05), sum(all$fdr_within_screen < 0.25)))
cat("\n== Drugs targeting a concordant dependency gene\n"); show(all[all$dep_targets != "", ])
d <- all[all$dep_targets != "", ]
cat(sprintf("\n   %d drug-screen pairs; p<0.05: %d, of which consistent with the dependency direction: %d\n",
            nrow(d), sum(d$p_welch < 0.05), sum(d$p_welch < 0.05 & d$consistent)))
cat("\n== Top 30 drugs overall (any target)\n"); show(head(all, 30))
# drugs significant in both screens, same direction
both <- merge(all[all$screen == "CTRP", c("drug", "targets", "diff_MYC_minus_mTOR", "p_welch")],
              all[all$screen == "GDSC2", c("drug", "diff_MYC_minus_mTOR", "p_welch")], by = "drug", suffixes = c("_CTRP", "_GDSC2"))
both$same_dir <- sign(both$diff_MYC_minus_mTOR_CTRP) == sign(both$diff_MYC_minus_mTOR_GDSC2)
cat(sprintf("\n== In both screens: %d drugs; Spearman of dAUC = %.2f; p<0.05 in both & same direction:\n", nrow(both),
            cor(both$diff_MYC_minus_mTOR_CTRP, both$diff_MYC_minus_mTOR_GDSC2, method = "spearman")))
bb <- both[both$p_welch_CTRP < 0.05 & both$p_welch_GDSC2 < 0.05 & both$same_dir, ]
print(data.frame(drug = bb$drug, target = substr(bb$targets, 1, 30), dCTRP = round(bb$diff_MYC_minus_mTOR_CTRP, 3),
                 pCTRP = signif(bb$p_welch_CTRP, 2), dGDSC2 = round(bb$diff_MYC_minus_mTOR_GDSC2, 3), pGDSC2 = signif(bb$p_welch_GDSC2, 2)), row.names = FALSE)

# ---------------- adjusted for general drug sensitivity -----------------------
# mTOR-High lines are more sensitive to most drugs, so raw differences mix a drug-specific effect
# with general sensitivity. General sensitivity = each line's mean per-drug z-scored AUC (leaving the
# drug itself out). Test = ANCOVA, AUC ~ group + general sensitivity; adj_diff = the group coefficient
# (MYC-High minus mTOR-High at equal general sensitivity).
adj_test <- function(key, s) {
  r <- results[[key]]; X <- attr(r, "X"); Z <- scale(X)
  g <- cl$grp[match(rownames(X), cl$cell_line)]
  zs <- rowSums(Z, na.rm = TRUE); zn <- rowSums(!is.na(Z))
  single <- !grepl(":", r$CompoundName) & !is.na(r$CompoundName) & r$CompoundName != ""
  do.call(rbind, lapply(which(single), function(i) {
    f <- r$feature[i]; v <- X[, f]; zi <- Z[, f]
    gen <- (zs - ifelse(is.na(zi), 0, zi)) / (zn - !is.na(zi))     # leave-one-drug-out general sensitivity
    ok <- !is.na(v) & !is.na(g) & is.finite(gen)
    if (sum(g[ok] == "MYC-High") < 4 || sum(g[ok] == "mTOR-High") < 4 || var(v[ok]) == 0) return(NULL)
    d <- data.frame(v = v[ok], g = factor(g[ok], levels = c("mTOR-High", "MYC-High")), gen = gen[ok])
    cf <- summary(lm(v ~ g + gen, d))$coefficients["gMYC-High", ]
    data.frame(screen = s, drug = r$CompoundName[i], targets = r$GeneSymbolOfTargets[i],
               adj_diff_MYC_minus_mTOR = cf[["Estimate"]], p_adj = cf[["Pr(>|t|)"]])
  }))
}
adj <- rbind(adj_test("drug_CTRP", "CTRP"), adj_test("drug_GDSC2", "GDSC2"))
all <- merge(all, adj[, c("screen", "drug", "adj_diff_MYC_minus_mTOR", "p_adj")], by = c("screen", "drug"), all.x = TRUE)
all$fdr_adj_within_screen <- ave(all$p_adj, all$screen, FUN = function(p) p.adjust(p, "BH"))
all$more_sensitive_adj <- ifelse(all$adj_diff_MYC_minus_mTOR < 0, "MYC-High", "mTOR-High")
all <- all[order(all$p_adj), ]
write.table(all, file.path(res_dir, "drugs_MYC_High_vs_mTOR_High_CTRP_GDSC2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

showa <- function(d) print(data.frame(d$screen, drug = substr(d$drug, 1, 22), target = substr(ifelse(is.na(d$targets), "", d$targets), 1, 26),
  sens = d$more_sensitive_adj, raw_dAUC = round(d$diff_MYC_minus_mTOR, 3), p_raw = signif(d$p_welch, 2),
  adj_dAUC = round(d$adj_diff_MYC_minus_mTOR, 3), p_adj = signif(d$p_adj, 2), fdr_adj = signif(d$fdr_adj_within_screen, 2), dep = d$dep_targets), row.names = FALSE)
cat(sprintf("\n== Adjusted for general sensitivity: p<0.05 %d (MYC-High %d / mTOR-High %d) | FDR<0.25: %d\n", sum(all$p_adj < 0.05, na.rm = TRUE),
            sum(all$p_adj < 0.05 & all$more_sensitive_adj == "MYC-High", na.rm = TRUE), sum(all$p_adj < 0.05 & all$more_sensitive_adj == "mTOR-High", na.rm = TRUE),
            sum(all$fdr_adj_within_screen < 0.25, na.rm = TRUE)))
cat("-- more active in MYC-High (adjusted), top 15\n"); showa(head(all[all$more_sensitive_adj == "MYC-High", ], 15))
cat("-- more active in mTOR-High (adjusted), top 15\n"); showa(head(all[all$more_sensitive_adj == "mTOR-High", ], 15))
cat("-- drugs targeting a concordant dependency gene (adjusted)\n"); showa(all[all$dep_targets != "", ])
