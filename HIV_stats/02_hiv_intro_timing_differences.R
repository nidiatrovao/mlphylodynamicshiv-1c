# ============================================================================
# Benchmarking Analysis for Introduction Timing
# HIV-1 Subtype C Phylodynamics Study
# Amanda Perofsky
# ============================================================================

# Use mixed effect models to test whether phylogenetic and temporal dating
# methods systematically shift introduction time estimates earlier or later.
# Dependent variable: introduction time (decimal year)
# Fixed effects: temporal dating method (dpt), phylogenetic method (mlt)
# Random effects: origin country, destination country, country pair, dataset
# ============================================================================
# 1. Load and prepare data
# ============================================================================
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(performance)
library(sjPlot)

d <- read.delim("input_data/intro_count_time.tsv")
head(d)

# ISO country code to name lookup
country_names <- c(
  AR = "Argentina", BE = "Belgium", BG = "Bulgaria", BR = "Brazil",
  BW = "Botswana", CN = "China", CY = "Cyprus", DE = "Germany",
  DK = "Denmark", ES = "Spain", ET = "Ethiopia", GB = "United Kingdom",
  GE = "Georgia", IL = "Israel", IN = "India", KE = "Kenya",
  MM = "Myanmar", MW = "Malawi", NP = "Nepal", NG = "Nigeria",
  PK = "Pakistan", SE = "Sweden", SN = "Senegal", SO = "Somalia",
  TH = "Thailand", TZ = "Tanzania", US = "United States",
  UY = "Uruguay", YE = "Yemen", ZA = "South Africa", ZM = "Zambia"
)

# Helper to convert pair codes to readable names
pair_label <- function(pair) {
  parts <- strsplit(pair, "->")[[1]]
  from_name <- ifelse(parts[1] %in% names(country_names), country_names[parts[1]], parts[1])
  to_name <- ifelse(parts[2] %in% names(country_names), country_names[parts[2]], parts[2])
  sprintf("%s -> %s", from_name, to_name)
}

d <- d %>%
  mutate(
    country_pair = paste0(from, "->", to),
    dataset_base = gsub("\\.\\d+$", "", dataset),
    replicate = as.integer(gsub(".*\\.", "", dataset))
  )

cat("Total introduction events:", nrow(d), "\n")

# ============================================================================
# 2. Identify shared country pairs (within-dataset filtering)
# ============================================================================
# For each dataset, keep country pairs present in all 16 method combos.
# A pair does not need to appear in every dataset — smaller subsamples
# naturally detect fewer routes. The mixed model handles the unbalanced
# structure through the random effects.
# ----------------------------------------------------------------------------

d <- d %>%
  mutate(method_combo = paste(dpt, mlt, sep = "_"))

total_combos <- n_distinct(d$method_combo)
total_datasets <- n_distinct(d$dataset)

# For each dataset, find pairs with all 16 method combos
complete_within_dataset <- d %>%
  group_by(dataset, country_pair) %>%
  summarise(n_combos = n_distinct(method_combo), .groups = "drop") %>%
  filter(n_combos == total_combos)

# Filter to keep only complete pair-dataset combinations
d_shared <- d %>%
  semi_join(complete_within_dataset, by = c("dataset", "country_pair"))

## FILTERING
# Keep pairs present in all 16 method combos
# Pairs don't need to appear in every dataset though
cat("Unique country pairs retained:", n_distinct(d_shared$country_pair), "\n")
cat("Rows after filtering:", nrow(d_shared), "\n")
cat("Datasets:", n_distinct(d_shared$dataset), "\n")
cat("Origin countries:", n_distinct(d_shared$from), "\n")
cat("Destination countries:", n_distinct(d_shared$to), "\n\n")

# How many datasets each pair appears in
pair_dataset_counts <- d_shared %>%
  distinct(dataset, country_pair) %>%
  count(country_pair, name = "n_datasets") %>%
  arrange(desc(n_datasets))

# Pairs by number of datasets (complete in all 16 combos)
pair_dataset_counts %>%
  count(n_datasets, name = "n_pairs") %>%
  arrange(desc(n_datasets)) %>%
  mutate(label = sprintf("  %2d/15 datasets: %d pairs", n_datasets, n_pairs)) %>%
  pull(label) %>%
  cat(sep = "\n")

# Pairs per dataset (shows effect of subsample size)
d_shared %>%
  distinct(dataset, country_pair) %>%
  count(dataset, name = "n_pairs") %>%
  arrange(dataset) %>%
  mutate(label = sprintf("  %s: %d pairs", dataset, n_pairs)) %>%
  pull(label) %>%
  cat(sep = "\n")

# ============================================================================
# 3. Descriptive statistics
# ============================================================================

# Introduction time by temporal method
d_shared %>%
  group_by(dpt) %>%
  summarise(
    n = n(),
    mean = round(mean(time), 2),
    median = round(median(time), 2),
    sd = round(sd(time), 2),
    .groups = "drop"
  ) %>%
  print()

# Introduction time by phylogenetic method
d_shared %>%
  group_by(mlt) %>%
  summarise(
    n = n(),
    mean = round(mean(time), 2),
    median = round(median(time), 2),
    sd = round(sd(time), 2),
    .groups = "drop"
  ) %>%
  print()

# Mean introduction time by method combination
d_shared %>%
  group_by(dpt, mlt) %>%
  summarise(mean_time = round(mean(time), 1), .groups = "drop") %>%
  pivot_wider(names_from = mlt, values_from = mean_time) %>%
  print()

# ============================================================================
# 4. Primary mixed effects model
# ============================================================================
# Full model: (1|from) + (1|to) + (1|country_pair) + (1|dataset)
# Separates variance due to origin country, destination country,
# route-specific effects, and subsampling replicate.
# ----------------------------------------------------------------------------

# Random effect levels
cat("  from (origin countries):", n_distinct(d_shared$from), "\n")
cat("  to (destination countries):", n_distinct(d_shared$to), "\n")
cat("  country_pair:", n_distinct(d_shared$country_pair), "\n")
cat("  dataset:", n_distinct(d_shared$dataset), "\n\n")

# Full model with phylo x temporal interaction
mod_full <- lmer(time ~ dpt * mlt + (1 | from) + (1 | to) + (1 | country_pair) + (1 | dataset),
  data = d_shared
)
print(summary(mod_full))

# Check convergence
isSingular(mod_full)

# print regression output as html table
tab_model(mod_full)

# Type III ANOVA: test phylo x temporal interaction
anova_full <- anova(mod_full, type = 3)
print(anova_full)
# Temporal dating method significantly influences the timing of introduction events (F=108.55, p<0.0001),
# while phylogenetic method has a small but statistically  significant effect (F=3.05, p=0.027)
# with no interaction between the two factors (p=0.18).

# Variance components
print(VarCorr(mod_full))

vc <- as.data.frame(VarCorr(mod_full))
total_var <- sum(vc$vcov)
for (i in 1:nrow(vc)) {
  cat(sprintf(
    "  %-15s %8.2f (%4.1f%%)\n",
    vc$grp[i], vc$vcov[i], 100 * vc$vcov[i] / total_var
  ))
}

# Additive model (no interaction between temporal x phylo method)
mod_add <- lmer(time ~ dpt + mlt + (1 | from) + (1 | to) + (1 | country_pair) + (1 | dataset),
  data = d_shared
)

# Additive model ANOVA
print(anova(mod_add, type = 3))

# Method interaction test (additive vs full)
print(anova(mod_add, mod_full, refit = TRUE))
# Interaction term doesn't improve the model.
# The additive model (dpt + mlt, no interaction) fits just as well as the full model (dpt * mlt),
# with a non-significant likelihood ratio test (p=0.18) and a lower AIC (214264 vs 214269)

# ============================================================================
# 5. Residual diagnostics for primary model
# ============================================================================
resid_vals <- residuals(mod_full)

# Skewness and kurtosis
n_resid <- length(resid_vals)
m_resid <- mean(resid_vals)
s_resid <- sd(resid_vals)
skew_resid <- sum((resid_vals - m_resid)^3) / (n_resid * s_resid^3)
kurt_resid <- sum((resid_vals - m_resid)^4) / (n_resid * s_resid^4)
cat(sprintf("Residual skewness: %.3f\n", skew_resid))
cat(sprintf("Residual kurtosis: %.3f\n", kurt_resid))

# Heteroscedasticity check
d_shared$resid <- resid_vals

# Residual SD by temporal method
d_shared %>%
  group_by(dpt) %>%
  summarise(resid_sd = round(sd(resid), 2), .groups = "drop") %>%
  arrange(-resid_sd) %>%
  print()

# Residual SD by phylogenetic method
d_shared %>%
  group_by(mlt) %>%
  summarise(resid_sd = round(sd(resid), 2), .groups = "drop") %>%
  arrange(-resid_sd) %>%
  print()

sds_dpt <- d_shared %>%
  group_by(dpt) %>%
  summarise(s = sd(resid), .groups = "drop") %>%
  pull(s)

sds_mlt <- d_shared %>%
  group_by(mlt) %>%
  summarise(s = sd(resid), .groups = "drop") %>%
  pull(s)

# Ratios < 2 are generally acceptable
cat(sprintf("\nMax/min SD ratio across temporal methods: %.2f\n", max(sds_dpt) / min(sds_dpt)))
cat(sprintf("Max/min SD ratio across phylo methods: %.2f\n", max(sds_mlt) / min(sds_mlt)))

# Residuals show moderate left-skew consistent with the observed temporal concentration of introduction events. # Fixed effect estimates are still robust to this departure given the large sample size.

# ============================================================================
# 6. Post-hoc comparisons of method effects
# ============================================================================
# Test whether each pair of temporal methods (or phylogenetic methods) differs
# significantly from each other, after accounting for the random effects.
# ----------------------------------------------------------------------------
# Note: Asymptotic df method is appropriate given large sample size (n=24,781);
# Kenward-Roger/Satterthwaite adjustments unnecessary

# Estimated marginal means for temporal methods
emm_dpt <- emmeans(mod_full, "dpt")
print(summary(emm_dpt))

# Pairwise comparisons for temporal methods (Bonferroni correction)
pairs_dpt <- pairs(emm_dpt, adjust = "bonferroni")
print(pairs_dpt) # almost all temporal method pairs are significantly different

# Estimated marginal means for phylogenetic methods
emm_mlt <- emmeans(mod_full, "mlt")
print(summary(emm_mlt))

# Pairwise comparisons for phylo methods (Bonferroni correction)
pairs_mlt <- pairs(emm_mlt, adjust = "bonferroni")
print(pairs_mlt) # Only FastTree and IQ-TREE are significantly different

# ============================================================================
# 7. Sensitivity check: Alternative random effects structures (simpler models)
# ============================================================================

# Model with country_pair + dataset only (no from/to)
mod_simple <- lmer(time ~ dpt * mlt + (1 | country_pair) + (1 | dataset),
  data = d_shared
)

# Singularity check
cat(sprintf("  Singular: %s\n", isSingular(mod_simple)))
# Variance components
print(VarCorr(mod_simple))

# Model with from + to + dataset only (no country_pair)
mod_no_pair <- lmer(time ~ dpt * mlt + (1 | from) + (1 | to) + (1 | dataset),
  data = d_shared
)

# Convergence check
isSingular(mod_no_pair)

# Variance components
print(VarCorr(mod_no_pair))

# Model performance comparison
print(compare_performance(mod_full, mod_simple, mod_no_pair))
# mod_full: best by AIC (0.975 weight), R2 cond = 0.371
# mod_simple: best by BIC (0.988 weight), R2 cond = 0.333
# mod_no_pair: worst on all criteria except R2 cond (0.466, inflated — see note below)

# Note: mod_no_pair has higher R2 cond because from/to absorb more variance without
# country_pair, but RMSE is worse (18.332 vs 18.120) and AIC is worst of all three.
# The higher R2 is misleading — the model actually predicts less well.

# Likelihood ratio tests
# full model (from + to + country_pair + dataset) vs simple model (country_pair + dataset only)
print(anova(mod_simple, mod_full, refit = TRUE))
# Does adding from/to directionality improve fit? yes

# full model (from + to + country_pair + dataset) vs no pair (from + to + dataset only)
print(anova(mod_no_pair, mod_full, refit = TRUE))
# Is country_pair needed beyond from/to? yes

# Compare fixed effect conclusions: primary vs simple vs no pair
anova_full <- anova(mod_full, type = 3)
print(anova_full)
anova_simple <- anova(mod_simple, type = 3)
print(anova_simple)
anova_no_pair <- anova(mod_no_pair, type = 3)
print(anova_no_pair)

# Conclusions identical across all three specifications:
# - Temporal method: always p<0.0001
# - Phylogenetic method: always significant (p=0.016-0.028)
# - Interaction: always non-significant (p~0.18)

# ============================================================================
# 6. Sensitivity check:
# dataset_base (subsample size) vs dataset (replicate identity)
# ============================================================================
# The primary model uses the specific dataset (e.g., locrisk260.1) as a random effect,
# accounting for replicate-specific variation. Here we check whether collapsing to
# dataset_base (subsample dataset size only, 5 levels instead of 15) changes any conclusions.
# ----------------------------------------------------------------------------

# Model with dataset_base (subsample size)
mod_base <- lmer(time ~ dpt * mlt + (1 | from) + (1 | to) + (1 | country_pair) + (1 | dataset_base),
  data = d_shared
)

# Convergence check
isSingular(mod_base)

# Type III ANOVA: test interaction between temporal x phylo methods
anova_base <- anova(mod_base, type = 3)
print(anova_base)

# Variance components
vc_base <- as.data.frame(VarCorr(mod_base))
total_var_base <- sum(vc_base$vcov)
for (i in 1:nrow(vc_base)) {
  cat(sprintf(
    "  %-15s %8.2f (%4.1f%%)\n",
    vc_base$grp[i], vc_base$vcov[i], 100 * vc_base$vcov[i] / total_var_base
  ))
}

# Side-by-side comparison: dataset vs dataset_base
compare_performance(mod_full, mod_base)
# Full model (considers specific dataset) fits better than collapsing to dataset size.

# Variance comparison: dataset vs data_base random effect
dataset_var_primary <- vc$vcov[vc$grp == "dataset"]
dataset_var_base <- vc_base$vcov[vc_base$grp == "dataset_base"]
cat(sprintf("  dataset:      %.4f (%.2f%%)\n", dataset_var_primary, 100 * dataset_var_primary / total_var))
cat(sprintf("  dataset_base: %.4f (%.2f%%)\n", dataset_var_base, 100 * dataset_var_base / total_var_base))
# Collapsing to (1|dataset_base) produces identical inferential conclusions.
# All fixed effect significance patterns and post-hoc contrasts are unchanged.

# ============================================================================
# 9. Sensitivity check: Rank-transformed mixed effect model
# ============================================================================
# Rank-transforming introduction times to eliminate distributional assumptions while
# preserving the random effects structure.
# ------------------------------------------------------------------------------

d_shared$rank_time <- rank(d_shared$time)

mod_rank <- lmer(rank_time ~ dpt * mlt + (1 | from) + (1 | to) + (1 | country_pair) + (1 | dataset),
  data = d_shared
)

# Type III ANOVA (rank-transformed)
anova_rank <- anova(mod_rank, type = 3)
print(anova_rank)

# Conclusions are consistent across original and rank-transformed models

# ============================================================================
# 8. Sensitivity check: Kruskal-Wallis tests (ignore random effects)
# ============================================================================
# By temporal method
kw_dpt <- kruskal.test(time ~ dpt, data = d_shared)
print(kw_dpt)

# By phylogenetic method
kw_mlt <- kruskal.test(time ~ mlt, data = d_shared)
print(kw_mlt)

# Non-parametric tests confirm the mixed effect model results:
# temporal method has bigger effect on introduction count estimates than phylogenetic method

# ============================================================================
# 9. Diving in to temporal method effects on intro counts
# ============================================================================
# Do temporal methods diverge equally across all time periods, or is the effect
# concentrated in early vs recent introduction events?
# ----------------------------------------------------------------------------

d_shared$period <- cut(d_shared$time,
  breaks = c(-Inf, 1980, 2000, Inf),
  labels = c("pre-1980", "1980-2000", "2000-2017")
)

# Density plots by method
# Temporal methods: shows TempEst shifting first wave earlier and filling the valley between waves,
# vs bimodal pattern in LSD/TreeTime/treedater
p_dpt_density <- ggplot(d_shared) +
  geom_density(aes(x = time)) +
  facet_wrap(~dpt) +
  labs(
    x = "Introduction time", y = "Density",
    title = "Distribution of introduction times by temporal method"
  ) +
  theme_minimal()
print(p_dpt_density)

# Phylogenetic methods: distributions nearly identical
p_mlt_density <- ggplot(d_shared) +
  geom_density(aes(x = time)) +
  facet_wrap(~mlt) +
  labs(
    x = "Introduction time", y = "Density",
    title = "Distribution of introduction times by phylogenetic method"
  ) +
  theme_minimal()
print(p_mlt_density)

# Event counts by period and temporal method
d_shared %>%
  group_by(dpt, period) %>%
  tally() %>%
  arrange(period, n) %>%
  print()
# TempEst places substantially more events in the 1980-2000 period (1283) and
# fewer in the 2000-2017 period (3326) compared to the other methods.
# LSD, TreeTime, and treedater are fairly consistent (~910-950 in 1980-2000, ~3700-3900 in 2000-2017).
# TempEst shifts events earlier and fills the trough between the two waves.

# Proportion of events by period and temporal method
period_props <- d_shared %>%
  count(dpt, period) %>%
  group_by(dpt) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  select(-n) %>%
  pivot_wider(names_from = period, values_from = pct)
print(period_props)

# Mean/median intro time by period and temporal method
period_stats <- d_shared %>%
  group_by(period, dpt) %>%
  summarise(
    n = n(),
    mean = round(mean(time), 1),
    median = round(median(time), 1),
    .groups = "drop"
  )

for (p in c("pre-1980", "1980-2000", "2000-2017")) {
  cat(sprintf("\n  %s:\n", p))
  sub <- period_stats %>% filter(period == p)
  for (i in 1:nrow(sub)) {
    cat(sprintf(
      "    %-10s n = %5d, mean = %.1f, median = %.1f\n",
      sub$dpt[i], sub$n[i], sub$mean[i], sub$median[i]
    ))
  }
}
# pre-1980: TempEst mean=1951.5 vs treedater mean=1966.3 (15-year gap)
# 1980-2000: methods agree within ~1 year
# 2000-2017: methods agree within ~0.5 years

# Methods largely agree on recent introductions (post-2000)
# The temporal method effect is driven by early divergence events (pre-1980),
# where TempEst estimates substantially earlier dates than other methods.

# Chi-squared test: do methods differ in the distribution of events across periods?
tbl <- table(d_shared$dpt, d_shared$period)
chi_result <- chisq.test(tbl)
print(chi_result)

# Take Home: Temporal dating methods differ not only in mean
# introduction time but also in the temporal distribution of events.
# LSD, TreeTime, and treedater infer two distinct waves of introduction
# separated by a clear valley (~1975-1985), suggesting discrete seeding
# events. TempEst shifts the first wave earlier and fills in the valley,
# producing a pattern of more continuous dispersal. This difference in
# inferred temporal structure has implications for interpreting epidemic
# history: the biological narrative (discrete waves vs sustained
# dispersal) depends on the temporal dating method used.

# ============================================================================
# 10. Per-pair (route specific) differences across methods
# ============================================================================
# Per-pair temporal method differences
# Median introduction time difference per pair:
pair_dpt_medians <- d_shared %>%
  group_by(country_pair, dpt) %>%
  summarise(median_time = median(time), .groups = "drop") %>%
  pivot_wider(names_from = dpt, values_from = median_time)

pair_n <- d_shared %>% count(country_pair, name = "total_events")

pair_diff <- pair_dpt_medians %>%
  mutate(diff_TempEst_treedater = TempEst - treedater) %>%
  left_join(pair_n, by = "country_pair") %>%
  arrange(diff_TempEst_treedater)

# Focus on TempEst vs treedater because these two methods show the largest overall contrast
# Pairs with > 5 year difference:
big_diff <- pair_diff %>% filter(abs(diff_TempEst_treedater) > 5)
for (i in 1:nrow(big_diff)) {
  cat(sprintf(
    "  %-35s TempEst: %.1f  treedater: %.1f  diff: %+.1f yrs  (n=%d)\n",
    pair_label(big_diff$country_pair[i]),
    big_diff$TempEst[i],
    big_diff$treedater[i],
    big_diff$diff_TempEst_treedater[i],
    big_diff$total_events[i]
  ))
}

n_earlier <- sum(pair_diff$diff_TempEst_treedater < -0.01, na.rm = TRUE)
n_same <- sum(abs(pair_diff$diff_TempEst_treedater) <= 0.01, na.rm = TRUE)
n_total <- nrow(pair_diff)
cat(sprintf("\nPairs where TempEst is earlier: %d/%d\n", n_earlier, n_total))
cat(sprintf("Pairs with no difference: %d/%d\n", n_same, n_total))
cat(sprintf(
  "Mean difference across all pairs: %.1f years\n",
  mean(pair_diff$diff_TempEst_treedater, na.rm = TRUE)
))

# The temporal method effect on specific routes is widespread but variable.
# Most pairs show TempEst estimating earlier introductions, but the magnitude ranges
# from negligible (0 years) to extreme (>30 years).
# The largest discrepancies involve routes with early inferred introduction times.

# Phylogenetic method effect
# Median introduction time difference per pair:
pair_mlt_medians <- d_shared %>%
  group_by(country_pair, mlt) %>%
  summarise(median_time = median(time), .groups = "drop") %>%
  pivot_wider(names_from = mlt, values_from = median_time)

pair_mlt_diff <- pair_mlt_medians %>%
  mutate(diff_FT_IQ = FastTree - `IQ-TREE`) %>%
  left_join(pair_n, by = "country_pair") %>%
  arrange(diff_FT_IQ)

# Focus on FastTree vs IQ-TREE because these two methods show the only significant pairwise contrast
big_mlt <- pair_mlt_diff %>% filter(abs(diff_FT_IQ) > 2)
if (nrow(big_mlt) > 0) {
  cat("Pairs with > 2 year difference:\n")
  for (i in 1:nrow(big_mlt)) {
    cat(sprintf(
      "  %-35s FastTree: %.1f  IQ-TREE: %.1f  diff: %+.1f yrs  (n=%d)\n",
      pair_label(big_mlt$country_pair[i]),
      big_mlt$FastTree[i],
      big_mlt$`IQ-TREE`[i],
      big_mlt$diff_FT_IQ[i],
      big_mlt$total_events[i]
    ))
  }
}

n_ft_earlier <- sum(pair_mlt_diff$diff_FT_IQ < -0.01, na.rm = TRUE)
cat(sprintf("\nPairs where FastTree is earlier: %d/%d\n", n_ft_earlier, n_total))
cat(sprintf(
  "Mean difference across all pairs: %.2f years\n",
  mean(pair_mlt_diff$diff_FT_IQ, na.rm = TRUE)
))
cat(sprintf(
  "Median difference across all pairs: %.2f years\n",
  median(pair_mlt_diff$diff_FT_IQ, na.rm = TRUE)
))

# The phylogenetic method effect shows no consistent directional bias across routes.
# The mean and median per-pair differences are near zero, and the statistically
# significant pairs are split in direction (some earlier, some later).
# This suggests the effect is driven by route-specific tree topology differences
# rather than a systematic methodological bias.

# ============================================================================
# 11. Summary
# ============================================================================
# Temporal dating method significantly affects introduction time estimates,
# driven primarily by disagreement on early divergence events (pre-1980) where
# TempEst estimates substantially earlier dates.
# Methods largely converge for recent introductions (post-2000).

# Notably, methods differ in the inferred temporal structure:
# LSD, TreeTime, and treedater show two distinct introduction waves, while
# TempEst infers more continuous dispersal (chi-squared, p<0.0001).

# Phylogenetic method reached statistical significance but shows no consistent
# directional bias — the effect is route-specific rather than systematic.

# No interaction between phylogenetic x temporal methods.
