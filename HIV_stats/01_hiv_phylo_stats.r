# ============================================================================
# Differences in TMRCA, ER, and introduction counts across phylo and temporal methods
# HIV-1 Subtype C Phylodynamics Study
# Amanda Perofsky
# ============================================================================

library(tidyverse)
library(rstatix)
library(car)
library(PMCMRplus)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(effectsize)
library(rcompanion) # For Scheirer-Ray-Hare test
library(psych) # For ICC calculations
library(moments)

# Set output directory
dir.create("figures", showWarnings = FALSE)

# ============================================================================
# PART 1: TMRCA AND EVOLUTIONARY RATE ANALYSIS
# ============================================================================

cat("Loading TMRCA and evolutionary rate data...\n")
tmrca_data <- read.csv("input_data/HIV_TMRCA_ER_by_method.csv")
tmrca_data$Method <- as.factor(tmrca_data$Method)
levels(tmrca_data$Method)[3] <- "treedater" # rename
tmrca_data$Method <- as.character(tmrca_data$Method)

# Clean column names
names(tmrca_data) <- str_replace_all(names(tmrca_data), "\\.", "_")
head(tmrca_data)

# Convert TMRCA dates to numeric years
tmrca_data <- tmrca_data %>%
  mutate(
    TMRCA_FastTree_year = as.numeric(format(as.Date(TMRCA_FastTree), "%Y")) +
      as.numeric(format(as.Date(TMRCA_FastTree), "%j")) / 365,
    TMRCA_IQ_TREE_year = as.numeric(format(as.Date(TMRCA_IQ_TREE), "%Y")) +
      as.numeric(format(as.Date(TMRCA_IQ_TREE), "%j")) / 365,
    TMRCA_PhyML_year = as.numeric(format(as.Date(TMRCA_PhyML), "%Y")) +
      as.numeric(format(as.Date(TMRCA_PhyML), "%j")) / 365,
    TMRCA_RAxML_NG_year = as.numeric(format(as.Date(TMRCA_RAxML_NG), "%Y")) +
      as.numeric(format(as.Date(TMRCA_RAxML_NG), "%j")) / 365
  )

# Reshape to long format

## TMRCA
tmrca_long <- tmrca_data %>%
  pivot_longer(
    cols = ends_with("_year"),
    names_to = "phylo_method",
    values_to = "TMRCA_year"
  ) %>%
  mutate(
    phylo_method = str_remove(phylo_method, "TMRCA_"),
    phylo_method = str_remove(phylo_method, "_year"),
    phylo_method = str_replace(phylo_method, "IQ_TREE", "IQ-TREE"),
    phylo_method = str_replace(phylo_method, "RAxML_NG", "RAxML-NG"),
    # Extract subsampling information
    replicate = as.numeric(str_extract(Dataset, "\\d+$")),
    dataset_base = str_extract(Dataset, "locrisk\\d+")
  ) %>%
  rename(temporal_method = Method, sample_size = N_genomes)

## ER
er_long <- tmrca_data %>%
  pivot_longer(
    cols = starts_with("ER_"),
    names_to = "phylo_method",
    values_to = "ER"
  ) %>%
  mutate(
    phylo_method = str_remove(phylo_method, "ER_"),
    phylo_method = str_replace(phylo_method, "IQ.TREE", "IQ-TREE"),
    phylo_method = str_replace(phylo_method, "RAxML.NG", "RAxML-NG"),
    # Extract subsampling information
    replicate = as.numeric(str_extract(Dataset, "\\d+$")),
    dataset_base = str_extract(Dataset, "locrisk\\d+")
  ) %>%
  rename(temporal_method = Method, sample_size = N_genomes)

# ------------------------------------------
# Analysis 1: Compare temporal methods for TMRCA
# ------------------------------------------
# Summary statistics

# Compare temporal methods
tmrca_summary <- tmrca_long %>%
  group_by(temporal_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(TMRCA_year), 2),
    SD = round(sd(TMRCA_year), 2),
    Median = round(median(TMRCA_year), 2),
    IQR = round(IQR(TMRCA_year), 2),
    Q1 = round(quantile(TMRCA_year, 0.25), 2),
    Q3 = round(quantile(TMRCA_year, 0.75), 2),
    Min = round(min(TMRCA_year), 2),
    Max = round(max(TMRCA_year), 2)
  )

print(tmrca_summary)

# Notes:
# 60 = 4 ML methods x 5 sequence datasets x 3 replicates
# tight ranges for TMRCA: IQR = 5-11 years

# 18 year spread for TMRCA:
# TreeTime gives earliest estimates (~1929), then TempEst (~1935), then LSD (~1941), then treedater (~1947)

# Kruskal-Wallis test
kw_tmrca <- kruskal.test(TMRCA_year ~ temporal_method, data = tmrca_long)
print(kw_tmrca)

# Result: significant differences between temporal methods

# Post-hoc Dunn's test - pairwise comparisons
dunn_tmrca <- dunn_test(
  TMRCA_year ~ temporal_method,
  data = tmrca_long,
  p.adjust.method = "bonferroni"
)
print(dunn_tmrca)

# Result: all pairwise comparisons are significant after Bonferroni correction

# Effect size (eta or epsilon squared for Kruskal-Wallis)
# Note: eta-squared is equivalent to adjusted R-squared for ANOVA on ranks (confusingly called epsilon)

## eta-squared
eta_tmrca <- kruskal_effsize(TMRCA_year ~ temporal_method, data = tmrca_long)
print(eta_tmrca)

# epsilon-squared
epsilonSquared(
  x = tmrca_long$TMRCA_year,
  g = tmrca_long$temporal_method,
  ci = T
)

# Result: both eta-squared and epsilon-squared indicate large effect size

# ------------------------------------------
# Analysis 2: Compare phylogenetic methods for TMRCA
# ------------------------------------------

tmrca_phylo_summary <- tmrca_long %>%
  group_by(phylo_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(TMRCA_year), 2),
    SD = round(sd(TMRCA_year), 2),
    Median = round(median(TMRCA_year), 2),
    IQR = round(IQR(TMRCA_year), 2)
  )
print(tmrca_phylo_summary)
# Result: all 4 methods give similar estimates: 1935-1939

kw_tmrca_phylo <- kruskal.test(TMRCA_year ~ phylo_method, data = tmrca_long)
print(kw_tmrca_phylo)
# Result: phylogenetic methods do NOT significantly affect TMRCA estimates

if (kw_tmrca_phylo$p.value < 0.05) {
  dunn_tmrca_phylo <- dunn_test(
    TMRCA_year ~ phylo_method,
    data = tmrca_long,
    p.adjust.method = "bonferroni"
  )
  print(dunn_tmrca_phylo)
}

# eta-squared statistic for KW test
eta_tmrca <- kruskal_effsize(TMRCA_year ~ phylo_method, data = tmrca_long)
print(eta_tmrca)
# Result: phylo method has minimal effect on TMRCA

# Takehomes for TMRCA:
# Temporal signal method matters (Analysis 1: p < 2.2e-16, large effect)
# Phylogenetic method does NOT matter (Analysis 2: p = 0.135, not significant)
# 18-year spread in TMRCA is driven by temporal signal methods, not phylogenetic methods

# ------------------------------------------
# Analysis 3: Compare temporal methods for Evolutionary Rate
# ------------------------------------------

er_summary <- er_long %>%
  group_by(temporal_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(ER), 4),
    SD = round(sd(ER), 4),
    Median = round(median(ER), 4),
    IQR = round(IQR(ER), 4),
    Q1 = round(quantile(ER, 0.25), 4),
    Q3 = round(quantile(ER, 0.75), 4),
    Min = round(min(ER), 4),
    Max = round(max(ER), 4)
  )
print(er_summary)

# Pattern of rates:
# Slowest to fastest: TreeTime, LSD, TempEst, treedater

kw_er <- kruskal.test(ER ~ temporal_method, data = er_long)
cat("\nKruskal-Wallis test:\n")
print(kw_er)

# Result: substantial ER differences across temporal methods (more pronounced than TMRCA)

dunn_er <- dunn_test(
  ER ~ temporal_method,
  data = er_long,
  p.adjust.method = "bonferroni"
)
print(dunn_er)

# Result: TempEst vs treedater is NOT significant (p = 0.628), but all other comparisons are highly significant

eta_er <- kruskal_effsize(ER ~ temporal_method, data = er_long)
print(eta_er)

# Result: temporal method has larger effect size for ER than for TMRCA

# ------------------------------------------
# Analysis 4: Compare phylogenetic methods for Evolutionary Rate
# ------------------------------------------

er_phylo_summary <- er_long %>%
  group_by(phylo_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(ER), 4),
    SD = round(sd(ER), 4),
    Median = round(median(ER), 4),
    IQR = round(IQR(ER), 4)
  )
print(er_phylo_summary)

# Result: FastTree produces significantly higher rates (median = 1.37) than all other methods
# IQ-TREE, PhyML, RAxML-NG all give similar rates (~1.16-1.18)

kw_er_phylo <- kruskal.test(ER ~ phylo_method, data = er_long)
print(kw_er_phylo)

# Result: highly significant differences in ER across phylo methods

if (kw_er_phylo$p.value < 0.05) {
  dunn_er_phylo <- dunn_test(
    ER ~ phylo_method,
    data = er_long,
    p.adjust.method = "bonferroni"
  )
  print(dunn_er_phylo)
}

# Result: FastTree produces significantly higher rates than all other methods
# No significant differences among IQ-TREE, PhyML, and RAxML-NG (all p > 0.6)

eta_er <- kruskal_effsize(ER ~ phylo_method, data = er_long)
print(eta_er)

# Result: phylo method has moderate effect size

# Summary:
# - Temporal method choice is critical for both TMRCA and ER
# - Phylogenetic method matters for ER but not TMRCA
# - FastTree is the "outlier" for ER (but not TMRCA)
# - For ER, IQ-TREE/PhyML/RAxML-NG are interchangeable

# ------------------------------------------
# Analysis 5: Two-way comparisons (Temporal × Phylogenetic methods)
# ------------------------------------------

# Two-way analysis for TMRCA
# Summary by method combination with confidence intervals
tmrca_combination_summary <- tmrca_long %>%
  group_by(temporal_method, phylo_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(TMRCA_year), 2),
    SD = round(sd(TMRCA_year), 2),
    SE = sd(TMRCA_year) / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    Median = round(median(TMRCA_year), 2),
    IQR = round(IQR(TMRCA_year), 2),
    .groups = "drop"
  ) %>%
  arrange(temporal_method, phylo_method)
print(tmrca_combination_summary)

# Since likely non-normal, use Scheirer-Ray-Hare test (non-parametric 2-way ANOVA)

# Scheirer-Ray-Hare test
srh_tmrca <- scheirerRayHare(
  TMRCA_year ~ temporal_method + phylo_method,
  data = tmrca_long
)
print(srh_tmrca)

# Intrepretation:
# - Temporal method has a very strong effect on TMRCA (p < 0.00001)
# - Phylogenetic method doesn't have a significant effect TMRCA (p = 0.13)
# - No interaction between the two (p = 0.62)

# Two-way analysis for Evolutionary Rate
er_combination_summary <- er_long %>%
  group_by(temporal_method, phylo_method) %>%
  summarise(
    N = n(),
    Mean = round(mean(ER), 4),
    SD = round(sd(ER), 4),
    SE = sd(ER) / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    Median = round(median(ER), 4),
    IQR = round(IQR(ER), 4),
    .groups = "drop"
  ) %>%
  arrange(temporal_method, phylo_method)
print(er_combination_summary)

# Scheirer-Ray-Hare test
srh_er <- scheirerRayHare(ER ~ temporal_method + phylo_method, data = er_long)
print(srh_er)

# Interpretation:
# - Temporal method has a very strong effect on ER (H = 143.7, p < 0.00001)
# - Phylogenetic method had a moderate effect on ER (H = 21.9, p = 0.00007)
# - No interaction between the two (p = 0.12); methods have independent effects

# Next: Identify method combination with most consistent estimates across replicates (lowest coefficient of variation, CV)
# CV = sd/mean x 100%

# Coefficient of Variation by method combination, all replicates combined
# Lower CV indicates more precise/consistent estimates
cv_by_combination <- tmrca_long %>%
  group_by(temporal_method, phylo_method) %>%
  summarise(
    mean_tmrca = mean(TMRCA_year),
    sd_tmrca = sd(TMRCA_year),
    cv = (sd_tmrca / mean_tmrca) * 100,
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(cv)
print(head(cv_by_combination, 10))

cat(
  "\n Most consistent method combination:",
  as.character(cv_by_combination$temporal_method[1]),
  "+",
  as.character(cv_by_combination$phylo_method[1]),
  "with CV =",
  round(cv_by_combination$cv[1], 2),
  "%\n"
)
# Most consistent method combination: LSD + FastTree with CV = 0.23 %, followed by treedater

# ------------------------------------------
# Part 1 Figures
# ------------------------------------------

# Interaction plot for TMRCA (mean + 95% CI)
p_interaction_tmrca <- ggplot(
  tmrca_combination_summary,
  aes(
    x = phylo_method,
    y = Mean,
    color = temporal_method,
    group = temporal_method
  )
) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.1,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Phylogenetic method",
    y = "Mean TMRCA (Date)",
    color = "Temporal method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_brewer(palette = "Set1")
p_interaction_tmrca

# ggsave(
#   "figures/fig_interaction_tmrca.pdf",
#   p_interaction_tmrca,
#   width = 10,
#   height = 6
# )
# ggsave(
#   "figures/fig_interaction_tmrca.png",
#   p_interaction_tmrca,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# Interaction plot for ER
p_interaction_er <- ggplot(
  er_combination_summary,
  aes(
    x = phylo_method,
    y = Mean,
    color = temporal_method,
    group = temporal_method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    linewidth = 0.8
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Phylogenetic method",
    y = "Mean evolutionary rate",
    color = "Temporal method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_brewer(palette = "Set1")
p_interaction_er

# ggsave(
#   "figures/fig_interaction_er.pdf",
#   p_interaction_er,
#   width = 10,
#   height = 6
# )
# ggsave(
#   "figures/fig_interaction_er.png",
#   p_interaction_er,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# Heatmap for TMRCA by method combination
tmrca_heatmap_data <- tmrca_combination_summary %>%
  select(temporal_method, phylo_method, Median) %>%
  pivot_wider(names_from = phylo_method, values_from = Median)

p_heatmap_tmrca <- tmrca_combination_summary %>%
  ggplot(aes(x = phylo_method, y = temporal_method, fill = Mean)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = round(Mean, 1)),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = mean(tmrca_combination_summary$Mean),
    limits = range(tmrca_combination_summary$Mean)
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Phylogenetic method", y = "Temporal method", fill = "Mean TMRCA") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p_heatmap_tmrca

# ggsave(
#   "figures/fig_heatmap_tmrca.pdf",
#   p_heatmap_tmrca,
#   width = 8,
#   height = 6
# )
# ggsave(
#   "figures/fig_heatmap_tmrca.png",
#   p_heatmap_tmrca,
#   width = 8,
#   height = 6,
#   dpi = 300
# )

# Heatmap for ER by method combination
er_heatmap_data <- er_combination_summary %>%
  select(temporal_method, phylo_method, Median) %>%
  pivot_wider(names_from = phylo_method, values_from = Median)

p_heatmap_er <- er_combination_summary %>%
  ggplot(aes(x = phylo_method, y = temporal_method, fill = Mean)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = round(Mean, 2)),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = mean(er_combination_summary$Mean),
    limits = range(er_combination_summary$Mean)
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Phylogenetic method", y = "Temporal method", fill = "Mean ER") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p_heatmap_er

# ggsave("figures/fig_heatmap_er.pdf",
#   p_heatmap_er,
#   width = 8,
#   height = 6
# )
# ggsave(
#   "figures/fig_heatmap_er.png",
#   p_heatmap_er,
#   width = 8,
#   height = 6,
#   dpi = 300
# )

# Take Home:
# The effect of temporal method is consistent across the 4 phylogenetic methods
# Can choose phylogenetic method based on other criteria (speed, accuracy) without worrying if it will affect the temporal estimates

# ------------------------------------------
# Analysis 6: Subsampling Strategy Consistency for TMRCA and ER
# ------------------------------------------

# TMRCA:

# Coefficient of Variation by Sample Size
# Calculate CV for each subsampling level (across replicates)
# Lower CV indicates more consistent estimates across 3 replicates
# Note: Only 3 replicates per sample size, so some SD/CV may be unstable

cv_by_sample_size <- tmrca_long %>%
  group_by(dataset_base, sample_size, temporal_method, phylo_method) %>% # Group replicates together
  summarise(
    mean_tmrca = mean(TMRCA_year),
    sd_tmrca = sd(TMRCA_year),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(cv = (sd_tmrca / mean_tmrca) * 100) %>%
  filter(!is.na(cv)) # Remove NA values if any
head(cv_by_sample_size)

# Average CV across methods for each sample size
cv_summary_by_size <- cv_by_sample_size %>%
  group_by(sample_size) %>%
  summarise(
    mean_cv = round(mean(cv, na.rm = TRUE), 2),
    median_cv = round(median(cv, na.rm = TRUE), 2),
    sd_cv = round(sd(cv, na.rm = TRUE), 2),
    min_cv = round(min(cv, na.rm = TRUE), 2),
    max_cv = round(max(cv, na.rm = TRUE), 2),
    n_method_combos = n(),
    n_na = sum(is.na(cv))
  ) %>%
  arrange(mean_cv)
print(cv_summary_by_size)

if (nrow(cv_summary_by_size) > 0) {
  cat(
    "\n Most consistent sample size:",
    cv_summary_by_size$sample_size[1],
    "with mean CV =",
    cv_summary_by_size$mean_cv[1],
    "%\n"
  )
} else {
  cat("\n Unable to calculate CV due to insufficient replicates\n")
}
# Most consistent sample size: 260 with mean CV = 0.1 %

# Also calculate CV by method combination (across all sample sizes)
cv_by_combination <- tmrca_long %>%
  group_by(dataset_base, temporal_method, phylo_method) %>% # Group replicates together
  summarise(
    mean_tmrca = mean(TMRCA_year),
    sd_tmrca = sd(TMRCA_year),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(cv = (sd_tmrca / mean_tmrca) * 100) %>%
  group_by(temporal_method, phylo_method) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    mean_tmrca = mean(mean_tmrca),
    .groups = "drop"
  ) %>%
  arrange(mean_cv)

cv_by_combination %>% arrange(mean_cv)

cat(
  "\n Most consistent method combination:",
  as.character(cv_by_combination$temporal_method[1]),
  "+",
  as.character(cv_by_combination$phylo_method[1]),
  "with mean CV =",
  round(cv_by_combination$mean_cv[1], 2),
  "%\n"
)
# Most consistent method combination: LSD + RAxML-NG with mean CV = 0.06 %

## ER

# Coefficient of Variation by Sample Size
cv_by_sample_size <- er_long %>%
  group_by(dataset_base, sample_size, temporal_method, phylo_method) %>% # Group replicates together
  summarise(
    mean_ER = mean(ER),
    sd_ER = sd(ER),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(cv = (sd_ER / mean_ER) * 100) %>%
  filter(!is.na(cv)) # Remove NA values if any

# Average CV across methods for each sample size
cv_summary_by_size <- cv_by_sample_size %>%
  group_by(sample_size) %>%
  summarise(
    mean_cv = round(mean(cv, na.rm = TRUE), 2),
    median_cv = round(median(cv, na.rm = TRUE), 2),
    sd_cv = round(sd(cv, na.rm = TRUE), 2),
    min_cv = round(min(cv, na.rm = TRUE), 2),
    max_cv = round(max(cv, na.rm = TRUE), 2),
    n_method_combos = n(),
    n_na = sum(is.na(cv))
  ) %>%
  arrange(mean_cv)
print(cv_summary_by_size)

# Result: ER shows more variability across replicates than TMRCA (higher CVs)

if (nrow(cv_summary_by_size) > 0) {
  cat(
    "\n Most consistent sample size:",
    cv_summary_by_size$sample_size[1],
    "with mean CV =",
    cv_summary_by_size$mean_cv[1],
    "%\n"
  )
} else {
  cat("\n Unable to calculate CV due to insufficient replicates\n")
}
#  Most consistent sample size: 260 with mean CV = 1.72 %

# Coefficient of Variation by method combination (across all sample sizes)
cv_by_combination <- er_long %>%
  group_by(dataset_base, temporal_method, phylo_method) %>% # Group replicates together
  summarise(
    mean_ER = mean(ER),
    sd_ER = sd(ER),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(cv = (sd_ER / mean_ER) * 100) %>%
  group_by(temporal_method, phylo_method) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    mean_ER = mean(mean_ER),
    .groups = "drop"
  ) %>%
  arrange(mean_cv)
cv_by_combination %>% arrange(mean_cv)

cat(
  "\n Most consistent method combination:",
  as.character(cv_by_combination$temporal_method[1]),
  "+",
  as.character(cv_by_combination$phylo_method[1]),
  "with mean CV =",
  round(cv_by_combination$mean_cv[1], 2),
  "%\n"
)
# Most consistent method combination: LSD + IQ-TREE with mean CV = 0.89 %

# Take Home:
# Smallest dataset (260 sequences) is most consistent for both TMRCA and ER
# LSD temporal method produces most consistent results across replicates
# treedater is least consistent (highest CVs for both TMRCA and ER)
# Overall, CVs are very low (< 0.25% for TMRCA, < 9% for ER), indicating good replicate reproducibility

# ============================================================================
# PART 2: INTRODUCTION COUNT ANALYSIS
# ============================================================================

intro_count <- read.delim("input_data/intro_count.tsv")
summary(intro_count)

# Create base dataset column (remove replicate suffix)
intro_count <- intro_count %>%
  mutate(
    replicate = as.numeric(str_extract(dataset, "\\d+$")),
    dataset_base = str_extract(dataset, "locrisk\\d+")
  )

cat("\nUnique base datasets:", length(unique(intro_count$dataset_base)), "\n")
cat(
  "Total datasets including replicates:",
  length(unique(intro_count$dataset)),
  "\n"
)

# ------------------------------------------
# Analysis 7: Introduction Count Distribution
# ------------------------------------------

# Overall summary statistics
overall_summary <- intro_count %>%
  summarise(
    N_observations = n(),
    N_routes = n_distinct(paste(from, to, dataset_base, dpt, mlt)),
    Mean_introductions = round(mean(n), 2),
    Median_introductions = median(n),
    SD = round(sd(n), 2),
    Min = min(n),
    Max = max(n),
    Q1 = quantile(n, 0.25),
    Q3 = quantile(n, 0.75),
    IQR = IQR(n)
  )
print(overall_summary)

n_samples <- nrow(intro_count)

cat(
  "\n Introduction counts: ",
  n_samples,
  " observations and discrete values (1-",
  overall_summary$Max,
  ")\n",
  sep = ""
)

# Frequency distribution of introduction counts
intro_freq <- intro_count %>%
  count(n, name = "frequency") %>%
  mutate(percent = round(frequency / sum(frequency) * 100, 2)) %>%
  arrange(n)
print(head(intro_freq, 10))

# Skewness and kurtosis for distribution shape
skew_val <- moments::skewness(intro_count$n)
kurt_val <- moments::kurtosis(intro_count$n)

cat(
  "  Skewness:",
  round(skew_val, 3),
  ifelse(
    abs(skew_val) < 0.5,
    "(approximately symmetric)",
    ifelse(skew_val > 0, "(right-skewed)", "(left-skewed)")
  ),
  "\n"
)
cat(
  "  Kurtosis:",
  round(kurt_val, 3),
  ifelse(
    abs(kurt_val - 3) < 0.5,
    "(approximately normal)",
    ifelse(kurt_val > 3, "(heavy-tailed)", "(light-tailed)")
  ),
  "\n"
)

# Check for overdispersion (variance > mean suggests negative binomial)
variance_mean_ratio <- overall_summary$SD^2 / overall_summary$Mean_introductions
cat("\nVariance-to-mean ratio:", round(variance_mean_ratio, 2), "\n")
if (variance_mean_ratio > 1.5) {
  cat("   → Data are overdispersed (variance > mean)\n")
  cat("   → Suggests negative binomial distribution rather than Poisson\n")
} else if (variance_mean_ratio < 0.7) {
  cat("   → Data are underdispersed (variance < mean)\n")
} else {
  cat("   → Variance ≈ mean, consistent with Poisson distribution\n")
}

# Recommendation: Use non-parametric tests, which don't make assumptions about the underlying distribution

# ------------------------------------------
# Analysis 8: Method Comparison for Introduction Counts
# ------------------------------------------

# Intro count summary by temporal method
intro_temporal_summary <- intro_count %>%
  group_by(dpt) %>%
  summarise(
    Total = sum(n),
    Mean = round(mean(n), 2),
    Median = median(n),
    SD = round(sd(n), 2),
    N_routes = n()
  ) %>%
  arrange(desc(Mean))
print(intro_temporal_summary)

# Summary by phylo method
intro_phylo_summary <- intro_count %>%
  group_by(mlt) %>%
  summarise(
    Total = sum(n),
    Mean = round(mean(n), 2),
    Median = median(n),
    SD = round(sd(n), 2),
    N_routes = n()
  ) %>%
  arrange(desc(Mean))
print(intro_phylo_summary)

# Kruskal-Wallis test for comparing counts across temporal methods
kw_temporal_intro <- kruskal.test(n ~ dpt, data = intro_count)
print(kw_temporal_intro)

if (kw_temporal_intro$p.value < 0.05) {
  cat("\nPost-hoc Dunn test:\n")
  dunn_temporal_intro <- dunn.test(
    intro_count$n,
    intro_count$dpt,
    method = "bonferroni"
  )
}

# KW test for phylo methods
kw_phylo_intro <- kruskal.test(n ~ mlt, data = intro_count)
print(kw_phylo_intro)

if (kw_phylo_intro$p.value < 0.05) {
  cat("\nPost-hoc Dunn test:\n")
  dunn_phylo_intro <- dunn.test(
    intro_count$n,
    intro_count$mlt,
    method = "bonferroni"
  )
}

## Take Home: Neither temporal signal method nor ML method significantly affects # introduction events
# - Introduction count estimates are robust across all method combinations
# - Phylogenetic or temporal method doesn't substantially bias the results
# - Conclusions about viral dispersal patterns are method-independent

# ------------------------------------------
# Analysis 9: Introduction Counts by Subsampling strategy
# ------------------------------------------

# Replicates grouped together by base dataset
intro_dataset_summary <- intro_count %>%
  group_by(dataset_base) %>%
  summarise(
    Total = sum(n),
    Mean = round(mean(n), 2),
    Median = median(n),
    SD = round(sd(n), 2),
    N_routes = n(),
    N_replicates = n_distinct(dataset)
  ) %>%
  arrange(desc(Mean))
print(intro_dataset_summary)

# Result: Larger datasets detect MORE introductions (locrisk574: mean = 2.81, locrisk260: mean = 2.02)
# As sample size increases, mean introductions increase
# More sequences enable more robust estimation of transmission routes

# Consistency across replicates (CV by base dataset)
intro_cv_by_dataset <- intro_count %>%
  group_by(dataset_base, from, to, dpt, mlt) %>%
  summarise(
    mean_n = mean(n),
    sd_n = sd(n),
    cv = (sd(n) / mean(n)) * 100,
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(cv))

intro_cv_summary <- intro_cv_by_dataset %>%
  group_by(dataset_base) %>%
  summarise(
    mean_cv = round(mean(cv, na.rm = TRUE), 2),
    median_cv = round(median(cv, na.rm = TRUE), 2),
    sd_cv = round(sd(cv, na.rm = TRUE), 2),
    min_cv = round(min(cv, na.rm = TRUE), 2),
    max_cv = round(max(cv, na.rm = TRUE), 2)
  ) %>%
  arrange(mean_cv)
print(intro_cv_summary)

# Result: Smaller datasets tend to be more consistent across replicates (locrisk378: CV = 11.4%)
# Larger datasets are less consistent (locrisk574: CV = 17.4%)
# All median CVs = 0, meaning many routes have perfect consistency (SD = 0) across replicates

# Take Home:
# Smaller datasets = fewer, more certain introduction events
# Larger datasets = more introduction events detected, but with more variability between replicates
# Interpretation: Larger datasets capture more fine-scale transmission dynamics but with increased stochastic variation

# Relationship to Part 1 findings:
#  The smallest dataset (locrisk260) is most consistent (CV = 13.1%)
#  This mirrors the finding from Part 1 where locrisk260 had the most consistent TMRCA/ER estimates
# Conclusion: locrisk260 provides the most reproducible estimates, but it may underestimate the total number of introductions compared to larger datasets

# ------------------------------------------
# Analysis 10: Two-way Analysis for Introduction Counts
# ------------------------------------------
intro_combination_summary <- intro_count %>%
  mutate(combination = paste(dpt, mlt, sep = "_")) %>%
  group_by(combination) %>%
  summarise(
    Mean = round(mean(n), 2),
    Median = median(n),
    SD = round(sd(n), 2),
    Min = min(n),
    Max = max(n),
    N = n()
  ) %>%
  arrange(desc(Mean))
print(intro_combination_summary)

# Scheirer-Ray-Hare test (non-parametric two-way ANOVA)
# Interaction for temporal method × phylogenetic method
intro_srh <- rcompanion::scheirerRayHare(
  n ~ dpt + mlt + dpt:mlt,
  data = intro_count
)
print(intro_srh)

# Effect sizes for main effects
# Temporal method
intro_es_temporal <- rstatix::kruskal_effsize(n ~ dpt, data = intro_count)
print(intro_es_temporal)

# Phylogenetic method
intro_es_phylo <- rstatix::kruskal_effsize(n ~ mlt, data = intro_count)
print(intro_es_phylo)

# Extract p-values using row names
p_dpt <- intro_srh["dpt", "p.value"]
p_mlt <- intro_srh["mlt", "p.value"]
p_interaction <- intro_srh["dpt:mlt", "p.value"]

# Temporal method does NOT have a significant effect (p = 0.7258 )

# Phylogenetic method does NOT have a significant effect (p = 0.9903 )

# No significant interaction between methods (p = 0.9998 ): methods have independent effects on introduction counts

# Take Home: confirms that introduction counts are robust across all method combinations

# ------------------------------------------
# Analysis 11: Subsampling Strategy Consistency for Introduction Counts
# ------------------------------------------

# Consistency by method combination (across replicates)
intro_cv_by_combination <- intro_count %>%
  mutate(combination = paste(dpt, mlt, sep = "_")) %>%
  group_by(dataset_base, combination, from, to) %>%
  summarise(
    mean_n = mean(n),
    sd_n = sd(n),
    cv = (sd(n) / mean(n)) * 100,
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(cv))

intro_cv_combination_summary <- intro_cv_by_combination %>%
  group_by(combination) %>%
  summarise(
    mean_cv = round(mean(cv, na.rm = TRUE), 2),
    median_cv = round(median(cv, na.rm = TRUE), 2),
    n_routes = n()
  ) %>%
  separate(
    combination,
    into = c("temporal_method", "phylo_method"),
    sep = "_",
    remove = FALSE
  ) %>%
  arrange(mean_cv)
print(intro_cv_combination_summary)

# Identify most and least consistent combinations
cat(
  "Most consistent combination:",
  intro_cv_combination_summary$combination[1],
  "with mean CV =",
  intro_cv_combination_summary$mean_cv[1],
  "%\n"
)
cat(
  "Least consistent combination:",
  intro_cv_combination_summary$combination[nrow(intro_cv_combination_summary)],
  "with mean CV =",
  intro_cv_combination_summary$mean_cv[nrow(intro_cv_combination_summary)],
  "%\n"
)

# Most consistent combination: TreeTime_IQ-TREE with mean CV = 11.77 %
# Least consistent combination: treedater_RAxML-NG with mean CV = 16.4 %

# Compare consistency across temporal methods
cv_by_temporal <- intro_cv_combination_summary %>%
  group_by(temporal_method) %>%
  summarise(
    mean_cv = round(mean(mean_cv), 2),
    median_cv = round(median(mean_cv), 2),
    min_cv = round(min(mean_cv), 2),
    max_cv = round(max(mean_cv), 2)
  ) %>%
  arrange(mean_cv)
print(cv_by_temporal)

# Compare consistency across phylogenetic methods
cv_by_phylo <- intro_cv_combination_summary %>%
  group_by(phylo_method) %>%
  summarise(
    mean_cv = round(mean(mean_cv), 2),
    median_cv = round(median(mean_cv), 2),
    min_cv = round(min(mean_cv), 2),
    max_cv = round(max(mean_cv), 2)
  ) %>%
  arrange(mean_cv)
print(cv_by_phylo)

cat(
  "Most consistent temporal method:",
  cv_by_temporal$temporal_method[1],
  "(mean CV =",
  cv_by_temporal$mean_cv[1],
  "%)\n"
)
cat(
  "Most consistent phylogenetic method:",
  cv_by_phylo$phylo_method[1],
  "(mean CV =",
  cv_by_phylo$mean_cv[1],
  "%)\n"
)
# Most consistent temporal method: TreeTime and LSD (both have mean CV = 14.2 %)
# Most consistent phylogenetic method: IQ-TREE (mean CV = 13.2 %)

# Take Home:
# - IQ-TREE produces the most consistent introduction counts across all temporal methods (CV = 13.2%)
# - TreeTime and LSD are the most consistent temporal methods (both CV = 14.2%)
# - Best combination overall: TreeTime + IQ-TREE (CV = 11.77%)
# - All CVs are relatively low (11.8-16.4%), indicating good reproducibility across all method combinations
# - PhyML detects slightly more introductions (mean = 2.61) but is less consistent than IQ-TREE

# ------------------------------------------
# PART 2 Figures
# ------------------------------------------

# Interaction Plot for Introduction Counts

# Calculate means and standard errors for interaction plot
intro_interaction_data <- intro_count %>%
  group_by(dpt, mlt) %>%
  summarise(
    mean_n = mean(n),
    se_n = sd(n) / sqrt(n()),
    .groups = "drop"
  )

# Interaction plot
p_intro_interaction <- ggplot(
  intro_interaction_data,
  aes(x = mlt, y = mean_n, color = dpt, group = dpt)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_n - se_n, ymax = mean_n + se_n),
    width = 0.1,
    linewidth = 0.8
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Phylogenetic method",
    y = "Mean introduction count",
    color = "Temporal method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_color_brewer(palette = "Set1")
p_intro_interaction

# ggsave(
#   "figures/fig_intro_interaction.pdf",
#   p_intro_interaction,
#   width = 10,
#   height = 6
# )
# ggsave(
#   "figures/fig_intro_interaction.png",
#   p_intro_interaction,
#   width = 10,
#   height = 6,
#   dpi = 300
# )

# Heatmap of Intro Count CV by Method Combination

p_intro_cv_heatmap <- ggplot(
  intro_cv_combination_summary,
  aes(x = phylo_method, y = temporal_method, fill = mean_cv)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = round(mean_cv, 2)),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "red",
    midpoint = mean(intro_cv_combination_summary$mean_cv),
    limits = range(intro_cv_combination_summary$mean_cv),
    name = "Mean count\nCV"
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Phylogenetic method", y = "Temporal method") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p_intro_cv_heatmap

# ggsave(
#   "figures/fig_intro_cv_heatmap.pdf",
#   p_intro_cv_heatmap,
#   width = 8,
#   height = 6
# )
# ggsave(
#   "figures/fig_intro_cv_heatmap.png",
#   p_intro_cv_heatmap,
#   width = 8,
#   height = 6,
#   dpi = 300
# )

# ============================================================================
# Manuscript Figures
# ============================================================================

p_legend <- get_legend(p_interaction_tmrca + theme(legend.position = "bottom"))

interact_grid <- plot_grid(
  p_interaction_tmrca + theme(legend.position = "none"),
  p_interaction_er + theme(legend.position = "none"),
  p_intro_interaction + theme(legend.position = "none"),
  nrow = 3,
  labels = "AUTO"
)

interaction_p <- plot_grid(
  interact_grid,
  p_legend,
  nrow = 2,
  rel_heights = c(3, 0.2)
)
interaction_p
ggsave(
  "figures/manuscript_fig_s18.png",
  interaction_p,
  width = 8,
  height = 12,
  dpi = 300
)
ggsave(
  "figures/manuscript_fig_s18.pdf",
  interaction_p,
  width = 8,
  height = 12,
  dpi = 300
)

heatmap_grid <- plot_grid(
  p_heatmap_tmrca,
  p_heatmap_er,
  p_intro_cv_heatmap,
  nrow = 3,
  labels = "AUTO", align = "hv"
)
heatmap_grid
ggsave(
  "figures/manuscript_fig_s19.png",
  heatmap_grid,
  width = 8,
  height = 14,
  dpi = 300
)
ggsave(
  "figures/manuscript_fig_s19.pdf",
  heatmap_grid,
  width = 8,
  height = 14,
  dpi = 300
)
