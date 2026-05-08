# ============================================================================
# Concordance Analysis of TreeTime Source-Sink Scores (SSS)
# HIV-1 Subtype C Phylodynamics Study
# Amanda Perofsky
# ============================================================================
# Tests whether source-sink rankings are robust across:
#   (1) ML tree methods (FastTree, IQ-TREE, PhyML, RAxML-NG)
#   (2) Temporal dating methods (LSD, TempEst, TreeTime, treedater)
#   (3) Subsampled datasets (locrisk260 to locrisk574)
#
# Approach: Kendall's W concordance coefficient
#   Do method combinations / datasets produce the same rankings of
#   locations/risk groups as sources vs sinks?
#
# Metrics tested: Export, Import, SSS
#   SSS = (Export - Import) / (Export + Import), ranges from -1 to +1
#
# Resolution levels: continent, region, country, risk6, risk5, risk4
#
# Inclusion criterion: for each trait × dataset, only locations present in
#   all 16 method combinations (4 temporal × 4 ML) are retained. A location
#   need not appear in every dataset — smaller subsamples naturally detect
#   fewer locations. This mirrors the within-dataset completeness filter
#   applied in the introduction timing LMM.
# ============================================================================

library(tidyverse)

dir.create("stat_results", showWarnings = FALSE)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

sss <- read.delim("input_data/All_treetime_metrics_SSS.tsv", sep = "\t")

cat(sprintf("Total rows: %d\n", nrow(sss)))
cat(sprintf("Columns: %s\n", paste(names(sss), collapse = ", ")))

# Parse dataset column
sss <- sss %>%
  mutate(
    dataset_base = str_extract(dataset, "locrisk\\d+"),
    replicate = as.numeric(str_extract(dataset, "(?<=\\.)\\d+$")),
    method_combo = paste(dpt, mlt, dataset, sep = "_"), # all 240 combos
    temporal_phylo = paste(dpt, mlt, sep = "_") # temporal x ML only
  )

cat(sprintf("Trait levels: %s\n", paste(sort(unique(sss$trait)), collapse = ", ")))
cat(sprintf("Temporal methods: %s\n", paste(sort(unique(sss$dpt)), collapse = ", ")))
cat(sprintf("ML methods: %s\n", paste(sort(unique(sss$mlt)), collapse = ", ")))
cat(sprintf("Datasets: %s\n", paste(sort(unique(sss$dataset_base)), collapse = ", ")))

total_combos <- n_distinct(sss$method_combo)
cat(sprintf(
  "Total unique method-dataset combos (dpt x mlt x dataset): %d\n\n",
  total_combos
))

# ============================================================================
# 2. APPLY COMPLETENESS FILTER
#
# For each trait × dataset combination, retain only locations present in all
# 16 method combinations (4 temporal × 4 ML). A location does not need to
# appear in every dataset — smaller subsamples naturally detect fewer
# locations. This mirrors the within-dataset completeness criterion used in
# the introduction timing LMM.
#
# Locations detected in only a subset of method combinations within a given
# dataset (e.g., CY and GB at country level) are excluded from that dataset.
# Imputing zeros for undetected locations would conflate true absence from a
# network with uninformative structural zeros, artificially inflating
# concordance for those nodes.
# ============================================================================

n_method_combos <- n_distinct(sss$temporal_phylo) # 16: 4 dpt x 4 mlt

complete_location_datasets <- sss %>%
  group_by(trait, dataset, location) %>%
  summarise(n_mc = n_distinct(temporal_phylo), .groups = "drop") %>%
  filter(n_mc == n_method_combos)

# Report excluded trait × dataset × location combinations
excluded <- sss %>%
  group_by(trait, dataset, location) %>%
  summarise(n_mc = n_distinct(temporal_phylo), .groups = "drop") %>%
  filter(n_mc < n_method_combos) %>%
  arrange(trait, location, dataset)

if (nrow(excluded) > 0) {
  cat("Excluded location × dataset combinations (not present in all 16 method combos):\n")
  excluded %>%
    mutate(label = sprintf(
      "  %-12s %-8s %-20s %d/16 combos",
      trait, dataset, location, n_mc
    )) %>%
    pull(label) %>%
    cat(sep = "\n")
  cat("\n\n")
} else {
  cat("No exclusions — all location × dataset combinations present in all 16 method combos.\n\n")
}

sss_filtered <- sss %>%
  semi_join(complete_location_datasets, by = c("trait", "dataset", "location"))

cat("Rows before filtering:", nrow(sss), "\n")
cat("Rows after filtering: ", nrow(sss_filtered), "\n")
cat("Rows removed:         ", nrow(sss) - nrow(sss_filtered), "\n\n")

# Report retained locations per trait (may vary across datasets)
sss_filtered %>%
  group_by(trait) %>%
  summarise(
    n_locations = n_distinct(location),
    locations = paste(sort(unique(location)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("  %-12s %d: %s", trait, n_locations, locations)) %>%
  pull(label) %>%
  cat(sep = "\n")

# Metrics and resolution levels
metrics <- c("Export", "Import", "SSS")
metric_labels <- c(
  "Export" = "Export (source activity)",
  "Import" = "Import (sink activity)",
  "SSS"    = "Source-Sink Score"
)
resolutions <- sort(unique(sss_filtered$trait))

# ============================================================================
# 3. KENDALL'S W FUNCTION
# ============================================================================

compute_kendall_w <- function(data, metric, rater_var, subject_var = "location") {
  avg_data <- data %>%
    group_by(across(all_of(c(subject_var, rater_var)))) %>%
    summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

  wide_data <- avg_data %>%
    pivot_wider(names_from = all_of(rater_var), values_from = value)

  mat <- wide_data %>%
    select(-all_of(subject_var)) %>%
    as.matrix()

  if (any(is.na(mat))) {
    return(NULL)
  }

  n_subjects <- nrow(mat)
  n_raters <- ncol(mat)
  if (n_subjects < 2 || n_raters < 2) {
    return(NULL)
  }

  rank_mat <- apply(mat, 2, rank)
  row_sums <- rowSums(rank_mat)
  mean_row_sum <- mean(row_sums)
  S <- sum((row_sums - mean_row_sum)^2)
  k <- n_raters
  n <- n_subjects
  W <- (12 * S) / (k^2 * (n^3 - n))

  chi_sq <- k * (n - 1) * W
  df <- n - 1
  p_value <- pchisq(chi_sq, df, lower.tail = FALSE)

  mean_ranks <- rowMeans(rank_mat)
  consensus <- tibble(
    location  = wide_data[[subject_var]],
    mean_rank = mean_ranks
  ) %>% arrange(desc(mean_rank))

  list(
    W = W, chi_sq = chi_sq, df = df, p = p_value,
    n_subjects = n, n_raters = k, consensus_ranking = consensus
  )
}

# ============================================================================
# 4. CONCORDANCE ANALYSIS
# ============================================================================
# W ranges from 0 (no agreement) to 1 (perfect agreement)
# W > 0.7 is generally considered strong concordance

all_concordance_results <- list()

for (res in resolutions) {
  cat(sprintf("\n%s\n", paste(rep("-", 60), collapse = "")))
  cat(sprintf("RESOLUTION: %s\n", toupper(res)))
  cat(sprintf("%s\n", paste(rep("-", 60), collapse = "")))

  res_data <- sss_filtered %>% filter(trait == res)

  cat(sprintf(
    "  Locations (%d): %s\n",
    n_distinct(res_data$location),
    paste(sort(unique(res_data$location)), collapse = ", ")
  ))
  cat(sprintf("  N rows: %d\n", nrow(res_data)))

  for (metric in metrics) {
    cat(sprintf("\n  %s:\n", metric_labels[metric]))

    # ---- 1. Concordance across ML methods (k = 4) ----
    w_mlt <- compute_kendall_w(res_data, metric, "mlt")
    if (!is.null(w_mlt)) {
      cat(sprintf(
        "    Across ML methods:       W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d)\n",
        w_mlt$W, w_mlt$chi_sq, w_mlt$df,
        ifelse(w_mlt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_mlt$p)),
        w_mlt$n_raters
      ))
      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res, metric = metric_labels[metric],
        rater = "ML tree method",
        W = round(w_mlt$W, 4), chi_sq = round(w_mlt$chi_sq, 2),
        df = w_mlt$df, p = w_mlt$p,
        p_formatted = ifelse(w_mlt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_mlt$p)),
        n_raters = w_mlt$n_raters, n_subjects = w_mlt$n_subjects
      )
    }

    # ---- 2. Concordance across temporal methods (k = 4) ----
    w_dpt <- compute_kendall_w(res_data, metric, "dpt")
    if (!is.null(w_dpt)) {
      cat(sprintf(
        "    Across temporal methods: W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d)\n",
        w_dpt$W, w_dpt$chi_sq, w_dpt$df,
        ifelse(w_dpt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_dpt$p)),
        w_dpt$n_raters
      ))
      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res, metric = metric_labels[metric],
        rater = "Temporal method",
        W = round(w_dpt$W, 4), chi_sq = round(w_dpt$chi_sq, 2),
        df = w_dpt$df, p = w_dpt$p,
        p_formatted = ifelse(w_dpt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_dpt$p)),
        n_raters = w_dpt$n_raters, n_subjects = w_dpt$n_subjects
      )
    }

    # ---- 3. Concordance across datasets (k = 5) ----
    w_ds <- compute_kendall_w(res_data, metric, "dataset_base")
    if (!is.null(w_ds)) {
      cat(sprintf(
        "    Across datasets:         W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d)\n",
        w_ds$W, w_ds$chi_sq, w_ds$df,
        ifelse(w_ds$p < 0.0001, "< 0.0001", sprintf("%.4f", w_ds$p)),
        w_ds$n_raters
      ))
      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res, metric = metric_labels[metric],
        rater = "Dataset",
        W = round(w_ds$W, 4), chi_sq = round(w_ds$chi_sq, 2),
        df = w_ds$df, p = w_ds$p,
        p_formatted = ifelse(w_ds$p < 0.0001, "< 0.0001", sprintf("%.4f", w_ds$p)),
        n_raters = w_ds$n_raters, n_subjects = w_ds$n_subjects
      )
    }

    # ---- 4. Concordance across all combinations (k = 240) ----
    w_all <- compute_kendall_w(res_data, metric, "method_combo")
    if (!is.null(w_all)) {
      cat(sprintf(
        "    Across all combos:       W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d)\n",
        w_all$W, w_all$chi_sq, w_all$df,
        ifelse(w_all$p < 0.0001, "< 0.0001", sprintf("%.4f", w_all$p)),
        w_all$n_raters
      ))
      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res, metric = metric_labels[metric],
        rater = "All combinations",
        W = round(w_all$W, 4), chi_sq = round(w_all$chi_sq, 2),
        df = w_all$df, p = w_all$p,
        p_formatted = ifelse(w_all$p < 0.0001, "< 0.0001", sprintf("%.4f", w_all$p)),
        n_raters = w_all$n_raters, n_subjects = w_all$n_subjects
      )

      cat("    Consensus ranking (highest = most source-like):\n")
      for (j in 1:nrow(w_all$consensus_ranking)) {
        cat(sprintf(
          "      %d. %s (mean rank = %.1f)\n",
          j,
          w_all$consensus_ranking$location[j],
          w_all$consensus_ranking$mean_rank[j]
        ))
      }
    }
  }
}

# ============================================================================
# 5. SUMMARY TABLES
# ============================================================================

concordance_summary <- bind_rows(all_concordance_results)

# Kendall's W Concordance Summary
concordance_summary %>%
  select(resolution, metric, rater, W, p_formatted, n_raters, n_subjects) %>%
  print(n = Inf)

# Concordance strength (all method-dataset combinations)
concordance_summary %>%
  filter(rater == "All combinations") %>%
  mutate(strength = case_when(
    W >= 0.9 ~ "very strong",
    W >= 0.7 ~ "strong",
    W >= 0.5 ~ "moderate",
    W >= 0.3 ~ "weak",
    TRUE ~ "very weak"
  )) %>%
  select(resolution, metric, W, strength, p_formatted) %>%
  print(n = Inf)

write.csv(concordance_summary, "stat_results/Supplementary_Table_S4_treetime_sss_concordance.csv",
  row.names = FALSE
)
# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

for (rater_level in c("ML tree method", "Temporal method", "Dataset", "All combinations")) {
  cat(sprintf("\n%s:\n", rater_level))
  concordance_summary %>%
    filter(rater == rater_level) %>%
    group_by(resolution) %>%
    summarise(
      mean_W = round(mean(W), 3),
      min_W = round(min(W), 3),
      max_W = round(max(W), 3),
      all_significant = all(p < 0.05),
      .groups = "drop"
    ) %>%
    print()
}
# SSS are consistent across ML tree methods, temporal methods, and datasets
