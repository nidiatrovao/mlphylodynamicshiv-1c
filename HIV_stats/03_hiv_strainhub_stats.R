# ============================================================================
# Concordance Analysis of StrainHub Centrality Metrics
# HIV-1 Subtype C Phylodynamics Study
# Amanda Perofsky
# ============================================================================
# Tests whether StrainHub network centrality rankings are robust across:
#   (1) ML tree methods (FastTree, IQ-TREE, PhyML, RAxML-NG)
#   (2) Subsampled datasets (locrisk260 – locrisk574)
#
# Approach: Kendall's W concordance coefficient
#   Do ML methods / datasets produce the same rankings of Metastates
#   (countries, regions, continents, risk groups) for each centrality metric?
#   W near 1 = perfect agreement in rankings across methods.
#
# Metrics tested:
#   Degree, In-degree, Out-degree, Betweenness, Closeness, Source.Hub.Ratio
#
# Resolution levels:
#   continent (Metastates: Africa, Americas, Asia, Europe)
#   country (Metastates: 29 countries)
#   region (Metastates: 13 regions)
#   risk6 (SM, PI, SH, MB, NR, OT — full 6-category classification)
#   risk5 (SM, PI, SH, MB, OT — NR merged into OT)
#   risk4 (SM, PI, SH, OT — NR and MB merged into OT)
# ============================================================================

library(tidyverse)

# Set output directory
dir.create("stat_results", showWarnings = FALSE)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

# Define file paths and resolution labels
file_info <- tibble(
  file = c(
    "input_data/strainhub/All_strainhub_metrics_continent.csv",
    "input_data/strainhub/All_strainhub_metrics_region.csv",
    "input_data/strainhub/All_strainhub_metrics_country.csv",
    "input_data/strainhub/All_strainhub_metrics_risk4.csv",
    "input_data/strainhub/All_strainhub_metrics_risk5.csv",
    "input_data/strainhub/All_strainhub_metrics_risk6.csv"
  ),
  resolution = c("continent", "region", "country", "risk4", "risk5", "risk6")
)

# Read all files and combine
all_data <- file_info %>%
  pmap_dfr(function(file, resolution) {
    read.csv(file) %>%
      mutate(resolution = resolution)
  })

# Parse dataset column into base dataset and replicate number
# dataset column format: "locrisk260.1" -> base = "locrisk260", replicate = 1
all_data <- all_data %>%
  mutate(
    dataset_base = str_extract(dataset, "locrisk\\d+"),
    replicate = as.numeric(str_extract(dataset, "(?<=\\.)\\d+$")),
    # Create a unique "rater" ID for each ML method x dataset combination
    method_dataset = paste(mlt, dataset, sep = "_")
  )

cat(sprintf("Total rows: %d\n", nrow(all_data)))
cat(sprintf("Resolutions: %s\n", paste(unique(all_data$resolution), collapse = ", ")))
cat(sprintf("Datasets: %s\n", paste(sort(unique(all_data$dataset_base)), collapse = ", ")))
cat(sprintf("ML methods: %s\n", paste(sort(unique(all_data$mlt)), collapse = ", ")))

# Define the 6 metrics to analyze
metrics <- c(
  "Degree.Centrality", "Indegree.Centrality", "Outdegree.Centrality",
  "Betweenness.Centrality", "Closeness.Centrality", "Source.Hub.Ratio"
)

# Nice labels for output
metric_labels <- c(
  "Degree.Centrality" = "Degree Centrality",
  "Indegree.Centrality" = "In-degree Centrality",
  "Outdegree.Centrality" = "Out-degree Centrality",
  "Betweenness.Centrality" = "Betweenness Centrality",
  "Closeness.Centrality" = "Closeness Centrality",
  "Source.Hub.Ratio" = "Source Hub Ratio"
)

# Note: we don't need to apply a completeness filter (as in SSS analysis)
# because all locations/risk groups are included across the 4 phylo ML methods.
# ============================================================================
# 2. KENDALL'S W FUNCTION
# ============================================================================

# Computes Kendall's W for ranking concordance
# Tests whether "raters" (ML methods or datasets) rank Metastates consistently
compute_kendall_w <- function(data, metric, rater_var, subject_var = "Metastates") {
  # Average metric values for each Metastate within each rater level
  # (averaging across the other grouping variables)
  avg_data <- data %>%
    group_by(across(all_of(c(subject_var, rater_var)))) %>%
    summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

  # Pivot to wide: rows = Metastates, columns = raters
  wide_data <- avg_data %>%
    pivot_wider(
      names_from = all_of(rater_var),
      values_from = value
    )

  # Check for complete cases
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

  # Rank each rater's scores (column-wise ranking)
  rank_mat <- apply(mat, 2, rank)

  # Kendall's W calculation
  # W = 12 * S / (k^2 * (n^3 - n))
  # where S = sum of squared deviations of row sums from their mean
  row_sums <- rowSums(rank_mat)
  mean_row_sum <- mean(row_sums)
  S <- sum((row_sums - mean_row_sum)^2)
  k <- n_raters
  n <- n_subjects
  W <- (12 * S) / (k^2 * (n^3 - n))

  # Friedman chi-squared approximation
  chi_sq <- k * (n - 1) * W
  df <- n - 1
  p_value <- pchisq(chi_sq, df, lower.tail = FALSE)

  # Get the consensus ranking (by mean rank across raters)
  mean_ranks <- rowMeans(rank_mat)
  consensus <- tibble(
    Metastate = wide_data[[subject_var]],
    mean_rank = mean_ranks
  ) %>%
    arrange(desc(mean_rank))

  list(
    W = W,
    chi_sq = chi_sq,
    df = df,
    p = p_value,
    n_subjects = n,
    n_raters = k,
    consensus_ranking = consensus
  )
}

# ============================================================================
# 3. CONCORDANCE ANALYSIS
# ============================================================================

# W ranges from 0 (no agreement) to 1 (perfect agreement)
# W > 0.7 is generally considered strong concordance

all_concordance_results <- list()

for (res in unique(all_data$resolution)) {
  cat(sprintf("\n%s\n", paste(rep("-", 60), collapse = "")))
  cat(sprintf("RESOLUTION: %s\n", toupper(res)))
  cat(sprintf("%s\n", paste(rep("-", 60), collapse = "")))

  res_data <- all_data %>% filter(resolution == res)

  for (metric in metrics) {
    cat(sprintf("\n  %s:\n", metric_labels[metric]))

    # ---- Concordance across ML methods ----
    # Each ML method is a "rater"; Metastates are the "subjects"
    # Average across datasets and replicates within each ML method
    w_mlt <- compute_kendall_w(res_data, metric, "mlt")

    if (!is.null(w_mlt)) {
      cat(sprintf(
        "    Across ML methods:  W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d raters, n = %d subjects)\n",
        w_mlt$W, w_mlt$chi_sq, w_mlt$df,
        ifelse(w_mlt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_mlt$p)),
        w_mlt$n_raters, w_mlt$n_subjects
      ))

      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res,
        metric = metric_labels[metric],
        rater = "ML tree method",
        W = round(w_mlt$W, 4),
        chi_sq = round(w_mlt$chi_sq, 2),
        df = w_mlt$df,
        p = w_mlt$p,
        p_formatted = ifelse(w_mlt$p < 0.0001, "< 0.0001", sprintf("%.4f", w_mlt$p)),
        n_raters = w_mlt$n_raters,
        n_subjects = w_mlt$n_subjects
      )
    } else {
      cat("    Across ML methods:  Could not compute (missing data)\n")
    }

    # ---- Concordance across datasets ----
    # Each dataset_base is a "rater"; Metastates are the "subjects"
    # Average across ML methods and replicates within each dataset
    w_ds <- compute_kendall_w(res_data, metric, "dataset_base")

    if (!is.null(w_ds)) {
      cat(sprintf(
        "    Across datasets:    W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d raters, n = %d subjects)\n",
        w_ds$W, w_ds$chi_sq, w_ds$df,
        ifelse(w_ds$p < 0.0001, "< 0.0001", sprintf("%.4f", w_ds$p)),
        w_ds$n_raters, w_ds$n_subjects
      ))

      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res,
        metric = metric_labels[metric],
        rater = "Dataset",
        W = round(w_ds$W, 4),
        chi_sq = round(w_ds$chi_sq, 2),
        df = w_ds$df,
        p = w_ds$p,
        p_formatted = ifelse(w_ds$p < 0.0001, "< 0.0001", sprintf("%.4f", w_ds$p)),
        n_raters = w_ds$n_raters,
        n_subjects = w_ds$n_subjects
      )
    } else {
      cat("    Across datasets:    Could not compute (missing data)\n")
    }

    # ---- Concordance across all method-dataset combinations ----
    # Each unique ML x dataset x replicate combo is a "rater"
    w_all <- compute_kendall_w(res_data, metric, "method_dataset")

    if (!is.null(w_all)) {
      cat(sprintf(
        "    Across all combos:  W = %.4f, chi2 = %.2f, df = %d, p = %s (k = %d raters, n = %d subjects)\n",
        w_all$W, w_all$chi_sq, w_all$df,
        ifelse(w_all$p < 0.0001, "< 0.0001", sprintf("%.4f", w_all$p)),
        w_all$n_raters, w_all$n_subjects
      ))

      all_concordance_results[[length(all_concordance_results) + 1]] <- tibble(
        resolution = res,
        metric = metric_labels[metric],
        rater = "All combinations",
        W = round(w_all$W, 4),
        chi_sq = round(w_all$chi_sq, 2),
        df = w_all$df,
        p = w_all$p,
        p_formatted = ifelse(w_all$p < 0.0001, "< 0.0001", sprintf("%.4f", w_all$p)),
        n_raters = w_all$n_raters,
        n_subjects = w_all$n_subjects
      )

      # Print consensus ranking for full combinations
      cat("    Consensus ranking (highest = most central):\n")
      for (j in 1:nrow(w_all$consensus_ranking)) {
        cat(sprintf(
          "      %d. %s (mean rank = %.1f)\n",
          j,
          w_all$consensus_ranking$Metastate[j],
          w_all$consensus_ranking$mean_rank[j]
        ))
      }
    }
  }
}

# ============================================================================
# 4. SUMMARY TABLES
# ============================================================================

# Combine all concordance results
concordance_summary <- bind_rows(all_concordance_results)

# Kendall's W Concordance Summary
concordance_summary %>%
  select(resolution, metric, rater, W, p_formatted, n_raters, n_subjects) %>%
  print(n = Inf)

# Concordance strength (all method-dataset combinations)
concordance_summary %>%
  filter(rater == "All combinations") %>%
  mutate(
    strength = case_when(
      W >= 0.9 ~ "very strong",
      W >= 0.7 ~ "strong",
      W >= 0.5 ~ "moderate",
      W >= 0.3 ~ "weak",
      TRUE ~ "very weak"
    )
  ) %>%
  select(resolution, metric, W, strength, p_formatted) %>%
  print(n = Inf)

# Save
write.csv(concordance_summary, "stat_results/Supplementary_Table_S3_strainhub_concordance.csv",
  row.names = FALSE
)

# ============================================================================
# 5. DESCRIPTIVE STATISTICS
# ============================================================================

for (res in unique(all_data$resolution)) {
  cat(sprintf("\n=== %s ===\n", toupper(res)))
  res_data <- all_data %>% filter(resolution == res)

  for (metric in metrics) {
    cat(sprintf("\n  %s:\n", metric_labels[metric]))

    # By ML method
    cat("    By ML method:\n")
    res_data %>%
      group_by(mlt) %>%
      summarise(
        N = n(),
        Mean = round(mean(.data[[metric]], na.rm = TRUE), 4),
        SD = round(sd(.data[[metric]], na.rm = TRUE), 4),
        Median = round(median(.data[[metric]], na.rm = TRUE), 4),
        .groups = "drop"
      ) %>%
      print()

    # By dataset
    cat("    By dataset:\n")
    res_data %>%
      group_by(dataset_base) %>%
      summarise(
        N = n(),
        Mean = round(mean(.data[[metric]], na.rm = TRUE), 4),
        SD = round(sd(.data[[metric]], na.rm = TRUE), 4),
        Median = round(median(.data[[metric]], na.rm = TRUE), 4),
        .groups = "drop"
      ) %>%
      print()
  }
}

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

# Are network statistics consistent across ML methods
# Yes, except for risk4 group
concordance_summary %>%
  filter(rater == "ML tree method") %>%
  group_by(resolution) %>%
  summarise(
    mean_W = round(mean(W), 3),
    min_W = round(min(W), 3),
    max_W = round(max(W), 3),
    all_significant = all(p < 0.05),
    .groups = "drop"
  ) %>%
  print()

concordance_summary %>%
  filter(rater == "ML tree method" & p > 0.05)
# In-degree centrality has W = 0.62 & p=0.06 for risk4

# Across datasets: yes
concordance_summary %>%
  filter(rater == "Dataset") %>%
  group_by(resolution) %>%
  summarise(
    mean_W = round(mean(W), 3),
    min_W = round(min(W), 3),
    max_W = round(max(W), 3),
    all_significant = all(p < 0.05),
    .groups = "drop"
  ) %>%
  print()

# Across all method-dataset combinations: yes
concordance_summary %>%
  filter(rater == "All combinations") %>%
  group_by(resolution) %>%
  summarise(
    mean_W = round(mean(W), 3),
    min_W = round(min(W), 3),
    max_W = round(max(W), 3),
    all_significant = all(p < 0.05),
    .groups = "drop"
  ) %>%
  print()
