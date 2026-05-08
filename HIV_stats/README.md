# Statistical Analyses (`HIV_stats/`)

This folder contains the R code and input data for statistical analyses comparing 16 method combinations (4 ML phylogenetic methods × 4 temporal dating methods) across subsampled HIV-1 Subtype C datasets. Maintained by Amanda Perofsky (a.perofsky@northeastern.edu).

## Summary

This folder covers statistical benchmarking of pre-computed phylogenetic and phylogeographic outputs. The upstream pipelines (sequence curation, alignment, subsampling, ML tree inference, temporal dating, TreeTime mugration, StrainHub centrality) live elsewhere in this repository.

Specifically, this folder produces:

- Non-parametric tests estimating differences in TMRCA, evolutionary rate, and introduction count estimates across methods combinations
- Linear mixed-effects models estimating the effects of temporal x ML methodological choice on introduction timing
- Kendall's W concordance analyses for StrainHub centrality and TreeTime Source-Sink Scores
- Manuscript Supplementary Figures S18, S19, S20 and Supplementary Tables S3, S4

## Folder structure

```
HIV_stats/
├── 01_hiv_phylo_stats.r                  # TMRCA, ER, introduction count statistics + figures
├── 02_hiv_intro_timing_differences.R     # Introduction timing LMM
├── 03_hiv_strainhub_stats.R              # StrainHub concordance (Kendall's W)
├── 04_hiv_treetime_sss_stats.R           # TreeTime Source-Sink Score concordance (Kendall's W)
├── make_wave_figure.R                    # Density plot of introduction times by temporal method
├── input_data/                           # Pre-computed phylogenetic / phylogeographic estimates
├── figures/                              # Figures produced by the scripts
├── stat_results/                         # Supplementary tables produced by scripts 03 and 04
├── HIV_stats.Rproj                       # RStudio project file (optional)
└── README.md
```

All scripts use paths relative to this folder. Set the working directory to `HIV_stats/` (or open `HIV_stats.Rproj` in RStudio) before sourcing.

## Software requirements

- Analyses were run using **R** version 4.5.2
- **`renv`** for environment management (installs the exact package versions used)

To reproduce the analysis environment:
```r
# from inside HIV_stats/
install.packages("renv")        # if not already installed
renv::restore()                  # installs all packages at the versions in renv.lock
```

Opening `HIV_stats.Rproj` in RStudio will activate the project-local library automatically via `.Rprofile`.

Packages used by the scripts (all pinned in `renv.lock`):

- Core: `tidyverse` (loads dplyr, tidyr, ggplot2, stringr, purrr); `dplyr`, `tidyr`, `ggplot2` also called individually in some scripts
- Statistics: `rstatix`, `car`, `PMCMRplus`, `effectsize`, `rcompanion`, `psych`, `moments`
- Mixed-effects models: `lme4`, `lmerTest`, `emmeans`, `performance`, `sjPlot`
- Plotting: `ggplot2`, `ggpubr`, `cowplot`

## Input data

All files in `input_data/` are pre-computed outputs from the upstream phylodynamic pipeline. Files are aggregated across 15 subsampled datasets (5 subsample sizes × 3 replicates: `locrisk260.{1,2,3}`, `locrisk325.{1,2,3}`, ..., `locrisk574.{1,2,3}`) and 16 method combinations (4 ML × 4 temporal dating).

| File | Description | Used by |
|---|---|---|
| `HIV_TMRCA_ER_by_method.csv` | Tree-level TMRCA (date) and substitution rate (subs/site/year × 1000) per dataset × method combination. | `01` |
| `intro_count.tsv` | Per-country pair introduction counts: `from`, `to` (ISO-2 country codes), `dataset`, `n` (count), `dpt` (temporal method), `mlt` (ML method). | `01` |
| `intro_count_time.tsv` | Individual introduction event timing (decimal year) for each `from→to` country pair × dataset × method combination. | `02`, `make_wave_figure.R` |
| `strainhub/All_strainhub_metrics_{continent,region,country,risk4,risk5,risk6}.csv` | StrainHub network centrality metrics (Degree, In-/Out-degree, Betweenness, Closeness, Source-Hub Ratio) per Metastate × dataset × ML method. Resolutions: continent, region, country, risk4–risk6. | `03` |
| `All_treetime_metrics_SSS.tsv` | TreeTime Source-Sink Scores per location: `Export`, `Import`, `SSS = (Export − Import) / (Export + Import)`. Resolutions: continent, region, country, risk4–risk6. | `04` |

Notes:
- The 16 combinations come from crossing `dpt ∈ {LSD, TempEst, TreeTime, treedater}` with `mlt ∈ {FastTree, IQ-TREE, PhyML, RAxML-NG}`. Across the 15 subsampled datasets, this yields 240 unique `dpt × mlt × dataset` combinations.
- Two-letter codes follow ISO 3166-1 alpha-2; risk-group codes follow the manuscript Methods.
- `Method` in `HIV_TMRCA_ER_by_method.csv` corresponds to the temporal dating method (`dpt`).

## Scripts

The four numbered scripts are independent; they do not need to be run in a specific order. Each reads from `input_data/` and writes to `figures/` or `stat_results/`, where applicable. Run each script by setting the working directory to this folder and sourcing the file (e.g., `source("01_hiv_phylo_stats.r")`).

### `01_hiv_phylo_stats.r` — TMRCA, ER, and introduction count statistics

Non-parametric tests (Kruskal–Wallis, Dunn's post-hoc, Scheirer–Ray–Hare for temporal × ML method interactions) and coefficient-of-variation (CV) analyses across temporal × ML method combinations and subsample sizes for:

- **TMRCA** (time to the most recent common ancestor)
- **Evolutionary rate (ER)** (substitutions/site/year)
- **Total introduction counts** per dataset × method combination

**Inputs:** `input_data/HIV_TMRCA_ER_by_method.csv`, `input_data/intro_count.tsv`

**Outputs (`figures/`):** `manuscript_fig_s18.{pdf,png}` (composite temporal x ML methods interaction figure → manuscript **Fig S18**), `manuscript_fig_s19.{pdf,png}` (composite CV-heatmap figure → manuscript **Fig S19**)

**Console output:** Kruskal–Wallis statistics, Dunn's post-hoc contrasts, Scheirer–Ray–Hare interaction tests, descriptive summaries.

### `02_hiv_intro_timing_differences.R` — Linear mixed-effects models of introduction timing

Linear mixed-effects models estimating the effect of temporal x ML method choice on introduction timing estimates:

- **Fixed effects:** temporal dating method (`dpt`), ML phylogenetic method (`mlt`), and their interaction
- **Random effects:** origin country, destination country, country pair, dataset
- Includes sensitivity analyses with alternative random effects structures (e.g., simplifying random effects, collapsing replicate-level random effects to dataset size level)

**Inputs:** `input_data/intro_count_time.tsv`

**Outputs:** Console only: ANOVA tables, EMMs, post-hoc contrasts, variance components. Model results are reported in the manuscript Results section.

### `03_hiv_strainhub_stats.R` — StrainHub concordance analysis

Kendall's coefficient of concordance (W) testing whether StrainHub centrality rankings are stable across ML methods within each dataset, computed for each centrality metric × resolution level (continent / region / country / risk4 / risk5 / risk6).

**Inputs:** `input_data/strainhub/All_strainhub_metrics_*.csv` (six files)

**Outputs (`stat_results/`):** `Supplementary_Table_S3_strainhub_concordance.csv`

### `04_hiv_treetime_sss_stats.R` — TreeTime SSS concordance analysis

Kendall's W concordance testing whether TreeTime Source-Sink Score rankings (Export, Import, SSS) are stable across the 16 temporal × ML method combinations within each dataset, for each resolution level.

**Inputs:** `input_data/All_treetime_metrics_SSS.tsv`

**Outputs (`stat_results/`):** `Supplementary_Table_S4_treetime_sss_concordance.csv`

### `make_wave_figure.R` — Density of introduction times by temporal method

Produces an overlaid density plot of introduction times stratified by temporal dating method, restricted to country-pair × dataset combinations observed across all 16 method combinations.

**Inputs:** `input_data/intro_count_time.tsv`

**Outputs (`figures/`):** `manuscript_fig_s20.{pdf,png}` → manuscript **Fig S20**

## Output directories

- **`figures/`** — All figures listed above in both PDF (vector) and PNG (raster, 300 dpi) format.
- **`stat_results/`** — Concordance analysis summary tables corresponding to manuscript Supplementary Tables S3 and S4.

The committed contents of `figures/` and `stat_results/` correspond to the versions used in the current manuscript. Re-running the scripts will overwrite these with newly generated copies.

## Reproducibility notes

- All analyses are deterministic; no random seeds are required.
- Run times are short (each script completes in well under a minute on a recent laptop).
- Some figures are produced by `cowplot::plot_grid()` after intermediate `ggplot` objects are built within the same script — they cannot be regenerated by sourcing partial sections.
