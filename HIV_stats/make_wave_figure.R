# ============================================================================
# Figure S20: Density plot to compare introduction timing across temporal methods
# HIV-1 Subtype C Phylodynamics Study
# Amanda Perofsky
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

# Load and filter data (same as main analysis)
d <- read.delim("input_data/intro_count_time.tsv")
d <- d %>%
  mutate(
    country_pair = paste0(from, "->", to),
    dataset_base = gsub("\\.\\d+$", "", dataset),
    method_combo = paste(dpt, mlt, sep = "_")
  )

total_combos <- n_distinct(d$method_combo)

complete_within_dataset <- d %>%
  group_by(dataset, country_pair) %>%
  summarise(n_combos = n_distinct(method_combo), .groups = "drop") %>%
  filter(n_combos == total_combos)

d_shared <- d %>%
  semi_join(complete_within_dataset, by = c("dataset", "country_pair"))

range(d_shared$time)
d_shared %>% arrange(time)
100 * nrow(d_shared %>% filter(time<1920))/nrow(d_shared) # < 0.01% of estimated intros occur before 1920

p_density <- ggplot(d_shared, aes(x = time, color = dpt, fill = dpt)) +
  geom_density(alpha = 0.12, adjust = 1.5) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(
    breaks = seq(1920, 2020, 10),
    limits = c(1920, 2020),
    expand = c(0.01,0.01)
  ) +
  scale_y_continuous(expand=c(0.001,0.001))+
  annotate("rect",
    xmin = 1975, xmax = 1985, ymin = -Inf, ymax = Inf,
    alpha = 0.08, fill = "gray40"
  ) +
  labs(
    x = "Introduction time",
    y = "Density",
    color = "Temporal method",
    fill = "Temporal method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.minor = element_blank()
  )
p_density
ggsave("figures/manuscript_fig_s20.png", p_density, width = 8, height = 5, dpi = 300)
ggsave("figures/manuscript_fig_s20.pdf", p_density, width = 8, height = 5)
