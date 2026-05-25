# =============================================================================
# Worked Example: Human–Mouse BTM Correlation (Influenza Fluad, Day 7)
# =============================================================================
#
# This script reproduces the cross-species Blood Transcription Module (BTM)
# correlation scatter plot for the Influenza Fluad condition at Day 7 —
# a key result of the manuscript (Figure 3).
#
# The central finding: while gene-level correlations between human and mouse
# are low, BTM-level (pathway-level) responses are highly conserved.
#
# Runtime: < 2 minutes
# Required files (already in tables/):
#   - tables/influenza_fluad_human_mouse_dge_limma_degs.rds
#   - tables/btm_annotation_genes.csv
#
# Output: Figures/example_btm_correlation_day7.png
# =============================================================================


# 1. Setup --------------------------------------------------------------------

library(here)        # project-root-relative paths
library(tidyverse)   # data manipulation + ggplot2
library(ggrepel)     # non-overlapping text labels

# Load shared theme and color palettes defined in required.R
source(here("scripts_notebooks", "required.R"))


# 2. Load data ----------------------------------------------------------------

# Differential expression results: mean log2FC per gene per condition
# Columns: human_symbol, organism, timepoint, condition, mean_l2fc, adj_p_val, ...
degs <- readRDS(here("tables", "influenza_fluad_human_mouse_dge_limma_degs.rds")) |>
  mutate(
    pathogen  = "Influenza",
    vaccine   = "Fluad",
    timepoint = as.numeric(timepoint)
  ) |>
  rename(mean_l2fc = mean, ci_lower = lower, ci_upper = upper)

# BTM annotation: maps gene symbols to Blood Transcription Modules
# Columns: symbol, process (module ID), group (immune cell type), subgroup, ...
btm_genes <- read_csv(here("tables", "btm_annotation_genes.csv"),
                      show_col_types = FALSE)


# 3. Compute mean log2FC per BTM module per organism -------------------------

# Join DEG results with BTM annotation (on human gene symbol)
btm_degs <- degs |>
  rename(symbol = human_symbol) |>
  inner_join(btm_genes, by = "symbol") |>
  filter(organism %in% c("Human", "Mouse"))

# Summarise: mean log2FC across all genes in each BTM module, per organism and
# timepoint. This is the "module-level" expression response.
btm_means <- btm_degs |>
  group_by(organism, timepoint, process, group, subgroup) |>
  summarise(
    mean_log2fc = mean(mean_l2fc, na.rm = TRUE),
    n_genes     = n(),
    .groups     = "drop"
  )


# 4. Reshape for scatter plot (Human vs. Mouse side-by-side) -----------------

# Pivot to wide format so each row is one BTM module at one timepoint,
# with separate columns for human and mouse mean log2FC.
btm_wide <- btm_means |>
  pivot_wider(
    names_from  = organism,
    values_from = mean_log2fc,
    names_prefix = "log2fc_"
  ) |>
  filter(!is.na(log2fc_Human), !is.na(log2fc_Mouse))


# 5. Filter to Day 7 (peak adaptive response) --------------------------------

btm_day7 <- btm_wide |>
  filter(timepoint == 7)


# 6. Compute Spearman correlation ---------------------------------------------

cor_result <- cor.test(
  btm_day7$log2fc_Human,
  btm_day7$log2fc_Mouse,
  method = "spearman"
)

# Format annotation label for the plot
cor_label <- sprintf(
  "Spearman r = %.2f\np = %.3f\nn = %d modules",
  cor_result$estimate,
  cor_result$p.value,
  nrow(btm_day7)
)


# 7. Plot: Human vs. Mouse BTM mean log2FC at Day 7 --------------------------

scatter_plot <- ggplot(btm_day7, aes(x = log2fc_Human, y = log2fc_Mouse)) +
  # Points colored by immune cell group
  geom_point(aes(color = group), size = 2.5, alpha = 0.8) +
  # Reference line (perfect agreement)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  # Horizontal and vertical reference lines at zero
  geom_hline(yintercept = 0, color = "gray80", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "gray80", linewidth = 0.4) +
  # Label the top modules by absolute human response
  geom_text_repel(
    data = btm_day7 |> slice_max(abs(log2fc_Human), n = 8),
    aes(label = process),
    size = 2.5, max.overlaps = 20
  ) +
  # Correlation annotation in the top-left corner
  annotate(
    "text", x = -Inf, y = Inf,
    label = cor_label,
    hjust = -0.1, vjust = 1.3,
    size = 3, color = "black"
  ) +
  # Color palette from required.R
  scale_color_manual(values = btm_immune_groups$group) +
  # Labels
  labs(
    title    = "Human vs. Mouse BTM responses — Influenza Fluad, Day 7",
    subtitle = "Each point = one Blood Transcription Module; dashed line = perfect agreement",
    x        = "Human mean log\u2082FC",
    y        = "Mouse mean log\u2082FC",
    color    = "BTM group"
  ) +
  theme_vaxgo() +
  theme(legend.position = "right")

scatter_plot


# 8. Save figure --------------------------------------------------------------

ggsave(
  filename = here("example", "example_btm_correlation_day7.png"),
  plot     = scatter_plot,
  width    = 9,
  height   = 6,
  dpi      = 300
)

message("Figure saved to: example/example_btm_correlation_day7.png")


## Worked Example

#The script [`example/example_btm_correlation.R`](example/example_btm_correlation.R) reproduces the cross-species BTM correlation scatter plot (manuscript Figure 3) using only pre-computed files already present in `tables/`. No GEO download required. Runtime < 2 minutes.

source(here::here("example", "example_btm_correlation.R"))


# What it does:
# 
# 1. Loads limma DEG results (`influenza_fluad_human_mouse_dge_limma_degs.rds`) and BTM annotations
# 2. Computes mean log₂FC per BTM module for human and mouse separately
# 3. Plots human vs. mouse module responses at Day 7 with Spearman *r* annotation
# 4. Saves the figure to `Figures/example_btm_correlation_day7.png`

------------------------------------------------------------------------

## Key Data Files

# For full variable definitions, table schemas, file naming conventions, and color palettes, see the [CODEBOOK.md](CODEBOOK.md).