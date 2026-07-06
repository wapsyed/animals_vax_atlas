# Skill: r-tidyverse-vaxgo

## Description
R tidyverse-oriented coding style used in the `animals_vax_atlas` project — a comparative transcriptomics pipeline comparing human and mouse vaccine/infection responses using Blood Transcription Modules (BTMs) and MSigDB Hallmarks.

---
## Code Documentation & Annotation Style
Annotations should be written in English as complete sentences, not as shorthand or bullet chains.

### Before each code chunk (in markdown)
Write a short paragraph describing the workflow of the upcoming chunk. Use full sentences, explain the *purpose* and *sequence of operations*, but keep it concise.

**Example:**
```
Load the pre-computed differential expression results and the Blood Transcription
Module annotation table. Join them on gene symbol to assign each gene to its BTM
module, then collapse genes to module-level mean log2FC for each organism and
timepoint. For each module, test whether its mean log2FC differs from zero using
a one-sample t-test with Benjamini-Hochberg correction. Save the resulting
module-level summary for downstream plotting.
```

### Inside code chunks
Every block of code **must** have an inline comment above it. No lines of code should be left uncommented.

Write one or two lines per logical block, explaining *why* a step exists (not just *what* it does). Use `# --------` as a visual separator between blocks.

**Always label save operations with `#Save` before `saveRDS`/`write_csv`/`ggsave`.**

**Example:**
```r
#Save
dge_btm_mean_process %>% 
  saveRDS(file = here("tables", paste0(filename, "_dge_btm_mean_process.rds")))
```

**Example:**
```r
# Summarise gene-level log2FC to module-level means --------
# For each organism x pathogen x timepoint x process:
#   - mean, median, sd, var of the constituent genes' log2FC
#   - n_genes: how many genes from the module were detected
#   - group/subgroup: preserved from BTM annotation (constant per process)
```

---

## Pipe
Always use `%>%` (magrittr), never `|>`.

## Paths
Always `here::here("folder", "file")`. Never use `setwd()` or relative paths without `here`.

## Data Import / Export
| Format | Read | Write |
|--------|------|-------|
| RDS (any R object) | `readRDS(here(...))` | `saveRDS(object, file = here(...))` |
| CSV | `read_csv(here(...))` | `write_csv(., file = here(...))` |
| Excel | `read_excel(here(...), sheet = "...")` | — |
| Delimited | `read_delim(here(...), delim = ";")` | `write_excel_csv2()` |

Typical pipeline: `readRDS` → transform → `saveRDS`, so each script can be re-run independently.

## Packages
- **Core**: `tidyverse` (dplyr, tidyr, ggplot2, purrr, stringr, forcats, readr)
- **Data cleaning**: `janitor` (`clean_names()`, `get_dupes()`)
- **Paths**: `here`
- **Strings**: `glue`
- **Bioinformatics**: `limma`, `GSEA` / `gseGO` (clusterProfiler), `enricher`, `fgsea`, `biomaRt`, `msigdbr`, `GSVA`, `ComplexHeatmap`, `edgeR`, `DESeq2`, `sva`, `pvca`, `variancePartition`
- **Stats**: `rstatix` (for `t_test`, `adjust_pvalue`, `add_significance`), `corrr`, `Hmisc::rcorr`
- **Plotting**: `ggplot2`, `ggrepel`, `ggnewscale`, `ggh4x` (`facet_nested_wrap`), `ggsci`, `ggridges`, `ggdist`, `ggprism`, `ggtext`, `ggbeeswarm`, `GGally`, `patchwork`, `corrr`, `ggcorrplot`
- **Machine learning**: `tidymodels`, `yardstick`
- **Misc**: `beepr` (sound alerts), `pracma`

All packages are listed in `scripts_notebooks/required.R`.

## Data Wrangling (dplyr / tidyr)
- `group_by(...)` → `summarise(...)` → `ungroup()` (always ungroup after summarise)
- `mutate()` with `case_when()` or `if_else()`
- `select()` for column subsetting, rename inside select: `select(new_name = old_name, ...)`
- `distinct()` for deduplication
- `drop_na()` from tidyr
- `fill(...)` for carrying values forward
- `pivot_longer(cols = ..., names_to = "...", values_to = "...")`
- `pivot_wider(names_from = "...", values_from = "...", names_prefix = "...", values_fn = ...)`
- `separate(col, into = c("a", "b"), sep = "...")` — often with `sep = " \\(|\\)"` for parenthetical splits
- `separate_rows(col, sep = "//")` for multi-value delimiters

## Joins
| Join | When |
|------|------|
| `inner_join(a, b, by = "...")` | Keep only matching rows |
| `left_join(a, b, by = "...")` | Keep all rows in a |
| `full_join(a, b, by = "...")` | Keep all rows in both |
| `anti_join(a, b, by = "...")` | Rows in a not in b |
| `bind_rows(a, b)` | Stack vertically |

Uses `join_by(...)` for complex joins (e.g., `join_by("condition", "process")`).

## Factor Handling (forcats)
- `fct_relevel(.f, "level1", "level2", ...)` — set factor level order
- `fct_reorder(.f, .x, .fun = median)` — reorder by another variable
- `fct_rev()` — reverse order
- `fct_inorder()` — order by first appearance
- `fct_relevel(as.factor(x), c(...))`

## String Handling (stringr)
- `str_detect(string, pattern)` — regex detection
- `str_remove(string, pattern)` / `str_remove_all()`
- `str_replace(string, pattern, replacement)` / `str_replace_all()`
- `str_to_lower()`, `str_to_title()`, `str_to_upper()`
- `str_c(...)` — string concatenation
- `str_remove(., "^\\d+\\.\\s")` — leading number removal pattern

## Vector Operations
- `case_when(...)` — multi-condition vectorized if/else (never nested `ifelse`)
- `if_else(condition, true, false)` — typed if/else (prefer over `ifelse`)
- `coalesce(x, y)` — first non-missing

## Namespace Management
Use explicit namespace when conflicts exist:
```r
dplyr::first(), dplyr::count(), dplyr::select()
rstatix::t_test(), rstatix::adjust_pvalue()
```
Common conflict: `first()` from `xts` masks `dplyr::first()`.

## Naming Conventions
- `snake_case` for everything
- Long, descriptive names: `dge_btm_mean_process`, `human_mouse_log2fc_avg_wide`, `degs_human_mouse_overlap_genes_summary`
- Plot objects: `name_plot` suffix (e.g., `volcano_plot`, `heatmap_plot_blood`)
- Contrast levels: `conditionFLUADD.Day.1` (from limma `makeContrasts`)
- Intermediate data: `_df`, `_long`, `_wide`, `_clean` suffixes

## ggplot2
### Pattern
```r
ggplot(data) +
  aes(x = ..., y = ...,
      fill = ..., color = ..., shape = ...) +
  geom_*(...) +
  scale_*_manual(values = named_vector, name = "...") +
  facet_nested_wrap(~var1+var2, scales = "free", nrow = ...,
                    nest_line = element_line(...)) +
  theme_vaxgo() +
  theme(...) +
  labs(x = "...", y = "...", title = "...")
```

### Custom Theme
`theme_vaxgo()` in `required.R`: modifies `theme_minimal()` — removes grid lines, adds axis lines/ticks, sets consistent text sizes.

### Key ggplot2 patterns
- `aes()` **outside** the geom call (not inside)
- `geom_text_repel()` from `ggrepel` for non-overlapping labels
- `new_scale_fill()` / `new_scale_color()` from `ggnewscale` per geom layer
- `stat_summary(fun.data = mean_se, geom = "errorbar")` for grouped distributions
- `geom_tile()` + `geom_text()` for correlation heatmaps
- `facet_nested_wrap(~var1+var2, ...)` from `ggh4x` for nested facet headers
- `scale_fill_gradientn(colors = c("low", "mid", "high"))` for continuous fills
- Color palettes stored as named vectors in `required.R` (e.g., `colors$organism`, `colors$comparison`, `btm_immune_groups$group`)

### Saving
```r
plot_object %>%
  ggsave(filename = here("Figures", "filename.png"),
         width = ..., height = ..., dpi = 300)
```

## Color Palettes (defined in `required.R`)
Key named vectors:
- `colors_vaccines` — vaccine-specific colors
- `colors$organism` — Human/Mouse/FIT
- `colors$comparison` — Shared/Not common
- `colors$treatment` — Vaccination/Infection/Injury
- `colors$timepoint` — time points (day/hour)
- `colors$pathogen` — pathogens
- `btm_immune_groups$group` — BTM immune cell groups
- `immune_colors` — shorter immune group palette
- `hallmarks_colors` — MSigDB Hallmark group colors

## Functions & Helpers (in `required.R`)
| Function | Purpose |
|----------|---------|
| `theme_vaxgo()` | Custom ggplot2 theme |
| `safe_cor_test(x, y, method)` | Safe `cor.test()` with NA handling and n≥3 guard |
| `autoGSEA(df, TERM2GENE, geneset_name)` | Run GSEA across all conditions |
| `gseGO_all(df, OrgDb, ...)` | Run `gseGO` across all conditions |
| `enricher_condition(df, term2gene_df)` | Run `enricher` across all conditions |
| `oraGO_condition(df, orgdb)` | Run `enrichGO` across all conditions |
| `cluster_by_day(df, day_column, condition_column)` | Hierarchical clustering by day |
| `overlap_genes(cond1, cond2, data)` | Overlap analysis between conditions |
| `get_sequences(genes, mart, symbol_attr)` | Retrieve cDNA sequences from biomaRt |

## Project Structure
```
animals_vax_atlas/
├── scripts_notebooks/    ← main Rmd scripts, numbered (0_, 1_, 2_, ...)
│   ├── required.R        ← central: packages, theme, colors, functions
│   ├── 0_Data_Curation.Rmd
│   ├── 1_QualityControl.Rmd
│   ├── 2_Preprocessing.Rmd
│   ├── 3.1_Comparing_Human_Mouse_ByCondition.Rmd
│   ├── 3.2_Comparing_Human_Mouse_UnifiedAnalyses.Rmd
│   ├── 4_Performance.Rmd
│   ├── 5_MachineLearning.Rmd
│   └── ...
├── tables/               ← intermediate RDS and CSV files
├── Figures/               ← output figures
├── DataCuration/          ← curated metadata
├── example/               ← reproducible example scripts
├── ai/                    ← opencode skill + config
└── renv/                  ← locked R environment
```

## R Notebook Structure (Rmd)
- YAML header with `title` and `output: html_document`
- First chunk: `library(here); source(here("scripts_notebooks", "required.R"))`
- Sections separated by markdown headers: `## Section`, `### Subsection`
- Code chunks use `# Comments -----` with dashed separators
- Intermediate saves: `saveRDS(...)` after each major processing step

## Bioinformatics Workflow
1. **Download**: `getGEO()` from GEOquery, ExpressionSet extraction
2. **Annotation**: biomaRt (`getBM`) for probe-to-gene mapping
3. **Normalization**: `normalizeBetweenArrays()` (quantile), `arrayWeights()`
4. **DE Analysis**: limma (`lmFit` → `contrasts.fit` → `eBayes` → `topTreat`)
   - Formula: `~ condition + sex + age + race + 0` (human), `~ timepoint + 0` (mouse)
   - Use `duplicateCorrelation()` with `block = participant_id` for repeated measures
5. **Functional**: `GSEA()`, `gseGO()`, `enricher()`, mean-per-module summaries
6. **Cross-species**: inner_join on human ortholog symbols, rcorr (Spearman)

## Error Handling
- `tryCatch({ ... }, error = function(e) NULL)` around GSEA, enrichment, t-test
- Guard conditions: `if (nrow(data) <= 4)` → skip
- `purrr::possibly()` as alternative pattern
