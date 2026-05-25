# Animals Vax Atlas

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![R version](https://img.shields.io/badge/R-%3E%3D4.5.2-276DC3?logo=r&logoColor=white)](https://cran.r-project.org/) [![renv](https://img.shields.io/badge/reproducibility-renv-blue)](https://rstudio.github.io/renv/) [![Journal](https://img.shields.io/badge/Genes%20%26%20Immunity-Under%20Review-orange)](https://www.nature.com/gi/)

> **Associated manuscript:**\
> *From Mice to Humans: Functional Modules Improve the Translatability of Transcriptomic Responses*\
> Wasim Aluísio Prates-Syed, Aline A. Lira, Nelson Cortes, Jaqueline D.Q. Silva, Bárbara Hamaguchi, Evelyn Carvalho, Adriana Castillo-Chávez, Ricardo Durães-Carvalho, Otavio Cabral-Marques, Ester Cerdeira Sabino, José Eduardo Krieger, Thomas Hagan, Gustavo Cabral-Miranda.\
> *Submitted to Genes and Immunity — currently under review.*

------------------------------------------------------------------------

## Overview

Mice are the dominant preclinical model in vaccine research, yet their translational value for human immune responses remains contested. This project systematically evaluates murine translatability across vaccination (Influenza, Hepatitis B), acute bacterial infection (*S. aureus*, *E. coli*), and sterile injury (burns and trauma) using publicly available blood transcriptome data from GEO/BioProject.

A key finding is that while individual orthologous gene correlations are often low, **higher-order pathway and module responses are highly conserved** between species. By shifting from a gene-level to a pathway-level analytical framework — using Blood Transcription Modules (BTMs) and rank-based statistics — murine models accurately predict human immune dynamics. Translational accuracy scales with stimulus intensity: acute infections and systemic injuries show the highest conservation, while milder vaccination stimuli reveal greater species-specific divergence. Divergent gene expression is linked to differences in *cis*-regulatory architecture, not protein sequence identity.

------------------------------------------------------------------------

## Conditions Covered

| Challenge                         | Vaccine / Agent    | Organisms    |
|-----------------------------------|--------------------|--------------|
| Influenza                         | Fluad (TIV + MF59) | Human, Mouse |
| Hepatitis B                       | Engerix B          | Human, Mouse |
| *Staphylococcus aureus* infection | —                  | Human, Mouse |
| *Escherichia coli* infection      | —                  | Human, Mouse |
| Burn                              | —                  | Human, Mouse |
| Trauma                            | —                  | Human, Mouse |

------------------------------------------------------------------------

## Analysis Workflow

![Flowchart](diagram_animal.png)

The notebooks are designed to be run in order:

1.  **`0_Data_Curation.Rmd`** — Searches BioProject for vaccination/immunization studies, filters by organism and design criteria, and produces a curated dataset list.

2.  **`1_QualityControl.Rmd`** — Runs array quality metrics (AQM) to flag low-quality samples and outliers before integration.

3.  **`2_Preprocessing.Rmd`** — Downloads ExpressionSets from GEO via `GEOquery`, applies quantile normalization, maps probes to gene symbols via `biomaRt`, and computes sample-level log2 fold-changes vs. Day 0 baseline.

4.  **`3_Comparing_Human_Mouse.Rmd`** — Core analysis:

    - Differential expression with `limma`
    - GSEA and ssGSEA using BTMs and MSigDB Hallmarks
    - Cross-species correlation of module NES scores and mean log2FCs
    - RRHO2 rank-rank hypergeometric overlap
    - Shared/unique DEG analysis (Fisher/Jaccard statistics)

5.  **`4_Performance.Rmd`** — Evaluates mouse-to-human predictive performance via ROC curves across all conditions, stratified by all genes vs. immune gene subsets.

6.  **`5_MachineLearning.Rmd`** — Builds machine learning classifiers for cross-species response prediction.

------------------------------------------------------------------------

## Repository Structure

```         
animals_vax_atlas/
├── scripts_notebooks/
│   ├── required.R                   # Packages, theme, utility functions, palettes
│   ├── 0_Data_Curation.Rmd
│   ├── 1_QualityControl.Rmd
│   ├── 2_Preprocessing.Rmd
│   ├── 3_Comparing_Human_Mouse.Rmd
│   ├── 4_Performance.Rmd
│   ├── 5_MachineLearning.Rmd
│   └── FIT_training_datasets.Rmd
├── DataCuration/                    # Raw BioProject exports and manual annotation files
├── tables/                          # Processed data (RDS, CSV) — main data store
├── VaxGO/                           # Gene set files (BTMs, ImmuneGO, VaxSigDB, Hallmarks)
├── Genomic/                         # Promoter sequences (FASTA)
├── Figures/                         # Generated plots
├── Figures_Article/                 # Publication-ready figures
├── ArrayQM/                         # AQM HTML reports
└── renv.lock                        # Package snapshot
```

------------------------------------------------------------------------

## Gene Sets Used

| Gene Set | Description |
|----|----|
| **BTMs** | Blood Transcription Modules (Li et al.) — immune cell-type and process modules |
| **MSigDB Hallmarks** | Broad hallmark gene sets from MSigDB |
| **ImmuneGO** | Custom mouse-adapted immune Gene Ontology annotations |
| **VaxSigDB** | Vaccination signature gene sets |

------------------------------------------------------------------------

## Software Requirements

| Component | Version |
|-----------|---------|
| R | 4.5.2 |
| Bioconductor | 3.22 |
| OS | Linux (tested on Ubuntu/Zorin); macOS and Windows (WSL2) expected to work |
| RAM | ≥ 16 GB recommended (large expression matrices in `tables/`) |
| Disk | ≥ 5 GB for `tables/` |

All R package versions are pinned via `renv.lock`. Key packages:

| Package | Version | Source |
|---------|---------|--------|
| tidyverse | 2.0.0 | CRAN |
| limma | 3.66.0 | Bioconductor 3.22 |
| GEOquery | 2.78.0 | Bioconductor 3.22 |
| clusterProfiler | 4.18.4 | Bioconductor 3.22 |
| fgsea | 1.36.2 | Bioconductor 3.22 |
| GSVA | 2.4.4 | Bioconductor 3.22 |
| biomaRt | 2.66.1 | Bioconductor 3.22 |
| ComplexHeatmap | 2.26.1 | Bioconductor 3.22 |
| RRHO2 | — | GitHub (RRHO2/RRHO2) |

> **Note:** A Docker image is planned for a future release. For now, `renv` provides full package-level reproducibility.

------------------------------------------------------------------------

## Data Acquisition

All **pre-processed intermediate files** are already provided in `tables/`, so for most analyses no raw data download is needed. Users who wish to re-run preprocessing from scratch can download the original datasets from GEO using the accessions below.

| Condition | Organism | GEO Accession | Platform |
|-----------|----------|--------------|----------|
| Influenza (Fluad) + Hepatitis B (Engerix B) | Mouse | [GSE120661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120661) | Agilent 8×60K |
| Influenza (Fluad) | Human | [GSE124689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124689) | Illumina HumanHT-12 |
| Hepatitis B (Engerix B) | Human | [GSE124533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124533) | Illumina HumanHT-12 |
| *S. aureus* infection | Human | [GSE19668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19668) | Affymetrix HuGene |
| *E. coli* infection | Human | [GSE33341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33341) | Affymetrix HuGene |
| Trauma | Human | [GSE36809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36809) | Affymetrix HuGene |
| Burn + Quadrivalent vaccine | Mouse | [GSE182858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182858) | Illumina MouseWG-6 |

Download example (run once; outputs saved to `tables/`):

```r
library(GEOquery)
library(here)

gse <- getGEO("GSE120661", GSEMatrix = TRUE)
eset <- gse[[1]]

# Save expression matrix and metadata
exprs(eset) |> as.data.frame() |> saveRDS(here("tables", "GSE120661_exprs.rds"))
Biobase::pData(eset) |> write.csv(here("tables", "GSE120661_metadata.csv"))
```

------------------------------------------------------------------------

## Reproducing the Analysis

### Prerequisites

- R ≥ 4.5.2 installed ([download](https://cran.r-project.org/))
- RStudio or any R IDE
- Internet access only required for GEO downloads (`2_Preprocessing.Rmd`) and `biomaRt` queries

### Step 0 — Clone and restore the environment

```bash
git clone https://github.com/wapsyed/animals_vax_atlas.git
cd animals_vax_atlas
```

```r
# In R, from the project root:
renv::restore()   # installs all packages at exact recorded versions
```

### Step 1–6 — Run notebooks in order

Each notebook sources `scripts_notebooks/required.R`, which loads all packages and defines the shared `theme_vaxgo` ggplot2 theme, color palettes, and utility functions.

> Chunks marked `eval=FALSE` in the notebooks correspond to one-time download or heavy computation steps. Pre-computed outputs are already in `tables/` and can be loaded directly.

| Step | Notebook | Key Inputs | Key Outputs | Est. time |
|------|----------|------------|-------------|-----------|
| 0 | `0_Data_Curation.Rmd` | `tables/animals_vaccines_bioproject_result.txt` | `tables/datacuration_step2.csv`, `DataCuration/AnimalVax_DataCuration_Annotated.csv` | ~10 min |
| 1 | `1_QualityControl.Rmd` | `tables/*_eset.rds`, `tables/*_metadata.rds` | `ArrayQM/` reports, QC plots in `Figures/` | ~20 min |
| 2 | `2_Preprocessing.Rmd` | GEO downloads or pre-saved `tables/*_exprs.rds` | `tables/*_eset.rds`, `tables/*_log2fc_sample_clean_long.rds` | ~60 min |
| 3 | `3_Comparing_Human_Mouse.Rmd` | `tables/*_dge_limma_degs.rds`, `tables/btm_annotation_genes.csv`, `tables/msigdb_hallmarks_grouped_genes.csv` | `tables/*_gsea_btm_results.rds`, `tables/*_samples_btm_correlation_all_timepoints_summary.rds`, figures | ~90 min |
| 4 | `4_Performance.Rmd` | `tables/all_human_mouse_metadata.rds`, `tables/all_fit_prediction_results.rds` | `tables/roc_curve_data_*.rds`, ROC figures | ~30 min |
| 5 | `5_MachineLearning.Rmd` | `tables/AllData_FIT_training_blood.rds` | Model objects, prediction outputs | ~20 min |

------------------------------------------------------------------------

## Worked Example

The script [`example/example_btm_correlation.R`](example/example_btm_correlation.R) reproduces the cross-species BTM correlation scatter plot (manuscript Figure 3) using only pre-computed files already present in `tables/`. No GEO download required. Runtime < 2 minutes.

```r
source(here::here("example", "example_btm_correlation.R"))
```

What it does:

1. Loads limma DEG results (`influenza_fluad_human_mouse_dge_limma_degs.rds`) and BTM annotations
2. Computes mean log₂FC per BTM module for human and mouse separately
3. Plots human vs. mouse module responses at Day 7 with Spearman *r* annotation
4. Saves the figure to `Figures/example_btm_correlation_day7.png`

------------------------------------------------------------------------

## Key Data Files

For full variable definitions, table schemas, file naming conventions, and color palettes, see the [CODEBOOK.md](CODEBOOK.md).

Raw and intermediate data live in `tables/`. Key files:

| File | Description |
|----|----|
| `all_human_mouse_metadata.rds` | Combined sample metadata across all conditions |
| `all_degs_human_mouse.rds` | Pooled differential expression results |
| `btm_annotation_genes.csv` | BTM gene set annotations |
| `msigdb_hallmarks_grouped_genes.csv` | Hallmark gene set annotations |
| `influenza_fluad_human_mouse_eset.rds` | Example integrated ExpressionSet (Fluad) |

Processed files follow the naming convention:\
`{pathogen}_{condition}_{organism}_{datatype}.rds`

------------------------------------------------------------------------

## Citation

If you use this code or data, please cite:

> Prates-Syed WA, Lira AA, Cortes N, Silva JDQ, Hamaguchi B, Carvalho E, Castillo-Chávez A, Durães-Carvalho R, Cabral-Marques O, Sabino EC, Krieger JE, Hagan T, Cabral-Miranda G. *From Mice to Humans: Functional Modules Improve the Translatability of Transcriptomic Responses.* Genes and Immunity (under review).

------------------------------------------------------------------------

## License

See [LICENSE](LICENSE).
