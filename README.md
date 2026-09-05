# Animals Vax Atlas

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![R version](https://img.shields.io/badge/R-%3E%3D4.5.2-276DC3?logo=r&logoColor=white)](https://cran.r-project.org/) [![renv](https://img.shields.io/badge/reproducibility-renv-blue)](https://rstudio.github.io/renv/) [![Journal](https://img.shields.io/badge/Genes%20%26%20Immunity-Under%20Review-orange)](https://www.nature.com/gi/)

> **Associated manuscript:**\
> *From Mice to Humans: Functional Modules Improve the Translatability of Transcriptomic Responses*\
> Wasim Aluísio Prates-Syed, Aline A. Lira, Nelson Cortes, Jaqueline D.Q. Silva, Bárbara Hamaguchi, Evelyn Carvalho, Adriana Castillo-Chávez, Ricardo Durães-Carvalho, Otavio Cabral-Marques, Ester Cerdeira Sabino, José Eduardo Krieger, Thomas Hagan, Gustavo Cabral-Miranda.\
> *Submitted to Genes and Immunity — currently under review.*

------------------------------------------------------------------------

## Overview

Mice are the dominant preclinical model in vaccine research, yet their translational value for human immune responses remains contested. This project systematically evaluates murine translatability across vaccination (Influenza, Hepatitis B), acute bacterial infection (*S. aureus*, *E. coli*), and sterile injury (burns and trauma) using publicly available blood transcriptome data from GEO and BioProject.

A central finding of this project is that while individual orthologous gene correlations are often weak, **higher-order pathway and module responses are highly conserved** between species. By shifting from a gene-centric to a pathway-level analytical framework—using Blood Transcription Modules (BTMs), MSigDB Hallmarks, and rank-based statistics—murine models accurately predict human immune dynamics. Translational concordance scales with stimulus intensity: acute infections and systemic injuries show the highest cross-species conservation, whereas milder vaccination stimuli exhibit greater species-specific divergence. Furthermore, divergent gene expression is primarily governed by divergence in *cis*-regulatory promoter architecture rather than protein-coding sequence identity.

------------------------------------------------------------------------

## Conditions Covered

| Challenge | Vaccine / Agent | Organisms | Human GEO | Mouse GEO | Platforms |
|:---|:---|:---|:---|:---|:---|
| Influenza | Fluad (TIV + MF59) | Human, Mouse | GSE124689 | GSE120661 | Illumina HumanHT-12, Agilent 8×60K |
| Hepatitis B | Engerix B | Human, Mouse | GSE124533 | GSE120661 | Illumina HumanHT-12, Agilent 8×60K |
| *Staphylococcus aureus* bacteremia | — | Human, Mouse | GSE19668 | GSE120661 | Affymetrix HuGene 1.0 ST, Agilent 8×60K |
| *Escherichia coli* sepsis | — | Human, Mouse | GSE33341 | GSE120661 | Affymetrix HuGene 1.0 ST, Agilent 8×60K |
| Burn injury | — | Human, Mouse | Clinical cohort | GSE182858 | Custom array, Illumina MouseWG-6 v2.0 |
| Trauma | — | Human | GSE36809 | — | Affymetrix HuGene 1.0 ST |
| Burn + Quadrivalent vaccine | — | Mouse | — | GSE182858 | Illumina MouseWG-6 v2.0 |

------------------------------------------------------------------------

## Analysis Workflow

![Flowchart](diagram_animal.png)

The computational pipeline is structured into 9 modular R Markdown notebooks designed to be executed sequentially:

1.  **`0_Data_Curation.Rmd`** — Programmatically scans and filters raw BioProject metadata from NCBI. Isolates time-course vaccination and infection studies, applying inclusion/exclusion criteria to remove oncology, autoimmune, or toxicology studies.
2.  **`1_QualityControl.Rmd`** — Evaluates data fidelity using Array Quality Metrics (`arrayQualityMetrics`) and Relative Log Expression (RLE) distributions to identify sample-level outliers and technical variation.
3.  **`2_Preprocessing_and_DGE.Rmd`** — Downloads ExpressionSets via `GEOquery`, normalizes array intensities (RMA for Affymetrix; Quantile normalization via `limma` for Illumina/Agilent), resolves probe redundancy by selecting the **probe with the maximum variance across samples**, and models differential expression with `limma` Empirical Bayes moderation (`adj. p-value <= 0.05`).
4.  **`3.1_Comparing_Human_Mouse_DGE_analyses.Rmd`** — Executes cross-species gene-level comparative analyses. Computes macroevolutionary effect size delta ($\Delta \text{log}_2\text{FC} = \text{log}_2\text{FC}_H - \text{log}_2\text{FC}_M$), coefficient of variation (CV), sampling stability from downsampling, and inverse-variance statistical weights ($1 / (SE_H^2 + SE_M^2)$).
5.  **`3.2_Comparing_Human_Mouse_GSEA.Rmd`** — Consolidates the multi-condition DGE data and runs unified pathway-level Gene Set Enrichment Analysis via `fgsea` on Blood Transcription Modules (BTMs) and MSigDB Hallmarks, using an exploratory threshold of $\text{padj} \le 0.25$, alongside single-sample GSEA (`GSVA/ssGSEA`).
6.  **`3.3_Comparing_Human_Mouse_Functional_Analyses.Rmd`** — Evaluates higher-order functional conservation. Generates module-level NES and mean log₂FC cross-species correlations over time (**Figure 43a**), quantifies shared vs. species-specific leading-edge genes (LEGs) (**Figure 43c**), and plots rank conservation for core modules such as "immune activation - generic cluster" (**Figure 43d**).
7.  **`4_Performance_EqualTImepoints.Rmd` & `4_Performance_DifferentTimepoints.rmd`** — Assesses murine predictive power for human module regulation. Generates ROC curves and computes Area Under the Curve (AUC) (**Figure 43b**) for matched (equal) timepoints and cross-temporal (different) timepoints, benchmarked against biological controls (e.g., Duchenne Muscular Dystrophy, DMD) and permutation null distributions.
8.  **`5.1_EvolutionaryAnalysis_Protein.Rmd` & `5.2_EvolutionaryAnalysis_Regulation.Rmd`** — Dissects evolutionary determinants. Retrieves Ensembl BioMart coding sequences (CDS) and amino acid identity %, computes codon-level pairwise alignment and **Kimura 2-Parameter (K80) genetic distances**, and integrates ENCODE candidate Cis-Regulatory Elements (cCREs: PLS, pELS, dELS, and CTCF-bound sites) across GRCh38 and mm10 to assess promoter conservation.
9.  **`6_Statistical_Modelling.Rmd`** — Builds multi-modal machine learning workflows using `tidymodels` (Random Forest via `ranger`, Elastic Net) combining coding sequence distance (`dist_k80`), amino acid identity, transcription factor networks, and promoter cCRE structures to model the genomic determinants of translatability.

------------------------------------------------------------------------

## Repository Structure

``` text
animals_vax_atlas/
├── scripts_notebooks/
│   ├── required.R                           # Global libraries, theme_vaxgo, palettes, utility functions
│   ├── 0_Data_Curation.Rmd                  # BioProject curation and filtering
│   ├── 1_QualityControl.Rmd                 # ArrayQM and RLE quality control
│   ├── 2_Preprocessing_and_DGE.Rmd          # GEO download, normalization, probe collapsing, limma DGE
│   ├── 2_Preprocessing_and_DGE_Simplified.Rmd # Streamlined preprocessing and DGE pipeline
│   ├── 3.1_Comparing_Human_Mouse_DGE_analyses.Rmd # Cross-species DGE comparison, noise & delta metrics
│   ├── 3.2_Comparing_Human_Mouse_GSEA.Rmd   # fgsea & ssGSEA unified pipeline (BTMs, Hallmarks)
│   ├── 3.3_Comparing_Human_Mouse_Functional_Analyses.Rmd # Functional correlations, LEGs, rank conservation (Fig 43a,c,d)
│   ├── 4_Performance_EqualTImepoints.Rmd    # Equal-timepoint ROC/AUC classification (Fig 43b)
│   ├── 4_Performance_DifferentTimepoints.rmd # Cross-temporal ROC/AUC benchmarking with controls
│   ├── 5.1_EvolutionaryAnalysis_Protein.Rmd # Protein sequence identity and Kimura K80 CDS distance
│   ├── 5.2_EvolutionaryAnalysis_Regulation.Rmd # ENCODE cCRE promoter/enhancer regulatory architecture
│   ├── 6_Statistical_Modelling.Rmd          # tidymodels predictive modeling of translatability drivers
│   └── FIT_training_datasets.Rmd            # Found In Translation (FIT) benchmarking
├── tables/                                  # Intermediate and processed RDS/CSV data files
│   ├── DataCuration/                        # BioProject search outputs and curation tables
│   ├── Genomic/                             # ENCODE cCRE BED files (PLS, pELS, dELS, CTCF-bound)
│   └── VaxGO/                               # Curated BTM, ImmuneGO, and Hallmark gene set definitions
├── example/
│   ├── example_btm_correlation.R            # Minimal reproducible script for cross-species BTM correlation
│   └── example_btm_correlation_day7.png     # Example output plot
├── Figures/                                 # Generated exploratory and diagnostic figures
├── Figures_Article/                         # High-resolution, publication-ready figures
├── ArrayQM/                                 # ArrayQualityMetrics HTML report directories
└── renv.lock                                # Pinned R dependency environment snapshot
```

------------------------------------------------------------------------

## Gene Sets Used

| Gene Set | Description | Source | Reference |
|:---|:---|:---|:---|
| **BTMs** | Blood Transcription Modules (346 consensus modules) | Li et al. | *Nat Immunol* 2014, 2021 |
| **MSigDB Hallmarks** | 50 well-defined hallmark biological processes | Broad Institute | Liberzon et al., *Cell Syst* 2015 |
| **ImmuneGO** | Mouse-adapted immune Gene Ontology annotations | VaxGO | Custom curated |
| **VaxSigDB** | Curated vaccination response signatures | VaxGO | Custom curated |

------------------------------------------------------------------------

## Software & Reproducibility Environment

All package dependencies are managed via `renv`. Pinned core specifications:

| Component | Version | Description |
|:---|:---|:---|
| **R** | ≥ 4.5.2 | Base language environment |
| **Bioconductor** | 3.22 | Genomic and microarray annotation suites |
| **Operating System** | Linux (Ubuntu/Zorin); macOS and Windows (WSL2) compatible | Tested on 64-bit Linux |
| **Hardware** | ≥ 16 GB RAM recommended | High-dimensional expression matrices |

### Key Packages Pinning

| Package | Version | Source | Key Usage |
|:---|:---|:---|:---|
| **tidyverse** | 2.0.0 | CRAN | Data wrangling, piping, and visualization |
| **limma** | 3.66.0 | Bioconductor | Linear modeling, empirical Bayes moderation, quantile normalization |
| **GEOquery** | 2.78.0 | Bioconductor | Programmatic retrieval of GEO datasets |
| **fgsea** | 1.36.2 | Bioconductor | Fast Gene Set Enrichment Analysis |
| **GSVA** | 2.4.4 | Bioconductor | Single-sample gene set enrichment (ssGSEA) |
| **biomaRt** | 2.66.1 | Bioconductor | Cross-species orthology and Ensembl sequence retrieval |
| **pwalign** | Bioconductor 3.22 | Bioconductor | Pairwise global sequence alignment |
| **ape** | 5.8 | CRAN | DNAbin conversion and Kimura K80 distance calculation |
| **tidymodels** | 1.2.0 | CRAN | Machine learning recipes, workflows, and evaluation |
| **ranger** | 0.16.0 | CRAN | High-performance Random Forest implementation |
| **ComplexHeatmap** | 2.26.1 | Bioconductor | High-dimensional heatmap visualizations |
| **pROC** | 1.18.5 | CRAN | ROC curve and AUC generation |

------------------------------------------------------------------------

## Data Acquisition

All pre-processed intermediate files are archived in `tables/`, allowing downstream analyses (steps 3–6) to run without re-downloading raw files. Users wishing to replicate preprocessing from scratch can query the original accessions:

| Condition | Organism | Accession | Platform | Platform ID |
|:---|:---|:---|:---|:---|
| Influenza (Fluad) + Hepatitis B | Mouse | [GSE120661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120661) | Agilent 8×60K | GPL21103 |
| Influenza (Fluad) | Human | [GSE124689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124689) | Illumina HumanHT-12 v4.0 | GPL10558 |
| Hepatitis B (Engerix B) | Human | [GSE124533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124533) | Illumina HumanHT-12 v4.0 | GPL10558 |
| *S. aureus* infection | Human | [GSE19668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19668) | Affymetrix Human Gene 1.0 ST | GPL6244 |
| *E. coli* infection | Human | [GSE33341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33341) | Affymetrix Human Gene 1.0 ST | GPL6244 |
| Trauma | Human | [GSE36809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36809) | Affymetrix Human Gene 1.0 ST | GPL6244 |
| Burn + Quadrivalent vaccine | Mouse | [GSE182858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182858) | Illumina MouseWG-6 v2.0 | GPL6887 |

------------------------------------------------------------------------

## Reproducing the Analysis

### Step 0 — Setup Environment

``` bash
git clone https://github.com/wapsyed/animals_vax_atlas.git
cd animals_vax_atlas
```

``` r
# In R console:
renv::restore()   # Restores exact package environment
```

### Notebook Execution Guide

Each notebook sources `scripts_notebooks/required.R`, initializing the shared workspace, custom ggplot2 themes (`theme_vaxgo`), palettes, and utility functions.

| Step | Notebook | Key Inputs | Key Outputs | Est. Time |
|:---|:---|:---|:---|:---|
| **0** | `0_Data_Curation.Rmd` | `tables/DataCuration/animals_vaccines_bioproject_result.csv` | `tables/DataCuration/datacuration_step2.csv` | \~10 min |
| **1** | `1_QualityControl.Rmd` | `tables/*_eset.rds`, `tables/*_metadata.rds` | `ArrayQM/` reports, RLE plots | \~20 min |
| **2** | `2_Preprocessing_and_DGE.Rmd` | Raw GEO ExpressionSets or `tables/*_exprs.rds` | `tables/*_dge_limma_degs.rds`, `tables/*_log2fc_sample_clean_long.rds` | \~60 min |
| **3.1** | `3.1_Comparing_Human_Mouse_DGE_analyses.Rmd` | `tables/all_human_mouse_dge_limma_degs.rds` | `tables/human_mouse_log2fc_avg_wide_all.rds`, divergence weights | \~30 min |
| **3.2** | `3.2_Comparing_Human_Mouse_GSEA.Rmd` | `all_human_mouse_dge_limma_degs_matched_control.rds`, BTM & Hallmark CSVs | `tables/all_human_mouse_gsea_btm_results.rds`, `tables/*_gsea_mean_wide.rds` | \~45 min |
| **3.3** | `3.3_Comparing_Human_Mouse_Functional_Analyses.Rmd` | `all_human_mouse_gsea_btm_results.rds`, `all_human_mouse_gsea_btm_legs.rds` | Module correlation over time (**Fig 43a**), LEG barplots (**Fig 43c**), Rank conservation (**Fig 43d**) | \~45 min |
| **4** | `4_Performance_EqualTImepoints.Rmd` & `4_Performance_DifferentTimepoints.rmd` | `all_human_mouse_dge_limma_degs_matched_filtered.rds`, BTM annotations | ROC curves, AUC summary tables (**Fig 43b**), PR curves | \~40 min |
| **5.1** | `5.1_EvolutionaryAnalysis_Protein.Rmd` | Ensembl BioMart CDS data, `all_alignments`, `all_human_mouse_gsea_btm_legs.rds` | `human_mouse_cds_distance.rds` (Kimura K80), protein identity vs $\Delta\text{log}_2\text{FC}$ | \~50 min |
| **5.2** | `5.2_EvolutionaryAnalysis_Regulation.Rmd` | ENCODE cCRE BED files (`tables/Genomic/*`), gene TSS coords | `cres_type_homology_comparison_wide.rds`, promoter conservation plots | \~40 min |
| **6** | `6_Statistical_Modelling.Rmd` | `human_mouse_statsmodelling_parameters_values.rds` | `tidymodels` Random Forest & Elastic Net models, VIP feature importance | \~30 min |

------------------------------------------------------------------------

## Minimal Worked Example

The standalone script [`example/example_btm_correlation.R`](example/example_btm_correlation.R) reproduces the cross-species BTM correlation scatter plot using pre-computed tables in `< 2 minutes`:

``` r
source(here::here("example", "example_btm_correlation.R"))
```

Outputs are saved directly to `Figures/example_btm_correlation_day7.png`.

------------------------------------------------------------------------

## Citation

If you use this code or data, please cite:

> Prates-Syed WA, Lira AA, Cortes N, Silva JDQ, Hamaguchi B, Carvalho E, Castillo-Chávez A, Durães-Carvalho R, Cabral-Marques O, Sabino EC, Krieger JE, Hagan T, Cabral-Miranda G. *From Mice to Humans: Functional Modules Improve the Translatability of Transcriptomic Responses.* Genes and Immunity (under review).

------------------------------------------------------------------------

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
