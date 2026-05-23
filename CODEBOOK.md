# Codebook — Animals Vax Atlas

This document describes the naming conventions, abbreviations, table schemas, and color palettes used throughout the project.

---

## Table of Contents

1. [Abbreviations & Acronyms](#1-abbreviations--acronyms)
2. [File Naming Conventions](#2-file-naming-conventions)
3. [Table Schemas](#3-table-schemas)
   - 3.1 [Sample metadata](#31-sample-metadata-_metadata)
   - 3.2 [Differential expression results](#32-differential-expression-results-_dge_limma_degs)
   - 3.3 [BTM module means](#33-btm-module-means-_dge_btm_mean_process)
   - 3.4 [GSEA results](#34-gsea-results-_gsea_btm_results--_gsea_hallmarks_results)
   - 3.5 [Human–mouse correlation](#35-humanmouse-correlation-_samples_btm_correlation_all_timepoints_summary)
   - 3.6 [Gene overlap (Fisher/Jaccard)](#36-gene-overlap-_sharedgenes_fisher_jaccard)
   - 3.7 [FIT predictor outputs](#37-fit-predictor-outputs-_fit_prediction_results)
   - 3.8 [BTM annotation](#38-btm-annotation-btm_annotation_genescsv)
   - 3.9 [Hallmark gene sets](#39-hallmark-gene-sets-msigdb_hallmarks_grouped_genescsv)
   - 3.10 [Human–mouse orthologs](#310-humanmouse-orthologs-homologous_mgi_hgncsv)
4. [BTM Immune Groups](#4-btm-immune-groups)
5. [Vaccine Types](#5-vaccine-types)
6. [Color Palettes](#6-color-palettes-defined-in-requiredr)

---

## 1. Abbreviations & Acronyms

| Abbreviation | Meaning |
|--------------|---------|
| **AUC** | Area Under the Curve |
| **BTM** | Blood Transcription Module |
| **CI** | Confidence Interval |
| **CONJ** | Conjugate vaccine |
| **DEG** | Differentially Expressed Gene |
| **DGE** | Differential Gene Expression |
| **FDR** | False Discovery Rate |
| **FIT** | Found In Translation (mouse-to-human predictor) |
| **GDS** | GEO DataSet |
| **GEO** | Gene Expression Omnibus |
| **GPL** | GEO Platform Accession |
| **GSE** | GEO Series Accession |
| **GSEA** | Gene Set Enrichment Analysis |
| **HGNC** | Human Gene Nomenclature Committee |
| **IM** | Intramuscular |
| **IN** | Intranasal |
| **LEG** | Leading-Edge Genes (from GSEA enrichment) |
| **MGI** | Mouse Genome Informatics |
| **NES** | Normalized Enrichment Score |
| **ORA** | Over-Representation Analysis |
| **PBMC** | Peripheral Blood Mononuclear Cell |
| **PCA** | Principal Component Analysis |
| **QN** | Quantile Normalized |
| **RLE** | Relative Log Expression |
| **ROC** | Receiver Operating Characteristic |
| **RRHO2** | Rank-Rank Hypergeometric Overlap (version 2) |
| **SE** | Standard Error |
| **ssGSEA** | Single-Sample Gene Set Enrichment Analysis |
| **VLP** | Virus-Like Particles |

---

## 2. File Naming Conventions

### General pattern

```
{pathogen/condition}_{vaccine}_{organism(s)}_{data_type}.{extension}
```

**Example:**
```
influenza_fluad_human_mouse_dge_limma_degs.rds
│          │     │            │              └─ R serialized format
│          │     │            └─ data type: DEG results (limma)
│          │     └─ organisms: human and mouse
│          └─ vaccine: Fluad (TIV + MF59)
└─ pathogen: influenza
```

---

### Condition/pathogen prefixes

| Prefix | Condition |
|--------|-----------|
| `influenza_fluad_` | Influenza vaccine — Fluad (TIV + MF59) |
| `hepatitisb_engerixb_` | Hepatitis B vaccine — Engerix-B |
| `ecoli_infection_` | Experimental *E. coli* infection |
| `saureus_infection_` | Experimental *S. aureus* infection |
| `burn_` | Burn injury (humans) |
| `trauma_` | Surgical/mechanical trauma |
| `all_human_mouse_` | Aggregated data across all conditions |
| `GSE######_` | Raw data from a specific GEO study |

---

### Data type suffixes

| Suffix | Description |
|--------|-------------|
| `_metadata` | Sample metadata |
| `_eset` | ExpressionSet (Bioconductor object) |
| `_exprs` | Expression matrix (genes × samples) |
| `_exprs_hgnc_symbol` | Expression matrix with HGNC gene symbols |
| `_log2fc_sample_clean_long` | Per-sample log₂FC in long (tidy) format |
| `_log2fc_paired_unpaired_ttest` | log₂FC with t-test results |
| `_dge_limma_degs` | Differential expression results via limma |
| `_dge_btm_process_genes` | Genes per BTM module with log₂FC |
| `_dge_btm_mean_process` | Mean log₂FC per BTM module |
| `_dge_btm_samples_log2fc_long` | Per-sample × BTM log₂FC, long format |
## Key Data Files

For full variable definitions, table schemas, file naming conventions, and color palettes, see the [CODEBOOK.md](CODEBOOK.md)._mean_process` | Mean log₂FC per MSigDB Hallmark |
| `_gsea_btm_results` | GSEA results using BTMs |
| `_gsea_hallmarks_results` | GSEA results using MSigDB Hallmarks |
| `_gsea_ssGSEA_btms_results` | ssGSEA results using BTMs |
| `_samples_btm_correlation_all_timepoints_summary` | Human–mouse correlation by BTM and timepoint |
| `_fit_prediction_results` | FIT predictor outputs |
| `_sharedgenes_Fisher_Jaccard` | Gene overlap analysis (Fisher + Jaccard) |

---

### File extensions

| Extension | Format |
|-----------|--------|
| `.rds` | R serialized object (`readRDS` / `saveRDS`) |
| `_long.rds` | Same object in long (tidy) format |
| `.csv` | Comma-separated values |
| `.xlsx` | Excel spreadsheet |
| `.gmt` | Gene Matrix Transposed (gene sets) |
| `.fasta` | Nucleotide sequences |

---

## 3. Table Schemas

### 3.1 Sample metadata (`*_metadata`)

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | chr | Unique sample identifier |
| `geo_sample` | chr | GEO accession (e.g., `GSM1234567`) |
| `subject_id` / `participant_id` | chr | Participant or animal identifier |
| `organism` | fct | Species: `"Human"` or `"Mouse"` |
| `timepoint` | num | Days post-treatment (0 = baseline) |
| `timepoints_bin` | chr | Binned time category (`"Day_0"`, `"Day_1"`, `"Day_7"`, etc.) |
| `treatment` / `condition` | chr | Vaccine or experimental condition |
| `pathogen` | chr | Target pathogen of the study |
| `vaccine` / `vaccine_type` | chr | Vaccine name (e.g., `"Fluad"`, `"Engerix B"`) |
| `tissue` | chr | Sample tissue (e.g., `"PBMC"`, `"Blood"`) |
| `sex` / `gender` | chr | Biological sex |
| `age` | num | Participant age (years) |
| `strain` | chr | Mouse strain (e.g., `"C57BL/6J"`) |
| `platform_id` | chr | Microarray platform accession (GPL) |
| `adjuvanted` | chr | Adjuvanted vaccine? (`"Yes"` / `"No"`) |

---

### 3.2 Differential expression results (`*_dge_limma_degs`)

| Column | Type | Description |
|--------|------|-------------|
| `symbol` / `human_symbol` | chr | Human gene symbol (HGNC) |
| `mgi_symbol` | chr | Mouse gene symbol (MGI) |
| `entrez` / `hs_entrez` | int | Human Entrez Gene ID |
| `mm_entrez` | int | Mouse Entrez Gene ID |
| `organism` | chr | Species (`"Human"` / `"Mouse"`) |
| `pathogen` | chr | Pathogen or condition |
| `vaccine` | chr | Vaccine |
| `timepoint` | num | Days post-treatment |
| `condition` | chr | Composite identifier (e.g., `"Human_Influenza_7"`) |
| `mean_l2fc` / `mean_log2fc` | num | Mean log₂ fold change |
| `ci_lower` / `ci_upper` | num | 95% confidence interval bounds |
| `pval` / `p_value` | num | Test p-value |
| `adj_p_val` / `padj` | num | FDR-adjusted p-value (BH method) |
| `t` | num | limma t-statistic |
| `sd` / `se` | num | Standard deviation / Standard error |
| `n` / `n_subjects` | int | Number of samples/subjects |

---

### 3.3 BTM module means (`*_dge_btm_mean_process`)

| Column | Type | Description |
|--------|------|-------------|
| `process` | chr | BTM module name/ID (e.g., `"M9"`) |
| `composite_name` | chr | Full BTM module name |
| `group` | chr | Immune cell group (e.g., `"B CELLS"`) |
| `subgroup` | chr | Functional subgroup |
| `organism` | chr | Species |
| `timepoint` | num | Days post-treatment |
| `condition` | chr | Condition identifier |
| `pathogen` | chr | Pathogen |
| `mean_log2fc` | num | Mean log₂FC of genes in the module |
| `median_log2fc` | num | Median log₂FC |
| `sd_log2fc` | num | Standard deviation of log₂FC |
| `var_log2fc` | num | Variance of log₂FC |
| `n_genes` | int | Number of module genes with available data |

---

### 3.4 GSEA results (`*_gsea_btm_results` / `*_gsea_hallmarks_results`)

| Column | Type | Description |
|--------|------|-------------|
| `ID` / `process` | chr | Gene set ID |
| `Description` | chr | Gene set description |
| `NES` / `nes` | num | Normalized Enrichment Score |
| `pvalue` | num | Nominal p-value |
| `qvalue` | num | Q-value (FDR) |
| `p.adjust` | num | BH-adjusted p-value |
| `setSize` | int | Number of genes in the gene set |
| `leading_edge` | chr | Leading-edge genes (separated by `"/"`) |
| `condition` | chr | Condition analyzed |
| `organism` | chr | Species |
| `gsea_enrichment` | chr | Gene set database used |
| `group` / `subgroup` | chr | Functional grouping (BTMs) |

---

### 3.5 Human–mouse correlation (`*_samples_btm_correlation_all_timepoints_summary`)

| Column | Type | Description |
|--------|------|-------------|
| `process` | chr | BTM module or pathway |
| `timepoint` | num | Days post-treatment |
| `pathogen` | chr | Pathogen or condition |
| `estimate` | num | Spearman or Pearson correlation coefficient |
| `p.value` | num | Correlation p-value |
| `p.adj` | num | FDR-adjusted p-value |
| `group` | chr | Functional group of the BTM module |
| `n_pairs` | int | Number of human–mouse sample pairs compared |

---

### 3.6 Gene overlap (`*_sharedgenes_Fisher_Jaccard`)

| Column | Type | Description |
|--------|------|-------------|
| `Cond1` | chr | Reference condition 1 (e.g., `"Mouse_24"`) |
| `Cond2` | chr | Reference condition 2 (e.g., `"Human_24"`) |
| `Shared` | int | Number of shared genes |
| `NotShared_cond1` | int | Genes exclusive to condition 1 |
| `NotShared_cond2` | int | Genes exclusive to condition 2 |
| `Total_Genes_Cond1` | int | Total genes in condition 1 |
| `Total_Genes_Cond2` | int | Total genes in condition 2 |
| `Percentage_Shared_Cond1` | num | % of shared genes (relative to Cond1) |
| `Percentage_Shared_Cond2` | num | % of shared genes (relative to Cond2) |
| `jaccard_distance` | num | Jaccard distance (0 = identical, 1 = no overlap) |
| `Genes_Names` | chr | Shared gene names (comma-separated) |
| `pvalue` / `p_adj` | num | Fisher's exact test p-value / adjusted p-value |

---

### 3.7 FIT predictor outputs (`*_fit_prediction_results`)

> **FIT** (Found In Translation) is a model that uses mouse log₂FCs to predict human log₂FCs for orthologous genes.

| Column | Type | Description |
|--------|------|-------------|
| `MM.Entrez` / `mm_entrez` | int | Mouse Entrez Gene ID (model input) |
| `HS.Entrez` / `hs_entrez` | int | Human ortholog Entrez Gene ID |
| `fc_mm` / `mean_log2fc` | num | Observed mouse log₂FC |
| `fc_hs` | num | Observed human log₂FC (for validation) |
| `predicted` | num | Human log₂FC predicted by FIT |
| `timepoint` | num | Time in hours or days |
| `condition` | chr | Condition or pathogen |

**FIT input files** (CSVs passed to the predictor):

| Column | Description |
|--------|-------------|
| `MM.Entrez` | Mouse Entrez Gene ID |
| `mean_log2fc` | Mean mouse log₂FC at that timepoint |

---

### 3.8 BTM annotation (`btm_annotation_genes.csv`)

| Column | Type | Description |
|--------|------|-------------|
| `process` | chr | Module name/ID (e.g., `"M9"`) |
| `composite_name` | chr | Full name (e.g., `"M9 (B cell development)"`) |
| `id` | chr | Numeric BTM ID |
| `symbol` | chr | Gene symbol (HGNC) |
| `module_size` | int | Number of genes in the module |
| `module_category` | chr | Category (e.g., `"immune"`) |
| `annotation_level` | chr | Annotation completeness (`"complete"`, `"partial"`) |
| `group` | chr | Immune cell group |
| `subgroup` | chr | Functional subgroup |
| `btm` | chr | Original BTM module reference |
| `transcription_factor` | chr | Associated transcription factors |
| `gene_ontology_terms` | chr | Associated GO terms |
| `jaccard_index` | num | Jaccard similarity with KEGG/Biocarta pathway |
| `SetSize` | int | Gene set size after filtering |

---

### 3.9 Hallmark gene sets (`msigdb_hallmarks_grouped_genes.csv`)

| Column | Type | Description |
|--------|------|-------------|
| `process` | chr | Hallmark name (e.g., `"APOPTOSIS"`, `"INTERFERON_RESPONSE"`) |
| `group` | chr | Functional biological grouping (e.g., `"Immune Response"`) |
| `symbol` | chr | Gene symbol (HGNC) |

---

### 3.10 Human–mouse orthologs (`homologous_mgi_hgnc.csv`)

| Column | Type | Description |
|--------|------|-------------|
| `MGI Marker Accession ID` | chr | MGI ID of the mouse gene |
| `Mouse Gene Symbol` | chr | Mouse gene symbol |
| `Mouse NCBI Gene ID` | int | Mouse Entrez Gene ID |
| `HGNC ID` | chr | Human HGNC ID |
| `Human Gene Symbol` | chr | Human gene symbol |
| `Human NCBI Gene ID` | int | Human Entrez Gene ID |

---

## 4. BTM Immune Groups

Defined in `required.R` (`btm_immune_groups`) and used to color and order BTM modules in figures.

| Group | Description |
|-------|-------------|
| `B CELLS` | B cells — development, activation, and function |
| `PLASMA CELLS` | Plasma cells — antibody secretion |
| `T CELLS` | T cells — CD4, CD8, regulatory |
| `NK CELLS` | Natural Killer cells |
| `MONOCYTES` | Monocytes — classical and non-classical |
| `NEUTROPHILS` | Neutrophils |
| `DC ACTIVATION` | Dendritic cell activation |
| `INNATE RESPONSE` | General innate immune response |
| `INFLAMMATORY/TLR/CHEMOKINES` | Inflammation, TLR signaling, chemokines |
| `INTERFERON/ANTIVIRAL SENSING` | Antiviral sensing and interferon response |
| `IFN` | Type I/II interferon response |
| `PLATELETS` | Platelets |
| `CELL CYCLE` | Cell cycle and proliferation |
| `ENERGY METABOLISM` | Energy metabolism |
| `ECM AND MIGRATION` | Extracellular matrix and cell migration |
| `SIGNAL TRANSDUCTION` | Signal transduction |

---

## 5. Vaccine Types

Used in the `type` column across metadata and study annotation files.

| Code | Vaccine type |
|------|--------------|
| `VLP` | Virus-Like Particles |
| `LA` | Live Attenuated |
| `CONJ` | Conjugate |
| `IN` / `IN/SU` | Inactivated / Subunit |
| `Inactivated` | Inactivated intramuscular |
| `VV` | Viral Vector (3rd generation) |
| `RNA` | mRNA vaccine |
| `SU` | Subunit (protein subunit) |
| `PS` | Polysaccharide |
| `I` | Inactivated (abbreviated) |
| `H` | Hybrid |
| `V-I` | Viral vector + Inactivated |

---

## 6. Color Palettes (defined in `required.R`)

All palettes are defined in the `colors` list object and in dedicated named vectors in `required.R`.

### Organisms (`colors$organism`)

| Category | Hex |
|----------|-----|
| `FIT` | `#4361ee` |
| `Mouse` | `gray50` |
| `Immune` | `#f72585` |
| `Human` | `black` |
| `Permutation` | `gray80` |

### Human–mouse comparison (`colors$comparison`)

| Category | Hex |
|----------|-----|
| `Shared` | `#4DBBD5FF` |
| `Mouse only` | `gray50` |
| `Human only` | `black` |
| `Not-shared` | `gray50` |
| `Immune BTMs` | `#4361ee` |
| `Non-Immune BTMs` | `#3a86ff` |

### Treatment type (`colors$treatment`)

| Category | Hex |
|----------|-----|
| `Vaccination` | `#006494` |
| `Infection` | `#3dccc7` |
| `Injury` | `#deaaff` |

### Temporal scale (`colors$timepoint`)

Gradient: gray (baseline) → light blue (early hours) → dark blue (late days).

| Timepoint | Hex |
|-----------|-----|
| `Day 0` / `0h` | `gray50` |
| `2h` | `#D0D8FB` |
| `6h` | `#4361ee` |
| `12h` | `#3a86ff` |
| `Day 1` | `#caf0f8` |
| `Day 3` | `#90e0ef` |
| `Day 7` | `#0096c7` |
| `Day 12` | `#03045e` |
| `Day 24` | `#184e77` |

### Pathogens (`colors$pathogen`)

| Pathogen | Hex |
|----------|-----|
| `Influenza` | `#4cc9f0` |
| `Hepatitis B` | `#a2d2ff` |
| `Saureus` | `#669bbc` |
| `Ecoli` | `#669bbc` |
| `Trauma` | `#3C5488FF` |
| `Burn` | `#184e77` |

---

*Auto-generated from the project source code — Animals Vax Atlas.*
