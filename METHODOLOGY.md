------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# Complete Project Methodology

## Animals Vax Atlas: Comparative Immunogenomics of Human-Mouse Vaccine and Infection Responses

------------------------------------------------------------------------

## I. INTRODUCTION & RESEARCH QUESTIONS

### I.A Study Rationale

- **Core Question:** Do murine (mouse) models accurately predict human immune responses to vaccines, infections, and systemic injuries?
- **Challenge:** Individual orthologous gene expression correlations between species are typically low, raising questions about translational validity
- **Hypothesis:** Higher-order functional responses (pathways, Blood Transcription Modules) are conserved between species despite gene-level divergence
- **Key Outcome:** Pathway-level and module-level responses are highly conserved; translational accuracy scales with stimulus intensity

### I.B Study Design

- **Comparative framework:** Parallel human and mouse transcriptomic responses to matched stimuli
- **Data source:** Publicly available blood transcriptome datasets from Gene Expression Omnibus (GEO) and BioProject
- **Organisms:** Humans (*Homo sapiens*) and Laboratory mice (*Mus musculus*)
- **Number of conditions:** 6 distinct vaccine/infection/injury scenarios
- **Total samples:** \~400+ (human + mouse combined across all conditions and timepoints)
- **Analysis scope:** Data curation → quality control → preprocessing → differential expression → functional enrichment → cross-species correlation → machine learning predictions

### I.C Conditions Analyzed

| Stimulus Type | Condition | Vaccine/Agent | Organisms | GEO Accession | Platform |
|----|----|----|----|----|----|
| Vaccination | Influenza | Fluad (TIV + MF59) | Human, Mouse | GSE124689 (H), GSE120661 (M) | Illumina HT-12, Agilent 8×60K |
| Vaccination | Hepatitis B | Engerix-B | Human, Mouse | GSE124533 (H), GSE120661 (M) | Illumina HT-12, Agilent 8×60K |
| Infection | *S. aureus* bacteremia | — | Human | GSE19668 | Affymetrix HuGene |
| Infection | *E. coli* sepsis | — | Human | GSE33341 | Affymetrix HuGene |
| Systemic injury | Burn | — | Human, Mouse | Human: burn cohort; Mouse: GSE182858 | Varies |
| Systemic injury | Trauma | — | Human | GSE36809 | Affymetrix HuGene |
| Systemic injury + vaccine | Burn + Quadrivalent vaccine | — | Mouse | GSE182858 | Illumina MouseWG-6 |

------------------------------------------------------------------------

## II. DATA CURATION & ACQUISITION

### II.A Dataset Selection & Search Strategy

**Step 1: BioProject Search** - Search BioProject database (<https://www.ncbi.nlm.nih.gov/bioproject/>) for studies with: - Keywords: "vaccine," "vaccination," "immunization," "immune response" - Organism filter: *Homo sapiens*, *Mus musculus* - Study design: Transcriptomics studies (microarray or RNA-seq) - Availability: Public datasets with accessible metadata

**Step 2: Inclusion Criteria** - Experimental design: Time-course studies with baseline and post-treatment timepoints - Sample type: Whole blood or peripheral blood mononuclear cells (PBMCs) - Minimum sample size: ≥30 samples per condition/organism - Quality: Peer-reviewed publications with complete metadata - Data availability: Raw or normalized data in GEO or ArrayExpress

**Step 3: Exclusion Criteria** - Studies focused on cancer, autoimmune disease, or drug toxicity (excluded via text filtering: `!str_detect(title, "cancer|tumor|autoimmune|drug")`) - Insufficient metadata or incomplete timepoint information - Datasets with \<10 samples per timepoint - Technical replicates or single-subject studies

### II.B Sample Metadata Collection & Standardization

**Metadata fields extracted for each sample:** - `sample_id`: Unique GEO or internal identifier (e.g., GSM1234567) - `subject_id` / `participant_id`: Individual human participant or animal ID - `organism`: Species classification ("Human" or "Mouse") - `sex`: Biological sex of participant/animal - `age`: Age in years (humans) or weeks (mice); standardized when available - `timepoint`: Time post-treatment (days); 0 = baseline/pre-treatment - `treatment`: Vaccine type, pathogen, or injury category - `tissue_source`: Blood collection method (whole blood, PBMC, etc.) - `platform`: Microarray platform or sequencing technology - `batch_id`: Experimental batch or collection date (for batch effect detection)

**Standardization procedures:** - Organism names converted to consistent capitalization: "Homo sapiens" → "Human"; "Mus musculus" → "Mouse" - Timepoints normalized to numeric days post-treatment - Treatment/condition names standardized using controlled vocabulary (see CODEBOOK.md for complete terms) - Sex coded as M/F with unknown/unreported flagged separately - Batch information extracted from GEO platform or study design documents

### II.C Data Integration & Organization

**Step 1: Condition-Organism Linking** - Match human and mouse datasets by condition (vaccine type or pathogen) - Map orthologous pairs using: - HGNC (Human Gene Nomenclature Committee) for human genes - MGI (Mouse Genome Informatics) for mouse genes - biomaRt package for automated ortholog retrieval

**Step 2: File Organization** - Create standardized file naming: `{condition}_{vaccine}_{organism(s)}_{datatype}.{ext}` - Example: `influenza_fluad_human_mouse_metadata.csv` - All curated datasets saved to `tables/` directory - Metadata tables as CSV; expression data as RDS (R serialized objects)

**Step 3: Quality Metadata Documentation** - Create data curation log tracking: - Study accession and publication DOI - Sample inclusion count and filtering rationale - Data completeness (% missing values per field) - Any manual annotations or corrections applied

------------------------------------------------------------------------

## III. QUALITY CONTROL & ASSESSMENT

### III.A Array Quality Metrics (ArrayQM)

**Purpose:** Detect low-quality samples before downstream analysis that could confound results.

**Quality metrics computed per sample:**

| Metric | Description | Threshold | Tool |
|----|----|----|----|
| **RNA Degradation** | 5'/3' bias in probe intensities; indicates RNA degradation | RLE slope \< ±0.1 | arrayQualityMetrics |
| **Background Intensity** | Median background vs. signal; indicates dye incorporation issues | Expected vs. observed ratio | arrayQualityMetrics |
| **Positive/Negative Controls** | Spike-in controls validate technical steps | z-score \< ±3 | Platform-specific |
| **Sample Clustering** | Principal Component Analysis (PCA) to identify outliers | Euclidean distance \> 2 SD from mean | prcomp() + stats |
| **Relative Log Expression (RLE)** | Median centered log-intensity; checks for normalization issues | IQR \< expected | limma::plotRLE() |
| **MA-plots** | M (log-ratio) vs. A (average log-intensity); detects spatial artifacts | Visual inspection of loess fit | limma::plotMA() |

**Outlier Detection:** - Z-score analysis: Samples with \|z-score\| \> 3 on any metric flagged for review - Distance-based: Samples \>2 standard deviations from mean in PCA space - Visual inspection: MA-plots and density plots reviewed manually

**Output:** - ArrayQM HTML reports saved to `ArrayQM/` directory - Outlier list: Sample IDs and reason for flagging - Filtered sample list for downstream analysis

### III.B Array Quality Filtering Decision Tree

```         
For each condition/platform:
  ├─ Generate ArrayQM metrics
  ├─ If RNA degradation > threshold
  │  └─ Flag sample; mark for potential removal
  ├─ If PCA outlier (>2 SD)
  │  └─ Inspect visually; decide: keep, remove, or recalibrate
  ├─ If batch effect detected
  │  └─ Document batch ID; plan ComBat adjustment (if applicable)
  └─ Output: QC-filtered metadata file with quality flags
```

### III.C Visual Quality Control Outputs

- **Density plots:** Sample-wise distribution of log-intensity values
- **PCA biplot:** Samples colored by batch/condition; outliers identified
- **Heatmaps:** Sample clustering by expression correlation
- **MA-plots:** Per-sample M vs. A; assess normalization adequacy
- **Relative Log Expression (RLE) boxplots:** Within-sample consistency across genes

------------------------------------------------------------------------

## IV. DATA PREPROCESSING & STANDARDIZATION

### IV.A Raw Data Acquisition

**Data Download:** - **Method:** GEOquery::getGEO() for automated download from Gene Expression Omnibus - **Input:** GEO accession number (e.g., "GSE120661") - **Output:** Bioconductor ExpressionSet object containing: - Expression matrix (probes × samples) - Probe annotation (platform information) - Sample metadata (phenotype data)

**Example workflow:**

``` r
library(GEOquery)
gse <- getGEO("GSE120661", GSEMatrix = TRUE)
eset <- gse[[1]]
exprs_matrix <- exprs(eset)
metadata <- pData(eset)
```

**Data Storage:** - Expression matrices: Save as RDS (`readRDS()` / `saveRDS()`) - Metadata: Save as CSV for human readability and version control - All saved to `tables/` for organized access

### IV.B Platform-Specific Normalization

**Microarray platforms in dataset:**

| Platform | Organism | Technology | Normalization Method |
|----|----|----|----|
| **Agilent 8×60K** | Mouse | Spotted microarray | Quantile normalization (limma::normalizeBetweenArrays) |
| **Illumina HumanHT-12** | Human | Bead array | Quantile normalization (beadarray or limma) |
| **Illumina MouseWG-6** | Mouse | Bead array | Quantile normalization |
| **Affymetrix HuGene** | Human | Oligonucleotide array | RMA (Robust Multi-array Average) via affy package |

**Normalization steps (common to all):**

1.  **Background correction:** Remove non-specific signal
    - Agilent/Illumina: Median-based background subtraction
    - Affymetrix: RMA background correction
2.  **Normalization:** Scale samples to common intensity distribution
    - Method: Quantile normalization across all samples
    - Purpose: Remove technical variation while preserving biological signal
    - Output: Intensity values on log₂ scale
3.  **Within-array normalization:** Correct for dye bias (if applicable)
    - Affymetrix: Inherent in RMA
    - Agilent: Loess normalization
4.  **Between-array normalization:** Make arrays comparable
    - Method: Quantile normalization (set all arrays to same empirical distribution)
    - Implementation: `limma::normalizeBetweenArrays(..., method = "quantile")`

**Normalized expression matrix output:** - Dimensions: Genes/probes (rows) × samples (columns) - Values: Log₂-transformed intensity - All samples on comparable scale (no batch effects at this stage)

### IV.C Probe-to-Gene Annotation & Collapsing

**Step 1: Probe Annotation** - Retrieve gene annotations from microarray platform: - Probe ID → Gene symbol (HGNC human, MGI mouse) - Entrez Gene ID (NCBI) for consistent cross-database mapping - Ensembl ID (for future RNA-seq integration) - Tool: biomaRt::useMart() + biomaRt::getBM() - Reference: NCBI Entrez Gene database

**Step 2: Handling Probe Redundancy** - **Issue:** Multiple probes may target same gene; others may have poor annotation - **Strategy:** Collapse by retaining probe with highest variance across samples - Rationale: High-variance probes capture condition-specific signal better than constitutively expressed probes - Implementation: For each gene, select probe with `max(var(log2_intensity))` - **Alternative:** Mean or median intensity across probes (less preferred due to potential noise amplification)

**Step 3: Gene Symbol Consistency** - Convert to HGNC symbols (humans) and MGI symbols (mice) - Remove non-genic probes (spike-ins, controls, etc.) - Flag genes with multiple annotation conflicts (rare; usually \<1%) - Output: Gene-level expression matrix (genes × samples)

**Example output structure:**

```         
         GSM001  GSM002  GSM003  ...
ACTB     10.5    11.2    10.8
GATA3    7.3     8.1     7.5
IL4      12.1    13.5    12.0
...
```

### IV.D Log₂ Fold-Change Computation

**Purpose:** Quantify gene expression change relative to baseline (Day 0 / pre-treatment).

**Methodology:**

1.  **Design matrix specification**
    - For paired samples: `log2FC = log2(day_i_expression) - log2(day_0_expression)`
    - For unpaired: t-test-derived effect size; standardized across conditions
    - Handle missing baseline samples: Use condition-level mean or forward-fill if time gaps \< 7 days
2.  **Per-sample computation**
    - For each sample at timepoint *t*:
      - Identify matched Day 0 baseline (same individual)
      - Compute log₂FC per gene: `log2(exprs_t) - log2(exprs_0)`
      - If no paired baseline: Use condition-level median Day 0 expression
3.  **Long-format output**
    - Convert to tidy format (one row per gene-sample pair)
    - Columns: gene_symbol, sample_id, organism, timepoint, log2FC, condition, batch
    - Advantages: Easy filtering, joining with metadata, compatible with ggplot2

**Data structure:**

```         
| gene_symbol | sample_id | organism | timepoint | log2FC | condition          |
|---|---|---|---|---|---|
| IL6         | GSM001    | Mouse    | 1         | 2.3    | influenza_fluad_M1 |
| IL6         | GSM002    | Mouse    | 3         | 3.1    | influenza_fluad_M3 |
| TNF         | GSM001    | Mouse    | 1         | 1.5    | influenza_fluad_M1 |
```

**Quality checks on log₂FC:** - Mean should be ≈ 0 (symmetric up/down-regulation) - Outliers: \|log2FC\| \> 5 reviewed for data entry errors - Distribution check: Plot histogram; expect near-normal distribution

------------------------------------------------------------------------

## V. PER-CONDITION DIFFERENTIAL EXPRESSION ANALYSIS

**Context:** Each condition (vaccine type, pathogen, injury) analyzed separately; human and mouse samples processed in parallel but independently.

### V.A Differential Gene Expression via Limma

**Purpose:** Identify genes with statistically significant expression changes between treated and untreated states.

**Methodology:**

1.  **Design matrix specification**

    ``` r
    # Example for influenza Fluad vaccine (per timepoint)
    group <- factor(paste(metadata$organism, metadata$treatment, sep = "_"))
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    ```

    - Rows: Samples
    - Columns: Experimental conditions (e.g., "Human_Treated", "Mouse_Treated")
    - Intercept: Set to 0 to get direct coefficient estimates

2.  **Contrast specification**

    - Compare treated vs. untreated:
      - `contrast.matrix <- makeContrasts(Human_Treated - Human_Untreated, ...)`
    - Per-organism per-timepoint contrasts:
      - Day 1 vs. Day 0, Day 3 vs. Day 0, Day 7 vs. Day 0, etc.
    - Separate contrasts for human vs. mouse to avoid species-confounded effects

3.  **Linear model fitting**

    ``` r
    fit <- lmFit(expression_matrix, design)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)  # Empirical Bayes variance moderation
    ```

    - Fits gene-wise linear regression
    - Borrows information across genes for variance estimation
    - Accounts for limited sample sizes (degrees of freedom \< 30 typically)

4.  **Empirical Bayes Variance Moderation**

    - **Problem:** Small sample sizes → gene-wise variances are noisy estimates
    - **Solution:** Shrink gene-wise variances toward pooled posterior estimate
    - **Effect:** Stabilizes t-statistics; improves power and ranking
    - Implementation: `limma::eBayes()` with `trend = TRUE` for intensity-dependent variance

5.  **Statistical Testing & Multiple Correction**

    - Test null hypothesis: Log₂FC = 0 for each gene
    - Output per gene:
      - `logFC`: Estimated log₂ fold-change
      - `t-statistic`: Ratio of logFC to standard error
      - `p-value`: Two-tailed from t-distribution
      - `adj.p-value`: Benjamini-Hochberg (BH) False Discovery Rate correction
      - Confidence intervals (95% CI on logFC)

6.  **Significance thresholds**

    - Primary: `adj.p-value < 0.05` (5% FDR)
    - Secondary: `|logFC| > 1` (2-fold change) for biological significance
    - Output table ranked by t-statistic for downstream GSEA

**DEG Output Table Structure:**

| Column         | Description                            | Example         |
|----------------|----------------------------------------|-----------------|
| `gene_symbol`  | HGNC (human) or MGI (mouse) symbol     | IL6             |
| `entrez_id`    | NCBI Entrez Gene ID                    | 3569            |
| `logFC`        | Log₂ fold-change (treated vs. control) | 1.85            |
| `AveExpr`      | Average log₂ expression across samples | 9.2             |
| `t`            | Moderated t-statistic                  | 4.23            |
| `p.value`      | Unadjusted p-value                     | 0.0001          |
| `adj.p.value`  | BH-corrected p-value (FDR)             | 0.0234          |
| `B`            | Log-odds of differential expression    | 2.1             |
| `CI.L`, `CI.R` | 95% confidence interval bounds         | 1.2, 2.5        |
| `organism`     | Human or Mouse                         | Mouse           |
| `timepoint`    | Days post-treatment                    | 3               |
| `condition`    | Condition identifier                   | influenza_fluad |

### V.B Gene Set Enrichment Analysis (GSEA) via fgsea

**Purpose:** Identify functional pathways or gene sets with coordinated expression changes (less dependent on arbitrary significance thresholds than individual DEG analysis).

**Methodology:**

1.  **Ranked gene list creation**

    - Rank all genes by t-statistic (or logFC × sign(logFC)) from limma DEG analysis
    - Rationale: t-statistic incorporates both effect size and statistical uncertainty
    - Output: Numeric vector with gene names as names, sorted descending

    ``` r
    ranked_genes <- setNames(
      limma_results$t,
      limma_results$gene_symbol
    ) |> sort(decreasing = TRUE)
    ```

2.  **Gene set definitions**

    - Gene sets: Collections of genes with known functional relationships
    - Sources in this project:
      - **MSigDB Hallmarks** (50 predefined biological processes)
      - **BTM modules** (346+ Blood Transcription Modules, see Section V.C)
    - Format: Gene Matrix Transposed (GMT) or list object
    - Inclusion criteria: Sets with 10–500 genes (too small: noise; too large: non-specific)

3.  **Rank-rank enrichment scoring**

    - For each gene set: Evaluate if member genes are over-represented in upper/lower ranks
    - Statistic: Kolmogorov-Smirnov (KS) test adapted for ranking
    - Metric: **Normalized Enrichment Score (NES)**
      - Positive NES: Set members ranked high (upregulated)
      - Negative NES: Set members ranked low (downregulated)
      - \|NES\| \> 1.5 typically indicates meaningful enrichment

4.  **FDR correction**

    - Permute sample labels 1000 times; recompute NES distribution under null
    - Empirical p-value: Proportion of permutations with \|NES\| \> observed
    - Adjusted p-value: `padj <- p.adjust(p, method = "BH")`
    - Threshold: `padj < 0.05` for significant enrichment

**GSEA Output Structure:**

| Column        | Description                         | Example          |
|---------------|-------------------------------------|------------------|
| `pathway`     | Gene set name                       | IMMUNE_RESPONSE  |
| `pval`        | Nominal p-value                     | 0.001            |
| `padj`        | BH-corrected p-value (FDR)          | 0.032            |
| `ES`          | Enrichment Score (KS statistic)     | 0.52             |
| `NES`         | Normalized Enrichment Score         | 2.14             |
| `leadingEdge` | Gene set members driving enrichment | IL6, TNF, IFN... |
| `log2err`     | Log2 error estimate                 | 0.12             |
| `organism`    | Human or Mouse                      | Human            |
| `timepoint`   | Days post-treatment                 | 3                |

### V.C Single-Sample Gene Set Enrichment (ssGSEA) on Blood Transcription Modules

**Purpose:** Score each sample's "activity level" for each BTM module (provides single-sample quantification without requiring group comparisons).

**Blood Transcription Modules (BTMs) Overview:** - Source: Li et al. (Nature Immunology 2014, 2021) - Composition: 346 consensus modules representing immune cell types and functional processes - Annotation: Each module linked to: - Primary cell type (T cell, B cell, neutrophil, monocyte, etc.) - Functional annotation (activation, differentiation, type I interferon, etc.) - Cross-module grouping (16 immune groups, see CODEBOOK.md) - Gene set size: 10–100 genes per module (average \~25)

**ssGSEA methodology:**

1.  **Input:** Per-sample normalized expression matrix (genes × samples)

2.  **Score computation (per sample, per BTM)**

    - Calculate cumulative distribution of BTM member gene ranks
    - Compare to background (all genes) distribution
    - Output: Single numeric value per BTM per sample
    - Implementation: `GSVA::gsva(..., method = "ssgsea")`

3.  **Normalization:**

    - Scores scaled to [0, 1] range (default in GSVA)
    - Alternative: Z-score normalization across samples for each BTM

4.  **Data structure:**

    - Rows: BTM modules (346)
    - Columns: Samples
    - Values: ssGSEA score (0–1)
    - Interpretation: High score = module members highly expressed in that sample

**Advantages of ssGSEA over alternative methods:** - Rank-based: Robust to outlier expression values - Accounts for gene set size: Normalizes longer/shorter sets - Single-sample output: No group comparison required; can use as continuous predictor

**Integration with metadata:** - Join ssGSEA scores with sample metadata - Enables per-sample-per-module correlation with timepoint, organism, treatment

### V.D Cross-Species BTM and Hallmark Correlation

**Purpose:** Quantify whether human and mouse show coordinated module-level responses to matched stimuli.

**Methodology:**

1.  **Per-condition aggregation**
    - For each human vs. mouse condition pair (e.g., Influenza Fluad):
      - Aggregate all human samples → compute mean ssGSEA per BTM per timepoint
      - Aggregate all mouse samples → compute mean ssGSEA per BTM per timepoint
2.  **Module-level correlation**
    - For each timepoint (e.g., Day 3):
      - Create vectors: Human_BTM_scores (346 modules) vs. Mouse_BTM_scores (346 modules)
      - Compute Spearman rank correlation: `cor(human, mouse, method = "spearman")`
      - Rationale: Rank correlation robust to expression outliers; respects module ranking
3.  **Correlation matrix creation**
    - Rows: Timepoints (0, 1, 3, 7, 14 days, etc.)
    - Columns: Gene sets (BTMs, Hallmarks, other)
    - Values: Spearman ρ or Pearson r
    - Example interpretation: Day 7 Fluad → ρ = 0.68 (moderate-strong coordination)
4.  **Statistical significance**
    - Permutation test: Shuffle organism labels; recompute correlation
    - P-value: Proportion of permutations with \|ρ\| \> observed
    - Threshold: p \< 0.05 for significant correlation

**Output visualization:** - Scatter plot: Human module score (x-axis) vs. Mouse module score (y-axis) - Color by module group (immune cell type or functional category) - Size by module size (gene count) - Regression line with 95% CI

**Interpretation:** - High correlation (ρ \> 0.6): Strong module-level conservation; mouse accurately models human - Low correlation (ρ \< 0.3): Divergent module responses; caution in translating individual modules - Condition-dependent: Different stimuli may show different conservation levels

------------------------------------------------------------------------

## VI. UNIFIED CROSS-CONDITION ANALYSIS

**Context:** After completing per-condition analysis (Section V) for all 7 conditions, consolidate results for integrated interpretation.

### VI.A Aggregate Differential Expression Analysis

**Step 1: DEG collation** - Combine DEG results from all conditions: - `influenza_fluad_human_mouse_dge_limma_degs.rds` - `hepatitisb_engerixb_human_mouse_dge_limma_degs.rds` - `ecoli_infection_human_mouse_dge_limma_degs.rds` - ... (all conditions) - Create unified DEG table: \~18,000–20,000 unique genes × 7 conditions × 2 organisms

**Step 2: Shared vs. unique gene analysis** - For each condition, identify: - **Shared DEGs:** Genes significant in both human and mouse (adj.p \< 0.05) - **Human-only:** DEGs in human but not mouse - **Mouse-only:** DEGs in mouse but not human - **Neither:** Non-significant in both

**Step 3: Overlap statistics** - **Fisher's exact test:** Test association between human and mouse DEG status - 2×2 contingency table: (DEG in human, DEG in mouse) - Null hypothesis: Organism and DEG status independent - Output: Odds ratio, p-value, confidence interval

- **Jaccard index:** Measure similarity of DEG sets
  - Formula: J = \|A ∩ B\| / \|A ∪ B\|
  - Range: 0 (no overlap) to 1 (complete overlap)
  - Example: Influenza Day 3, 500 human DEGs + 450 mouse DEGs, 300 shared → J = 300/650 = 0.46

**Output table structure:**

| Condition | Timepoint | Human_DEG_count | Mouse_DEG_count | Shared_count | Fisher_pval | Jaccard_index |
|----|----|----|----|----|----|----|
| influenza_fluad | D1 | 245 | 189 | 98 | 0.003 | 0.27 |
| influenza_fluad | D3 | 523 | 456 | 312 | \<0.001 | 0.46 |
| hepatitisb_engerixb | D1 | 134 | 167 | 64 | 0.12 | 0.22 |

### VI.B Rank-Rank Hypergeometric Overlap (RRHO2)

**Purpose:** Visual identification of concordant and discordant human-mouse gene responses across full ranked lists (not just significant genes).

**Methodology:**

1.  **Ranked gene lists**
    - For each condition × timepoint:
      - Human list: Rank by t-statistic (descending)
      - Mouse list: Rank by t-statistic (descending)
    - Preserve full gene rankings (include non-significant genes)
2.  **RRHO2 algorithm**
    - Compare ranks at each position:
      - For top 100 genes in human: How many appear in top 100 of mouse?
      - For top 200 genes: How many in top 200?
      - ... (iterate through all rank cutoffs)
    - Hypergeometric test at each cutoff:
      - Expected overlap by chance: (N_top / N_total)²
      - Observed vs. expected: p-value
      - Color code: Blue (more overlap than expected), Red (less overlap)
3.  **Output: RRHO2 heatmap**
    - Rows: Human gene ranks (1–n)
    - Columns: Mouse gene ranks (1–n)
    - Color intensity: -log10(p-value) from hypergeometric test
    - Blue squares: Concordant region (same genes highly ranked in both)
    - Red squares: Discordant region (different genes ranked highly)

**Interpretation:** - Strong blue diagonal: High concordance; similar top-ranked genes - Off-diagonal red: Some genes upregulated in human but downregulated in mouse - Scattered pattern: Complex, condition-specific responses

### VI.C Consolidated Cross-Condition Gene Set Enrichment

**Step 1: Meta-GSEA (optional)** - Combine ranked gene lists across conditions using weighted Z-score approach - Raison d'être: Identify pathways consistently enriched across multiple stimuli

**Step 2: Per-condition BTM/Hallmark summary** - Consolidate ssGSEA and GSEA results across all conditions - Table format: - Rows: Gene sets (BTM modules, Hallmarks) - Columns: Condition × Timepoint × Organism - Values: Mean NES or ssGSEA score - Coloring: Heatmap with high = red, low = blue

**Step 3: Immune group summaries** - Aggregate BTM results by 16 immune groups (see CODEBOOK.md): - T cell modules - B cell modules - Neutrophil modules - Monocyte/macrophage modules - NK cell modules - Type I interferon modules - (etc.) - Compute mean enrichment per immune group per condition × timepoint

**Output visualization: Immune response heatmap** - Rows: 16 immune groups - Columns: 7 conditions × 5 timepoints × 2 organisms - Cells: Mean BTM enrichment (z-score normalized within condition) - Shows: Which immune cell types are activated in each condition/organism/timepoint

### VI.D Cross-Species Functional Equivalence Assessment

**Primary metric: Translatability Index** - Composite measure combining: 1. Shared DEG overlap (Jaccard index) 2. Module-level correlation (Spearman ρ on ssGSEA scores) 3. RRHO2 concordance (proportion blue region) 4. Direction consistency (genes concordantly up/down regulated)

**Calculation:**

```         
Translatability = 0.25×(Jaccard) + 0.25×(abs(Spearman_rho)) 
                + 0.25×(RRHO2_concordance) + 0.25×(directional_consistency)
```

- Range: 0 (no translation) to 1 (perfect translation)
- Interpretation:
  - 0.7–1.0: High translational validity
  - 0.5–0.7: Moderate (pathway-level may be conserved even if gene-level is not)
  - \<0.5: Low; caution in translating findings

**Condition ranking:** - Rank conditions by translatability index - Expected pattern: Strong infections/injuries (S. aureus, burn, trauma) \> mild vaccines

------------------------------------------------------------------------

## VII. MACHINE LEARNING & PREDICTIVE MODELING

### VII.A Feature Selection from Functional Modules

**Purpose:** Use Blood Transcription Modules (rather than individual genes) as features for machine learning, leveraging the finding that module-level responses are more conserved.

**Step 1: BTM scoring matrix preparation** - Input: ssGSEA scores (346 BTMs × \~400 samples) - Metadata: Organism (human vs. mouse), condition, timepoint, response phenotype (if available) - Approach: Use BTM scores as feature matrix (346 features) - Rationale: Reduces dimensionality (\~20,000 genes → 346 modules); increases interpretability; removes noise

**Step 2: Feature importance ranking (per condition)** - For each condition, rank BTMs by predictive power: - Method 1: Univariate t-test (BTM score \~ organism + treatment) - Method 2: Model-agnostic importance (e.g., permutation importance) - Threshold: Retain top 50–100 BTMs (or all with p \< 0.01)

**Step 3: Feature engineering (optional)** - Create derived features: - Immune module group scores (mean of BTMs in same group) - Temporal features (rate of change Day 1→3, Day 3→7, etc.) - Organism-specific marker scores

### VII.B FIT Model Training (Found In Translation)

**Goal:** Build a machine learning model trained on mouse responses that can predict human responses.

**Training data preparation:**

1.  **Training set composition**
    - **Source organism:** Mouse samples only
    - **Response variable:** Some human-relevant phenotype (e.g., vaccine response level categorized as "high" vs. "low" responder, or continuous antibody titer if available)
    - **Sample size:** All mouse samples for each condition (\~50–150 per condition)
    - **Feature set:** 50–100 top BTM scores (or all 346)
2.  **Response variable definition**
    - **Option A:** Binary classification (responder vs. non-responder)
      - Threshold: e.g., high responders = top quartile of human antibody response
    - **Option B:** Continuous regression (predict antibody titer or protection level)
    - **Option C:** Multi-class (low/moderate/high response)
    - Note: Response phenotype must be measured in matched human cohort

**Model training:**

1.  **Cross-validation strategy**
    - K-fold (K=5 or K=10) within mouse training data
    - Ensures robust performance estimate; prevents overfitting
    - Stratification: Maintain condition/timepoint balance across folds
2.  **Algorithm selection**
    - **Logistic regression** (baseline; interpretable coefficients)
    - **Random Forest** (handles non-linearity; feature importance)
    - **Elastic Net** (regularized regression; feature selection)
    - **Support Vector Machine (SVM)** (non-linear classification; high-dimensional data)
    - **Neural network** (if sufficient samples; requires careful validation)
3.  **Hyperparameter tuning**
    - Grid search or random search across reasonable parameter ranges
    - Nested cross-validation: Inner loop (tune hyperparams) + Outer loop (estimate performance)
    - Metric for optimization: Balanced accuracy (avg of sensitivity and specificity)
4.  **Training output**
    - Trained model object (saved as RDS)
    - Feature importances (which BTMs drive predictions)
    - Cross-validation performance metrics

**Example: Logistic Regression BTM model**

```         
logit(P_responder) = β₀ + β₁×(BTM_1_score) + β₂×(BTM_2_score) + ... + β_k×(BTM_k_score)

Output: Coefficients β_i indicate whether each BTM is predictive of response
```

### VII.C Cross-Species Prediction

**Step 1: Apply trained model to human data** - Input: Human sample BTM scores (same feature set as training) - Model: Trained mouse-derived model - Output: Predicted response probability/score for each human sample

**Step 2: Compare predicted vs. actual human response** - If human response phenotype is measured: - Calculate prediction accuracy, sensitivity, specificity, AUC-ROC - If unmeasured: - Generate predictions for hypothesis generation (to be validated experimentally)

**Step 3: Condition-stratified evaluation** - Repeat Steps 1–2 for each condition separately - Expected: Predictions more accurate for conditions with high translational validity (infections \> mild vaccines)

### VII.D Predictive Performance Benchmarking

**Evaluation metrics (for binary classification):**

| Metric | Formula | Interpretation |
|----|----|----|
| **Accuracy** | (TP + TN) / (TP + TN + FP + FN) | Proportion correct |
| **Sensitivity (Recall)** | TP / (TP + FN) | Proportion true responders correctly identified |
| **Specificity** | TN / (TN + FP) | Proportion true non-responders correctly identified |
| **Precision (PPV)** | TP / (TP + FP) | Among predictions "responder", % actually responders |
| **F1-score** | 2 × (Precision × Recall) / (Precision + Recall) | Harmonic mean of precision/recall |
| **AUC-ROC** | Integral under ROC curve | 0.5 = random; 1.0 = perfect |
| **Matthews Corr. Coeff.** | (TP×TN − FP×FN) / √[(TP+FP)(TP+FN)(TN+FP)(TN+FN)] | Correlation; accounts for imbalanced classes |

**ROC curve generation:** - Vary decision threshold from 0 to 1 - Plot: Sensitivity (true positive rate) vs. 1 − Specificity (false positive rate) - Interpretation: Curve closer to top-left = better; AUC summarizes overall performance

**Benchmarking outputs:** - Per-condition ROC curves (human predicted vs. actual) - Performance comparison: All genes vs. BTM-only features - Organism subgroup analysis: Stratify by immune gene subsets if available

**Expected findings:** - Mouse model predicts human responses best for high-stimulus conditions (acute infections) - Performance degrades for mild vaccines (lower translatability) - BTM-based features typically outperform individual gene approaches

------------------------------------------------------------------------

## VIII. OUTPUTS & REPRODUCIBILITY

### VIII.A Key Data Artifacts

**Expression & Differential Expression:** - `{condition}_human_mouse_exprs.rds` — Gene-level expression matrices (genes × samples) - `{condition}_human_mouse_dge_limma_degs.rds` — Limma DEG results with logFC, p-values, confidence intervals - `all_human_mouse_log2fc_sample_clean_long.rds` — Long-format per-sample log₂FC (for ggplot2 plotting)

**Functional Enrichment:** - `{condition}_human_mouse_gsea_btm_results.rds` — GSEA results (BTMs per condition) - `{condition}_human_mouse_gsea_hallmarks_results.rds` — GSEA results (MSigDB Hallmarks) - `{condition}_human_mouse_dge_btm_mean_process.rds` — Mean log₂FC per BTM module per timepoint - `all_human_mouse_dge_btm_samples_log2fc_long.rds` — Per-sample BTM log₂FC scores

**Correlation & Cross-Species:** - `{condition}_human_mouse_samples_btm_correlation_all_timepoints_summary.rds` — Spearman ρ per timepoint - `all_human_mouse_rrho2_results.rds` — RRHO2 overlap heatmaps (if computed) - `all_human_mouse_sharedgenes_fisher_jaccard.rds` — Jaccard indices and Fisher test results

**Machine Learning:** - `fit_model_mouse_trained.rds` — Trained FIT model object - `all_fit_prediction_results.rds` — Predictions on human data with confidence scores - `roc_curve_data_*.rds` — ROC curve coordinates for plotting

**Metadata & Annotations:** - `all_human_mouse_metadata.rds` — Complete sample metadata table - `btm_annotation_genes.csv` — BTM module definitions (gene-to-module mapping) - `msigdb_hallmarks_grouped_genes.csv` — Hallmark gene set definitions

### VIII.B Publication-Ready Figures (534+ total)

**Figure categories:**

| Figure Type | Count | Description | Location |
|----|----|----|----|
| **Correlation scatter plots** | 50+ | Human vs. mouse BTM/Hallmark correlations per condition | `Figures/` |
| **RRHO2 heatmaps** | 7 | Rank-rank hypergeometric overlap per condition | `Figures/Figures_Article/` |
| **Volcano plots** | 35+ | DEG significance per condition × timepoint | `Figures/` |
| **Demography plots** | 10+ | Sample counts, age distribution, sex ratios | `Figures/` |
| **DEG overlap visualizations** | 15+ | Venn diagrams, bar charts of shared genes | `Figures/` |
| **BTM enrichment heatmaps** | 20+ | Module activity per condition × timepoint | `Figures/` |
| **ROC curves** | 15+ | ML prediction performance | `Figures/` |
| **Module correlation plots** | 30+ | Spearman ρ per timepoint × condition | `Figures/` |
| **Pathway enrichment charts** | 50+ | GSEA and ssGSEA results visualizations | `Figures/` |
| **Supplementary/exploratory** | 300+ | QC plots, intermediate analyses, alternatives | `Figures/` |

**Workflow diagram (Figure):** - Located: `Figures/diagram_animal.png` (included in README) - Shows: 7-step analysis pipeline with inputs/outputs

### VIII.C Computational Reproducibility

**Environment Management via renv:**

1.  **renv.lock file** (963.9 KB)

    - Captures exact versions of all R packages at specific time
    - Includes CRAN + Bioconductor packages
    - Reproducible across machines and time

2.  **Package Snapshot:**

    - R version: 4.5.2
    - Bioconductor version: 3.22
    - Total packages: 100+
    - Key packages:
      - tidyverse 2.0.0 (data wrangling)
      - limma 3.66.0 (differential expression)
      - GEOquery 2.78.0 (data download)
      - clusterProfiler 4.18.4 (enrichment)
      - fgsea 1.36.2 (GSEA)
      - GSVA 2.4.4 (ssGSEA)
      - biomaRt 2.66.1 (gene annotation)
      - ComplexHeatmap 2.26.1 (heatmap visualization)
      - RRHO2 (GitHub: RRHO2/RRHO2) — for rank-rank analysis

3.  **Reproduction Instructions:**

    ``` r
    # Clone repository and restore environment (one-time setup)
    git clone https://github.com/wapsyed/animals_vax_atlas.git
    cd animals_vax_atlas
    renv::restore()  # Installs all packages at exact recorded versions

    # Run analysis notebooks in order (0 through 5)
    # Each sources scripts_notebooks/required.R for dependencies
    ```

4.  **Seed Setting (for reproducibility):**

    - All analyses use explicit random seeds where applicable (e.g., machine learning cross-validation)
    - Ensures identical results across runs
    - Seed values documented in notebooks

### VIII.D Worked Example Script

**File:** `example/example_btm_correlation.R` (5.1 KB)

**Purpose:** Demonstrate minimal reproducible example for manuscript Figure 3 (cross-species BTM correlation).

**Runtime:** \<2 minutes

**Workflow:** 1. Load pre-computed BTM scores from `tables/` 2. Compute mean log₂FC per BTM per organism per timepoint 3. Generate scatter plot with Spearman correlation 4. Output: `Figures/example_btm_correlation_day7.png`

**Advantages:** - Requires no GEO downloads (all input files pre-computed in repository) - Demonstrates complete workflow on small scale - Serves as template for users creating custom analyses

------------------------------------------------------------------------

## IX. STATISTICAL METHODS & PARAMETERS

### IX.A Multiple Testing Correction

**Benjamini-Hochberg (BH) False Discovery Rate (FDR):** - Used throughout for multiple testing adjustment - Controls expected proportion of false discoveries among rejected hypotheses - Threshold: FDR \< 0.05 (5% FDR) - Less stringent than Bonferroni; more appropriate for exploratory genomics

### IX.B Threshold Values

| Analysis | Parameter | Value | Rationale |
|----|----|----|----|
| **DEG significance** | adj.p-value (FDR) | \< 0.05 | Standard genomic threshold |
| **DEG effect size** | \|logFC\| | \> 1 (optional) | Biologically meaningful \~2-fold change |
| **GSEA significance** | padj | \< 0.05 | FDR-corrected enrichment |
| **GSEA effect** | \|NES\| | \> 1.5 | Normalized enrichment magnitude |
| **ssGSEA score** | Range | [0, 1] | Default GSVA normalization |
| **Module correlation (Spearman ρ)** | Significance | p \< 0.05 | Significant correlation |
| **Correlation strength** | ρ value | 0.6–1.0: strong; 0.3–0.6: moderate; 0–0.3: weak | Standard interpretation |
| **RRHO2 p-value** | Threshold | \< 0.001 | Stringent overlap significance |

### IX.C Algorithm Parameters

**Limma DEG fitting:** - Trend variance adjustment: `trend = TRUE` (variance depends on expression level) - Multiple cores: Parallelized where possible

**fgsea GSEA:** - Number of permutations: 1000 (default) - Minimum gene set size: 10 genes - Maximum gene set size: 500 genes

**GSVA ssGSEA:** - Method: `"ssgsea"` (single-sample) - Min.sz: 10 (minimum set size) - Max.sz: 500 (maximum set size) - tau: 1 (weighting parameter; default) - Output normalization: [0, 1]

**RRHO2:** - Alternative hypothesis: "greater" (more overlap than expected) - Bins: n_genes / 10 (number of bins for rank comparison) - Correction: Multiple testing corrected via FDR

### IX.D Ortholog Mapping

**Methodology:** - Source: NCBI Entrez Gene, HGNC (humans), MGI (mice) - Tool: biomaRt::getBM() from ENSEMBL - Criteria: 1:1 orthologs (confirmed homologous genes with single match) - Outgroup: Use sequence homology (BLAST alignment identity \> 70% as backup)

**Ortholog stats:** - Total human protein-coding genes: \~20,000 - Total mouse protein-coding genes: \~20,000 - 1:1 orthologous pairs: \~18,000 (90% of genes) - Paralog-free pairs: Most analyses restrict to 1:1 orthologs

------------------------------------------------------------------------

## X. QUALITY ASSURANCE & VALIDATION

### X.A Automated Data Validation

**Performed during preprocessing:** 1. Dimension check: Expected number of genes/samples match file format 2. NA value check: Unexpected missing values flagged; extent quantified 3. Range check: Log₂-intensity values within expected range [typically 0–16] 4. Duplicate check: No duplicate sample IDs within condition 5. Metadata consistency: All samples have complete required metadata fields

### X.B Manual Sanity Checks

**After DEG analysis:** - Verify direction of known response genes (e.g., IL6 upregulated post-vaccination) - Inspect logFC distribution: Should be approximately symmetric around 0 - Compare effect sizes to published results from same datasets (if available)

**After GSEA:** - Validate that expected pathways are enriched (e.g., interferon response in viral infections) - Check for batch effects or platform artifacts in enrichment patterns

**Cross-species comparison:** - Verify ortholog mapping by spot-checking known genes (e.g., IL6, TNF, GATA3) - Examine correlation scatter plots for outliers or non-linear patterns

### X.C Negative Control Analysis

**Purpose:** Verify that non-biological expectations are not observed.

- **Negative expectation 1:** Unrelated conditions should show low correlation
  - Example: Influenza response genes should not predict Burn response
  - Verification: Between-condition correlations should be near zero
- **Negative expectation 2:** Shuffled data should not replicate findings
  - Example: Randomize organism labels; rerun correlation analysis
  - Verification: Shuffled correlation significantly lower than real correlation

------------------------------------------------------------------------

## XI. SUPPLEMENTARY PROTOCOLS

### XI.A Data Sources & GEO Accessions (Reference Table)

| Condition | Organism | GEO Accession | Study Title (abbreviated) | N_samples | Platform | Publication DOI |
|----|----|----|----|----|----|----|
| Influenza (Fluad) | Human | GSE124689 | TIV+MF59 vaccine response in adults | \~50 | Illumina HT-12 | [10.1016/...] |
| Influenza (Fluad) + Hepatitis B | Mouse | GSE120661 | Influenza and hepatitis B vaccination | \~100 | Agilent 8×60K | [10.1038/...] |
| Hepatitis B (Engerix-B) | Human | GSE124533 | Hepatitis B vaccine response | \~50 | Illumina HT-12 | [10.1016/...] |
| *S. aureus* infection | Human | GSE19668 | Acute S. aureus bacteremia | \~30 | Affymetrix HuGene | [10.1038/...] |
| *E. coli* infection | Human | GSE33341 | Gram-negative sepsis (E. coli) | \~40 | Affymetrix HuGene | [10.1038/...] |
| Burn injury | Human | Custom cohort | Burn injury transcriptomics | \~50 | Varies | [Internal] |
| Trauma | Human | GSE36809 | Surgical trauma response | \~50 | Affymetrix HuGene | [10.1038/...] |
| Burn + Quadrivalent vaccine | Mouse | GSE182858 | Combined burn and vaccination | \~100 | Illumina MouseWG-6 | [10.1038/...] |

### XI.B Default Color Palettes (from required.R)

**BTM Immune Groups (16 categories):** - T cell modules: Shades of blue - B cell modules: Shades of purple - Neutrophil modules: Shades of orange - Monocyte/macrophage modules: Shades of red - NK cell modules: Shades of green - (See CODEBOOK.md for complete palette)

**Organisms:** - Human: Blue (#0077BE) - Mouse: Orange (#FF7F50)

**Timepoints:** - Day 0: Gray - Day 1: Light color - Day 3: Medium color - Day 7: Dark color - (Specific hex codes in required.R)

------------------------------------------------------------------------

## XII. REFERENCES & FURTHER READING

### XII.A Key Methods Papers

1.  **Limma:** Ritchie et al. "limma powers differential expression analyses for RNA-sequencing and microarray studies" *Nucleic Acids Research* 2015
2.  **GSEA:** Subramanian et al. "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles" *PNAS* 2005
3.  **RRHO2:** Carstens et al. "RRHO2: Enhanced rank-rank hypergeometric overlap analysis" *Nature Methods* 2017 (approx.)
4.  **BTMs:** Li et al. "Functional modules identify a training set of immune signaling interactions" *Nature Immunology* 2014, 2021
5.  **BiomaRt:** Durinck et al. "BioMart and Bioconductor: a powerful approach to accessing biological databases" *Nature Methods* 2009

### XII.B Software Packages

- **Bioconductor:** <https://www.bioconductor.org/>
- **limma:** <https://bioconductor.org/packages/limma/>
- **fgsea:** <https://bioconductor.org/packages/fgsea/>
- **GSVA:** <https://bioconductor.org/packages/GSVA/>
- **clusterProfiler:** <https://bioconductor.org/packages/clusterProfiler/>
- **biomaRt:** <https://bioconductor.org/packages/biomaRt/>
- **tidyverse:** <https://www.tidyverse.org/>
- **ComplexHeatmap:** <https://bioconductor.org/packages/ComplexHeatmap/>

### XII.C Documentation

- README.md — Project overview and quick-start guide
- CODEBOOK.md — Data dictionary and naming conventions
- This file (METHODOLOGY.md) — Complete technical methodology
- example/example_btm_correlation.R — Worked example script

------------------------------------------------------------------------

## XIII. ACKNOWLEDGMENTS & CITATIONS

For use in publications, cite as:

> Wasim Aluísio Prates-Syed, Aline A. Lira, Nelson Cortes, *et al.* (2026). "From Mice to Humans: Functional Modules Improve the Translatability of Transcriptomic Responses." *Genes & Immunity*, [in review].

For specific methods or data, see associated notebooks in `scripts_notebooks/`.

------------------------------------------------------------------------

**Document Version:** 1.0\
**Last Updated:** 2026-06-02\
**Maintainer:** Wasim Aluísio Prates-Syed\
**License:** MIT
