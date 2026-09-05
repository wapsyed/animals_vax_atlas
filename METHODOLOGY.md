# Complete Study Methodology

## Animals Vax Atlas: Comparative Immunogenomics of Human-Mouse Vaccine and Infection Responses

---

## I. INTRODUCTION & RESEARCH FRAMEWORK

### I.A Study Rationale & Objectives
- **Core Scientific Question:** To what degree do murine models accurately mirror human blood transcriptomic dynamics during vaccination, acute bacterial infection, and systemic sterile injury?
- **The Translational Conundrum:** Individual orthologous gene-level correlations between mice and humans are frequently poor or inconsistent, casting doubt on preclinical murine translatability.
- **Central Hypothesis:** Higher-order biological structures—specifically functional pathways, Blood Transcription Modules (BTMs), and coordinated gene networks—are evolutionarily conserved across species, retaining high predictive fidelity even when individual gene effect sizes diverge.
- **Key Findings:** 
  1. Shifting from gene-centric to pathway-level metrics substantially elevates cross-species concordance.
  2. Translational fidelity scales monotonically with stimulus intensity: systemic injuries and acute bacterial infections show robust conservation, while milder vaccination stimuli reveal pronounced species-specific divergence.
  3. Expression divergence between orthologs is driven by divergent *cis*-regulatory promoter architecture (ENCODE cCRE rewiring) rather than protein-coding sequence identity or codon evolution.

### I.B Study Design
- **Comparative Design:** Parallel time-course blood transcriptomic profiling of matched biological perturbations across humans (*Homo sapiens*) and laboratory mice (*Mus musculus*).
- **Perturbations Analyzed:** Six matched immune and injury conditions spanning viral-antigen vaccination, Gram-positive bacteremia, Gram-negative sepsis, and major sterile tissue trauma.
- **Data Repositories:** Publicly available transcriptomic cohorts retrieved from the NCBI Gene Expression Omnibus (GEO) and NCBI BioProject databases.
- **Analytical Trajectory:** Automated data curation $\to$ multi-level quality control $\to$ platform-specific normalization $\to$ maximum variance probe collapsing $\to$ linear modeling (limma DGE) $\to$ pathway enrichment (fgsea & ssGSEA) $\to$ macroevolutionary divergence modeling $\to$ predictive classification (ROC/AUC) $\to$ structural/regulatory evolutionary genomics $\to$ multi-modal statistical modeling (`tidymodels`).

---

## II. DATA ACQUISITION & STUDY COHORTS

### II.A Public Repository Search Strategy
1. **BioProject Systematic Curation:** NCBI BioProject was queried for functional genomics studies matching keywords (`"vaccine"`, `"vaccination"`, `"immunization"`, `"immune response"`, `"infection"`), restricted to *Homo sapiens* (Taxonomy ID 9606) and *Mus musculus* (Taxonomy ID 10090).
2. **Inclusion Criteria:**
   - Longitudinal, time-course experimental designs incorporating pre-treatment baseline (Day 0) and post-challenge timepoints.
   - Peripheral blood sampling (whole blood or isolated peripheral blood mononuclear cells [PBMCs]).
   - Minimum sample size of $\ge 3$ biological replicates per timepoint/group.
   - Fully accessible raw or normalized microarray expression matrices and complete sample metadata.
3. **Exclusion Criteria:**
   - Oncology, tumor immunology, chronic autoimmune pathogenesis, or drug-toxicity screens (filtered via regex: `!str_detect(title, "cancer|tumor|autoimmune|drug")`).
   - Single-cell RNA-seq studies without bulk equivalent resolution.
   - Incomplete timepoint or phenotypic annotation.

### II.B Dataset Cohort Matrix

| Stimulus Category | Condition / Agent | Organism | GEO Accession | Platform ID | Array Technology | Tissue Source |
|:---|:---|:---|:---|:---|:---|:---|
| **Vaccination** | Influenza (Fluad: TIV + MF59) | Mouse | GSE120661 | GPL21103 | Agilent-074309 8×60K | Whole Blood |
| **Vaccination** | Influenza (Fluad: TIV + MF59) | Human | GSE124689 | GPL10558 | Illumina HumanHT-12 v4.0 | Whole Blood |
| **Vaccination** | Hepatitis B (Engerix-B) | Mouse | GSE120661 | GPL21103 | Agilent-074309 8×60K | Whole Blood |
| **Vaccination** | Hepatitis B (Engerix-B) | Human | GSE124533 | GPL10558 | Illumina HumanHT-12 v4.0 | Whole Blood |
| **Acute Infection**| *S. aureus* bacteremia | Human | GSE19668 | GPL6244 | Affymetrix Human Gene 1.0 ST | Whole Blood |
| **Acute Infection**| *S. aureus* systemic challenge | Mouse | GSE120661 | GPL21103 | Agilent-074309 8×60K | Whole Blood |
| **Acute Infection**| *E. coli* sepsis | Human | GSE33341 | GPL6244 | Affymetrix Human Gene 1.0 ST | Whole Blood |
| **Acute Infection**| *E. coli* endotoxemia/challenge| Mouse | GSE120661 | GPL21103 | Agilent-074309 8×60K | Whole Blood |
| **Sterile Injury** | Severe Burn Injury | Human | Clinical cohort | Custom Array | Microarray | Whole Blood |
| **Sterile Injury** | Severe Burn Injury | Mouse | GSE182858 | GPL6887 | Illumina MouseWG-6 v2.0 | Whole Blood |
| **Sterile Injury** | Severe Blunt Trauma | Human | GSE36809 | GPL6244 | Affymetrix Human Gene 1.0 ST | Whole Blood |
| **Injury + Vaccine**| Burn + Quadrivalent Vaccine | Mouse | GSE182858 | GPL6887 | Illumina MouseWG-6 v2.0 | Whole Blood |

---

## III. QUALITY CONTROL & PREPROCESSING

### III.A Array Quality Metrics & Outlier Detection
Prior to integration, all expression datasets were evaluated for technical fidelity via `1_QualityControl.Rmd`:
1. **ArrayQualityMetrics (`arrayQualityMetrics`):** Assessed distance matrices, boxplots of signal intensities, and pooled RNA degradation gradients.
2. **Relative Log Expression (RLE):** 
   $$\text{RLE}_{gi} = \log_2(E_{gi}) - \operatorname{median}_{j}(\log_2(E_{gj}))$$
   Samples exhibiting anomalous interquartile range (IQR) divergence or median shifts $>2.5$ standard deviations from the cohort centroid were flagged and excluded.
3. **Principal Component Analysis (PCA):** Outliers exceeding 2 standard deviations across PC1 and PC2 were visually inspected and removed if technical confounding (e.g., severe hybridisation failure) was confirmed.

### III.B Platform-Specific Normalization
Microarray intensity matrices were processed using technology-specific algorithms in `2_Preprocessing_and_DGE.Rmd`:
- **Affymetrix Oligonucleotide Arrays (Human Gene 1.0 ST):** Raw probe cell intensity files (.CEL) were preprocessed with the **Robust Multi-array Average (RMA)** algorithm via `affy`/`oligo`, executing background correction, quantile normalization, and median-polish probe set summarization.
- **Illumina BeadChips & Agilent Arrays:** Raw expression intensities were $\log_2$-transformed and subjected to **between-array quantile normalization** using `limma::normalizeBetweenArrays(method = "quantile")` to enforce identical empirical distributions across arrays while preserving relative biological rank orders.

### III.C Probe-to-Gene Annotation & Maximum Variance Collapsing
To map platform-specific probe identifiers to standardized gene nomenclatures (HGNC symbols for human; MGI symbols for mouse), annotation tables were retrieved via Bioconductor platform packages (`illuminaHumanv4.db`, `hugene10sttranscriptcluster.db`, `biomaRt`).

**Probe Collapsing Strategy:** Microarray platforms frequently feature multiple independent probe sets interrogating the same gene locus. To prevent signal dampening or dilution caused by unweighted averaging across non-responsive or poorly hybridising probes, we employed a **maximum variance probe selection strategy**:
$$\text{Representative Probe}(g) = \arg\max_{p \in \mathcal{P}_g} \left( \operatorname{Var}_{s \in \mathcal{S}} \left( \log_2 I_{p,s} \right) \right)$$
where $\mathcal{P}_g$ is the set of all probes mapping to gene $g$, and $\mathcal{S}$ is the set of all biological samples. For each gene, only the single probe displaying the greatest expression variance across the experiment was retained.

---

## IV. DIFFERENTIAL GENE EXPRESSION & STATISTICAL THRESHOLDS

### IV.A Linear Modeling with Limma
Differential gene expression was modeled independently for each condition and species using linear models with Empirical Bayes variance moderation via `limma`:
1. **Model Formulation:** For each condition, an unintercepted design matrix was constructed:
   $$\mathbf{Y}_{g} = \mathbf{X}\boldsymbol{\beta}_g + \boldsymbol{\varepsilon}_g, \quad \boldsymbol{\varepsilon}_g \sim \mathcal{N}(0, \sigma_g^2 \mathbf{I})$$
   where $\mathbf{X}$ represents experimental groups parameterized by organism, stimulus, and timepoint.
2. **Contrasts:** Post-treatment timepoints were contrasted directly against their matched pre-treatment Day 0 baseline (e.g., $\text{Day 1} - \text{Day 0}$, $\text{Day 3} - \text{Day 0}$, $\text{Day 7} - \text{Day 0}$).
3. **Variance Shrinkage (eBayes):** Squeeze gene-wise sample variances toward a global intensity-dependent prior using `limma::eBayes(fit, trend = TRUE)`:
   $$\tilde{s}_g^2 = \frac{d_0 s_0^2 + d_g s_g^2}{d_0 + d_g}$$
   moderating the test statistics against small sample size instability.

### IV.B Statistical Significance Limits
- **DGE Threshold:** Genes were classified as significantly differentially expressed genes (DEGs) if they satisfied:
  $$\text{adj. } P\text{-value (BH FDR)} \le 0.05$$
  Effect size magnitude ($\log_2\text{FC}$) was preserved continuously without arbitrary threshold truncation to power downstream rank-based pathway analyses.
- **GSEA Significance Threshold:** When evaluating functional modules via `fgsea`, the standard False Discovery Rate limit was established at:
  $$\text{padj} \le 0.25$$
  This $\le 0.25$ FDR threshold is standard for GSEA, capturing broader coordinated pathway trends (especially relevant for weaker stimuli like vaccination).

---

## V. FUNCTIONAL PATHWAY ENRICHMENT (BTMs & HALLMARKS)

### V.A Blood Transcription Modules (BTMs)
The primary analytical unit for immune response evaluation comprised 346 Blood Transcription Modules (Li et al., *Nat Immunol* 2014), capturing specific leukocyte subsets (T cells, B cells, NK cells, monocytes, neutrophils, dendritic cells) and intrinsic functional states (interferon response, inflammatory chemokines, cell cycle). Modules are grouped hierarchically into 16 broader physiological domains.

### V.B Fast Gene Set Enrichment Analysis (fgsea)
For each contrast, orthologous genes were ranked by their moderated $t$-statistic ($t_g$):
$$\text{Rank}(g) = t_g$$
Enrichment scores ($ES$) across BTMs and MSigDB Hallmarks (50 gene sets) were computed using `fgsea::fgsea()`, testing for non-random distribution of module members within the ranked transcriptome via an adaptive multi-level split Monte Carlo permutation scheme (1,000 to 10,000 permutations). Enrichment scores were normalized for gene set size ($NES$).

### V.C Single-Sample Gene Set Enrichment Analysis (ssGSEA)
To quantify module activation at individual sample resolution without requiring group-level contrast specification, single-sample GSEA was conducted using `GSVA::gsva(method = "ssgsea", ssgsea.norm = TRUE)`. For sample $s$ and module $m$, the ssGSEA score reflects the degree to which members of $m$ are coordinately up- or down-regulated relative to all other genes within that individual's transcriptome.

---

## VI. UNIFIED CROSS-SPECIES COMPARATIVE ANALYSES & FIGURE 43 ARCHITECTURE

Notebooks `3.1_Comparing_Human_Mouse_DGE_analyses.Rmd`, `3.2_Comparing_Human_Mouse_GSEA.Rmd`, and `3.3_Comparing_Human_Mouse_Functional_Analyses.Rmd` consolidate all 6 conditions into unified comparative frameworks.

### VI.A Macroevolutionary Expression Divergence & Inverse-Variance Weighting
For each 1:1 orthologous gene pair across matched experimental conditions:
1. **Directional Effect Size Divergence:**
   $$\Delta \log_2\text{FC} = \log_2\text{FC}_{\text{Human}} - \log_2\text{FC}_{\text{Mouse}}$$
2. **Coefficient of Variation (CV):**
   $$CV = \frac{SD}{|\text{mean } \log_2\text{FC}|}$$
3. **Inverse-Variance Statistical Weighting:**
   To ensure that downstream cross-species correlation and distance evaluations were not biased by noisy low-expression probes, joint standard error-based inverse weights were formulated:
   $$W_{SE} = \frac{1}{SE_{\text{Human}}^2 + SE_{\text{Mouse}}^2}$$
   where $SE$ represents the standard error from the linear model fit.

### VI.B Figure 43: Detailed Panel Methodologies

#### Figure 43a | Functional Module NES Correlation Dynamics
- **Objective:** Quantify the temporal concordance of higher-order immune programs across human and murine systems.
- **Computation:** For each condition and timepoint pair, the vector of module Normalized Enrichment Scores ($NES$) was extracted for humans ($\mathbf{NES}_H$) and mice ($\mathbf{NES}_M$).
- **Correlation Metric:** Spearman rank correlation ($\rho$) and Pearson linear correlation ($r$) were computed between species across functional modules over time.
- **Stratified Filtering:** To evaluate the effect of statistical stringency on translatability, correlations were computed and compared across three nested module subsets:
  1. **All Enriched Modules:** Unfiltered module sets.
  2. **Mouse-Significant Modules:** Modules meeting $\text{padj} \le 0.25$ in the mouse model.
  3. **Dual-Significant Modules:** Modules meeting $\text{padj} \le 0.25$ concordantly in both human and mouse cohorts.

#### Figure 43b | Predictive Classification Performance (ROC & AUC)
- **Objective:** Determine whether murine pathway signatures can accurately predict the regulatory status of human biological modules.
- **Classification Task:** Murine effect metrics ($NES$ or mean $\log_2\text{FC}$) were deployed as continuous decision scores to classify human modules as significantly up- or down-regulated.
- **Timepoint Regimes:**
  - **Equal Timepoints:** Matched chronological intervals (e.g., Human Day 1 vs. Mouse Day 1).
  - **Different Timepoints:** Cross-temporal evaluations accommodating species differences in kinetic velocity (e.g., murine Day 1 modeling human Day 3 or Day 7).
- **Validation Controls:** Murine model performance was benchmarked against:
  - **Permutation Controls:** Null distributions generated by randomly shuffling module labels ($N = 1,000$ iterations).
  - **Biological Controls:** Unrelated baseline pathological cohorts (e.g., Duchenne Muscular Dystrophy, DMD) to establish empirical specificity baselines.
- **Performance Metrics:** Receiver Operating Characteristic (ROC) curves and Area Under the Curve (AUC) were generated using `pROC` and `yardstick`.

#### Figure 43c | Leading-Edge Gene (LEG) Conservation
- **Objective:** Dissect whether module-level conservation is driven by identical core genes or alternative pathway components.
- **Extraction:** For each enriched BTM ($\text{padj} \le 0.25$), the Leading-Edge Genes (LEGs)—genes accounting for the core enrichment signal prior to the peak running enrichment score—were extracted for both species.
- **Categorization:** Genes within each module were partitioned into:
  1. **Shared LEGs:** Present in the leading edge of both human and mouse.
  2. **Human-Specific LEGs:** Driving enrichment exclusively in humans.
  3. **Mouse-Specific LEGs:** Driving enrichment exclusively in mice.
  4. **Non-LEGs:** Module members not contributing to the core leading-edge signal.
- **Visualization:** Proportions and counts are visualized via stacked and faceted barplots across conditions.

#### Figure 43d | Gene Rank Conservation in Core Immune Modules
- **Objective:** Evaluate the conservation of intra-module gene prioritization between human and mouse.
- **Focal Module:** Illustrated using the *"immune activation - generic cluster"* (BTM M37.0; 347 genes), representing innate inflammatory, Toll-like receptor, and chemokine signaling programs.
- **Composite Rank Metric:** For each gene $g$, a rank score combining effect size and statistical significance was formulated:
  $$\text{Rank Score}_g = \log_2(\text{FC}_g) \times \left( -\log_{10}(P\text{-value}_g) \right)$$
- **Ranking Vector:** All orthologous module genes were ranked in descending order separately within human ($\text{Rank}_H$) and mouse ($\text{Rank}_M$). Spearman rank correlation was calculated between the two vectors.
- **Visual Mapping:** A parallel-axis rank trajectory plot connects gene positions between human (left axis, H) and mouse (right axis, M):
  - **Blue Lines:** Shared LEGs (conserved core drivers).
  - **Black Lines:** Human-specific LEGs.
  - **Dark Gray Lines:** Mouse-specific LEGs.
  - **Light Gray Lines:** Non-LEGs (background module members).

---

## VII. EVOLUTIONARY GENOMICS: PROTEIN CODING VS. CIS-REGULATORY ARCHITECTURE

Notebooks `5.1_EvolutionaryAnalysis_Protein.Rmd` and `5.2_EvolutionaryAnalysis_Regulation.Rmd` assess whether transcriptomic divergence is governed by structural coding evolution or regulatory rewiring.

### VII.A Coding Sequence (CDS) Acquisition & Longest Isoform Selection
1. **Source:** High-confidence 1:1 orthologous gene pairs were identified via Ensembl BioMart (`hsapiens_gene_ensembl` and `mmusculus_gene_ensembl`) using `biomaRt`.
2. **Partitioned Retrieval:** Due to query volume and server timeout constraints, mouse CDS records were programmatically split into 3 balanced partitions (`target_map_mouse_part1`, `part2`, `part3`) and queried in chunks of 100 genes.
3. **Isoform Resolution (Longest Canonical CDS):** Because alternative splicing produces multiple transcripts per Ensembl Gene ID, a filtering step retained only the longest coding sequence:
   $$\text{Canonical CDS}_g = \arg\max_{t \in \mathcal{T}_g} \left( \operatorname{nchar}(\text{coding}_{g,t}) \right)$$
   Records with `"Sequence unavailable"` or non-coding annotations were purged.

### VII.B Pairwise Nucleotide Global Alignment & Distance Computation
1. **Sequence Sanitization:** Leading/trailing whitespace and newline delimiters were stripped, and nucleotide sequences were converted to uppercase.
2. **Global Needleman-Wunsch Alignment:** Pairwise global alignment was performed using `pwalign::pairwiseAlignment(type = "global")` between human and mouse CDS strings:
   ```r
   dna_align <- pwalign::pairwiseAlignment(
     Biostrings::DNAString(human_cds),
     Biostrings::DNAString(mouse_cds),
     type = "global"
   )
   ```
3. **Conversion to Alignment Object:** Aligned pattern and subject strings—preserving gap characters (`'-'`)—were converted to character matrices and converted into `DNAbin` objects via `ape::as.alignment()` and `ape::as.DNAbin()`.
4. **Kimura 2-Parameter (K80) Distance Calculation:** Molecular evolutionary distances were computed using `ape::dist.dna(bin_dna, model = "K80")`. The Kimura K80 model corrects for multiple hits and differential rates between transitions ($P$) and transversions ($Q$):
   $$d_{\text{K80}} = -\frac{1}{2}\ln(1 - 2P - Q) - \frac{1}{4}\ln(1 - 2Q)$$
   Jukes-Cantor (JC69) distances were computed in parallel as a sensitivity control.
5. **Protein Sequence Identity:** Amino acid sequence identity percentages (`identity_human2mouse` and `identity_mouse2human`) were extracted from Ensembl BioMart homology tables.

### VII.C Cis-Regulatory Element (cCRE) Architecture from ENCODE
To interrogate promoter and enhancer rewiring, Candidate Cis-Regulatory Elements (cCREs) were integrated from the ENCODE project:
- **Human Reference:** GRCh38 cCRE catalogue.
- **Mouse Reference:** mm10 cCRE catalogue.
- **Functional Classes Analyzed:**
  - **PLS:** Promoter-Like Signatures (high DNase, high H3K4me3, centered within 200 bp of TSS).
  - **pELS / dELS:** Proximal and Distal Enhancer-Like Signatures (high DNase, high H3K27ac).
  - **CTCF-bound:** Elements demonstrating high CTCF occupancy (chromatin architectural boundaries).
- **Homology Mapping:** cCREs were assigned to target orthologs based on genomic distance to the primary Transcription Start Site (TSS) ($\pm 2\text{ kb}$ for promoters, $\pm 50\text{ kb}$ for enhancers). Regulatory conservation scores were calculated based on reciprocal element presence and transcription factor binding site conservation.

---

## VIII. MULTI-MODAL STATISTICAL MODELING (`tidymodels`)

Script `6_Statistical_Modelling.Rmd` integrates heterogeneous evolutionary, regulatory, and transcriptomic metrics into a unified predictive modeling framework.

### VIII.A Feature Matrix Construction
A master feature matrix (`human_mouse_statsmodelling_parameters_values.rds`) was assembled by merging:
1. **Structural Coding Evolution:** Kimura K80 CDS distance (`dist_k80`), amino acid identity %.
2. **Transcriptional Concordance:** $\Delta \log_2\text{FC}$, sampling stability $SD$, linear model standard error $SE$.
3. **Cis-Regulatory Features:** Promoter cCRE class presence (PLS, pELS, dELS, CTCF) in human and mouse, cCRE distance metrics.
4. **Transcription Factor (TF) Networks:** Shared vs. species-unique TF binding site counts within orthologous promoters.
5. **Functional Domain Annotations:** BTM module category, immune vs. non-immune designation.

### VIII.B Model Training & Variable Importance
- **Algorithms:** Random Forest classifiers/regressors (via the `ranger` engine) and regularized Elastic Net regression (via `glmnet`) implemented within the `tidymodels` ecosystem.
- **Cross-Validation:** 10-fold cross-validation repeated 5 times, stratified by condition and immune module category.
- **Variable Importance in Projection (VIP):** Permutation-based variable importance metrics were computed via `vip::vip()` to rank the genomic and epigenetic predictors governing cross-species translatability.

---

## IX. COMPUTATIONAL REPRODUCIBILITY & REPOSITORIES

### IX.A Environment Management via `renv`
All analyses were executed under R version 4.5.2 (Bioconductor 3.22). Complete computational environments are pinned in `renv.lock`. The project environment can be restored via:
```r
renv::restore()
```

### IX.B Code Availability & Reproducibility Scripts
- **Primary GitHub Repository:** `https://github.com/wapsyed/animals_vax_atlas`
- **Minimal Worked Example:** `example/example_btm_correlation.R` demonstrates end-to-end BTM correlation computation from cached tables in $< 2$ minutes.

---

## X. REFERENCES

1. **BTMs:** Li S, et al. Molecular signatures of antibody responses derived from a systems biology approach. *Nature Immunology*. 2014;15(2):195-204.
2. **limma:** Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43(7):e47.
3. **fgsea:** Korotkevich G, et al. Fast gene set enrichment analysis. *bioRxiv*. 2021; doi:10.1101/060012.
4. **GSVA:** Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7.
5. **Kimura Distance (K80):** Kimura M. A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. *J Mol Evol*. 1980;16(2):111-120.
6. **ENCODE cCREs:** The ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature*. 2020;583:699-710.
7. **tidymodels:** Kuhn M, Wickham H. Tidymodels: a collection of packages for modeling and machine learning using tidyverse principles. *https://www.tidymodels.org*. 2020.
