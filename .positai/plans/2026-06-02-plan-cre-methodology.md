# Plan: Comprehensive Project Methodology & Visual Flowchart

**Date:** 2026-06-02  
**Goal:** Create a complete sequential methodology document describing ALL analyses in the Animals Vax Atlas project (data curation through validation) + generate a visual flowchart showing how all scripts connect.  
**Status:** Planning

---

## Understanding the Request

The user wants:

1. **Consolidated Methodology Document (Bullet Points)**
   - Describe EVERYTHING done in the project sequentially
   - From data importation & standardization → differential expression → enrichment → ML → validation
   - Format as hierarchical bullet points (inspired by the CRE example provided)
   - Include all 6 notebooks: 0→1→2→3.1→3.2→4→5
   - Make it publication/dissertation ready

2. **Visual Flowchart Diagram**
   - Similar style to `Posdoc_Projeto_Metodologia.png` (the example image)
   - Main stages on left (boxes)
   - Detailed steps/tools on right (bullet lists)
   - Show data flow and connections between scripts/analyses
   - Easy to understand the complete workflow

---

## Project Overview from Exploration

### Project Scope
- **Title:** Animals Vax Atlas / Mouse2Human
- **Goal:** Assess translational value of murine models for human immune responses
- **Conditions:** 7 vaccine/infection scenarios (Influenza, Hepatitis B, S. aureus, E. coli, Trauma, Burn, Quadrivalent)
- **Organisms:** Human + Mouse (comparative immunogenomics)
- **Total Samples:** ~400+ across GEO datasets
- **Key Finding:** Pathway/BTM-level responses conserved between species (even if individual genes aren't)

### Sequential Analysis Pipeline
1. **0_Data_Curation.Rmd** → Collect & curate datasets from BioProject
2. **1_QualityControl.Rmd** → ArrayQM, quality metrics, outlier detection
3. **2_Preprocessing.Rmd** → Download GEO, normalize, probe-to-gene mapping, compute log2FC
4. **3.1_Comparing_Human_Mouse_ByCondition.Rmd** → Per-condition DEG, GSEA, BTM correlation
5. **3.2_Comparing_Human_Mouse_UnifiedAnalyses.Rmd** → Cross-condition unified analysis, RRHO2
6. **4_Performance.Rmd** → ROC curves, predictive performance
7. **5_MachineLearning.Rmd** → Build ML classifiers (FIT model)

### Data Organization
- **tables/**: 426 processed files (expression matrices, DEG results, GSEA outputs)
- **VaxGO/**: Gene set annotations (BTMs, MSigDB Hallmarks, VaxSigDB)
- **Genomic/**: cCRE BED files, ortholog registries, gene-CRE associations
- **Figures/**: 534+ publication-ready visualizations
- **DataCuration/**: Raw metadata & manual annotations

---

## Plan: Create Comprehensive Methodology Document

### Document Structure: METHODOLOGY.md

```
I. INTRODUCTION
   • Research questions
   • Study design overview
   • Comparative framework (human vs mouse)

II. DATA CURATION & ACQUISITION
   A. Dataset Selection
      • BioProject search strategy
      • Inclusion/exclusion criteria
      • Data sources (7 GEO accessions)
      • Sample metadata collection
   
   B. Sample Filtering & Annotation
      • Organism filtering (human, mouse)
      • Vaccine/treatment classification
      • Timepoint standardization
      • Sex & age annotations

III. QUALITY CONTROL & ASSESSMENT
   A. Array Quality Metrics (ArrayQM)
      • RNA degradation analysis
      • Background intensity distribution
      • Positive/negative controls
      • Sample clustering assessment
      • Outlier detection and flagging
   
   B. Quality Filtering
      • Sample removal criteria
      • Platform-specific QC thresholds
      • Batch effect assessment

IV. DATA PREPROCESSING & STANDARDIZATION
   A. Raw Data Acquisition
      • GEO download (GEOquery)
      • Format standardization
      • Metadata parsing
   
   B. Normalization
      • Platform-specific methods (RMA, quantile, etc.)
      • ExpressionSet object creation
      • Sample/probe filtering
   
   C. Probe-to-Gene Mapping
      • Gene annotation (Entrez, HGNC, ENSEMBL)
      • Collapsing duplicate probes
      • Gene-level expression matrices
   
   D. Log2 Fold-Change Computation
      • Per-sample normalization
      • Long-format expression tables
      • Metadata integration

V. PER-CONDITION DIFFERENTIAL EXPRESSION ANALYSIS
   A. Limma-based DEG Analysis
      • Design matrix specification
      • Contrast definitions (treated vs control, by timepoint)
      • Linear model fitting
      • Empirical Bayes shrinkage
      • Multiple testing correction (Benjamini-Hochberg)
      • Output: DEG tables with logFC, t-stat, p-values
   
   B. Gene Set Enrichment Analysis (GSEA)
      • Ranked gene lists (by t-statistic or logFC)
      • msigDB Hallmarks (6 categories)
      • Fgsea implementation
      • NES (normalized enrichment score) computation
      • FDR correction across gene sets
   
   C. Blood Transcription Modules (BTM) Enrichment
      • ssGSEA scoring (single-sample GSEA)
      • BTM module annotation (346+ modules)
      • 16 immune cell/process groups
      • Per-sample module expression scores
   
   D. Human-Mouse Expression Correlation
      • Ortholog mapping (biomaRt)
      • Per-gene Spearman/Pearson correlation
      • Module-level correlation (mean logFC per BTM)
      • Correlation matrices and scatter plots

VI. UNIFIED CROSS-CONDITION ANALYSIS
   A. Integrated DEG Analysis
      • Aggregate DEG results across conditions
      • Shared vs condition-specific genes
      • Overlap statistics (Fisher exact, Jaccard index)
      • Ranking by effect size and consistency
   
   B. RRHO2 (Rank-Rank Hypergeometric Overlap)
      • Pairwise condition comparisons (human vs mouse)
      • Rank transformation and statistical scoring
      • Heatmap visualizations
      • Identification of consistent/discordant responses
   
   C. Consolidated Gene Set Results
      • Meta-GSEA (combining conditions)
      • BTM module patterns across conditions
      • Hallmark enrichment summaries
   
   D. Cross-Species Functional Equivalence
      • Pathway-level conservation assessment
      • Module response coordination (human ≈ mouse)
      • Quantification of translational validity

VII. MACHINE LEARNING & PREDICTIVE MODELING
   A. FIT Model Training (Found In Translation)
      • Feature selection from BTM/hallmark scores
      • Training data: mouse response signatures
      • Algorithm: [specify which ML method]
      • Cross-validation scheme
   
   B. Cross-Species Prediction
      • Apply trained model to human data
      • Accuracy metrics (precision, recall, F1-score)
      • ROC curves and AUC scores
      • Per-condition performance breakdown
   
   C. Performance Benchmarking
      • Compare mouse-derived predictions vs human ground truth
      • Condition-by-condition accuracy
      • Translational reliability quantification

VIII. VALIDATION & EXPERIMENTAL DESIGN
   A. Immunization Experiments
      • YF17D immunization in camels (if conducted)
      • Sample collection protocol (days 0,1,3,7,14)
      • Blood collection and processing
   
   B. Omics Profiling
      • RNAseq vs microarray methodologies
      • Multi-omics integration (if applicable)
   
   C. Model Verification
      • Comparison with experimental validation data
      • Testing predictions from mouse model on camel data
      • Statistical validation of pathway predictions

IX. OUTPUTS & REPRODUCIBILITY
   A. Data Artifacts
      • Processed expression matrices (tables/)
      • DEG & enrichment results
      • Correlation tables
      • Model predictions
   
   B. Publication Figures (534+ total)
      • Correlation scatter plots
      • RRHO2 heatmaps
      • Volcano plots (DEGs)
      • BTM/hallmark enrichment visualizations
      • ROC curves
      • Demography plots
   
   C. Computational Reproducibility
      • renv lock file (package versioning)
      • Exact R/Bioconductor versions
      • Seed setting for randomization
      • Worked example script (example_btm_correlation.R)
```

### Format Details
- **Hierarchical structure:** I → II.A.1 → II.A.1.a etc.
- **Bullet style:** Main steps + sub-steps + implementation details
- **Cross-references:** Link to specific notebooks (e.g., "See 2_Preprocessing.Rmd line 45")
- **Tool/Package names:** Include all software with versions (limma 3.66.0, etc.)
- **Data flow:** Show input files → processing → output files
- **Parameters:** Key threshold values, statistical cutoffs, algorithm choices

---

## Plan: Create Visual Flowchart Diagram

### Flowchart Structure (Inspired by provided image)

**Left Column: Main Analysis Stages**
```
Compilação de dados (Data Curation)
           ↓
Padronização e engenharia de dados (QC & Preprocessing)
           ↓
Análise exploratória (Per-Condition DEG & Enrichment)
           ↓
Análise integrada (Unified Cross-Condition Analysis)
           ↓
Construção de modelos (ML Training & Prediction)
           ↓
Validação experimental (Experimental Verification)
           ↓
Disponibilizar modelo (Publication & Deployment)
```

**Right Column: Detailed Steps per Stage**

For each main stage box, list 3-5 key sub-steps:

1. **Compilação de dados**
   - BioProject search & filtering
   - Sample metadata collection
   - Dataset standardization
   - Linking human & mouse datasets

2. **Padronização e engenharia de dados**
   - Array Quality Metrics (ArrayQM)
   - Outlier detection & removal
   - Platform-specific normalization
   - Probe-to-gene mapping
   - Log2FC computation

3. **Análise exploratória (Per-Condition)**
   - Differential expression (limma)
   - Gene Set Enrichment (GSEA)
   - Blood Transcription Modules (BTM)
   - ssGSEA scoring
   - Human-mouse correlation

4. **Análise integrada**
   - Aggregate DEG analysis
   - RRHO2 overlaps
   - Consolidated gene sets
   - Cross-species functional equivalence

5. **Construção de modelos**
   - FIT model training (mouse features)
   - Cross-species prediction
   - ROC curves & AUC
   - Performance benchmarking

6. **Validação experimental**
   - Immunization experiments
   - RNAseq/microarray profiling
   - Model verification
   - Experimental comparison

7. **Disponibilizar modelo**
   - Publication figures (534+)
   - Processed data tables
   - Reproducibility infrastructure (renv)
   - Worked examples
   - Code sharing (GitHub)

### Data Types & Tools Annotations
- Add file formats: GEO → RDS → CSV/TSV
- Include tool names: GEOquery, limma, fgsea, RRHO2, etc.
- Color-code by analysis type:
  - **Gray boxes:** Data ingestion
  - **Brown boxes:** Analysis
  - **Blue boxes:** Validation
  - **Green boxes:** Output/deployment

### Technology Implementation
- **Format:** Mermaid flowchart syntax (works in Quarto/R Markdown)
- **Output:** Export as SVG/PNG for presentations
- **Style:** Match provided image aesthetic (clean, professional, publication-ready)

---

## Detailed Workplan

### Step 1: Extract & Synthesize Existing Documentation
- [ ] Read README.md (project overview)
- [ ] Read CODEBOOK.md (data dictionary & naming conventions)
- [ ] Review all 6 notebooks (0-5) for methods, parameters, outputs
- [ ] Review required.R for custom functions & dependencies
- [ ] Compile list of tools, versions, packages

### Step 2: Write METHODOLOGY.md (Main Document)
- [ ] Section I: Introduction (copy/adapt from README)
- [ ] Section II: Data Curation (from 0_Data_Curation.Rmd)
- [ ] Section III: QC (from 1_QualityControl.Rmd)
- [ ] Section IV: Preprocessing (from 2_Preprocessing.Rmd)
- [ ] Section V: Per-Condition Analysis (from 3.1)
- [ ] Section VI: Unified Analysis (from 3.2)
- [ ] Section VII: ML (from 4_Performance.Rmd + 5_MachineLearning.Rmd)
- [ ] Section VIII: Validation (design section + experimental context)
- [ ] Section IX: Outputs & Reproducibility
- [ ] Format as hierarchical bullet points
- [ ] Add cross-references to notebooks & tools

### Step 3: Create Mermaid Flowchart Diagram
- [ ] Design main 7-stage pipeline structure
- [ ] Define sub-steps for each stage (3-5 per stage)
- [ ] Add data type annotations (GEO, RDS, CSV, TSV)
- [ ] Add tool/package names
- [ ] Create color coding scheme
- [ ] Write Mermaid syntax
- [ ] Generate PNG/SVG export
- [ ] Save as PROJECT_FLOWCHART.md

### Step 4: Create Supporting Documents
- [ ] PIPELINE_QUICKREF.md (table format: Notebook → Purpose → I/O → Time)
- [ ] SOFTWARE_VERSIONS.md (all packages with versions from renv.lock)
- [ ] GEO_ACCESSIONS.md (data sources table)

### Step 5: Integration & Updates
- [ ] Update README.md with link to METHODOLOGY.md
- [ ] Create index file linking all documentation
- [ ] Verify all cross-references work
- [ ] Check formatting for consistency

---

## Deliverables

### Primary
1. **METHODOLOGY.md** (10-15 pages)
   - Complete project methodology (Sections I-IX)
   - Hierarchical bullet points with sub-details
   - Cross-references to notebooks, tools, data files
   - ~100-150 bullet points with proper nesting
   - Publication-quality formatting

2. **PROJECT_FLOWCHART.md** (Mermaid diagram + visual)
   - 7-stage sequential pipeline
   - Sub-steps with tool/data type labels
   - PNG/SVG export for presentations
   - Matches style of provided example image

### Secondary
3. **PIPELINE_QUICKREF.md** (Quick reference table)
   - Notebook name → Purpose → Key inputs → Key outputs → Est. runtime
   - Dependency information (which notebooks must run before others)

4. **SOFTWARE_VERSIONS.md** (Complete package list)
   - R version, Bioconductor version
   - All packages with exact versions
   - Key dependencies highlighted

5. **Updated README.md**
   - Link to METHODOLOGY.md
   - Link to PROJECT_FLOWCHART.md

### Optional
6. **GEO_ACCESSIONS.md** (Data sources reference)
   - Table of all 7 conditions with GEO accession codes
   - Sample counts, platform info, download links

---

## File Locations

```
animals_vax_atlas/
├── METHODOLOGY.md                    [NEW - Main deliverable]
├── PROJECT_FLOWCHART.md              [NEW - Mermaid + diagram]
├── PIPELINE_QUICKREF.md              [NEW - Quick reference]
├── SOFTWARE_VERSIONS.md              [NEW - Package versions]
├── README.md                         [UPDATED - Add links]
├── CODEBOOK.md                       [REFERENCE]
├── scripts_notebooks/
│   ├── 0_Data_Curation.Rmd
│   ├── 1_QualityControl.Rmd
│   ├── 2_Preprocessing.Rmd
│   ├── 3.1_Comparing_Human_Mouse_ByCondition.Rmd
│   ├── 3.2_Comparing_Human_Mouse_UnifiedAnalyses.Rmd
│   ├── 4_Performance.Rmd
│   ├── 5_MachineLearning.Rmd
│   ├── required.R
│   └── example_btm_correlation.R
└── .positai/plans/
    └── 2026-06-02-plan-project-methodology.md [THIS FILE]
```

---

## Success Criteria

✅ Complete project methodology from data curation to validation documented  
✅ METHODOLOGY.md describes all 5 main analysis notebooks sequentially  
✅ All tools, packages, and versions specified  
✅ Visual flowchart clearly shows pipeline stages & connections  
✅ Format suitable for dissertation appendix or supplementary methods  
✅ Cross-references between all documentation files  
✅ Hierarchical bullet point structure (readable & scannable)  
✅ Matches professional style of provided example image  

---

## User Approval Responses (CONFIRMED ✓)

1. **Language:** ✅ **English** (document primary language)
2. **Flowchart Detail:** ✅ **Detailed** (7 stages + sub-steps + tools)
3. **Experimental Validation:** ✅ **Skip** (focus on computational analysis)
4. **Appendices:** ✅ **All** (data table + thresholds + parameters)

---

## Timeline

- **Phase 1 (Documentation extraction & synthesis):** ~45 min
- **Phase 2 (Write METHODOLOGY.md):** ~90 min
- **Phase 3 (Create flowchart diagram):** ~45 min
- **Phase 4 (Supporting docs + integration):** ~30 min

**Total Estimated Time:** ~3-4 hours

---

## Next Steps (upon user approval)

1. **ExitMode** → Present plan to user
2. User approves/adjusts questions above
3. Implementation begins:
   - Extract notebook content
   - Write hierarchical METHODOLOGY.md
   - Create Mermaid flowchart
   - Generate PNG/SVG visualization
   - Update documentation links
   - Create quick-reference tables

