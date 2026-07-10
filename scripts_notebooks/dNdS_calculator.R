################################################################################
# PIPELINE: De Novo Pairwise dN/dS Calculation via biomaRt and Sequence Alignment
# Target Species: Homo sapiens (Human) vs Mus musculus (Mouse)
# Methodology: Codon-guided alignment & Li (1993) method for Ka/Ks
################################################################################

# ---- 1. Prerequisites & Installation ----
# Ensure all required packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(biomaRt)
library(Biostrings)
library(seqinr)
library(tidyverse)

# ---- 2. Helper Functions for Codon Alignment & dN/dS ----

#' Codon-Guided Back-Translation
#' Maps a protein alignment back to the original nucleotide sequence
back_translate_codon <- function(dna_sequence, aligned_protein_str) {
  dna_str <- as.character(dna_sequence)
  prot_chars <- strsplit(aligned_protein_str, "")[[1]]
  
  aligned_dna <- ""
  dna_idx <- 1
  
  for (aa in prot_chars) {
    if (aa == "-") {
      # Insert a 3-nucleotide gap for a protein deletion/insertion
      aligned_dna <- paste0(aligned_dna, "---")
    } else {
      # Extract the corresponding triplet codon from the original DNA
      codon <- substr(dna_str, dna_idx, dna_idx + 2)
      aligned_dna <- paste0(aligned_dna, codon)
      dna_idx <- dna_idx + 3
    }
  }
  return(aligned_dna)
}

#' Pairwise dN/dS Calculator
#' Performs global protein alignment, guides cDNA, and calculates stats
calculate_pairwise_dnds <- function(human_cds, mouse_cds) {
  # Clean sequences (remove any spaces/newlines and enforce uppercase)
  human_cds <- toupper(gsub("\\s+", "", human_cds))
  mouse_cds <- toupper(gsub("\\s+", "", mouse_cds))
  
  # Validation: Ensure sequences are valid CDS (multiples of 3)
  if (nchar(human_cds) %% 3 != 0 || nchar(mouse_cds) %% 3 != 0) {
    return(data.frame(dN = NA, dS = NA, dN_dS_ratio = NA, note = "Invalid CDS length"))
  }
  
  # Convert to Biostrings objects
  seq_h_dna <- DNAString(human_cds)
  seq_m_dna <- DNAString(mouse_cds)
  
  # Translate Nucleotides to Amino Acids
  prot_h <- Biostrings::translate(seq_h_dna, if.fuzzy.codon = "solve")
  prot_m <- Biostrings::translate(seq_m_dna, if.fuzzy.codon = "solve")
  
  # Global Peptide Alignment (Needleman-Wunsch with BLOSUM62)
  prot_align <- pairwiseAlignment(prot_h, prot_m, substitutionMatrix = "BLOSUM62", 
                                  gapOpening = 10, gapExtension = 0.5)
  
  aligned_prot_h_str <- as.character(pattern(prot_align))
  aligned_prot_m_str <- as.character(subject(prot_align))
  
  # Perform back-translation to obtain codon-aligned cDNA
  codon_aligned_h <- back_translate_codon(seq_h_dna, aligned_prot_h_str)
  codon_aligned_m <- back_translate_codon(seq_m_dna, aligned_prot_m_str)
  
  # Format alignment list for 'seqinr'
  align_list <- list(
    nb = 2,
    nam = c("Human", "Mouse"),
    seq = c(tolower(codon_aligned_h), tolower(codon_aligned_m)),
    com = NA
  )
  class(align_list) <- "alignment"
  
  # Calculate Ka (dN) and Ks (dS) using Li (1993) method
  kaks_res <- tryCatch({
    seqinr::kaks(align_list)
  }, error = function(e) { NULL })
  
  if (is.null(kaks_res)) {
    return(data.frame(dN = NA, dS = NA, dN_dS_ratio = NA, note = "seqinr calculation failed"))
  }
  
  dn <- kaks_res$ka[1, 1]
  ds <- kaks_res$ks[1, 1]
  omega <- if (is.na(ds) || ds == 0) NA else dn / ds
  
  return(data.frame(
    dN = dn,
    dS = ds,
    dN_dS_ratio = omega,
    note = "Success"
  ))
}

# ---- 3. Connect to Ensembl BioMart ----
message("Connecting to Ensembl BioMart...")
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# ---- 4. Retrieve Human-Mouse 1:1 Orthologs Map ----
message("Fetching orthology mapping...")
ortholog_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", 
                 "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"),
  mart = human_mart
) %>% 
  as_tibble() %>%
  filter(mmusculus_homolog_orthology_type == "ortholog_one2one") %>%
  filter(mmusculus_homolog_ensembl_gene != "" & ensembl_gene_id != "")

# NOTE: For demonstration purposes, we sample 10 random genes to prevent 
# BioMart from timing out during heavy sequence downloads. 
# Remove the 'slice_sample' line to run the pipeline for your entire dataset.
target_map <- ortholog_map %>% 
  slice_sample(n = 10) 

message(paste("Processing dN/dS for", nrow(target_map), "sampled genes..."))

# ---- 5. Fetch Coding Sequences (CDS) from BioMart ----
message("Downloading Human coding sequences...")
human_cds_raw <- getSequence(
  id = target_map$ensembl_gene_id,
  type = "ensembl_gene_id",
  seqType = "coding",
  mart = human_mart
) %>% as_tibble()

message("Downloading Mouse coding sequences...")
mouse_cds_raw <- getSequence(
  id = target_map$mmusculus_homolog_ensembl_gene,
  type = "ensembl_gene_id",
  seqType = "coding",
  mart = mouse_mart
) %>% as_tibble()

# ---- 6. Filter for Longest Canonical/Isoform CDS per Gene ----
# BioMart returns multiple rows per gene due to alternative splicing transcripts.
# We pick the longest sequence per stable gene ID to avoid truncation artifacts.
human_cds_clean <- human_cds_raw %>%
  filter(coding != "Sequence unavailable") %>%
  mutate(seq_len = nchar(coding)) %>%
  group_by(ensembl_gene_id) %>%
  slice_max(seq_len, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ensembl_gene_id, human_cds = coding)

mouse_cds_clean <- mouse_cds_raw %>%
  filter(coding != "Sequence unavailable") %>%
  mutate(seq_len = nchar(coding)) %>%
  group_by(ensembl_gene_id) %>%
  slice_max(seq_len, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(mmusculus_homolog_ensembl_gene = ensembl_gene_id, mouse_cds = coding)

# ---- 7. Build the Master Sequence Matrix ----
pipeline_matrix <- target_map %>%
  inner_join(human_cds_clean, by = "ensembl_gene_id") %>%
  inner_join(mouse_cds_clean, by = "mmusculus_homolog_ensembl_gene")

# ---- 8. Execute Pairwise dN/dS Computations ----
message("Running sequence alignments and evolutionary rate calculations...")

results_df <- pipeline_matrix %>%
  mutate(dnds_metrics = map2(human_cds, mouse_cds, ~ calculate_pairwise_dnds(.x, .y))) %>%
  unnest(dnds_metrics)

# ---- 9. Final Cleanup and Output ----
final_dnds_table <- results_df %>%
  select(
    gene_symbol = external_gene_name,
    human_ensembl = ensembl_gene_id,
    mouse_ensembl = mmusculus_homolog_ensembl_gene,
    dN,
    dS,
    dN_dS_ratio,
    status_note = note
  )

print(head(final_dnds_table))

# Export data table for Mixed-Effects Modeling (lme4)
# saveRDS(final_dnds_table, "human_mouse_custom_dnds.rds")
