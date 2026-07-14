#Required
##Packages --------

# 1. Setup: Ensure pacman is available
# if (!require("pacman")) install.packages("pacman")

# 2. Define your package collections
cran_pkgs <- c(
  "tidyverse", "rentrez", "yardstick", "shadowtext", "here", "glue", "ggsci", 
  "NGLVieweR", "janitor", "readr", "maditr", "ggmsa", "ggdist", "ggridges", 
  "see", "Matrix", "RColorBrewer", "ggrepel", "plotly", "corrr", "ggcorrplot", 
  "beepr", "FactoMineR", "factoextra", "esquisse", "gghighlight", "ggh4x", 
  "readxl", "gt", "pracma", "ggnewscale", "ggprism", "ggtext", "devtools", 
  "ggpp", "corto", "Hmisc", "patchwork", "tidymodels", "TidyDensity", "forcats", 
  "GGally", "ggbeeswarm", "geomtextpath", "ggfx", "rstatix", "matrixTests",
  "pkgconfig", "wCorr"
)

bioc_pkgs <- c(
  "pvca", "msa", "Biostrings", "biomaRt", "GEOquery", "circlize", "celldex", 
  "org.Hs.eg.db", "DESeq2", "msigdbr", "org.Mm.eg.db", "ape", "variancePartition", 
  "IRanges", "AnnotationHub", "GenomeInfoDb", "GenomicRanges", "pwalign", 
  "rtracklayer", "GSVA", "sva", "clusterProfiler", "ComplexHeatmap", "edgeR", 
  "limma", "fgsea"
)



# 3. Safe, Non-Destructive Loading Pipeline
# Combining lists to load sequentially via base R
all_packages <- c(cran_pkgs, bioc_pkgs)

for (pkg in all_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste0(
      "\n[ERRO] O pacote '", pkg, "' nao esta instalado no ambiente renv.\n",
      "Por favor, execute no console: renv::install('", pkg, "')\n",
      "Depois disso, reinicie a sessao e rode o script novamente."
    ))
  }
}

# github_pkgs <- c("RRHO2/RRHO2", "YuLab-SMU/ggmsa")
#Download FIT 
# pak::pak('shenorrLabTRDF/FIT.mouse2man')



######### RENV
# # Synchronize your physical project directories with your lockfile
# renv::restore()

# Validate that your sandbox library is fully healthy and synchronized
# renv::status()

# Capture this clean, functional state to lock it down permanently
# renv::snapshot()

# Aesthetics -----
#Custom theme
theme_vaxgo = function(){
  theme_minimal() +
    theme(
      legend.position = "right",
      plot.caption = element_text(hjust = 0, size = 5),
      axis.text.x = element_text(size = 10,
                                 color = "black",
                                 angle = 0),
      axis.text.y = element_text(size = 10,
                                 color = "black",
                                 angle = 0), 
      strip.text = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      axis.title.y = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major.y  = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.key.width = unit(0.4, 'cm'),
      legend.key.height = unit(0.4, 'cm'),
      axis.line.x = element_line(size = 0.5,
                                 colour = "black",
                                 linetype = 1),
      axis.line.y = element_line(size = 0.5 ,
                                 colour = "black",
                                 linetype = 1),
      axis.ticks.x = element_line(size = 0.5, color = "black"),
      axis.ticks.y = element_line(size = 0.5, color = "black"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
      ggh4x.facet.nestline = element_line(colour = "black", linetype = 1)
    )}



#Functions ------

#Function for correlation
safe_cor_test <- function(x, y, method) {
  # Remove NAs and infinite values first
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  # cor.test requires n >= 3 to run without error
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  } else {
    tryCatch({
      ct <- cor.test(x, y, method = method)
      return(list(estimate = unname(ct$estimate), p.value = ct$p.value))
    }, error = function(e) {
      return(list(estimate = NA_real_, p.value = NA_real_))
    })
  }
}


# Function for overlapping -------
overlap_genes <- function(cond1, cond2, data) {
  genes_cond1 <- data$genes[data$process == cond1]
  genes_cond2 <- data$genes[data$process == cond2]
  
  genes_shared <- intersect(genes_cond1, genes_cond2)
  
  genes_notshared_cond1 <- setdiff(genes_cond1, genes_cond2)
  genes_notshared_cond2 <- setdiff(genes_cond2, genes_cond1)
  
  total_genes_cond1 <- length(genes_cond1)
  total_genes_cond2 <- length(genes_cond2)
  
  percentage_shared_cond1 <- length(genes_shared) / total_genes_cond1 * 100
  percentage_shared_cond2 <- length(genes_shared) / total_genes_cond2 * 100
  
  shared_genes <- data.frame(
    Cond1 = cond1,
    Cond2 = cond2,
    Shared = length(genes_shared),
    NotShared_cond1 = length(genes_notshared_cond1),
    NotShared_cond2 = length(genes_notshared_cond2),
    Total_Genes_Cond1 = total_genes_cond1,
    Total_Genes_Cond2 = total_genes_cond2,
    Genes_Names = paste(genes_shared, collapse = ", "),
    Percentage_Shared_Cond1 = percentage_shared_cond1,
    Percentage_Shared_Cond2 = percentage_shared_cond2
  )
  
  return(shared_genes)
}

#Run GSEA for a given contrast
autoGSEA <- function(df, TERM2GENE, geneset_name) {
  
  gsea_results <- list()
  
  conditions <- df$condition %>% unique() %>% as.character()
  
  for (condition_i in conditions) {
    
    degs_condition <- df %>%
      filter(condition == condition_i) %>%
      select(genes, rank) %>%
      distinct() %>%
      arrange(desc(rank)) %>% 
      deframe()   
    
    auto_gsea <- tryCatch({
      
      GSEA(
        geneList = degs_condition,
        TERM2GENE = TERM2GENE,
        minGSSize = 1,
        maxGSSize = 1000,
        pvalueCutoff = 1,
        pAdjustMethod = "BH"
      ) %>%
        as.data.frame() %>%
        arrange(qvalue) %>%
        mutate(
          condition = condition_i,
          gsea_enrichment = geneset_name
        )
      
    }, error = function(e) NULL)
    
    if (!is.null(auto_gsea)) {
      key <- paste(condition_i, geneset_name, sep = "_")
      gsea_results[[key]] <- auto_gsea
    }
  }
  
  return(list(gsea = bind_rows(gsea_results)))
}

#Function for GSEA (Gene Ontology database)
gseGO_all <- function(df, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ontology = "BP") {
  all_conditions <- df %>%
    distinct(condition) %>%
    pull(condition)
  
  gsea_results <- lapply(all_conditions, function(cond) {
    message("Running for condition: ", cond)
    
    ranked_df <- df %>%
      filter(condition == cond) %>%
      arrange(value)
    
    gene_list <- ranked_df %>%
      arrange(-value) %>%
      pull(value, name = genes)
    
    gsea_result <- tryCatch({
      gseGO(
        geneList = gene_list,
        ont = ontology,
        keyType = keyType,
        pvalueCutoff = 1,
        verbose = TRUE,
        OrgDb = OrgDb,
        pAdjustMethod = "BH"
      ) %>%
        as.data.frame() %>%
        arrange(qvalue) %>%
        mutate(
          condition = cond,
          geneset = "Gene ontology, BP"
        ) %>%
        clean_names()
    }, error = function(e) {
      warning("Error in GSEA for condition: ", cond)
      return(NULL)
    })
    
    return(gsea_result)
  })
  
  bind_rows(gsea_results)
}


#Function for ORA/ENRICHER custom
enricher_condition <- function(df, term2gene_df) {
  df %>%
    distinct(condition) %>%
    pull(condition) %>%
    set_names() %>%
    map_dfr(function(cond) {
      genes <- df %>%
        filter(condition == cond) %>%
        pull(genes)
      
      result <- tryCatch({
        enricher(
          gene = genes,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          TERM2GENE = term2gene_df, 
        ) %>%
          as.data.frame() %>%
          clean_names() %>%
          mutate(condition = cond)
      }, error = function(e) {
        message("Failed on condition: ", cond)
        return(NULL)
      })
      
      return(result)
    })
}

# Function for ORA Gene Ontology
oraGO_condition <- function(df, orgdb = org.Hs.eg.db) {
  df %>%
    distinct(condition) %>%
    pull(condition) %>%
    set_names() %>%
    map_dfr(function(cond) {
      genes <- df %>%
        filter(condition == cond) %>%
        pull(genes)
      
      result <- tryCatch({
        enrichGO(
          gene = genes,
          OrgDb = orgdb,
          keyType = 'SYMBOL',
          readable = T,
          ont = "BP",
          pAdjustMethod = "BH",
          qvalueCutoff = 1) %>%
          as.data.frame() %>%
          clean_names() %>%
          mutate(condition = cond)
      }, error = function(e) {
        message("Failed on condition: ", cond)
        return(NULL)
      })
      
      return(result)
    })
}


# Function for clustering by day
cluster_by_day <- function(df, day_column = "day", condition_column = "condition") {
  
  days <- levels(df[[day_column]])
  orders <- list()
  
  # Loop por cada nível (dia) da coluna fator day
  for (day_value in days) {
    # Filtrar os dados para o dia atual
    day_data <- df %>%
      filter(!!sym(day_column) == day_value) %>%
      dplyr::select(!!sym(condition_column)) %>%
      distinct()
    
    # Criar uma matriz binária onde 1 representa a presença da condição
    if (nrow(day_data) > 1) {
      # Matriz com condições como rótulos e 1 para presença
      presence_matrix <- table(day_data[[condition_column]]) > 0
      dist_matrix <- as.matrix(presence_matrix)
      
      # Clustering
      col_clustered <- dist(dist_matrix) %>% hclust()
      
      # Armazenar a ordem dos clusters
      orders[[day_value]] <- col_clustered$labels[col_clustered$order]
    } else {
      # Se houver apenas uma condição, armazenar o rótulo
      orders[[day_value]] <- as.character(day_data[[condition_column]])
    }
  }
  
  # Converter a lista de ordens para um data frame
  order_df <- data.frame(
    day = rep(names(orders), sapply(orders, length)),
    condition = unlist(orders)
  )
  
  return(order_df)
}

#Safe correlation 
safe_cor_test <- function(x, y, method) {
  # Remove NAs and infinite values first
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  # cor.test requires n >= 3 to run without error
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  } else {
    tryCatch({
      ct <- cor.test(x, y, method = method)
      return(list(estimate = unname(ct$estimate), p.value = ct$p.value))
    }, error = function(e) {
      return(list(estimate = NA_real_, p.value = NA_real_))
    })
  }
}


# Function to retrieve sequences
get_sequences <- function(genes, mart, symbol_attr) {
  getBM(
    attributes = c(symbol_attr, "cdna", "transcript_length"), 
    filters = symbol_attr,
    values = genes,
    mart = mart
  )
}


#Colors --------



colors_vaccines = c("Agrippal" = "#669bbc",
                    "Fluad" = "#219ebc",
                    "Quadri 2019-2020" = "#4DBBD5FF",
                    "Quadri 19-20" = "#4DBBD5FF",
                    "Engerix" = "#4361ee")



colors = list(organism = c("FIT" = "#4361ee",
                           "Mouse" = "gray50",
                           "Immune" = "#f72585",
                           "Human" = "black",
                           "Permutation" = "gray80"),
              #Fill
              comparison = c("Mouse only" = "gray50",
                             "Mouse\nonly" = "gray50",
                             "Human only" = "black",
                             "Human\nonly" = "black",
                             "Shared" = "#4DBBD5FF",
                             "Not-shared" = "gray50",
                             "Not core" = "gray90",
                             "Not\ncore" = "gray90",
                             "BTMs" = "#4361ee",
                             "Immune BTMs" = "#4361ee",
                             "Immune\nBTMs" = "#4361ee",
                             "Non-Immune BTMs" = "#3a86ff",
                             "Non-Immune\nBTMs" = "#3a86ff",
                             "All" = "gray90"),
              #Color
              comparison_color = c("Mouse only" = "black",
                                   "Mouse\nonly" = "black",
                                   "Human only" = "white",
                                   "Human\nonly" = "white",
                                   "Not-shared" = "white",
                                   "Shared" = "black",
                                   "BTMs" = "black",
                                   "Immune BTMs" = "black",
                                   "Immune\nBTMs" = "black",
                                   "Non-Immune BTMs" = "black",
                                   "Non-Immune\nBTMs" = "black"),
              comparison_versus = c("Human vs. Mouse" = "#4DBBD5FF",
                                    "Mouse vs. Mouse" = "gray50",
                                    "Human vs. Human" = "black"),
              shared = c("All genes" = "black",
                        "All DEGs" = "#3dccc7",
                        "Shared DEGs" = "purple",
                        "Immune" ="#f72585"),
              direction = c("Up-Up" = "#4DBBD5FF", 
                            "Up-Down" = "#4DBBD5FF", 
                            "Down-Down" = "#4361ee", 
                            "Down-Up" = "#4361ee"),
              treatment = c("Vaccination" = "#006494",
                            "Infection" = "#3dccc7",
                            "Injury" = "#deaaff"),
                     timepoint = c("0h" = "gray50",
                                   "Day 0" = "gray50",
                                   "Early hours" = "#D0D8FB",
                                   "2h" = "#D0D8FB", 
                                   "4h" = "#A1B0F7",
                                   "6h" = "#4361ee",
                                   "12h" = "#3a86ff",
                                   "24h" = "purple",
                                   "Day 1" = "#caf0f8",
                                   "Day 2" = "#ade8f4",
                                   "Day 3" = "#90e0ef",
                                   "Day 4" = "#6CD5EA",
                                   "Day 6" = "#00b4d8",
                                   "Day 7" = "#0096c7",
                                   "Day 12" = "#03045e",
                                   "Day 24" = "#184e77"),
                     sex = c("Male" = "#4DBBD5FF",
                             "Female" = "#4361ee"),
                     race = c("White" = "gray80",
                              "Other" = "#4DBBD5FF",
                              "Black" = "#4361ee", 
                              "Asian" = "#df65b0", 
                              "CB6F1" = "#90A4AEFF"),
                     colors_vaccines = c("Agrippal" = "#669bbc",
                                         "Fluad" = "#219ebc",
                                         "Quadri 2019-2020" = "#4DBBD5FF",
                                         "Quadri 19-20" = "#4DBBD5FF",
                                         "Engerix" = "#4361ee"),
                     type = c("VLP" = "#7E6148FF",
                              "LA" =  "#EFC000FF",
                              "CONJ" = "#F39B7FFF", 
                              'IN'= "#00A087FF", 
                              'Inactivated'= "#00A087FF", 
                              'VV' = "#3C5488FF", #3rd Gen vaccines
                              'RNA' = "#4DBBD5FF",
                              'SU' = "#8491B4FF",
                              "IN/SU" = "#8491B4FF",
                              "PS" = "#F39B7FFF",
                              'I'= "#DC0000FF",
                              "H" = "grey95",
                              "V-I" = "#ffadc7"),
              pathogen = c("Influenza" = "#4cc9f0",
                           "Hepatitis B" = "#a2d2ff",
                           "Saureus" = "#669bbc",
                           "Ecoli" = "#669bbc",
                           "Trauma" = "#3C5488FF",
                           "Burn" = "#184e77"))


#Define colors

btm_immune_groups = list(group = c("SIGNAL TRANSDUCTION" = "#708d81",
                                   "CELL CYCLE" = "#06d6a0",
                                   "ECM AND MIGRATION" = "#52b788",
                                   "ENERGY METABOLISM" = "#95d5b2", 
                                   "INNATE RESPONSE" = "#184e77",
                                   "INFLAMMATORY/TLR/CHEMOKINES"  = "#3a86ff",
                                   "INTERFERON/ANTIVIRAL SENSING" = "#669bbc",
                                   "NEUTROPHILS" = "#219ebc",
                                   "NK CELLS" = "#a2d2ff",
                                   "IFN"= "#4361ee",
                                   "MONOCYTES" = "#4cc9f0",
                                   "DC ACTIVATION" = "#9f86c0",
                                   "PLATELETS" = "#7209b7",
                                   "B CELLS" = "#e5383b",
                                   "T CELLS" = "#f72585",
                                   "PLASMA CELLS" = "#ffafcc"),
                         subgroup = c("SIGNAL TRANSDUCTION" = "#708d81",
                                      "CELL CYCLE" = "#06d6a0",
                                      "ECM AND MIGRATION" = "#52b788",
                                      "ENERGY METABOLISM" = "#95d5b2", 
                                      "INNATE RESPONSE" = "#184e77",
                                      "INFLAMMATORY/TLR/CHEMOKINES"  = "#3a86ff",
                                      "INTERFERON/ANTIVIRAL SENSING" = "#669bbc",
                                      "NEUTROPHILS" = "#219ebc",
                                      "NK CELLS" = "#a2d2ff",
                                      "IFN"= "#4361ee",
                                      "MONOCYTES" = "#4cc9f0",
                                      "ANTIGEN PRESENTATION"= "#7209b7", 
                                      "DC ACTIVATION" = "#9f86c0",
                                      "PLATELETS" = "#7209b7",
                                      "B CELLS" = "#e5383b",
                                      "T CELLS" = "#f72585",
                                      "PLASMA CELLS" = "#ffafcc"))


immune_colors = c("SIGNAL TRANSDUCTION" = "#CA6702",
                  "CELL CYCLE" = "#9B2226",
                  "ECM AND MIGRATION" = "#d08c60",
                  "ENERGY METABOLISM" = "#ffc300", 
                  "PLATELETS" = "#a4b75c",
                  "INNATE RESPONSE" = "#3D405B",
                  "NEUTROPHILS" = "#90A4AEFF",
                  "NK CELLS" = "#001219",
                  "IFN"= "gray80",
                  "MONOCYTES" = "#8491B4FF",
                  "B CELLS" = "#4cc9f0",
                  "T CELLS" = "#4361ee")

immune_order = c("SIGNAL TRANSDUCTION",
                 "CELL CYCLE",
                 "ECM AND MIGRATION",
                 "ENERGY METABOLISM", 
                 "INNATE RESPONSE",
                 "NEUTROPHILS",
                 "NK CELLS",
                 "IFN",
                 "MONOCYTES",
                 "PLATELETS",
                 "B CELLS",
                 "T CELLS")


hallmarks_colors = c("Immune Response" = "#4cc9f0",
                     "Apoptosis and Hormonal Response" = "#06d6a0",
                     "Differentiation and Cell Structure" = "gray80",
                     "Metabolism"  = "#4361ee", 
                     "Proliferation and Repair" = "#001219",
                     "Signaling and Stress Response"   = "#0A9396")

colors_all <- c(
  "#4361ee",
  "#4DBBD5FF",
  "#90A4AEFF",
  "#06d6a0",
  "#46eec6",
  "#175289",
  "#caf0f8",
  "#4cc9f0",
  "#669bbc",
  "purple",
  "#7014cc",
  "#cd61a1",
  "gray80",
  #Neutral harmony bliss
  "#F4F1DE",
  "#E07A5F",
  "#3D405B",
  "#81B29A",
  "#F2CC8F",
  
  #Ocean sunset
  "#001219",
  "#005F73",
  "#0A9396",
  "#94D2BD",
  "#CA6702",
  "#AE2012",
  "#9B2226"
)


