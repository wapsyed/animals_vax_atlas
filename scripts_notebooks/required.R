#Required
##Packages --------

# 1. Setup: Ensure pacman is available
# if (!require("pacman")) install.packages("pacman")

# 2. Define your package collections
cran_pkgs <- c(
  "tidyverse", "rentrez", "yardstick", "shadowtext", "here", "glue", "ggsci", 
  "NGLVieweR", "janitor", "readr", "maditr", "ggmsa", "ggdist", "ggridges", 
  "see", "Matrix", "RColorBrewer", "ggrepel", "plotly", "corrr", "ggcorrplot", 
  "beepr", "FactoMineR", "factoextra", "esquisse", "gghighlight", "ggh4x", "ggExtra",
  "readxl", "gt", "pracma", "ggnewscale", "ggprism", "ggtext", "devtools", 
  "ggpp", "corto", "Hmisc", "patchwork", "tidymodels", "TidyDensity", "forcats", "deeptime",
  "GGally", "ggbeeswarm", "geomtextpath", "ggfx", "rstatix", "matrixTests",  "seqinr", "ape",
  "pkgconfig", "wCorr", "weights", "ggpubr", "rstatix", "emmeans", "effectsize", "coin", "ggtext",
  "ranger", "gghalves", "pROC"
)

bioc_pkgs <- c(
  "pvca", "msa", "Biostrings", "biomaRt", "GEOquery", "circlize", "celldex", 
  "org.Hs.eg.db", "DESeq2", "msigdbr", "org.Mm.eg.db", "ape", "variancePartition", 
  "IRanges", "AnnotationHub", "GenomeInfoDb", "GenomicRanges", "pwalign", 
  "rtracklayer", "GSVA", "sva", "clusterProfiler", "ComplexHeatmap", "edgeR", 
  "limma", "fgsea"
)

# devtools::install_github('erocoar/gghalves')
# install.packages("https://cran.r-project.org/src/contrib/Archive/notifier/notifier_1.0.0.tar.gz")
# install.packages("vip", repos = c("https://bgreenwell.r-universe.dev", "https://cloud.r-project.org"))
library(vip)

# 3. Safe, Non-Destructive Loading Pipeline
# Combining lists to load sequentially via base R
all_packages <- c(cran_pkgs, "notifier", bioc_pkgs)

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



#Colors --------

colors_all <- list(
  #Blues
  blues = c("#4361ee", "#4DBBD5FF", "#3D405B", "#90A4AEFF", "#06d6a0", "#46eec6", 
            "#175289", "#caf0f8", "#4cc9f0", "#669bbc", "purple", "#7014cc", "#cd61a1", "gray80"),
  #Neutral harmony bliss
  neutral_harmony = c("#F4F1DE", "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F"),
  
  #Ocean sunset
  ocean_sunset = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#ffc300", "#CA6702", "#AE2012", "#9B2226")
)



colors_vaccines = c("Agrippal" = "#669bbc",
                    "Fluad" = "#219ebc",
                    "Quadri 2019-2020" = "#4DBBD5FF",
                    "Quadri 19-20" = "#4DBBD5FF",
                    "Engerix" = "#4361ee")



colors = list(organism = c("FIT" = "#4361ee",
                           "Mouse" = "gray50",
                           "Immune" = "#f72585",
                           "Human" = "black",
                           "Control" = "#00A087FF",
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
                             "Not\nLEG" = "gray90",
                             "Not LEG" = "gray90",
                             "BTMs" = "#4361ee",
                             "Immune BTMs" = "#4361ee",
                             "Immune\nBTMs" = "#4361ee",
                             "Non-Immune BTMs" = "#3a86ff",
                             "Non-Immune\nBTMs" = "#3a86ff",
                             "All" = "black",
                             "All genes" = "black"),
              #Color
              comparison_color = c("Mouse only" = "black",
                                   "Mouse\nonly" = "black",
                                   "Human only" = "white",
                                   "Human\nonly" = "white",
                                   "Not-shared" = "white",
                                   "Shared" = "black",
                                   "BTMs" = "black",
                                   "All" = "white",
                                   "All genes" = "white",
                                   "Immune BTMs" = "black",
                                   "Immune\nBTMs" = "black",
                                   "Non-Immune BTMs" = "black",
                                   "Non-Immune\nBTMs" = "black"),
              comparison_versus = c("Human vs. Mouse" = "#4DBBD5FF",
                                    "Mouse vs. Mouse" = "gray50",
                                    "Human vs. Human" = "black"),
              degs_comparison = c("DEG" = "#4DBBD5FF",
                                  "DEG, both" = "#4DBBD5FF", 
                                  "DEG, human" = "black", 
                                  "DEG, mouse" = "gray50",
                                  "not DEG" = NA,
                                  "not DEG, both" = NA), 
              degs_comparison_line = c("DEG" = "#4DBBD5FF",
                                  "DEG, both" = "#4DBBD5FF", 
                                  "DEG, human" = "black", 
                                  "DEG, mouse" = "gray25",
                                  "not DEG" = "gray75",
                                  "not DEG, both" = "gray90"), 
              shared = c("All genes" = "gray75",
                         "All DEGs" = "gray25",
                         "Shared DEGs" = "#3dccc7",
                         "Immune" ="#4361ee"),
              direction = c("Up-Up" = "#4DBBD5FF", 
                            "Up-Down" = "#4DBBD5FF", 
                            "Down-Down" = "#4361ee", 
                            "Down-Up" = "#4361ee",
                            "Same" = "#4DBBD5FF",
                            "Different" = "#4361ee"),
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
                           "HepB" = "#a2d2ff",
                           "S. aureus" = "#669bbc",
                           "E. coli" = "#669bbc",
                           "Trauma" = "#3C5488FF",
                           "Burn" = "#184e77"),
              treatment = c("Vaccination" = "#b9dbf4",
                            "Infection" = "#317ec2",
                            "Injury" = "#155289"),
              feature_sets = c(
                "DGE Baseline"      = "black",
                "+ Sequence"        = "gray25",
                "+ Sequence + TFs"      = "gray50",
                "+ Sequence + TFs + CRE type"      = "#a2d2ff",
                "+ Sequence + TFs + CRE type + CTCF"        = "#4cc9f0",
                "+ Sequence + TFs + CRE type + CTCF + Homology"        = "#184e77",
                "Full + BTM"        = "#4361ee"
              )
              )



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

colors_genesets = list(immune_colors = c("SIGNAL TRANSDUCTION" = "#CA6702",
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
                                         "T CELLS" = "#4361ee"),
                       hallmarks_colors = c("Immune Response" = "#4cc9f0",
                                            "Apoptosis and Hormonal Response" = "#06d6a0",
                                            "Differentiation and Cell Structure" = "gray80",
                                            "Metabolism"  = "#4361ee", 
                                            "Proliferation and Repair" = "#001219",
                                            "Signaling and Stress Response"   = "#0A9396")
                       )


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














# Aesthetics -----
#Custom theme
theme_vaxgo <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      # Geral
      ggh4x.facet.nestline = ggplot2::element_line(colour = "black", linetype = 1),
      
      # Axis
      axis.text            = ggplot2::element_text(size = 10),
      axis.text.x          = ggplot2::element_text(size = 10, color = "black", angle = 0),
      axis.text.y          = ggplot2::element_text(size = 10, color = "black", angle = 0),
      axis.title.x         = ggplot2::element_text(color = "black", face = "bold"),
      axis.title.y         = ggplot2::element_text(color = "black", face = "bold"),
      axis.line.x          = ggplot2::element_line(linewidth = 0.5, colour = "black", linetype = 1),
      axis.line.y          = ggplot2::element_line(linewidth = 0.5, colour = "black", linetype = 1),
      axis.ticks.x         = ggplot2::element_line(linewidth = 0.5, color = "black"),
      axis.ticks.y         = ggplot2::element_line(linewidth = 0.5, color = "black"),
      
      # Legend
      legend.position      = "right",
      legend.location      = "plot",
      legend.text          = ggplot2::element_text(size = 10),
      legend.title         = ggplot2::element_text(size = 10),
      legend.key.width     = grid::unit(0.4, "cm"),
      legend.key.height    = grid::unit(0.4, "cm"),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,-10,0),
      
      # Panel
      panel.border         = ggplot2::element_blank(),
      panel.grid.major.x   = ggplot2::element_blank(),
      panel.grid.major.y   = ggplot2::element_blank(),
      panel.grid.minor     = ggplot2::element_blank(),
      panel.spacing        = grid::unit(0.2, "cm"),
      
      # Strip
      strip.text           = ggplot2::element_text(size = 10, color = "black", margin = margin(1, 0, 1, 0)),
      
      # Plot
      plot.title           = ggplot2::element_text(size = 10),
      plot.subtitle        = ggplot2::element_text(size = 10),
      plot.caption         = ggplot2::element_text(hjust = 0, size = 5),
      plot.margin          = ggplot2::margin(0.1, 0.2, 0.1, 0.1, "cm")
    )
}



#Functions ------

#Function for correlation
# safe_cor_test <- function(x, y, method) {
#   # Remove NAs and infinite values first
#   ok <- is.finite(x) & is.finite(y)
#   x <- x[ok]
#   y <- y[ok]
#   
#   # cor.test requires n >= 3 to run without error
#   if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
#     return(list(estimate = NA_real_, p.value = NA_real_))
#   } else {
#     tryCatch({
#       ct <- cor.test(x, y, method = method)
#       return(list(estimate = unname(ct$estimate), p.value = ct$p.value))
#     }, error = function(e) {
#       return(list(estimate = NA_real_, p.value = NA_real_))
#     })
#   }
# }


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


#Weighted correlation ====
calculate_weighted_correlation_spearman <- function(data,
                                           x_col = mean_log2fc_Human,
                                           y_col = mean_log2fc_Mouse,
                                           weight_col = weight_joint_gsea_fdr,
                                           group_cols = c("pathogen", "timepoint_comparison")) {
  
  data %>%
    # Select target variables and rename them for internal calculation
    dplyr::select(
      dplyr::all_of(group_cols),
      x_val = {{ x_col }},
      y_val = {{ y_col }},
      weight_val = {{ weight_col }}
    ) %>%
    
    # Compute weighted Spearman correlation per specified group with safety checks
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      {
        # Clean vectors by filtering out NAs and non-positive weights
        valid_idx <- !is.na(x_val) & !is.na(y_val) & !is.na(weight_val) & weight_val > 0
        x_clean <- x_val[valid_idx]
        y_clean <- y_val[valid_idx]
        w_clean <- weight_val[valid_idx]
        
        # Guard clause: Require at least 3 complete cases with non-zero variance
        if (length(x_clean) < 3 || stats::var(x_clean) == 0 || stats::var(y_clean) == 0) {
          data.frame(
            cor_value = NA_real_,
            se        = NA_real_,
            t         = NA_real_,
            p_value   = NA_real_,
            n_obs     = length(x_clean)
          )
        } else {
          # Safely execute weighted correlation on ranks (Spearman transformation)
          res_try <- tryCatch(
            weights::wtd.cor(
              x = rank(x_clean),
              y = rank(y_clean),
              weight = w_clean
            ),
            error = function(e) NULL
          )
          
          if (is.null(res_try)) {
            data.frame(
              cor_value = NA_real_,
              se        = NA_real_,
              t         = NA_real_,
              p_value   = NA_real_,
              n_obs     = length(x_clean)
            )
          } else {
            data.frame(
              cor_value = as.numeric(res_try[1, "correlation"]),
              se        = as.numeric(res_try[1, "std.err"]),
              t         = as.numeric(res_try[1, "t.value"]),
              p_value   = as.numeric(res_try[1, "p.value"]),
              n_obs     = length(x_clean)
            )
          }
        }
      },
      .groups = "drop"
    ) %>%
    
    # Calculate R-squared, FDR-adjusted p-values, and formatted label strings
    dplyr::mutate(
      # Raw p-value significance stars
      p_label = dplyr::case_when(
        is.na(p_value)  ~ "ns",
        p_value <= 0.01 ~ "***",
        p_value <= 0.05 ~ "**",
        p_value <= 0.10 ~ "*",
        TRUE            ~ "ns"
      ),
      
      # BH adjusted p-values and significance stars
      p_adj = stats::p.adjust(p_value, method = "BH"),
      p_adj_label = dplyr::case_when(
        is.na(p_adj)    ~ "ns",
        p_adj <= 0.01   ~ "***",
        p_adj <= 0.05   ~ "**",
        p_adj <= 0.10   ~ "*",
        TRUE            ~ "ns"
      )
    ) %>% 
    mutate(
      # Consolidated text label for plot annotations
      cor_label = dplyr::if_else(
        !is.na(cor_value),
        stringr::str_c(
          "r=", round(cor_value, 2),
          p_label, ", n=", n_obs
        ),
        "NA"
      )
    )
}

calculate_weighted_correlation_pearson <- function(data,
                                           x_col = mean_log2fc_Human,
                                           y_col = mean_log2fc_Mouse,
                                           weight_col = weight_joint_gsea_fdr,
                                           group_cols = c("pathogen", "timepoint_comparison")) {
  
  data %>%
    # Select target variables and rename them for internal calculation
    dplyr::select(
      dplyr::all_of(group_cols),
      x_val = {{ x_col }},
      y_val = {{ y_col }},
      weight_val = {{ weight_col }}
    ) %>%
    
    # Compute weighted Spearman correlation per specified group with safety checks
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      {
        # Clean vectors by filtering out NAs and non-positive weights
        valid_idx <- !is.na(x_val) & !is.na(y_val) & !is.na(weight_val) & weight_val > 0
        x_clean <- x_val[valid_idx]
        y_clean <- y_val[valid_idx]
        w_clean <- weight_val[valid_idx]
        
        # Guard clause: Require at least 3 complete cases with non-zero variance
        if (length(x_clean) < 3 || stats::var(x_clean) == 0 || stats::var(y_clean) == 0) {
          data.frame(
            cor_value = NA_real_,
            se        = NA_real_,
            t         = NA_real_,
            p_value   = NA_real_,
            n_obs     = length(x_clean)
          )
        } else {
          # Safely execute weighted correlation on ranks (Spearman transformation)
          res_try <- tryCatch(
            weights::wtd.cor(
              x = x_clean,
              y = y_clean,
              weight = w_clean
            ),
            error = function(e) NULL
          )
          
          if (is.null(res_try)) {
            data.frame(
              cor_value = NA_real_,
              se        = NA_real_,
              t         = NA_real_,
              p_value   = NA_real_,
              n_obs     = length(x_clean)
            )
          } else {
            data.frame(
              cor_value = as.numeric(res_try[1, "correlation"]),
              se        = as.numeric(res_try[1, "std.err"]),
              t         = as.numeric(res_try[1, "t.value"]),
              p_value   = as.numeric(res_try[1, "p.value"]),
              n_obs     = length(x_clean)
            )
          }
        }
      },
      .groups = "drop"
    ) %>%
    
    # Calculate R-squared, FDR-adjusted p-values, and formatted label strings
    dplyr::mutate(
      r_squared = dplyr::if_else(!is.na(cor_value), cor_value^2, NA_real_),
      
      # Raw p-value significance stars
      p_label = dplyr::case_when(
        is.na(p_value)  ~ "ns",
        p_value <= 0.01 ~ "***",
        p_value <= 0.05 ~ "**",
        p_value <= 0.10 ~ "*",
        TRUE            ~ "ns"
      ),
      
      # BH adjusted p-values and significance stars
      p_adj = stats::p.adjust(p_value, method = "BH"),
      p_adj_label = dplyr::case_when(
        is.na(p_adj)    ~ "ns",
        p_adj <= 0.01   ~ "***",
        p_adj <= 0.05   ~ "**",
        p_adj <= 0.10   ~ "*",
        TRUE            ~ "ns"
      )
    ) %>% 
    mutate(
      # Consolidated text label for plot annotations
      cor_label = dplyr::if_else(
        !is.na(cor_value),
        stringr::str_c(
          "r=", round(cor_value, 2),
          ", R²=", round(r_squared, 2),
          p_label, ", n=", n_obs
        ),
        "NA"
      )
    )
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


# CLASSIFICATION PERFORMANCE
# Compute classification performance for GENES and DEGs ==============
classification_performance_compute_genes_degs = function(df = genes_df_input,
                                              pval_filter = 0.05) {
  
  #Compute probabilities
  dge_long_genes_up =
    genes_df_input %>%
    mutate(
      df = n - 1,
      pval_one_sided = pt(t, df, lower.tail = FALSE),
      pval_one_sided_adj = p.adjust(pval_one_sided, method = "BH"),
      direction = if_else(pval_one_sided_adj <= pval_filter & t > 0, "up", "not_up"),
      prob = 1 - pval_one_sided_adj
    ) %>%
    select(-adj_p_val) %>%
    rename(adj_p_val = pval_one_sided_adj)
  
  # Down
  dge_long_genes_down = genes_df_input %>%
    mutate(
      df = n - 1,
      pval_one_sided = pt(t, df, lower.tail = TRUE),
      pval_one_sided_adj = p.adjust(pval_one_sided, method = "BH"),
      direction = if_else(pval_one_sided_adj <= pval_filter & t < 0, "down", "not_down"),
      prob = 1 - pval_one_sided_adj
    ) %>%
    select(-adj_p_val) %>%
    rename(adj_p_val = pval_one_sided_adj)
  
  #Check
  dge_long_genes_up %>% 
    dplyr::count(pathogen, timepoint_comparison, organism, direction)
  #Check
  dge_long_genes_down %>% 
    dplyr::count(pathogen, timepoint_comparison, organism, direction)
  
  # Assign human truth==========
  truth_human_up = dge_long_genes_up %>% 
    select(human_symbol, pathogen, organism, timepoint_comparison, adj_p_val, direction) %>% 
    filter(organism == "Human") %>% 
    distinct() %>% 
    pivot_wider(names_from = organism, values_from = direction) %>% 
    rename(truth = Human) %>% 
    mutate(organism = "Human")
  
  
  truth_human_down = dge_long_genes_down %>% 
    select(human_symbol, pathogen, organism, timepoint_comparison, adj_p_val, direction) %>% 
    filter(organism == "Human") %>% 
    distinct() %>% 
    pivot_wider(names_from = organism, values_from = direction) %>% 
    rename(truth = Human) %>% 
    mutate(organism = "Human")
  
  #Generate reference distributions ========
  perfect_prediction_up = truth_human_up %>%
    mutate(prob = if_else(truth == "up", 1, 0),
           model = "Human") %>%
    rename(pred = truth)
  
  perfect_prediction_down = truth_human_down %>%
    mutate(prob = if_else(truth == "down", 1, 0),
           model = "Human") %>%
    rename(pred = truth)
  
  # Combine models (Mouse + Immune + Permutation + Perfect) ========
  models_up = dge_long_genes_up %>% 
    filter(organism != "Human") %>% 
    mutate(pred = direction,
           organism,
           model = organism) %>%
    bind_rows(perfect_prediction_up)
  
  models_down = dge_long_genes_down %>% 
    filter(organism != "Human") %>% 
    mutate(pred = direction,
           model = organism) %>%
    bind_rows(perfect_prediction_down)
  
  
  # Classify probabilities  ========
  #Up
  classified_df_probs_up = truth_human_up %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_up %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("up", "not_up"))
    ) %>% 
    distinct()
  
  #Down
  classified_df_probs_down = truth_human_down %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_down %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("down", "not_down"))
    ) %>% 
    distinct() 
  
  #Up DEGs
  classified_df_probs_up_DEGs = truth_human_up %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_up %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("up", "not_up"))
    ) %>% 
    distinct() %>% 
    inner_join(degs_genes_human %>% 
                 select(-organism), 
               by = join_by(human_symbol, pathogen, timepoint_comparison))
  
  #Down DEGs
  classified_df_probs_down_DEGs = truth_human_down %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_down %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("down", "not_down"))
    ) %>% 
    distinct()  %>% 
    inner_join(degs_genes_human %>% 
                 select(-organism), 
               by = join_by(human_symbol, pathogen, timepoint_comparison))
  
  #Check
  classified_df_probs_up  %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  classified_df_probs_up_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  classified_df_probs_down %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  classified_df_probs_down_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  #Compute ROC  ========

  #UP-notUP
  roc_curve_up_notup = classified_df_probs_up %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  roc_curve_up_notup_DEGs = classified_df_probs_up_DEGs %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  roc_curve_up_notup %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  roc_curve_up_notup_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #DOWN-notDOWN
  roc_curve_down_notdown = classified_df_probs_down %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation",  "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  roc_curve_down_notdown_DEGs = classified_df_probs_down_DEGs %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  roc_curve_down_notdown %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  roc_curve_down_notdown_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #Compute AUC  ========
  auc_results_up = classified_df_probs_up %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  auc_results_up_DEGs = classified_df_probs_up_DEGs %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  auc_results_down = classified_df_probs_down %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  auc_results_down_DEGs = classified_df_probs_down_DEGs %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  
  # Unite all results  ========
  roc_curve_auc_df = bind_rows(auc_results_up,
                               auc_results_down) %>% 
    inner_join(
      roc_curve_up_notup %>% 
        bind_rows(roc_curve_down_notdown %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison)) 
  
  
  roc_curve_auc_df_DEGs = bind_rows(auc_results_up_DEGs,
                                    auc_results_down_DEGs) %>% 
    inner_join(
      roc_curve_up_notup_DEGs %>% 
        bind_rows(roc_curve_down_notdown_DEGs %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison)) 
  #Check
  roc_curve_auc_df %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, auc)
  
  roc_curve_auc_df_DEGs %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, auc)
  
  #Get n genes
  ngenes = bind_rows(
    classified_df_probs_up %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  
  # Combine AUC
  roc_auc_df <- bind_rows(auc_results_up,
                          auc_results_down) %>%
    inner_join(ngenes, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  #Get n genes
  ngenes_DEGs = bind_rows(
    classified_df_probs_up_DEGs %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down_DEGs %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  # Visualize AUC
  roc_auc_df_DEGs <- bind_rows(auc_results_up_DEGs,
                               auc_results_down_DEGs) %>%
    inner_join(ngenes_DEGs, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  #Compute PR  ========

  pr_curve_up_notup = classified_df_probs_up %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  pr_curve_up_notup_DEGs = classified_df_probs_up_DEGs %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  pr_curve_up_notup %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  pr_curve_up_notup_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  pr_curve_down_notdown = classified_df_probs_down %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  pr_curve_down_notdown_DEGs = classified_df_probs_down_DEGs %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  pr_curve_down_notdown %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  pr_curve_down_notdown_DEGs %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #Compute PR-AUC  ========
  pr_auc_results_up = classified_df_probs_up %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc, comparison)
  
  pr_auc_results_up_DEGs = classified_df_probs_up_DEGs %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc,  comparison)
  
  pr_auc_results_down = classified_df_probs_down %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc,  comparison)
  
  pr_auc_results_down_DEGs = classified_df_probs_down_DEGs %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc,  comparison)
  
  
  #Unite all results  ========
  pr_curve_auc_df = bind_rows(pr_auc_results_up,
                              pr_auc_results_down) %>% 
    inner_join(
      pr_curve_up_notup %>% 
        bind_rows(pr_curve_down_notdown %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison)) 
  
  
  
  pr_curve_auc_df_DEGs = bind_rows(pr_auc_results_up_DEGs,
                                   pr_auc_results_down_DEGs) %>% 
    inner_join(
      pr_curve_up_notup_DEGs %>% 
        bind_rows(pr_curve_down_notdown_DEGs %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison))
  
  #Check
  pr_curve_auc_df %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, pr_auc)
  
  pr_curve_auc_df_DEGs %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, pr_auc)
  
  #Get n genes
  ngenes = bind_rows(
    classified_df_probs_up %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  
  # Combine AUC
  pr_auc_df <- bind_rows(pr_auc_results_up,
                         pr_auc_results_down) %>%
    inner_join(ngenes, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  
  #Get n genes
  ngenes_DEGs = bind_rows(
    classified_df_probs_up_DEGs %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down_DEGs %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  
  # Combine and add N genes AUC
  pr_auc_DEGs_df <- bind_rows(pr_auc_results_up_DEGs,
                              pr_auc_results_down_DEGs) %>%
    inner_join(ngenes_DEGs, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  
  
  # RETURN all computed tables --------
  return(list(
    roc_auc_df = roc_auc_df,
    roc_auc_df_DEGs = roc_auc_df_DEGs,
    roc_curve_auc_df = roc_curve_auc_df,
    roc_curve_auc_df_DEGs = roc_curve_auc_df_DEGs,
    pr_curve_auc_df = pr_curve_auc_df,
    pr_curve_auc_df_DEGs = pr_curve_auc_df_DEGs,
    pr_auc_df = pr_auc_df,
    pr_auc_DEGs_df = pr_auc_DEGs_df
  ))
}










# Compute classification performance MODULES ==============

classification_performance_compute_modules = function(df = btms_df_input,
                                                            effect_var = mean_log2fc) {
  
  #Compute probabilities  ===============
  dge_long_genes_up =
    btms_df_input %>%
    mutate(
      df = n_genes - 1,
      t = mean_log2fc/(sd_log2fc/sqrt(n_genes)), 
      pval_one_sided = pt(t, df, lower.tail = FALSE),
      pval_one_sided_adj = p.adjust(pval_one_sided, method = "BH"),
      direction = if_else(pval_one_sided_adj <= 0.05 & t > 0, "up", "not_up"),
      prob = 1 - pval_one_sided_adj
    ) %>%
    rename(adj_p_val = pval_one_sided_adj)
  
  # Down
  dge_long_genes_down = btms_df_input %>% 
    mutate(
      df = n_genes - 1,
      t = mean_log2fc/(sd_log2fc/sqrt(n_genes)), 
      pval_one_sided = pt(t, df, lower.tail = TRUE),
      pval_one_sided_adj = p.adjust(pval_one_sided, method = "BH"),
      direction = if_else(pval_one_sided_adj <= 0.05 & t < 0, "down", "not_down"),
      prob = 1 - pval_one_sided_adj
    ) %>%
    rename(adj_p_val = pval_one_sided_adj) %>% 
    drop_na(adj_p_val)
  
  #Check
  dge_long_genes_up %>% 
    dplyr::count(pathogen, timepoint_comparison, organism, direction)
  #Check
  dge_long_genes_down %>% 
    dplyr::count(pathogen, timepoint_comparison, organism, direction)
  
  # Assign human truth ===============
  truth_human_up = dge_long_genes_up %>% 
    select(pathogen, human_symbol, organism, timepoint_comparison, adj_p_val, direction) %>% 
    filter(organism == "Human") %>% 
    distinct() %>% 
    pivot_wider(names_from = organism, values_from = direction) %>% 
    rename(truth = Human) %>% 
    mutate(organism = "Human") 
  
  truth_human_down = dge_long_genes_down %>% 
    select(pathogen, human_symbol, organism, timepoint_comparison, adj_p_val, direction) %>% 
    filter(organism == "Human") %>% 
    distinct() %>% 
    pivot_wider(names_from = organism, values_from = direction) %>% 
    rename(truth = Human) %>% 
    mutate(organism = "Human")
  
  #Generate reference distributions ===============
  perfect_prediction_up = truth_human_up %>%
    mutate(prob = if_else(truth == "up", 1, 0),
           model = "Human") %>%
    rename(pred = truth)
  
  perfect_prediction_down = truth_human_down %>%
    mutate(prob = if_else(truth == "down", 1, 0),
           model = "Human") %>%
    rename(pred = truth)
  
  # Combine models   ===============
  models_up = dge_long_genes_up %>% 
    filter(organism != "Human") %>% 
    mutate(pred = direction,
           organism,
           model = organism) %>%
    bind_rows(perfect_prediction_up)
  
  models_down = dge_long_genes_down %>% 
    filter(organism != "Human") %>% 
    mutate(pred = direction,
           model = organism) %>%
    bind_rows(perfect_prediction_down)
  
  
  # Classify probabilities  ===============
  #Up
  classified_df_probs_up = truth_human_up %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_up %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("up", "not_up"))
    ) %>% 
    distinct()
  
  #Down
  classified_df_probs_down = truth_human_down %>% 
    select(human_symbol, pathogen, timepoint_comparison, truth) %>%
    inner_join(models_down %>% 
                 select(human_symbol, organism, pathogen, timepoint_comparison, pred, prob, model),
               by = join_by(human_symbol, timepoint_comparison, pathogen)) %>% 
    mutate(
      truth = factor(truth, levels = c("down", "not_down"))
    ) %>% 
    distinct() 
  
  
  #Check
  classified_df_probs_up  %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  classified_df_probs_down %>% 
    dplyr::count(pathogen, model, timepoint_comparison, pred, truth)
  
  #Compute ROC  ===============

  #UP-notUP
  roc_curve_up_notup = classified_df_probs_up %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  roc_curve_up_notup %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #DOWN-notDOWN
  roc_curve_down_notdown = classified_df_probs_down %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation",  "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  roc_curve_down_notdown %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #Compute AUC  ===============
  auc_results_up = classified_df_probs_up %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  auc_results_down = classified_df_probs_down %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Control", "Permutation", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    roc_auc(truth, prob) %>%
    rename(auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, auc, comparison)
  
  
  ##### Unite all results
  roc_curve_auc_df = bind_rows(auc_results_up,
                               auc_results_down) %>% 
    inner_join(
      roc_curve_up_notup %>% 
        bind_rows(roc_curve_down_notdown %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison)) 

  #Check
  roc_curve_auc_df %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, auc)
 
  #Get n genes
  ngenes = bind_rows(
    classified_df_probs_up %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  # Unite 
  roc_auc_df <- bind_rows(auc_results_up,
                          auc_results_down) %>%
    inner_join(ngenes, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  #Compute PR  ===============
  pr_curve_up_notup = classified_df_probs_up %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Up")
  
  pr_curve_up_notup %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  pr_curve_down_notdown = classified_df_probs_down %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Control", "Human", "FIT")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_curve(truth, prob) %>%
    mutate(comparison = "Down")
  
  pr_curve_down_notdown %>% 
    dplyr::count(pathogen, model, timepoint_comparison)
  
  #Compute PR-AUC  ===============
  pr_auc_results_up = classified_df_probs_up %>% 
    mutate(truth = fct_relevel(truth, "up", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Up") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc, comparison)
  
  pr_auc_results_down = classified_df_probs_down %>% 
    mutate(truth = fct_relevel(truth, "down", after = 0)) %>% 
    filter(model %in% c("Mouse", "Immune", "Permutation", "Human", "FIT", "Control")) %>%
    group_by(pathogen, model, timepoint_comparison) %>%
    pr_auc(truth, prob) %>%
    rename(pr_auc = .estimate) %>%
    mutate(comparison = "Down") %>% 
    ungroup() %>% 
    select(pathogen, model, timepoint_comparison, pr_auc,  comparison)
  
  ##### Unite all results
  pr_curve_auc_df = bind_rows(pr_auc_results_up,
                              pr_auc_results_down) %>% 
    inner_join(
      pr_curve_up_notup %>% 
        bind_rows(pr_curve_down_notdown %>% 
                    filter(model != "Human")),
      by = join_by(pathogen, model, timepoint_comparison, comparison)) 
  
  #Check
  pr_curve_auc_df %>% 
    dplyr::count(pathogen, model, comparison, timepoint_comparison, pr_auc)
  
  #Get n genes
  ngenes = bind_rows(
    classified_df_probs_up %>%
      filter(pred == "up") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Up") %>% 
      ungroup(),
    classified_df_probs_down %>%
      filter(pred == "down") %>% 
      group_by(pathogen, timepoint_comparison, model) %>% 
      summarise(n_genes = n(),
                comparison = "Down") %>% 
      ungroup(),
  )
  
  
  # Unite results
  pr_auc_df <- bind_rows(pr_auc_results_up,
                         pr_auc_results_down) %>%
    inner_join(ngenes, by = join_by(comparison, pathogen, timepoint_comparison, model)) 
  
  
  
  # RETURN all computed tables --------
  return(list(
    roc_auc_df = roc_auc_df,
    roc_curve_auc_df = roc_curve_auc_df,
    pr_curve_auc_df = pr_curve_auc_df,
    pr_auc_df = pr_auc_df
  ))
}



# Codon Alignment ----

# Molecular distance
calculate_pairwise_dna_distance <- function(human_cds, mouse_cds) {
  # Clean sequence inputs (remove whitespace/newlines and convert to uppercase)
  human_cds <- toupper(gsub("\\s+", "", human_cds))
  mouse_cds <- toupper(gsub("\\s+", "", mouse_cds))
  
  # Validation: Ensure non-empty sequences
  if (nchar(human_cds) == 0 || nchar(mouse_cds) == 0) {
    return(tibble::tibble(
      dist_jc69    = NA_real_,
      dist_k80     = NA_real_,
      status_dist  = "Empty sequence"
    ))
  }
  
  # Step 1: Global nucleotide alignment to guarantee matching length and gap positions
  dna_align <- tryCatch({
    pwalign::pairwiseAlignment(
      Biostrings::DNAString(human_cds),
      Biostrings::DNAString(mouse_cds),
      type = "global"
    )
  }, error = function(e) { NULL })
  
  if (is.null(dna_align)) {
    return(tibble::tibble(
      dist_jc69    = NA_real_,
      dist_k80     = NA_real_,
      status_dist  = "Alignment failed"
    ))
  }
  
  # Step 2: Extract aligned sequence strings with gap characters ('-')
  aligned_h <- as.character(pwalign::alignedPattern(dna_align))
  aligned_m <- as.character(pwalign::alignedSubject(dna_align))
  
  # Step 3: Construct character matrix and convert via ape::as.alignment & ape::as.DNAbin
  ab <- rbind(
    unlist(strsplit(aligned_h, "")),
    unlist(strsplit(aligned_m, ""))
  )
  
  bin_dna <- tryCatch({
    ape::as.DNAbin(ape::as.alignment(ab))
  }, error = function(e) { NULL })
  
  if (is.null(bin_dna)) {
    return(tibble::tibble(
      dist_jc69    = NA_real_,
      dist_k80     = NA_real_,
      status_dist  = "DNAbin conversion failed"
    ))
  }
  
  # Step 4: Calculate distances using Jukes-Cantor (JC69) and Kimura (K80)
  jc69_val <- tryCatch({
    as.numeric(ape::dist.dna(bin_dna, model = "JC69")[1])
  }, error = function(e) { NA_real_ })
  
  k80_val <- tryCatch({
    as.numeric(ape::dist.dna(bin_dna, model = "K80")[1])
  }, error = function(e) { NA_real_ })
  
  return(tibble::tibble(
    dist_jc69   = jc69_val,
    dist_k80    = k80_val,
    status_dist = if (!is.na(jc69_val)) "Success" else "Distance calculation failed"
  ))
}



