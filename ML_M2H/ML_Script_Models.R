

#Library
library(tidyverse)
library(here)
library(ggsci)
library(janitor)
library(readr)
library(explore)
library(maditr)
library(GEOquery)
library(Matrix)
library(circlize)
library(RColorBrewer)
library(celldex)
library(biomaRt)
library(org.Hs.eg.db)
library(plotly)
library(DESeq2)
library(msigdbr)
library(ape)
library(GSVA)
library(sva)
library(clusterProfiler)
library(pheatmap)
library(EnhancedVolcano)
library(ggbeeswarm)
library(tidyverse)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(esquisse)
library(corto)
library(factoextra)
library(reshape2)
library(ComplexHeatmap)
library(gghighlight)
library(janitor)
library(paletteer)
library(edgeR)
library(limma)
library(mogene10sttranscriptcluster.db)
library(ggh4x)
library(readxl)
library(patchwork)
library(tidymodels)
# devtools::install_github('shenorrLabTRDF/FIT.mouse2man')
library(FIT.mouse2man)


theme_vaxgo = function(){
  theme_minimal() +
    theme(
      legend.position = "right",
      plot.caption = element_text(hjust = 0, size = 5),
      axis.text.x = element_text(size = 10,
                                 color = "black",
                                 angle = 90),
      axis.text.y = element_text(size = 10,
                                 color = "black",
                                 angle = 0),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      axis.title.y = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(size = 0.25, 
                                        linetype = 2,
                                        color = "gray75"),
      panel.grid.major.x = element_line(size = 0.25, 
                                        linetype = 2,
                                        color = "gray75"),
      axis.line.x = element_line(size = 0.5,
                                 colour = "black",
                                 linetype = 1),
      axis.line.y = element_line(size = 0.5 ,
                                 colour = "black",
                                 linetype = 1),
      axis.ticks.x = element_line(size = 0.5, color = "black"),
      axis.ticks.y = element_line(size = 0.5, color = "black"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0.5, 0.2, 0.2, 0.2, "cm")
    )}




#Import files
HS_MM_Symbol_Entrez = FIT.mouse2man::HS_MM_Symbol_Entrez
AllData_V2.0 = FIT.mouse2man::AllData_V2.0 #Training data
AllData_V2.0_blood = read_rds(here("new_models", "AllData_FIT_training_blood.rds"))

