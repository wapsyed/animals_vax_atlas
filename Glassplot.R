library(ggplot2)
library(ggfx)
library(dplyr)
library(MASS)
library(here) # Para salvar

# ---------------------------------------------------------
# 1. GERAR DADOS
# ---------------------------------------------------------
set.seed(101)
n_genes <- 1000

# Correlação alvo (r = 0.9)
target_r <- 0.9
sigma_matrix <- matrix(c(1, target_r, target_r, 1), nrow = 2)

dados_brutos <- mvrnorm(n = n_genes, mu = c(0, 0), Sigma = sigma_matrix)

df_log2fc <- data.frame(
  human_logfc = dados_brutos[, 1] * 2.5, 
  mouse_logfc = dados_brutos[, 2] * 2.5 
) %>%
  mutate(
    # --- LÓGICA DE RANKING (Para as CAMADAS do Blur) ---
    distancia = abs(human_logfc) + abs(mouse_logfc),
    fake_pval = 1 / (distancia + runif(n(), 0, 2)), 
    pval_rank = rank(fake_pval), # 1 a 1000
    
    # --- LÓGICA DE COR (Para o FILL/COLOR) ---
    # Aqui criamos grupos que NÃO dependem do ranking.
    # Simula vias biológicas diferentes misturadas nos dados.
    pathway_group = sample(c("Inflammation", "Cell Cycle"), 
                           n(), replace = TRUE)
  )

# Calcular R² para o texto
r_squared <- round(cor(df_log2fc$human_logfc, df_log2fc$mouse_logfc)^2, 2)
label_text <- paste0("R² ≈ ", r_squared)

# ---------------------------------------------------------
# 2. DEFINIÇÃO DE CORES
# ---------------------------------------------------------
# Cores para os grupos (e não para o rank)
my_colors <- c(
  "Inflammation" = "#607d8b", # Vermelho
  "Cell Cycle"   = "#455a64"  # Azul
)

# ---------------------------------------------------------
# 3. PLOT COM GGFX
# ---------------------------------------------------------
glass_plot <- ggplot(df_log2fc, aes(x = human_logfc, y = mouse_logfc)) +
  
  # Guias de fundo
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  
  # --- CAMADA 1: BACKGROUND (Todos os 1000 pontos) ---
  # Blur forte. Usamos shape 16 (sem borda) para um desfoque mais suave.
  with_blur(
    geom_point(
      data = . %>% filter(pval_rank > 500),
      aes(color = pathway_group), # Colore pelo grupo
      size = 3,
      alpha = 0.5,
      shape = 16
    ),
    sigma = 5 # Blur alto
  ) +
  
  # --- CAMADA 2: MIDGROUND (Top 500 Rank) ---
  # Blur médio. Filtrado pelo Rank, mas colorido pelo Grupo.
  with_blur(
    geom_point(
      data = . %>% filter(pval_rank <= 500),
      aes(color = pathway_group),
      size = 3,
      alpha = 0.8,
      shape = 16
    ),
    sigma = 2 # Blur médio
  ) +
  
  # --- CAMADA 3: FOREGROUND (Top 100 Rank - Foco) ---
  # Nítido. Usamos shape 21 para ter borda preta e preenchimento (fill) colorido.
  with_shadow(
    geom_point(
    data = . %>% filter(pval_rank <= 300),
    aes(fill = pathway_group,
        color = pathway_group), # Note que aqui usamos FILL, não color
    size = 3,
    # color = "black",
    shape = 21,                # Círculo preenchido
    stroke = 0.2,                # Espessura da borda preta
    alpha = 1
  ),
  x_offset = 2,
  y_offset = 2,
  sigma = 3
  ) +
  
  geom_smooth(method = "lm", se = FALSE, color = "#c59d45", size = 1.5) +
  
  # --- TEMA E ESCALAS ---
  scale_color_manual(values = my_colors) + # Para as camadas 1 e 2 (blur)
  scale_fill_manual(values = my_colors) +  # Para a camada 3 (foco)
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 7)) +
  theme_classic() +
  labs(
    x = "Human Log2FC", 
    y = "Mouse Log2FC", 
    fill = "Pathway", 
    color = "Pathway"
  ) +
  theme(
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "top"
  )

# Exibir o plot
print(glass_plot)

# Salvar
ggsave(filename = here("Figures", "glassplot_test.png"), plot = glass_plot, width = 6, height = 6)

