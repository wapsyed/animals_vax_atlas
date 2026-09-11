library(here)
source(here("scripts_notebooks", "required.R"))
library(tidyverse)
library(patchwork)

fig_dir <- here("tests", "v2", "Figures")
tab_dir <- here("tests", "v2", "Tables")

task1 <- read_csv(file.path(tab_dir, "task1_convrel_loco.csv"), show_col_types = FALSE)
task2 <- read_csv(file.path(tab_dir, "task2_shared_loco.csv"), show_col_types = FALSE)
task3a <- read_csv(file.path(tab_dir, "task3a_rankhuman_loco.csv"), show_col_types = FALSE)
task3b <- read_csv(file.path(tab_dir, "task3b_signconc_loco.csv"), show_col_types = FALSE)
score_tbl <- readRDS(here("tests", "v2", "Models", "score_table_v2.rds"))
pooled <- score_tbl %>% filter(treatment == "Pooled")

# ============================== (b) classification ==============================
b_dat <- task2 %>%
  pivot_longer(c(roc_auc, pr_auc), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, roc_auc = "ROC-AUC", pr_auc = "PR-AUC"),
         framing = factor(recode(framing, A_predictive = "Predictive (A)", B_explanatory = "Explanatory (B)"),
                          levels = c("Predictive (A)", "Explanatory (B)")))

p_b <- ggplot(b_dat, aes(x = feature_set, y = value, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  facet_grid(metric ~ framing, scales = "free_y") +
  scale_fill_manual(values = c("Random Forest" = "#4DBBD5FF", "Logistic Regression" = "gray70")) +
  theme_vaxgo() +
  labs(x = "", y = "LOCO CV estimate", fill = "Algorithm",
       title = "b  Shared LEG classification")

# ============================== (c) regression + direction ==============================
c_conv <- task1 %>% filter(model == "Random Forest") %>%
  transmute(metric = "R2", task = "Convergence", framing = recode(framing, A_predictive = "Predictive (A)", B_explanatory = "Explanatory (B)"),
            step = feature_set, value = rsq)
c_rank <- task3a %>% filter(model == "Random Forest") %>%
  transmute(metric = "R2", task = "Human rank transfer", framing = "Transfer", step = feature_set, value = rsq)
c_dir <- task3b %>% filter(model == "Random Forest") %>%
  transmute(metric = "ROC-AUC", task = "Directional concordance", framing = "Transfer", step = feature_set, value = roc_auc)
c_dat <- bind_rows(c_conv, c_rank, c_dir) %>%
  mutate(task = factor(task, levels = c("Convergence", "Human rank transfer", "Directional concordance")),
         metric = factor(metric, levels = c("R2", "ROC-AUC")))

p_c <- ggplot(c_dat, aes(x = step, y = value, fill = framing)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  facet_grid(metric ~ task, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Predictive (A)" = "#4361ee",
                               "Explanatory (B)" = "#4cc9f0",
                               "Transfer" = "gray50")) +
  theme_vaxgo() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7)) +
  labs(x = "", y = "LOCO CV estimate", fill = "Framing",
       title = "c  Convergence, transfer and direction")

# ============================== (d) score_shared vs rank ==============================
p_d <- ggplot(pooled, aes(x = score_shared, y = rank_mouse, fill = status_original)) +
  geom_point(shape = 21, size = 2, alpha = 0.3, stroke = 0) +
  scale_fill_manual(values = colors$comparison) +
  theme_vaxgo() +
  labs(x = "score_shared", y = "Rank (Mouse)", fill = "LEG status",
       title = "d  Shared LEG score")

# ============================== (e) score_conv vs rank_diff ==============================
cor_lab <- function(x, y) {
  tibble(rho = cor(x, y, method = "spearman", use = "complete.obs"),
         r2  = summary(lm(y ~ x))$r.squared,
         p   = cor.test(x, y, method = "spearman", exact = FALSE)$p.value) %>%
    mutate(padj = p.adjust(p, method = "BH"),
           label = sprintf("rho=%.2f\nR2=%.2f\npadj=%.1e", rho, r2, padj))
}
lab_e <- cor_lab(pooled$score_conv, pooled$rank_diff)
p_e <- ggplot(pooled, aes(x = score_conv, y = rank_diff, fill = score_conv)) +
  geom_point(shape = 21, size = 2, alpha = 0.7, stroke = 0) +
  scale_fill_gradientn(colours = c("white", "#4cc9f0", "#4361ee")) +
  geom_text(data = lab_e, aes(x = Inf, y = Inf, label = label), hjust = 1.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  theme_vaxgo() +
  labs(x = "score_conv", y = "rank_diff (convergence)", fill = "Score",
       title = "e  Convergence score")

# ============================== (f) rank_mouse vs rank_human ==============================
p_f <- ggplot(pooled, aes(x = rank_mouse, y = rank_human, fill = score_conv)) +
  geom_point(data = . %>% filter(score_conv < 90), shape = 21, size = 2, alpha = 0.2, stroke = 0) +
  geom_point(data = . %>% filter(score_conv >= 90), shape = 21, size = 2, alpha = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_gradientn(colours = c("white", "#4cc9f0", "#4361ee")) +
  theme_vaxgo() +
  labs(x = "Rank (Mouse)", y = "Rank (Human)", fill = "Score",
       title = "f  Dual-species rank")

# ============================== (g) heatmap of top 20 ==============================
priority <- read_csv(file.path(tab_dir, "priority_relevant_convergent_top20.csv"), show_col_types = FALSE)
gal <- readRDS(here("tables", "human_mouse_statsmodelling_gene_annotated_layers.rds")) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

heat_dat <- priority %>%
  left_join(gal, by = "hgnc_symbol") %>%
  transmute(
    hgnc_symbol,
    `Rank mouse` = rank_mouse,
    `Rank human` = rank_human,
    `Rank diff`  = rank_diff,
    `Dual rank`  = dual_rank,
    `conv_rel`   = conv_rel,
    `score_shared` = score_shared,
    `score_conv`   = score_conv,
    `%Protein identity` = identity_human2mouse,
    `Kimura` = dist_k80,
    `nTotal TFs` = n_tf_total,
    `%Shared TFs` = pct_tf_shared,
    `nTotal cCREs` = n_total_cres_gene,
    `CTCF match` = `pct_match_ctcf_dELS_CTCF-bound`
  )

heat_long <- heat_dat %>%
  pivot_longer(-hgnc_symbol, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup() %>%
  mutate(hgnc_symbol = factor(hgnc_symbol, levels = rev(priority$hgnc_symbol)),
         feature = factor(feature, levels = colnames(heat_dat)[-1]))

p_g <- ggplot(heat_long, aes(x = feature, y = hgnc_symbol, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("#4361ee", "white", "#e5383b"), limits = c(-2.5, 2.5), oob = scales::squish) +
  theme_vaxgo() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7)) +
  labs(x = "", y = "", fill = "z-score",
       title = "g  Prioritized relevant-convergent genes")

# ============================== assemble ==============================
fig7 <- (p_b | p_c) / (p_d | p_e) / (p_f | p_g) +
  plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(title = "Fig. 7 | Multilayer modelling of cross-species translatability")

ggsave(file.path(fig_dir, "Fig7_multilayer_modelling.png"), fig7, width = 13, height = 15, dpi = 300)
ggsave(file.path(fig_dir, "Fig7_panel_b_classification.png"), p_b, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Fig7_panel_c_regression.png"), p_c, width = 6, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "Fig7_panel_g_heatmap.png"), p_g, width = 10, height = 6, dpi = 300)

fig7
