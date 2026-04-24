# =============================================================================
# Script 03 : Analyse statistique des résultats et visualisations
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Chargement et préparation des données ------------------------------------
results_df  <- readRDS("data/results_df.rds")
div_avant   <- readRDS("data/diversity_avant.rds")

H_ref    <- div_avant$H_obs
poly_ref <- div_avant$prop_poly

# Correction/Recalcul du taux de récupération
results_df$recovery_H_obs <- results_df$H_obs / H_ref * 100

summary_df <- results_df %>%
  group_by(scenario_id, scenario, n_males, description) %>%
  summarise(
    H_obs_mean         = mean(H_obs),
    H_obs_sd           = sd(H_obs),
    kinship_mean_mean  = mean(kinship_mean),
    prop_poly_mean     = mean(prop_poly),
    recovery_H_mean    = mean(recovery_H_obs, na.rm = TRUE),
    recovery_poly_mean = mean(recovery_poly, na.rm = TRUE),
    .groups = "drop"
  )

# Ordonner les facteurs pour les graphiques
results_df$scenario_id <- factor(results_df$scenario_id, levels = c("S0","S1","S2","S3","S4"))
summary_df$scenario_id <- factor(summary_df$scenario_id, levels = c("S0","S1","S2","S3","S4"))

# Palette et labels
palette_scenarios <- c("S0"="#d73027", "S1"="#fc8d59", "S2"="#fee090", "S3"="#91bfdb", "S4"="#4575b4")
labels_sc <- c(S0="S0 — Témoin", S1="S1 — 5 mâles", S2="S2 — 15 mâles", S3="S3 — 30 mâles", S4="S4 — Optimisé")

# 2. Figure 1 : Hétérozygotie observée ----------------------------------------
fig1 <- ggplot(results_df, aes(x = scenario_id, y = H_obs, fill = scenario_id)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21) +
  geom_hline(yintercept = H_ref, linetype = "dashed", alpha = 0.7) +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  labs(title = "Hétérozygotie observée après reconstruction",
       x = "Scénario", y = "H_obs") +
  theme_bw()

ggsave("outputs/fig1_heterozygotie.png", fig1, width = 9, height = 6)

# 3. Figure 2 : Taux de récupération ------------------------------------------
fig2 <- ggplot(summary_df, aes(x = scenario_id, y = recovery_H_mean, fill = scenario_id)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", recovery_H_mean), y = recovery_H_mean + 2), size = 3.5) +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  labs(title = "Taux de récupération de la diversité",
       x = "Scénario", y = "Récupération (%)") +
  theme_bw()

ggsave("outputs/fig2_recuperation.png", fig2, width = 9, height = 6)

# 4. Figure 3 : Relation mâles/diversité --------------------------------------
fig3 <- ggplot(summary_df %>% filter(n_males > 0), 
               aes(x = n_males, y = recovery_H_mean, color = scenario_id)) +
  geom_point(size = 4) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, linetype = "dotted") +
  scale_color_manual(values = palette_scenarios, guide = "none") +
  labs(title = "Effet du nombre de mâles sur la récupération",
       x = "Nombre de mâles cryoconservés", y = "Récupération H_obs (%)") +
  theme_bw()

ggsave("outputs/fig3_nmales_vs_diversite.png", fig3, width = 8, height = 6)

# 5. Figure 4 : Consanguinité (Kinship) ---------------------------------------
fig4 <- ggplot(results_df, aes(x = scenario_id, y = kinship_mean, fill = scenario_id)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  labs(title = "Consanguinité moyenne par scénario",
       x = "Scénario", y = "Kinship") +
  theme_bw()

ggsave("outputs/fig4_kinship.png", fig4, width = 9, height = 6)

# 6. Analyse statistique ------------------------------------------------------
model_anova <- aov(H_obs ~ scenario_id, data = results_df)
tukey_res   <- TukeyHSD(model_anova)

# Export des tests
sink("outputs/tests_statistiques.txt")
print(summary(model_anova))
print(tukey_res)
sink()

# 7. Tableau de synthèse ------------------------------------------------------
table_finale <- summary_df %>%
  mutate(
    `H_obs (moy ± sd)` = sprintf("%.4f ± %.4f", H_obs_mean, H_obs_sd),
    `Récupération H`   = sprintf("%.1f%%", recovery_H_mean)
  ) %>%
  select(scenario, n_males, description, `H_obs (moy ± sd)`, `Récupération H`)

write.csv(table_finale, "outputs/tableau_synthese.csv", row.names = FALSE)
saveRDS(results_df, "data/results_df_final.rds")
saveRDS(summary_df, "data/summary_df_final.rds")
