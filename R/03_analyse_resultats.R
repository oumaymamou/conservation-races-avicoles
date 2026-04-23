# =============================================================================
# Script 03 : Analyse statistique des résultats et visualisations
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
# Correction du taux de récupération
H_ref <- readRDS("data/diversity_avant.rds")$H_obs
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

saveRDS(results_df, "data/results_df.rds")
saveRDS(summary_df, "data/summary_df.rds")

cat("=== Analyse statistique des résultats ===\n\n")

# -----------------------------------------------------------------------------
# 1. Chargement des données
# -----------------------------------------------------------------------------
results_df  <- readRDS("data/results_df.rds")
summary_df  <- readRDS("data/summary_df.rds")
div_avant   <- readRDS("data/diversity_avant.rds")

H_ref       <- div_avant$H_obs
poly_ref    <- div_avant$prop_poly

# Palette de couleurs cohérente pour les scénarios
palette_scenarios <- c(
  "S0" = "#d73027",   # rouge — perte totale
  "S1" = "#fc8d59",   # orange — minimal
  "S2" = "#fee090",   # jaune — intermédiaire
  "S3" = "#91bfdb",   # bleu clair — extensif
  "S4" = "#4575b4"    # bleu foncé — optimisé
)

# Labels pour les figures
labels_sc <- c(
  S0 = "S0 — Aucune\nsauvegarde",
  S1 = "S1 — Minimal\n(5 mâles)",
  S2 = "S2 — Intermédiaire\n(15 mâles)",
  S3 = "S3 — Extensif\n(30 mâles)",
  S4 = "S4 — Optimisé\n(15 mâles, OGC)"
)

# Ordonner les scénarios
results_df$scenario_id <- factor(results_df$scenario_id,
                                  levels = c("S0","S1","S2","S3","S4"))
summary_df$scenario_id <- factor(summary_df$scenario_id,
                                  levels = c("S0","S1","S2","S3","S4"))

# -----------------------------------------------------------------------------
# 2. Figure 1 : Hétérozygotie observée — boxplot comparatif
# -----------------------------------------------------------------------------
fig1 <- ggplot(results_df, aes(x = scenario_id, y = H_obs, fill = scenario_id)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  geom_hline(yintercept = H_ref, linetype = "dashed", color = "black",
             linewidth = 0.8, alpha = 0.7) +
  annotate("text", x = 0.6, y = H_ref + 0.005,
           label = sprintf("Référence pré-catastrophe\n(H_obs = %.3f)", H_ref),
           hjust = 0, size = 3, color = "black") +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  labs(
    title    = "Hétérozygotie observée après reconstruction",
    subtitle = "Comparaison des scénarios de sauvegarde vs. référence pré-catastrophe",
    x        = "Scénario de sauvegarde",
    y        = "Hétérozygotie observée (H_obs)",
    caption  = sprintf("n = %d répétitions par scénario. Ligne pointillée = référence avant épidémie.",
                       max(results_df$rep))
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    axis.text.x   = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave("figures/fig1_heterozygotie_boxplot.png", fig1, width = 9, height = 6, dpi = 150)
cat("✓ Figure 1 sauvegardée : fig1_heterozygotie_boxplot.png\n")

# -----------------------------------------------------------------------------
# 3. Figure 2 : Taux de récupération de la diversité (%)
# -----------------------------------------------------------------------------
# Figure 2
summary_rep <- results_df %>%
  group_by(scenario_id) %>%
  summarise(
    mean_rec = mean(recovery_H_obs, na.rm = TRUE),
    sd_rec   = sd(recovery_H_obs,   na.rm = TRUE)
  )

fig2 <- ggplot(summary_df, aes(x = scenario_id, y = recovery_H_mean, fill = scenario_id)) +
  geom_col(alpha = 0.85, color = "white") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = sprintf("%.1f%%", recovery_H_mean),
                y = recovery_H_mean + 2),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 20)) +
  labs(
    title   = "Taux de récupération de l'hétérozygotie",
    subtitle = "% de la diversité pré-catastrophe récupérée après reconstruction",
    x = "Scénario de sauvegarde",
    y = "Taux de récupération (%)",
    caption = "100% = récupération complète"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave("figures/fig2_taux_recuperation.png", fig2, width = 9, height = 6, dpi = 150)
cat("✓ Figure 2 sauvegardée : fig2_taux_recuperation.png\n")

# -----------------------------------------------------------------------------
# 4. Figure 3 : Relation entre nombre de mâles sauvegardés et diversité récupérée
# -----------------------------------------------------------------------------
fig3 <- ggplot(summary_df %>% filter(n_males > 0),
               aes(x = n_males, y = recovery_H_mean, color = scenario_id)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "gray50", linetype = "dotted", linewidth = 0.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40", alpha = 0.6) +
  geom_text(aes(label = scenario_id), nudge_y = 2.5, size = 3.5, fontface = "bold") +
  scale_color_manual(values = palette_scenarios, guide = "none") +
  scale_x_continuous(breaks = c(5, 15, 30)) +
  labs(
    title    = "Effet du nombre de mâles sauvegardés sur la récupération génétique",
    subtitle = "S4 = même effectif que S2 mais sélection génétiquement optimisée",
    x        = "Nombre de mâles dont la semence est cryoconservée",
    y        = "Taux de récupération de H_obs (%)",
    caption  = "Droite = tendance linéaire (intervalle de confiance 95%)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("figures/fig3_nmales_vs_recuperation.png", fig3, width = 8, height = 6, dpi = 150)
cat("✓ Figure 3 sauvegardée : fig3_nmales_vs_recuperation.png\n")

# -----------------------------------------------------------------------------
# 5. Figure 4 : Consanguinité moyenne par scénario
# -----------------------------------------------------------------------------
fig4 <- ggplot(results_df, aes(x = scenario_id, y = kinship_mean, fill = scenario_id)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 1.5) +
  scale_fill_manual(values = palette_scenarios, guide = "none") +
  scale_x_discrete(labels = labels_sc) +
  labs(
    title    = "Consanguinité moyenne (kinship) après reconstruction",
    subtitle = "Une kinship élevée indique une perte de diversité génétique",
    x        = "Scénario",
    y        = "Kinship moyenne entre individus",
    caption  = "Violin plot + boxplot. Kinship élevée = effectifs reproducteurs réduits."
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave("figures/fig4_consanguinite_violin.png", fig4, width = 9, height = 6, dpi = 150)
cat("✓ Figure 4 sauvegardée : fig4_consanguinite_violin.png\n")

# -----------------------------------------------------------------------------
# 6. Analyse statistique formelle (ANOVA + tests post-hoc)
# -----------------------------------------------------------------------------
cat("\n--- Tests statistiques ---\n")

# ANOVA à un facteur sur H_obs
model_anova <- aov(H_obs ~ scenario_id, data = results_df)
cat("\nANOVA — H_obs ~ scénario :\n")
print(summary(model_anova))

# Test de Tukey pour comparaisons par paires
tukey_res <- TukeyHSD(model_anova)
cat("\nTest de Tukey (HSD) — comparaisons par paires :\n")
print(tukey_res$scenario_id)

# Sauvegarder les résultats statistiques
sink("outputs/tests_statistiques.txt")
cat("=== ANOVA — H_obs ~ scénario ===\n")
print(summary(model_anova))
cat("\n=== Test de Tukey (HSD) ===\n")
print(tukey_res$scenario_id)
sink()
cat("\n✓ Résultats statistiques sauvegardés dans outputs/tests_statistiques.txt\n")

# -----------------------------------------------------------------------------
# 7. Tableau de synthèse final
# -----------------------------------------------------------------------------
table_finale <- summary_df %>%
  mutate(
    `H_obs (moy ± sd)` = sprintf("%.4f ± %.4f", H_obs_mean, H_obs_sd),
    `Kinship moy.`     = sprintf("%.4f", kinship_mean_mean),
    `% Poly.`          = sprintf("%.1f%%", prop_poly_mean * 100),
    `Récupération H`   = sprintf("%.1f%%", recovery_H_mean)
  ) %>%
  select(scenario, n_males, description, `H_obs (moy ± sd)`,
         `Kinship moy.`, `% Poly.`, `Récupération H`)

cat("\n=== Tableau de synthèse ===\n")
print(table_finale, n = Inf)

write.csv(table_finale, "outputs/tableau_synthese.csv", row.names = FALSE)

cat("\n✓ Tableau de synthèse sauvegardé : outputs/tableau_synthese.csv\n")
cat("=== Script 03 terminé ===\n")
cat(" Compiler maintenant le rapport : rmarkdown::render('rapport/rapport_conservation.Rmd')\n")
