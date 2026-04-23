# =============================================================================
# Script 02 : Simulation des scénarios de sauvegarde
# =============================================================================
# Logique : On part de la population avant catastrophe, on "simule" la perte
# du noyau femelle, et on reconstruit la population avec différents scénarios
# de sauvegarde (nombre de semences de mâles conservées en cryobanque).
# On compare ensuite la diversité récupérée.
# =============================================================================

library(MoBPS)
library(dplyr)

source("R/functions_utils.R")

cat("=== Simulation des scénarios de sauvegarde ===\n\n")

# -----------------------------------------------------------------------------
# 1. Chargement de la population de référence
# -----------------------------------------------------------------------------
population_ref  <- readRDS("data/population_avant_catastrophe.rds")
diversity_avant <- readRDS("data/diversity_avant.rds")

cat(sprintf("Population de référence chargée : %d individus, génération %d\n",
            get.nindi(population_ref), length(population_ref$breeding)))

# -----------------------------------------------------------------------------
# 2. Définition des scénarios
# -----------------------------------------------------------------------------
# Paramètre clé : nombre de mâles dont la semence a été cryoconservée.
# Ces mâles sont les seuls contributeurs génétiques masculins disponibles
# après l'épidémie. Les femelles sont reconstituées par croisement.

scenarios <- list(
  S0 = list(name = "S0_controle",       n_males_cryo = 0,  description = "Aucune sauvegarde (extinction)"),
  S1 = list(name = "S1_minimal",        n_males_cryo = 5,  description = "Sauvegarde minimale (5 mâles)"),
  S2 = list(name = "S2_intermediaire",  n_males_cryo = 15, description = "Sauvegarde intermédiaire (15 mâles)"),
  S3 = list(name = "S3_extensive",      n_males_cryo = 30, description = "Sauvegarde extensive (30 mâles)"),
  S4 = list(name = "S4_optimise",       n_males_cryo = 15, description = "Sauvegarde optimisée — sélection génétiquement diverse (15 mâles)")
)

N_SIM    <- 5     # Nombre de répétitions par scénario (augmenter pour résultats robustes)
N_RECON  <- 100   # Effectif cible après reconstruction (5 générations)
N_GEN_RECON <- 5  # Générations de reconstruction

cat("Scénarios définis :\n")
for (sc in scenarios) {
  cat(sprintf("  [%s] %s\n", sc$name, sc$description))
}
cat(sprintf("\nRépétitions par scénario : %d\n\n", N_SIM))

# -----------------------------------------------------------------------------
# 3. Fonction principale de simulation d'un scénario
# -----------------------------------------------------------------------------

simulate_rescue <- function(pop_ref, n_males_cryo, n_gen = 5, n_recon = 100,
                            optimized = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Cas S0 : aucune semence disponible
  if (n_males_cryo == 0) {
    return(list(H_obs = 0, kinship_mean = 0.5, prop_poly = 0, n_poly = 0))
  }
  
  # Reconstruction : on repart de la population de référence
  # en limitant l'effectif reproducteur aux mâles sauvegardés
  pop_recon <- pop_ref
  
  for (g in 1:n_gen) {
    target_size <- min(n_recon, 20 * g)
    
    # On contrôle la diversité via selection.size
    # n_males_cryo = nombre de mâles sélectionnés
    n_f <- max(5, target_size - n_males_cryo)
    
    pop_recon <- breeding.diploid(
      pop_recon,
      breeding.size  = target_size,
      selection.size = c(n_males_cryo, n_f),
      verbose        = FALSE
    )
  }
  
  # Métriques finales
  last_gen_recon  <- length(pop_recon$breeding)
  geno_matrix     <- get.geno(pop_recon, gen = last_gen_recon)
  
  H_obs       <- mean(apply(geno_matrix, 1, function(x) mean(x == 1)))
  n_poly      <- sum(apply(geno_matrix, 1, function(x) length(unique(x)) > 1))
  prop_poly   <- n_poly / nrow(geno_matrix)
  
  k_mat        <- kinship.emp(pop_recon, gen = last_gen_recon)
  kinship_mean <- mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)
  
  return(list(H_obs = H_obs, kinship_mean = kinship_mean,
              prop_poly = prop_poly, n_poly = n_poly))
}

# -----------------------------------------------------------------------------
# 4. Lancement des simulations
# -----------------------------------------------------------------------------

all_results <- list()

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  cat(sprintf("--- Scénario %s : %s ---\n", sc$name, sc$description))

  sc_results <- list()
  for (rep in 1:N_SIM) {
    cat(sprintf("  Répétition %d/%d...\r", rep, N_SIM))
    flush.console()

    res <- simulate_rescue(
      pop_ref      = population_ref,
      n_males_cryo = sc$n_males_cryo,
      n_gen        = N_GEN_RECON,
      n_recon      = N_RECON,
      optimized    = (sc_name == "S4"),
      seed         = rep * 100 + which(names(scenarios) == sc_name)
    )
    res$scenario   <- sc$name
    res$scenario_id <- sc_name
    res$rep        <- rep
    res$n_males    <- sc$n_males_cryo
    res$description <- sc$description
    sc_results[[rep]] <- as.data.frame(res)
  }

  all_results[[sc_name]] <- bind_rows(sc_results)
  cat(sprintf("  ✓ Scénario %s terminé. H_obs moyen = %.4f\n",
              sc$name, mean(all_results[[sc_name]]$H_obs)))
}

# -----------------------------------------------------------------------------
# 5. Compilation et sauvegarde des résultats
# -----------------------------------------------------------------------------
results_df <- bind_rows(all_results)

# Ajout des valeurs de référence (avant catastrophe)
H_obs_ref    <- diversity_avant$H_obs
prop_poly_ref <- diversity_avant$prop_poly

# Calcul du taux de récupération
results_df <- results_df %>%
  mutate(
    recovery_H_obs   = ifelse(H_obs_ref > 0, H_obs / H_obs_ref * 100, NA),
    recovery_poly    = ifelse(prop_poly_ref > 0, prop_poly / prop_poly_ref * 100, NA)
  )

# Résumé par scénario
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
  ) %>%
  arrange(n_males)

# Affichage du résumé
cat("\n=== Résumé des résultats ===\n")
print(summary_df %>% select(scenario, n_males, H_obs_mean, H_obs_sd,
                              recovery_H_mean, recovery_poly_mean))

# Sauvegarde
write.csv(results_df,  "outputs/results_all_reps.csv",  row.names = FALSE)
write.csv(summary_df,  "outputs/results_summary.csv",   row.names = FALSE)
saveRDS(results_df,    "data/results_df.rds")
saveRDS(summary_df,    "data/summary_df.rds")

cat("\n✓ Résultats sauvegardés dans outputs/\n")
cat("=== Script 02 terminé ===\n")
cat("Lancer ensuite : source('R/03_analyse_resultats.R')\n")
