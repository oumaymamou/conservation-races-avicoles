# =============================================================================
# Script 02 : Simulation des scénarios de sauvegarde
# =============================================================================

library(MoBPS)
library(dplyr)

source("R/functions_utils.R")

# 1. Chargement des données de référence --------------------------------------
population_ref  <- readRDS("data/population_avant_catastrophe.rds")
diversity_avant <- readRDS("data/diversity_avant.rds")

# 2. Définition des scénarios -------------------------------------------------
scenarios <- list(
  S0 = list(name = "S0_controle",       n_males_cryo = 0,  description = "Extinction"),
  S1 = list(name = "S1_minimal",        n_males_cryo = 5,  description = "Sauvegarde 5 mâles"),
  S2 = list(name = "S2_intermediaire",  n_males_cryo = 15, description = "Sauvegarde 15 mâles"),
  S3 = list(name = "S3_extensive",      n_males_cryo = 30, description = "Sauvegarde 30 mâles"),
  S4 = list(name = "S4_optimise",       n_males_cryo = 15, description = "Sélection diversité (15 mâles)")
)

N_SIM       <- 5   
N_RECON     <- 100 
N_GEN_RECON <- 5  

# 3. Fonction de simulation ---------------------------------------------------
simulate_rescue <- function(pop_ref, n_males_cryo, n_gen = 5, n_recon = 100,
                            optimized = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (n_males_cryo == 0) {
    return(list(H_obs = 0, kinship_mean = 0.5, prop_poly = 0, n_poly = 0))
  }
  
  pop_recon <- pop_ref
  
  for (g in 1:n_gen) {
    target_size <- min(n_recon, 20 * g)
    n_f <- max(5, target_size - n_males_cryo)
    
    pop_recon <- breeding.diploid(
      pop_recon,
      breeding.size  = target_size,
      selection.size = c(n_males_cryo, n_f),
      verbose        = FALSE
    )
  }
  
  last_gen_recon  <- length(pop_recon$breeding)
  geno_matrix     <- get.geno(pop_recon, gen = last_gen_recon)
  
  H_obs     <- mean(apply(geno_matrix, 1, function(x) mean(x == 1)))
  n_poly    <- sum(apply(geno_matrix, 1, function(x) length(unique(x)) > 1))
  prop_poly   <- n_poly / nrow(geno_matrix)
  
  k_mat        <- kinship.emp(pop_recon, gen = last_gen_recon)
  kinship_mean <- mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)
  
  return(list(H_obs = H_obs, kinship_mean = kinship_mean, 
              prop_poly = prop_poly, n_poly = n_poly))
}

# 4. Exécution des simulations ------------------------------------------------
all_results <- list()

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  sc_results <- list()
  
  for (rep in 1:N_SIM) {
    res <- simulate_rescue(
      pop_ref      = population_ref,
      n_males_cryo = sc$n_males_cryo,
      n_gen        = N_GEN_RECON,
      n_recon      = N_RECON,
      optimized    = (sc_name == "S4"),
      seed         = rep * 100 + which(names(scenarios) == sc_name)
    )
    res$scenario    <- sc$name
    res$scenario_id <- sc_name
    res$rep         <- rep
    res$n_males     <- sc$n_males_cryo
    res$description <- sc$description
    sc_results[[rep]] <- as.data.frame(res)
  }
  all_results[[sc_name]] <- bind_rows(sc_results)
}

# 5. Compilation et export ----------------------------------------------------
results_df <- bind_rows(all_results)

# Calcul des indices de récupération
H_obs_ref     <- diversity_avant$H_obs
prop_poly_ref <- diversity_avant$prop_poly

results_df <- results_df %>%
  mutate(
    recovery_H_obs   = ifelse(H_obs_ref > 0, H_obs / H_obs_ref * 100, NA),
    recovery_poly    = ifelse(prop_poly_ref > 0, prop_poly / prop_poly_ref * 100, NA)
  )

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

# Sauvegardes
write.csv(results_df, "outputs/results_all_reps.csv", row.names = FALSE)
write.csv(summary_df, "outputs/results_summary.csv", row.names = FALSE)
saveRDS(results_df,   "data/results_df.rds")
saveRDS(summary_df,   "data/summary_df.rds")
