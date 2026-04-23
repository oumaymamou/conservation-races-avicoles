# =============================================================================
# Script 01 : Initialisation de la population de base
# =============================================================================

library(MoBPS)
library(ggplot2)
library(dplyr)

cat("=== Initialisation de la population de base ===\n")

N_SNP     <- 1000
N_MALES   <- 20
N_FEMALES <- 80
N_TOTAL   <- N_MALES + N_FEMALES
N_QTL_ADD <- 30

cat(sprintf("Population : %d mâles + %d femelles = %d individus\n",
            N_MALES, N_FEMALES, N_TOTAL))

set.seed(42)

population_base <- creating.diploid(
  nsnp       = N_SNP,
  nindi      = N_TOTAL,
  n.additive = N_QTL_ADD,
  verbose    = FALSE
)

cat(sprintf(" Population créée : %d individus\n",
            get.nindi(population_base)))

cat("\n--- 5 générations de conservation ---\n")
population <- population_base

for (gen in 1:5) {
  population <- breeding.diploid(
    population,
    breeding.size  = N_TOTAL,
    selection.size = c(N_MALES, N_FEMALES),
    verbose        = FALSE
  )
  cat(sprintf("  Génération %d simulée\n", gen))
}

cat("\n--- Calcul de la diversité ---\n")
last_gen    <- length(population$breeding)
geno_matrix <- get.geno(population, gen = last_gen)

H_obs_avant     <- mean(apply(geno_matrix, 1, function(x) mean(x == 1)))
n_poly_avant    <- sum(apply(geno_matrix, 1, function(x) length(unique(x)) > 1))
prop_poly_avant <- n_poly_avant / N_SNP

kinship_matrix <- kinship.emp(population, gen = last_gen)
kinship_vals   <- kinship_matrix[upper.tri(kinship_matrix)]
F_mean_avant   <- mean(kinship_vals, na.rm = TRUE)
F_sd_avant     <- sd(kinship_vals,   na.rm = TRUE)

cat(sprintf("  H_obs     : %.4f\n", H_obs_avant))
cat(sprintf("  Kinship   : %.4f ± %.4f\n", F_mean_avant, F_sd_avant))
cat(sprintf("  Loci poly : %d / %d (%.1f%%)\n",
            n_poly_avant, N_SNP, prop_poly_avant * 100))

saveRDS(population,     "data/population_avant_catastrophe.rds")
saveRDS(kinship_matrix, "data/kinship_avant.rds")

diversity_avant <- data.frame(
  phase        = "avant_catastrophe",
  H_obs        = H_obs_avant,
  kinship_mean = F_mean_avant,
  kinship_sd   = F_sd_avant,
  prop_poly    = prop_poly_avant,
  n_poly       = n_poly_avant
)
saveRDS(diversity_avant, "data/diversity_avant.rds")

cat("\n Données sauvegardées dans data/\n")
cat("=== Script 01 terminé ===\n")
cat("Lancer ensuite : source('R/02_simulation_scenarios.R')\n")