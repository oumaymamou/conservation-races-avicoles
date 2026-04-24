# =============================================================================
# Script 01 : Initialisation de la population de base
# Auteur : oumaymamou - Avril 2026
# =============================================================================

library(MoBPS)
library(ggplot2)
library(dplyr)

# Paramètres de simulation
N_SNP      <- 1000
N_MALES    <- 20
N_FEMALES  <- 80
N_TOTAL    <- N_MALES + N_FEMALES
N_QTL_ADD  <- 30

set.seed(42)

# Création de la population fondatrice
population_base <- creating.diploid(
  nsnp       = N_SNP,
  nindi      = N_TOTAL,
  n.additive = N_QTL_ADD,
  verbose    = FALSE
)

# Simulation de 5 générations de conservation (stabilité)
population <- population_base
for (gen in 1:5) {
  population <- breeding.diploid(
    population,
    breeding.size  = N_TOTAL,
    selection.size = c(N_MALES, N_FEMALES),
    verbose        = FALSE
  )
}

# Analyse de la diversité génétique initiale
last_gen    <- length(population$breeding)
geno_matrix <- get.geno(population, gen = last_gen)

H_obs_avant     <- mean(apply(geno_matrix, 1, function(x) mean(x == 1)))
n_poly_avant    <- sum(apply(geno_matrix, 1, function(x) length(unique(x)) > 1))
prop_poly_avant <- n_poly_avant / N_SNP

kinship_matrix <- kinship.emp(population, gen = last_gen)
kinship_vals   <- kinship_matrix[upper.tri(kinship_matrix)]
F_mean_avant   <- mean(kinship_vals, na.rm = TRUE)
F_sd_avant     <- sd(kinship_vals,   na.rm = TRUE)

# Sauvegarde des objets de données
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
