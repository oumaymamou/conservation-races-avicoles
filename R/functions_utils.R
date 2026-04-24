# =============================================================================
# functions_utils.R : Fonctions de calcul de diversité génétique
# =============================================================================

# Calcule l'hétérozygotie observée
calc_H_obs <- function(pop, gen = NULL) {
  if (is.null(gen)) gen <- length(pop$breeding)
  geno <- get.geno(pop, gen = gen)
  mean(apply(geno, 1, function(x) mean(x == 1)))
}

# Calcule la kinship moyenne (hors diagonale)
calc_kinship_mean <- function(k_mat) {
  mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)
}

# Calcule la proportion de loci polymorphes
calc_prop_poly <- function(pop, gen = NULL) {
  if (is.null(gen)) gen <- length(pop$breeding)
  geno <- get.geno(pop, gen = gen)
  n_poly <- sum(apply(geno, 1, function(x) length(unique(x)) > 1))
  return(n_poly / nrow(geno))
}

# Résumé complet de la diversité pour une génération donnée
get_diversity_summary <- function(pop, gen = NULL, label = "population") {
  if (is.null(gen)) gen <- length(pop$breeding)

  geno      <- get.geno(pop, gen = gen)
  H_obs     <- mean(apply(geno, 1, function(x) mean(x == 1)))
  n_poly    <- sum(apply(geno, 1, function(x) length(unique(x)) > 1))
  prop_poly <- n_poly / nrow(geno)
  
  k_mat  <- kinship.emp(pop, gen = gen)
  k_mean <- mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)

  data.frame(
    label        = label,
    generation   = gen,
    n_individus  = get.nindi(pop, gen = gen),
    H_obs        = round(H_obs, 5),
    kinship_mean = round(k_mean, 5),
    prop_poly    = round(prop_poly, 4),
    n_poly       = n_poly,
    stringsAsFactors = FALSE
  )
}
