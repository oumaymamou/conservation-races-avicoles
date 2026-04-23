# =============================================================================
# functions_utils.R — Fonctions utilitaires pour mon projet sur la conservation
# =============================================================================

#' Calcule l'hétérozygotie observée depuis une matrice de génotypes MoBPS
#' @param pop     objet MoBPS population
#' @param gen     numéro de génération (défaut : dernière)
#' @return H_obs moyenne
calc_H_obs <- function(pop, gen = NULL) {
  if (is.null(gen)) gen <- length(pop$breeding)
  geno <- get.geno(pop, gen = gen)
  mean(apply(geno, 1, function(x) mean(x == 1)))
}

#' Calcule la kinship moyenne (hors diagonale) depuis une matrice de kinship
#' @param k_mat   matrice de kinship (symétrique)
#' @return moyenne des valeurs hors diagonale
calc_kinship_mean <- function(k_mat) {
  mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)
}

#' Calcule le taux de loci polymorphes
#' @param pop     objet MoBPS population
#' @param gen     numéro de génération
#' @return proportion de loci polymorphes
calc_prop_poly <- function(pop, gen = NULL) {
  if (is.null(gen)) gen <- length(pop$breeding)
  geno <- get.geno(pop, gen = gen)
  n_poly <- sum(apply(geno, 1, function(x) length(unique(x)) > 1))
  return(n_poly / nrow(geno))
}

#' Extrait un résumé complet de diversité pour une génération
#' @param pop       objet MoBPS
#' @param gen       numéro de génération
#' @param label     étiquette pour l'output
#' @return data.frame avec les métriques de diversité
get_diversity_summary <- function(pop, gen = NULL, label = "population") {
  if (is.null(gen)) gen <- length(pop$breeding)

  geno     <- get.geno(pop, gen = gen)
  H_obs    <- mean(apply(geno, 1, function(x) mean(x == 1)))
  n_poly   <- sum(apply(geno, 1, function(x) length(unique(x)) > 1))
  prop_poly <- n_poly / nrow(geno)
  k_mat    <- kinship.emp(pop, gen = gen)
  k_mean   <- mean(k_mat[upper.tri(k_mat)], na.rm = TRUE)

  data.frame(
    label        = label,
    generation   = gen,
    n_individus  = get.nindi(pop, gen = gen),
    H_obs        = round(H_obs, 5),
    kinship_mean = round(k_mean, 5),
    prop_poly    = round(prop_poly, 4),
    n_poly       = n_poly
  )
}

#' Formate un nombre en pourcentage pour l'affichage
pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), x * 100)

#' Affiche un message de progression formaté
progress_msg <- function(msg, width = 60) {
  cat(paste0("\n", strrep("─", width), "\n", msg, "\n", strrep("─", width), "\n"))
}
