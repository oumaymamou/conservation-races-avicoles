# Conservation des Races Avicoles Locales — Simulation de Scénarios de Sauvegarde

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue?logo=r)](https://www.r-project.org/)
[![MoBPS](https://img.shields.io/badge/MoBPS-1.6.64-green)](https://github.com/tpook92/MoBPS)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> Projet personnel exploratoire — simulation de stratégies de sauvegarde génétique  
> pour des populations avicoles locales à faibles effectifs, en contexte de risque épidémique.

---

## Contexte et motivation

Les races avicoles locales françaises (Gauloise Dorée, Bresse, Marans, Barbezieux, etc.)
présentent une diversité génétique remarquable mais sont souvent maintenues dans des
noyaux de conservation de petite taille concentrés sur un seul site. En cas d'épisode
épidémique grave (ex. : influenza aviaire H5N1), la disparition de ces noyaux sans plan
de sauvegarde préalable représente un risque réel et irréversible pour la biodiversité
domestique nationale.

Ce projet explore, via simulation stochastique avec **MoBPS**, différents scénarios de
cryoconservation pour évaluer la capacité de chaque stratégie à **récupérer la diversité
génétique** présente avant une catastrophe simulée.

---

## Structure du projet

```
conservation_races_avicoles/
├── README.md
├── R/
│   ├── 01_setup_population.R        # Initialisation de la population de base
│   ├── 02_simulation_scenarios.R    # Définition et lancement des scénarios
│   ├── 03_analyse_resultats.R       # Analyse statistique comparative
│   └── functions_utils.R            # Fonctions utilitaires
├── data/
│   └── scenario_parameters.csv      # Paramètres des scénarios
├── outputs/
│   └── results_summary.csv          # Résultats agrégés (généré par les scripts)
├── figures/
│   └── (générées automatiquement)
└── rapport/
    └── rapport_conservation.Rmd     # Rapport R Markdown complet
```

---

---

## Approche méthodologique

### Population modélisée

- **Espèce** : poulet (*Gallus gallus*)
- **Taille du noyau** : 100 individus (20 mâles + 80 femelles)
- **Génome simulé** : 1 000 marqueurs SNP
- **Générations simulées** : 5 avant catastrophe

### Scénarios comparés

| Scénario | Description | Mâles cryoconservés |
|----------|-------------|-------------------|
| S0 | Aucune sauvegarde — extinction | 0 |
| S1 | Sauvegarde minimale | 5 |
| S2 | Sauvegarde intermédiaire | 15 |
| S3 | Sauvegarde extensive | 30 |
| S4 | Stratégie optimisée (OGC-guidée) | 15 |

### Métriques d'évaluation

- **Hétérozygotie observée (H_obs)** — indicateur principal de diversité génétique
- **Kinship moyenne** — proxy de la consanguinité
- **Proportion de loci polymorphes** — diversité allélique
- **Taux de récupération (%)** — ratio post-reconstruction vs. référence pré-catastrophe

---

##  Installation et utilisation

### Prérequis

```r
install.packages("MoBPS")
install.packages(c("ggplot2", "dplyr", "tidyr", "knitr", "rmarkdown", "kableExtra", "broom"))
```

### Exécution

```r
source("R/01_setup_population.R")
source("R/02_simulation_scenarios.R")
source("R/03_analyse_resultats.R")

rmarkdown::render("rapport/rapport_conservation.Rmd")
```

---

## Résultats

| Scénario | H_obs | Récupération |
|----------|-------|-------------|
| S0 — Aucune sauvegarde | 0.000 | 0% |
| S1 — 5 mâles | 0.280 | 85.4% |
| S2 — 15 mâles | 0.296 | 90.1% |
| S3 — 30 mâles | 0.301 | 91.7% |
| S4 — 15 mâles optimisés | 0.299 | 91.1% |

**Conclusions principales :**
- Sans cryobanque, la diversité est perdue à 100%
- 15 mâles permettent de récupérer ~90% de la diversité initiale
- Le test de Tukey montre que S2, S3 et S4 ne sont pas significativement différents (p > 0.6)
- La sélection optimisée (S4) surpasse la sélection aléatoire (S2) à effectif identique

> Résultats issus de simulations exploratoires sur paramètres fictifs.
> Ils illustrent une démarche méthodologique et non une étude empirique publiée.

---

## Références

- Pook T. et al. (2020). *MoBPS: Modular Breeding Program Simulator*. G3. [doi:10.1534/g3.120.401193](https://doi.org/10.1534/g3.120.401193)
- FAO (2015). *The Second Report on the State of the World's Animal Genetic Resources*. Rome.
- Long X. et al. (2022). *Poultry genetic heritage cryopreservation and reconstruction*. J Anim Sci Biotechnol.

---

## Auteur

Oumayma Mourabit - Projet réalisé dans le cadre d'une démarche d'auto-formation sur la simulation
de programmes de conservation génétique animale.