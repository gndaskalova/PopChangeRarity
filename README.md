# PopChangeRarity

_This repository contains the code for data extraction, data integration and statistical analyses in "All is not decline across global vertebrate populations"._

Authors: Gergana N. Daskalova, Isla H. Myers-Smith, John L. Godlee

### Disclaimer
__The aim of this repository is to give reviewers of the manuscript private access to the code for all analyses in the manuscript. Copyright for the code remains with the authors. Please do not distribute or share the code with anyone.__

#### Contact: Gergana Daskalova gndaskalova@gmail.com

# Aim
_Our aim was to determine how the trends and fluctuations of vertebrate populations vary with biogeography, taxa, phylogenetic relationships and across species’ rarity metrics and IUCN conservation and threat categories._

# Preprint
The manuscript is available as a preprint on bioRxiv.

__Daskalova, G.N., Myers-Smith, I.H. & Godlee, J.L. (2018) Rarity and conservation status do not predict vertebrate population trends. bioRxiv. https://doi.org/10.1101/272898 __

# Workflow
For a detailed explanation of the methods, please see the preprint.

We analyzed vertebrate population time-series from the Living Planet Database (133,092 records) covering the period between 1970 and 2014. These time-series represent repeated monitoring surveys of the number of individuals in a given area (species’ abundance over time), hereafter called “populations”. We focus on two aspects of population change – overall changes in abundance over time (population trend) and abundance variability over time (population fluctuations). 

In the first stage of our analyses, we quantified trends and fluctuations for each population using state-space models that account for observation error and random fluctuations36 (Figure S1). 

<p align="center">
  <img src="/img/Si_workflow1.png" width=800 align="middle">
</p> 

In the second stage, we modelled the trend and fluctuation estimates from the first stage across latitude, realm, biome, taxa, rarity metrics, phylogenetic relatedness, species’ conservation status and threat type using a Bayesian modelling framework (Figure S2).

<p align="center">
  <img src="/img/SI_workflow2.png" width=800 align="middle">
</p> 

# Data

### Living Planet Database
The population time-series analysed here come from the Living Planet Database, with the raw data publicly available from http://www.livingplanetindex.org

### Key data objects within the scripts and what they contain

__`global_mus_scaled.csv`__
Population change (mu values for overall trend), taxa and biome classification and coordinates of the Living Planet Database time series. Calculated using state-space models in `06-population-models.R`.



# Scripts

_Please note that the majority of the code requires high computing power to successfully run._

### Calculate habitat specificity and geographic range

__`01-calculate-habspec.R`__
Calculate habitat specificity for each species based on the number of habitats they occupy, as per their IUCN species profile.

__`02-calculate-range.R`__
Calculate geographic range size using GBIF occurrence data for species monitored in the UK.

### Quantify population change (overall trends in population abundance over time and fluctuations over time)

__`03-derive-pop-change-state-space.R`__
Calculate population trends (mu values for overall trend) and population fluctuations (sigma values for process error) across the duration of each Living Planet Database time-series. The code uses state-space models.

__`04-derive-pop-change-linear.R`__
Calculate population trends (slope values for overall trend) across the duration of each Living Planet Database time-series. The code uses general linear models of abundance versus time.

### Test biogeographic, taxonomic, phylogenetic, rarity, conservation status and threat patterns in population trends and fluctuations

__`05-manuscript-analysis.R`__
Format data and run all statistical models. The script includes

### Visualise findings

__`06-figures.R`__
Generate figures included in the manuscript maintext and supplementary information.

### Sensitivity analyses

__`07-sensitivity-left-truncation.R`__


__`08-sensitivity-right-truncation.R`__

# Requirements

### Software
R version 3.5.1 or greater

### Packages
library(readr)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(ggthemes)
library(viridis)
library(png)
library(mapdata)
library(maps)
library(gridExtra)
library(broom)
library(MCMCglmm)
library(stargazer)
library(diptest)
library(plotrix)
library(scales)
library(rredlist)
library(stringr)
library(corrplot)
library(ggtree)
library(proj4)
library(ggalt)
library(RColorBrewer)
library(ggridges)
library(ape)
library(forcats)

