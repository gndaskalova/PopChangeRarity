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

__Daskalova, G.N., Myers-Smith, I.H. & Godlee, J.L. (2018) <a href="https://www.biorxiv.org/content/10.1101/272898v3" target="_blank">Rarity and conservation status do not predict vertebrate population trends.</a> bioRxiv. https://doi.org/10.1101/272898 __

# Workflow
For a detailed explanation of the methods, please see the preprint.

We analyzed vertebrate population time-series from the Living Planet Database (133,092 records) covering the period between 1970 and 2014. These time-series represent repeated monitoring surveys of the number of individuals in a given area (species’ abundance over time), hereafter called “populations”. We focus on two aspects of population change – overall changes in abundance over time (population trend) and abundance variability over time (population fluctuations). 

In the first stage of our analyses, we quantified trends and fluctuations for each population using state-space models that account for observation error and random fluctuations36 (Figure S1). 

<p align="center">
  <img src="/img/SI_workflow1.png" width=800 align="middle">
</p> 
_We analyzed vertebrate population time-series from the Living Planet Database (133,092 records) covering the period between 1970 and 2014. These time-series represent repeated monitoring surveys of the number of individuals in a given area (species’ abundance over time), to which we refer as “populations”. Diagram shows one sample population of Red squirrel (Sciurus vulgaris). We quantified two aspects of population change – overall changes in abundance over time (population trend) and abundance variability over time (population fluctuations). We used state-space models that account for observation error and random fluctuations1. See methods for additional details._

In the second stage, we modelled the trend and fluctuation estimates from the first stage across latitude, realm, biome, taxa, rarity metrics, phylogenetic relatedness, species’ conservation status and threat type using a Bayesian modelling framework.

<p align="center">
  <img src="/img/SI_workflow2.png" width=800 align="middle">
</p>
_We modelled the trend and fluctuation estimates from the first stage (Figure S1) across latitude, realm, biome, taxa, rarity metrics, phylogenetic relatedness, species’ conservation status and threat type using a Bayesian modelling framework2. Each model included a species random intercept effect to account for the possible correlation between the trends of populations from the same species. The prior structure (weakly informative priors) was identical across all models except the phylogeny models from the taxonomic patterns section, where the prior structure allowed for an additional phylogeny random effect. See methods for additional details._

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
Following Fournier et al. 2019, we tested the time-series we analyzed for site-selection bias. Removing the first five survey points reduces the bias stemming from starting population surveys at points when individual density is high, whereas removing the last five years reduces the bias of starting surveys when species are very rare. The distribution of population trend values across time-series was not sensitive to the omission of the first five (left-truncation) or the last five years (right-truncation) of population records (Figure S8). Overall, our sensitivity analyses confirmed that our findings were robust to the potential confounding effects of differences in monitoring duration, sampling method and site-selection.

Reference:
Fournier, A. M. V., White, E. R. & Heard, S. B. <a href="https://peerj.com/preprints/27507/" target="_blank">Site-selection bias can drive apparent population declines in long-term studies.</a> doi:10.7287/peerj.preprints.27507v1

# Requirements

### Software
_R version 3.5.1 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .
If you would like to reproduce our analyses but are not familiar with using `R`, you can find an introduction to `R` and running code from `RStudio` on the <a href="https://ourcodingclub.github.io/2016/11/13/intro-to-r.html" target=_blank">Coding Club website</a>.

### Packages
The following `R` packages are required for our analyses. If they are not already installed on your computer, you can install them using the function `install.packages("package-name")`, where `package-name` is the name of the specific package, e.g. `install.packages("readr")`. The specific versions of the packages we used are outlined in <a href="https://github.com/gndaskalova/PopChangeRarity/blob/master/package-versions.md" target="_blank">`package-versions.md`</a>. Note that installing a package using `install.packages()` will automatically install the latest version of the package, which might be different from what we used for our analyses. For code that shows how to install older versions of packages, check out this <a href="https://stackoverflow.com/questions/17082341/installing-older-version-of-r-package" target="_blank">thread on Stackoverflow</a>. 

`readr, data.table, tidyr, dplyr, ggplot2, ggExtra, ggthemes, viridis, png, mapdata, maps, gridExtra, broom, MCMCglmm, stargazer, diptest, plotrix, scales, rredlist, stringr, corrplot, ggtree, proj4, ggalt, RColorBrewer, ggridges, ape, forcats, CoordinateCleaner, geosphere, parallel, scrubr`
