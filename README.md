# PopChangeRarity

_This repository contains the code for data extraction, data integration and statistical analyses in "Rare and common vertebrates span a wide spectrum of population trends"._

Authors: Gergana N. Daskalova, Isla H. Myers-Smith, John L. Godlee

### Disclaimer
__The aim of this repository is to give reviewers of the manuscript private access to the code for all analyses in the manuscript. Copyright for the code remains with the authors. Please do not distribute or share the code with anyone.__

#### Contact: Gergana Daskalova gndaskalova@gmail.com

# Aim
_Our aim was to determine how the trends and fluctuations of vertebrate populations vary with biogeography, taxa, phylogenetic relationships and across species’ rarity metrics and IUCN conservation and threat categories._

# Preprint
The manuscript is available as a preprint on bioRxiv.

__Daskalova, G.N., Myers-Smith, I.H. & Godlee, J.L. (2018) Rare and common vertebrates span a wide spectrum of population trends. bioRxiv. https://doi.org/10.1101/272898__

# Workflow
For a detailed explanation of the methods, please see the preprint.

We analyzed vertebrate population time-series from the Living Planet Database (133,092 records) covering the period between 1970 and 2014. These time-series represent repeated monitoring surveys of the number of individuals in a given area (species’ abundance over time), hereafter called “populations”. We focus on two aspects of population change – overall changes in abundance over time (population trend) and abundance variability over time (population fluctuations). 

In the first stage of our analyses, we quantified trends and fluctuations for each population using state-space models that account for observation error and random fluctuations. 

<p align="center">
  <img src="/img/SI_workflow1.png" width=800 align="middle">
</p> 

_We analyzed vertebrate population time-series from the Living Planet Database (133,092 records) covering the period between 1970 and 2014. These time-series represent repeated monitoring surveys of the number of individuals in a given area (species’ abundance over time), to which we refer as “populations”. Diagram shows one sample population of Red squirrel (Sciurus vulgaris). We quantified two aspects of population change – overall changes in abundance over time (population trend) and abundance variability over time (population fluctuations). We used state-space models that account for observation error and random fluctuations1. See methods for additional details._


<p>In the second stage, we modelled the trend and fluctuation estimates from the first stage across latitude, realm, biome, taxa, rarity metrics, phylogenetic relatedness, species’ conservation status and threat type using a Bayesian modelling framework.</p.

<p align="center">
  <img src="/img/SI_workflow2.png" width=800 align="middle">
</p>

_We modelled the trend and fluctuation estimates from the first stage (Figure S1) across latitude, realm, biome, taxa, rarity metrics, phylogenetic relatedness, species’ conservation status and threat type using a Bayesian modelling framework2. Each model included a species random intercept effect to account for the possible correlation between the trends of populations from the same species. The prior structure (weakly informative priors) was identical across all models except the phylogeny models from the taxonomic patterns section, where the prior structure allowed for an additional phylogeny random effect. See methods for additional details._

# Data

### Living Planet Database
The population time-series analysed here come from the Living Planet Database, with the raw data publicly available from http://www.livingplanetindex.org

### Key data objects within the scripts and what they contain

__`LPIdata_Feb2016.csv`__
Raw Living Planet Database population time-series. Following the data use regulations, we cannot upload the raw data here, but they can be downloaded freely from http://www.livingplanetindex.org .

__`global_mus_scaled.csv`__
Population change (mu values for overall trend and sigma values for population fluctuations), taxa and biome classification and coordinates of the Living Planet Database time-series. Calculated using state-space models in `03-derive-pop-change-state-space.R`.

__`global_slopes.csv`__
Population change (model slopes of population abundance versus time), taxa and biome classification and coordinates of the Living Planet Database time series. Calculated using general linear models in `04-derive-pop-change-linear.R`.

__`IUCNall.csv`__
Species IUCN conservation status (Red List level), extracted from https://www.iucnredlist.org .

__`iucn_sp_habitats_count_global.csv`__
Species habitat specificity (number of distinct habitat in which each species occurs), calculated using the `rredlist` package in `01-calculate-habspec.R`.

__`habspec_profiles.csv`__
Habitat specificity for species monitored in the UK, derived by surveying species' profiles on the IUCN website https://www.iucnredlist.org .

__`gbif_ranges_clean.csv`__
Geographic ranges for species monitored in the UK, quantified based on occurrence records from GBIF in `01-calculate-habspec.R`.

__`bird_ranges.csv`__
Geographic ranges for bird species, extracted from BirdLife http://datazone.birdlife.org/home .

__`PanTHERIA_1-0_WR05_Aug2008.txt`__
Geographic ranges for mammal species, downloaded from the PanTHERIA database https://esajournals.onlinelibrary.wiley.com/doi/10.1890/08-1494.1 .

__`threats_lpi.RData`__
Species' threats, extracted from species' IUCN Red List classification https://www.iucnredlist.org .

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
Format data and run all statistical models.

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
If you would like to reproduce our analyses but are not familiar with using `R`, you can find an introduction to `R` and running code from `RStudio` on the <a href="https://ourcodingclub.github.io/2016/11/13/intro-to-r.html" target="_blank">Coding Club website</a>.

### Packages
The following `R` packages are required for our analyses. If they are not already installed on your computer, you can install them using the function `install.packages("package-name")`, where `package-name` is the name of the specific package, e.g. `install.packages("readr")`. The specific versions of the packages we used are outlined in <a href="https://github.com/gndaskalova/PopChangeRarity/blob/master/package-versions.md" target="_blank">`package-versions.md`</a>. Note that installing a package using `install.packages()` will automatically install the latest version of the package, which might be different from what we used for our analyses. For code that shows how to install older versions of packages, check out this <a href="https://stackoverflow.com/questions/17082341/installing-older-version-of-r-package" target="_blank">thread on Stackoverflow</a>. 

`readr, data.table, tidyr, dplyr, ggplot2, ggExtra, ggthemes, viridis, png, mapdata, maps, gridExtra, broom, MCMCglmm, stargazer, diptest, plotrix, scales, rredlist, stringr, corrplot, ggtree, proj4, ggalt, RColorBrewer, ggridges, ape, forcats, CoordinateCleaner, geosphere, parallel, scrubr`
