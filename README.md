# PopChangeRarity

_This repository contains the code for data extraction, data integration and statistical analyses in "All is not decline across global vertebrate populations"._

Authors: Gergana N. Daskalova, Isla H. Myers-Smith, John L. Godlee

### Disclaimer
__The aim of this repository is to give reviewers of the manuscript private access to the code for all analyses in the manuscript. Copyright for the code remains with the authors. Please do not distribute or share the code with anyone.__

#### Contact: Gergana Daskalova gndaskalova@gmail.com

# Aim
_Our aim was to quantify how populations (trends in numerical abundance) and biodiversity (trends in species richness and community composition) across vertebrate, invertebrate and plant taxa vary according to the timing and magnitude of forest cover change and habitat conversions._

# Preprint
The manuscript is available as a preprint on bioRxiv.

__Daskalova, G.N., Myers-Smith, I.H. & Godlee, J.L. (2018) Rarity and conservation status do not predict vertebrate population trends. bioRxiv. https://doi.org/10.1101/272898 __

## will update below tomorrow

# Workflow
For a detailed explanation on the methods, please see the preprint.

<p align="center">
  <img src="/img/Si_methods_diagram.png" width=800 align="middle">
</p> 

The workflow for one example time series is outlined below. Those steps were then repeated for all population and biodiversity time series included in the analysis.

<p align="center">
  <img src="/img/SI_workflow_time_series.png" width=800 align="middle">
</p> 

# Data

### Living Planet Database
The population time series analysed here came from the Living Planet Database, publicly available from http://www.livingplanetindex.org

### Key data objects within the scripts and what they contain
_Note that because we cannot make the whole BioTIME database public (see data section above) and because of the size of the files, we can only make our scripts, and not the data outputs they create, publicly available. However, we note that 90% of the raw BioTIME data and all other raw data are already publicly available (see above)._

__`rarefied_medians.RData`__
Yearly records for species richness and turnover, taxa and biome classification and coordinates of the BioTIME biodiversity time series. Output of Blowes and Supp et al.

__`global_mus_scaled.csv`__
Population change (mu values for overall trend), taxa and biome classification and coordinates of the Living Planet Database time series. Calculated using state-space models in `06-population-models.R`.

__`slopes_richness.RData`__
Slopes of richness change over time for the BioTIME biodiversity time series. Output of running `06-richness-models.R` with `rarefied_medians.RData`.

__`mus_luh.RData`__
Integrated data object containing the amount of population change and the amount of forest loss (over the same time period, using the Land Use Harmonisation Database) for each studied population. Created using the scripts from sections 1-7 below.

__`lpi_forest_change.RData`__
Integrated data object containing the amount of population change and the amount of forest cover gain and loss (over the same time period, using the Global Forest Change Database) for each studied population. Created using the scripts from sections 1-7 below.

__`biotime_luh_polys_change.RData`__
Integrated data object containing the amount of richness change and turnover and the amount of forest loss (over the same time period, using the Land Use Harmonisation Database) for each studied population. Created using the scripts from sections 1-7 below.

__`biotime_forest_change.RData`__
Integrated data object containing the amount of richness change and turnover and the amount of forest cover gain and loss (over the same time period, using the Global Forest Change Database) for each studied population. Created using the scripts from sections 1-7 below.

# Scripts

_Please note that the majority of the code was written to run on a HPC cluster and thus requires high computing power to successfully run._

### 1-4. Extract data

__`01-cell-size-sensitivity-luh.R`__
This script extracts land cover types from the Land Use Harmonisation Database for cells of different sizes to test the sensitity of our findings to observational scale. The amount of forest loss scaled proportionately with cell size, and we conducted all of our analyses using a cell size of approximately 96 squared kilometers.

__`01-cell-size-sensitivity-prep-hansen.R`__
Same as above, but with the Hansen Forest Cover Change Database.

__`01-prep-for-earth-engine-bt.R`__
This script prepares the coordinates of the cells of the BioTIME biodiversity time series in the right format so that they can imported as a fusion table in the Google Earth Engine.

__`01-prep-for-earth-engine-lpd.R`__
Same as above, but for the coordinates of the cells of the Living Planet Database population time series.

__`02-luh-land-cover-bt.R`__
This script extracts the land cover composition of the BioTIME cells over time (from 850 to 2015).

__`02-luh-land-cover-lpd.R`__
Same as above, but for the cells of the Living Planet Database.

__`02-luh-primf-bt-lpd.R`__
This script extracts the amount of primary forest cover of the BioTIME and Living Planet Database cells over time (from 850 to 2015).

__`02-gfc.js`__
This script extracts forest cover gain and loss over time using the Google Earth Engine and the Global Forest Change Database. See 01-prep-for-earth-engine-bt.R and 01-prep-for-earth-engine-lpd.R for details on how the cell coordinates were prepared for import in the GEE as fusion tables.

__`02-modis.js`__
This script extracts land cover composition over time using the Google Earth Engine and the MODIS Landcover Database. See 01-prep-for-earth-engine-bt.R and 01-prep-for-earth-engine-lpd.R for details on how the cell coordinates were prepared for import in the GEE as fusion tables.

__`03-extract-even-periods-peak-forest-loss-bt-lpd.R`__
This script extracts forest cover change in the same periods as the population and biodiversity monitoring data that were available before and after peak contemporary forest loss. The script contains the calculations for both the population and biodiversity time series.

__`03-luh-habitat-transitions-bt.R`__
This script integrates the data for all land cover types in the BioTIME cells over time and calculates the dominant land cover at the start and end of the monitoring period for each time series.

__`03-luh-habitat-transitions-lpd.R`__
Same as above, but for the population time series part of the Living Planet Database.

__`03-modis-habitat-transitions-bt.R`__
Combine MODIS land cover data from all years in which MODIS data were available (2000 - 2013) and determine dominant habitat type at the start and end of monitoring of the BioTIME biodiversity time series.

__`03-modis-habitat-transitions-lpd.R`__
Combine MODIS land cover data from all years in which MODIS data were available (2000 - 2013) and determine dominant habitat type at the start and end of monitoring of the Living Planet Database population time series.

__`04-extract-threats-IUCN.R`__
Extract the types of threats associated with the species part of the Living PLanet and BioTIME databases from their IUCN classifications.

### 5-6. Calculate population and richness change

Note that turnover was calculated by Blowes and Supp _et al._ The code for the rarefaction of the BioTIME database, as well as the calculation of turnover, is <a href="https://doi.org/10.5281/zenodo.1475218" target="_blank">archived on Zenodo</a>. 

__`05-population-models-period1.R`__
Calculate population change (_mu_ values for overall trend) in the first period - before contemporary peak forest loss. The code uses state-space models.

__`05-population-models-period2.R`__
Calculate population change (_mu_ values for overall trend) in the second period - after contemporary peak forest loss. The code uses state-space models.

__`05-richness-models-periods1-2.R`__
Calculate richness change (slopes of species richness over time) in the periods before and after contemporary peak forest loss.

__`06-biotime-popchange.R`__
Calculate population change (_mu_ values for overall trend) across the duration of each BioTIME population time series. The code uses state-space models.

__`06-population-models-modis.R`__
Calculate population change (_mu_ values for overall trend) across the period matching the duration of the MODIS database (2000-2013). The code uses state-space models.

__`06-population-models.R`__
Calculate population change (_mu_ values for overall trend) across the duration of each Living Planet Database time series. The code uses state-space models.

__`06-richness-models-modis.R`__
Calculate richness change (slopes of species richness over time) across the duration of each BioTIME biodiversity time series.


__`06-richness-models.R`__
Calculate richness change (slopes of species richness over time) across the duration of each BioTIME biodiversity time series.

### 7. Integrate land-use change, population change and biodiversity change data

__`07-time-matching-hansen-bt-cell-sizes.R`__
This script calculates the amount of forest cover gain and loss at different cell sizes (sensitivity analysis).

__`07-time-matching-hansen-bt.R`__
This script calculates the amount of forest cover gain and loss across the period of biodiversity monitoring (i.e., matching the time scales of the Hansen Global Forest Change Dataset (available from 2000 to 2016) and the biodiversity time series after 2000). For example, for a biodiversity time series ranging from 2000 to 2009, we calculated forest cover gain and loss also from 2000 to 2009.

__`07-time-matching-hansen-lpd.R`__
This script calculates the amount of forest cover gain and loss across the period of population monitoring (i.e., matching the time scales of the Hansen Global Forest Change Dataset (available from 2000 to 2016) and the population and biodiversity time series after 2000). For example, for a population time series ranging from 2000 to 2009, we calculated forest cover gain and loss also from 2000 to 2009.

__`07-time-matching-luh-bt-pop-change.R`__
Same as above, but for population change based on the BioTIME time series.

__`07-time-matching-luh-bt.R`__
This script calculates the amount of forest cover gain and loss across the period of biodiversity monitoring (i.e., matching the time scales of the Land Use Harmonisation Database and the full duration of the biodiversity time series). For example, for a biodiversity time series ranging from 1990 to 2009, we calculated forest loss also from 1990 to 2009.

__`07-time-matching-luh-lpd.R`__
This script calculates the amount of forest cover gain and loss across the period of population monitoring (i.e., matching the time scales of the Land Use Harmonisation Database and the full duration of the population time series). For example, for a population time series ranging from 1990 to 2009, we calculated forest cover gain and loss also from 1990 to 2009.

### Categorise time series and calculate lags

__`08-calculate-lags.R`__
This script calculates our estimates for temporal lag - the time period between when contemporary peak forest loss and maximum population/biodiversity change occurred.

__`08-categorise-before-after-during.R`__
This script categorises population time series based on whether the population monitoring started before, during or after the all time historic peak forest loss period. Forest loss during the period between 850 and 2015 was estimated using the Land Use Harmonisation Database (see `02-luh-primf-bt-lpd.R`).

### Statistical analyses

__`09-before-after-peak-forest-loss-models.R`__
Tests if population and biodiversity change differ before and after contemporary peak forest loss.

__`09-forest-cover-change-continuous-models.R`__
Test the relationships between forest loss, forest gain and population and biodiversity change.

__`09-models-lags.R`__
Tests if population and biodiversity change lags following contemporary peak forest loss are longer for species with longer generation times.

# Requirements

### Software
R version 3.5.1 or greater

### Packages
tidyverse, brms, ncdf4, raster, rgdal, splitstackshape, data.table, ggthemes, gridExtra, ggalt, viridis, ggridges, ggstatsplot, rredlist, tidybayes, modelr, bayesplot, sjstats, maps, mapdata, scales, forcats, proj4, hrbrthemes
