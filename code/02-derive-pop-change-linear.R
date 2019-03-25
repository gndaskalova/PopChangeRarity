# How does rarity influence rate of population change in UK vetebrates?
# How has global vertebrate abundance changed across biomes and taxa?

# Gergana Daskalova (gndaskalova@gmail.com)
# Isla Myers-Smith (isla.myers-smith@ed.ac.uk)
# John Godlee (johngodlee@gmail.com)

# Packages ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(broom)

# Load data ----
LPI <- read.csv("data/input/LPIdata_Feb2016.csv")  # raw LPI data

# Transform from wide to long format
LPI.long <- gather(data = LPI, key = "year", value = "pop", select = 26:70)

# Get rid of the X in front of years
LPI.long$year <- parse_number(LPI.long$year)

# Create new column with genus and species together
LPI.long$species <- paste(LPI.long$Genus, LPI.long$Species)

# ** Data manipulation ----

# Keep species with >5 years worth of data 
# Calculate length of monitoring and scale population trend data
LPI.long <- LPI.long %>%
  drop_na(pop) %>%
  group_by(id) %>%   # group rows so that each group is one population
  mutate(scalepop = rescale(pop, to = c(-1, 1))) %>%
  filter(length(unique(year)) > 5) %>%
  drop_na(scalepop) %>%
  mutate(meanpop = mean(pop),  # Create column for mean population
         minyear = min(year),
         maxyear = max(year),
         lengthyear = maxyear - minyear) %>%
  ungroup()

# Number of species = 2074
length(unique(LPI.long$species))

# Number of populations = 9288
length(unique(LPI.long$id))

# ** Models ----
# Run linear models of abundance trends over time for each population and extract model coefficients
LPI.models <- LPI.long %>%
  group_by(biome, system, Country.list, Class, species, lengthyear, meanpop, Decimal.Latitude, Decimal.Longitude, id) %>% 
  do(mod = lm(scalepop ~ year, data = .)) %>%  # Create a linear model for each group
  mutate(., n = df.residual(mod),  # Create columns: degrees of freedom
         intercept = summary(mod)$coeff[1],  # intercept coefficient
         slope = summary(mod)$coeff[2],  # slope coefficient
         intercept_se = summary(mod)$coeff[3],  # standard error of intercept
         slope_se = summary(mod)$coeff[4],  # standard error of slope
         intercept_p = summary(mod)$coeff[7],  # p value of intercept
         slope_p = summary(mod)$coeff[8]) %>%  # p value of slope
  ungroup() %>%
  mutate(id = id,
         biome = biome,
         system = system,
         Country.list = Country.list,
         Class = Class,
         species = species,
         lengthyear = lengthyear,
         meanpop = meanpop,
         Decimal.Latitude = Decimal.Latitude,
         Decimal.Longitude = Decimal.Longitude)

# Number of species = 2074
length(unique(LPI.models$species))

# Number of populations = 9284
length(unique(LPI.models$id))

# Merge model outputs with IUCN conservation status
IUCNall <- read.csv("data/input/IUCNall.csv")  # species IUCN categories

IUCN <- IUCNall %>% mutate(species = paste(Genus, Species, sep = ' ')) %>%
  dplyr::select(species, Red.List.status)

LPI.models.IUCN <- merge(LPI.models, IUCN, by = "species", all.x = T)
LPI.models.IUCN$Red.List.status <- recode(LPI.models.IUCN$Red.List.status, "LR/lc" = "LC")
LPI.models.IUCN$Red.List.status <- recode(LPI.models.IUCN$Red.List.status, "LR/nt" = "NT")
LPI.models.IUCN$Red.List.status <- recode(LPI.models.IUCN$Red.List.status, "LR/cd" = "CD")

global_slopes <- LPI.models.IUCN %>% dplyr::select(-mod)

write.csv(global_slopes, "data/input/global_slopes.csv")
