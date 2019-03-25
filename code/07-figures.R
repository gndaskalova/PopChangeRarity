# How does vertebrate population change differ across biomes and taxa?
# How does rarity influence population change in UK vertebrates?

# Gergana Daskalova (gndaskalova@gmail.com)

# 19th Feb 2019
# Figures

# Start code with an empty working directory

packrat::init()

# Packages ----
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
library(ape)
library(proj4)
library(ggalt)
library(RColorBrewer)
library(ggridges)
library(forcats)

# Define ggplot2 theme functions ----
theme_LPI <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 20, hjust = -2),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

# With angled x axis labels
theme_LPI2 <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 20),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12, face = "italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

# Without x axis labels
theme_LPI3 <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title = element_text(size = 20),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12, face = "italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))


# Define function to extract MCMC outputs (no random effects) ----
# Function by Gabriela Hajduk
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # adding the model name

# Load data ----
LPI <- read.csv("data/input/LPIdata_Feb2016.csv")  # raw LPI data
mus <- read.csv("data/input/global_mus_scaled.csv")  # overall population change using mu
slopes <- read.csv("data/input/global_slopes.csv")  # overall pop change using slopes
IUCNall <- read.csv("data/input/IUCNall.csv")  # species IUCN categories
habspec <- read.csv("data/input/iucn_sp_habitats_count_global.csv")  # hab spec from rredlist package
habspec_prof <- read.csv("data/input/habspec_profiles.csv")  # hab spec from Red List profiles for UK species
ranges <- read.csv("data/input/gbif_ranges_clean.csv") # range data from GBIF
bird_ranges <- read.csv("data/input/bird_ranges.csv")  # data from BirdLife
mammal_ranges <- read_delim("data/input/PanTHERIA_1-0_WR05_Aug2008.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
mammal_ranges <- mammal_ranges %>% dplyr::select(MSW05_Binomial, `26-1_GR_Area_km2`)
load("data/input/threats_lpi.RData")  # threat data for each species from IUCN

# Add duration to the data frame from the state-space models
durations <- slopes %>% dplyr::select(id, lengthyear)
colnames(durations) <- c("id", "duration")
mus <- left_join(mus, durations, by = "id")

# Rename column names and make a species column
habspec <- habspec[, c(2,3)]
colnames(habspec) <- c("species", "hab_spec")
mus <- mus %>% mutate(species = paste(Genus, Species, sep = ' '))

# Checking data ----
str(mus)
# Number of populations = 9284
length(unique(mus$id))

# Number of populations with the same lat and long
mus$coord <- paste(mus$Decimal.Latitude, mus$Decimal.Longitude, sep = " ")
length(mus$coord)  # 9286
length(unique(mus$coord))  # 3033
3033/9286
# 6253 geographic locations have more than one population

# Number of species with more than one population
pops <- mus %>% group_by(species) %>% tally() %>% 
  filter(n > 3)

# 567 species
673/2074

# Number of species = 2074
length(unique(mus$species))

sp.geo <- mus %>% group_by(coord, species) %>%
  summarise(loc.count = n()) %>% filter(loc.count > 1) 
# 294 species have more than one geographic location
294/2074 # 14%

# Creating data frames for the diferent analyses ----
# Remove "unknown" biome
# Excluding Unknown and Mangroves (very little data for mangroves)
mus.biome <- filter(mus, biome != "Unknown" & biome != "Mangroves" & 
                      biome != "Oceanic islands")

# Renaming some of the habitats, so that the names fit on the plot
print(unique(mus.biome$biome))

mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical floodplain rivers and wetland complexes" = 
                            "Tropical wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical coastal rivers" = 
                            "Tropical wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical coniferous forests" = 
                            "Tropical and subtropical forests")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical moist broadleaf forests" = 
                            "Tropical and subtropical forests")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical dry broadleaf forests" = 
                            "Tropical and subtropical forests")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical upland rivers" = 
                            "Tropical wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Tropical and subtropical grasslands savannas and shrublands" = 
                            "Trop. and subtrop. grasslands savannas and shrublands")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Flooded grasslands and savannas" = 
                            "Trop. and subtrop. grasslands savannas and shrublands")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate coastal rivers" = 
                            "Temperate wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate floodplain rivers and wetlands" = 
                            "Temperate wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate broadleaf and mixed forests" = 
                            "Temperate forests")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate upland rivers" = 
                            "Temperate wetlands and rivers")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate coniferous forests" = 
                            "Temperate forests")
mus.biome$biome <- recode(mus.biome$biome, 
                          "Temperate upwelling" = 
                            "Temperate wetlands and rivers")

# Get data for birds, mammals, reptiles, amphibians and fish
# The taxa for which there are most data
mu.class <- filter(mus, Class == "Aves" | Class == "Mammalia" | 
                     Class == "Reptilia" | Class == "Amphibia" |
                     Class == "Actinopterygii" | Class == "Elasmobranchii")

# Merge mus with IUCN conservation status
IUCN <- IUCNall %>% mutate(species = paste(Genus, Species, sep = ' ')) %>%
  dplyr::select(species, Red.List.status) %>% distinct()
mus.IUCN <- inner_join(mus, IUCN, by = "species")

# For some species there are no IUCN categories
mus.IUCN$Red.List.status <- recode(mus.IUCN$Red.List.status, "LR/lc" = "LC")
mus.IUCN$Red.List.status <- recode(mus.IUCN$Red.List.status, "LR/nt" = "NT")
mus.IUCN$Red.List.status <- recode(mus.IUCN$Red.List.status, "LR/cd" = "CD")

# Exclude DD (Data deficient species) and EW (Extinct in the wild)
mus.IUCN <- filter(mus.IUCN, Red.List.status != "DD" & Red.List.status != "EW")
# Reorder levels to match increasing risk
mus.IUCN$Red.List.status <- factor(mus.IUCN$Red.List.status, 
                                   levels = c("LC", "NT", "VU", "EN", "CR"),
                                   labels=c("Least concern", "Near threatened", 
                                            "Vulnerable", "Endangered",
                                            "Critically endangered"))

# Merge mus and hab spec data
mus.habspec <- inner_join(mus, habspec, by = "species")
# Species = 1709
length(unique(mus.habspec$species))
# Populations = 7901
length(unique(mus.habspec$id))

# Calculate mean pop size for all populations that have count data

# Transform from wide to long format
LPI.long <- gather(data = LPI, key = "year", value = "pop", select = 26:70)
# Get rid of the X in front of years
LPI.long$year <- parse_number(LPI.long$year)

# Create new column with genus and species together
LPI.long$species <- paste(LPI.long$Genus, LPI.long$Species)

# Remove duplicate rows
LPI.long <- LPI.long %>% distinct() %>% filter(is.finite(pop))

# Get just the populatations with count data
LPI.long.pop <- LPI.long

load("data/input/LPI_mean_pop.Rdata")

LPI.long <- LPI.long %>%
  drop_na(pop) %>%
  group_by(id) %>%   # group rows so that each group is one population
  mutate(scalepop = scales::rescale(pop, to = c(-1, 1))) %>%
  filter(length(unique(year)) > 5) %>%
  drop_na(scalepop)

mus.popsize <- inner_join(mus, LPI.mean.pop, by = "id") %>%
  dplyr::select(mu, id, meanpop, lCI, uCI, tau.2, sigma.2, species)
mus.popsize <- na.omit(mus.popsize)
mus.popsize$log.meanpop <- log(mus.popsize$meanpop)
mus.popsize <- mus.popsize %>% filter(is.finite(log.meanpop))

# Ranges for bird species
birds_mu <- mus %>% filter(Class == "Aves")
bird_mu_ranges <- inner_join(birds_mu, bird_ranges, by = "species")
bird_mu_ranges <- bird_mu_ranges %>% filter(range != "EXTINCT") %>%
  drop_na(range)
bird_mu_ranges$range <- as.numeric(as.character(bird_mu_ranges$range))

# Ranges for mammal species
mammals_mu <- mus %>% filter(Class == "Mammalia")
colnames(mammal_ranges) <- c("species", "range")
mammal_mu_ranges <- left_join(mammals_mu, mammal_ranges, by = "species")

# Combine birds and mammals
bird_mu_ranges <- bird_mu_ranges %>% dplyr::select(-person, -notes)
b_m_ranges <- rbind(bird_mu_ranges, mammal_mu_ranges)

# Combine population change and threat data
threats_sum <- threats_lpi %>% group_by(title) %>%
  tally() %>% arrange(desc(n))

threats_sum <- filter(threats_sum, !title %in% c("Unspecified species",
                                                 "Type Unknown/Unrecorded",
                                                 "Scale Unknown/Unrecorded",
                                                 "Motivation Unknown/Unrecorded"))
threats_sum$title <- as.factor(threats_sum$title)
top_threats <- threats_sum[1:10,]
top_threats <- arrange(top_threats, desc(n))

# Threats and population change
species_mus <- mus %>% dplyr::select(mu, sigma.2, species, id)
threats_mus <- inner_join(species_mus, threats_lpi, by = "species")


top_threats_mus <- threats_mus %>% filter( title %in% c("Fishing & harvesting aquatic resources",
                                                        "Hunting & trapping terrestrial animals",
                                                        "Annual & perennial non-timber crops",
                                                        "Intentional use (species is the target)",
                                                        "Habitat shifting & alteration",
                                                        "Housing & urban areas",
                                                        "Unintentional effects: (large scale) [harvest]",
                                                        "Intentional use: (subsistence/small scale) [harvest]",
                                                        "Agricultural & forestry effluents",
                                                        "Intentional use: (large scale) [harvest]"))

top_threats_mus$title <- factor(top_threats_mus$title,
                                levels = c("Fishing & harvesting aquatic resources",
                                           "Hunting & trapping terrestrial animals",
                                           "Annual & perennial non-timber crops",
                                           "Intentional use (species is the target)",
                                           "Habitat shifting & alteration",
                                           "Housing & urban areas",
                                           "Unintentional effects: (large scale) [harvest]",
                                           "Intentional use: (subsistence/small scale) [harvest]",
                                           "Agricultural & forestry effluents",
                                           "Intentional use: (large scale) [harvest]"),
                                labels = c("Fishing & harvesting aquatic resources",
                                           "Hunting & trapping terrestrial animals",
                                           "Annual & perennial non-timber crops",
                                           "Intentional use (species is the target)",
                                           "Habitat shifting & alteration",
                                           "Housing & urban areas",
                                           "Unintentional effects: (large scale) [harvest]",
                                           "Intentional use: (subsistence/small scale) [harvest]",
                                           "Agricultural & forestry effluents",
                                           "Intentional use: (large scale) [harvest]"))
# Summary by the numbers
threat_mu_sum <- top_threats_mus %>% group_by(title) %>% tally()

# Data frames for slope models
# Remove "unknown" biome
# Excluding Unknown and Mangroves (no data for mangroves)
slopes.biome <- filter(slopes, biome != "Unknown" & biome != "Mangroves" & 
                         biome != "Oceanic islands")

slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical floodplain rivers and wetland complexes" = 
                               "Tropical wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical coastal rivers" = 
                               "Tropical wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical coniferous forests" = 
                               "Tropical and subtropical forests")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical moist broadleaf forests" = 
                               "Tropical and subtropical forests")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical dry broadleaf forests" = 
                               "Tropical and subtropical forests")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical upland rivers" = 
                               "Tropical wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Tropical and subtropical grasslands savannas and shrublands" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Flooded grasslands and savannas" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate coastal rivers" = 
                               "Temperate wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate floodplain rivers and wetlands" = 
                               "Temperate wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate broadleaf and mixed forests" = 
                               "Temperate forests")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate upland rivers" = 
                               "Temperate wetlands and rivers")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate coniferous forests" = 
                               "Temperate forests")
slopes.biome$biome <- recode(slopes.biome$biome, 
                             "Temperate upwelling" = 
                               "Temperate wetlands and rivers")

# Get data for birds, mammals, reptiles, amphibians and fish
# The taxa for which there are most data
slope.class <- filter(slopes, Class == "Aves" | Class == "Mammalia" | 
                        Class == "Reptilia" | Class == "Amphibia" |
                        Class == "Actinopterygii" | Class == "Elasmobranchii")

# Exclude DD (Data deficient species) and EW (Extinct in the wild)
slopes.IUCN <- filter(slopes, Red.List.status != "EW" & Red.List.status !="DD" &
                        Red.List.status != "EX" & Red.List.status != "CD")

# Reorder levels to match increasing risk
slopes.IUCN$Red.List.status <- factor(slopes.IUCN$Red.List.status, 
                                      levels = c("LC", "NT", "VU", "EN", "CR"),
                                      labels=c("Least concern", "Near threatened", 
                                               "Vulnerable", "Endangered",
                                               "Critically endangered"))

# Merge mus and hab spec data
slopes.habspec <- inner_join(slopes, habspec, by = "species")

slopes.popsize <- inner_join(slopes, LPI.mean.pop, by = "id") %>%
  dplyr::select(slope, id, meanpop, Country.list, slope_se, species)

# Global birds
birds_slopes <- slopes %>% filter(Class == "Aves")
bird_slopes_ranges <- inner_join(birds_slopes, bird_ranges, by = "species")
bird_slopes_ranges <- bird_slopes_ranges %>% filter(range != "EXTINCT")
bird_slopes_ranges$range <- as.numeric(as.character(bird_slopes_ranges$range))
bird_slopes_ranges <- bird_slopes_ranges %>%
  mutate(log.range = log(range)) %>% filter(is.finite(log.range))

# Ranges for mammal species
mammals_slopes <- slopes %>% filter(Class == "Mammalia")
mammal_slopes_ranges <- left_join(mammals_slopes, mammal_ranges, by = "species")
mammal_slopes_ranges <- mammal_slopes_ranges %>%
  mutate(log.range = log(range)) %>% drop_na(log.range)

# Combine birds and mammals
bird_slopes_ranges <- bird_slopes_ranges %>% dplyr::select(-person, -notes)
b_m_slopes_ranges <- rbind(bird_slopes_ranges, mammal_slopes_ranges)

# Number of threats and pop trend
threats_lpi2 <- filter(threats_lpi, !title %in% c("Unspecified species",
                                                  "Type Unknown/Unrecorded",
                                                  "Scale Unknown/Unrecorded",
                                                  "Motivation Unknown/Unrecorded"))

threats_sum_species <- threats_lpi2 %>% group_by(species) %>%
  tally() %>% arrange(desc(n))

threats_sum_species_mus <- inner_join(species_mus, threats_sum_species,
                                      by = "species")

# Birds and mammals
b_m_ranges$log.range <- log(b_m_ranges$range)
b_m_ranges <- b_m_ranges %>% drop_na(log.range)

# Global mean pop size model
# Add taxa back to the data frame

taxa <- mus %>% dplyr::select(id, Class) %>% distinct()
mus.popsize <- left_join(mus.popsize, taxa, by = "id")

# Global mean pop size
slopes.popsize <- na.omit(slopes.popsize)
slopes.popsize$log.meanpop <- log(slopes.popsize$meanpop)
slopes.popsize <- slopes.popsize %>% filter(is.finite(log.meanpop))

# System
mus$CI <- (mus$uCI - mus$lCI)/2
mus$CI <- scales::rescale(mus$CI, to = c(0, 1))
system.ci <- mus %>% group_by(system) %>%
  mutate(sd = sd(CI)) %>% dplyr::select(system, sd) %>% distinct()

system.ci$Model <- "CI"

mus$system <- factor(mus$system, 
                     levels = c("Terrestrial", "Marine", "Freshwater"),
                     labels = c("Terrestrial", "Marine", "Freshwater"))

# Biome
mus.biome$CI <- (mus.biome$uCI - mus.biome$lCI)/2
mus.biome$CI <- scales::rescale(mus.biome$CI, to = c(0, 1))

# Taxa
mu.class$CI <- (mu.class$uCI - mu.class$lCI)/2
mu.class$CI <- scales::rescale(mu.class$CI, to = c(0, 1))

class.ci <- mu.class %>% group_by(Class) %>%
  mutate(sd = sd(CI)) %>% dplyr::select(Class, sd) %>% distinct()

class.ci$Model <- "CI"

# IUCN
mus.IUCN$CI <- (mus.IUCN$uCI - mus.IUCN$lCI)/2
mus.IUCN$CI <- scales::rescale(mus.IUCN$CI, to = c(0, 1))

IUCN.ci <- mus.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(CI)) %>% dplyr::select(Red.List.status, sd) %>% distinct()

IUCN.ci$Model <- "CI"

# Global hab spec
mus.habspec$CI <- (mus.habspec$uCI - mus.habspec$lCI)/2
mus.habspec$CI <- scales::rescale(mus.habspec$CI, to = c(0, 1))

# Global mean pop size
mus.popsize$CI <- (mus.popsize$uCI - mus.popsize$lCI)/2
mus.popsize$CI <- scales::rescale(mus.popsize$CI, to = c(0, 1))

# Global bird and mammal range
b_m_ranges$CI <- (b_m_ranges$uCI - b_m_ranges$lCI)/2
b_m_ranges$CI <- scales::rescale(b_m_ranges$CI, to = c(0, 1))

# Data frame for the sd of the input data
system.ci <- mus %>% group_by(system) %>%
  mutate(sd = sd(CI)) %>% dplyr::select(system, sd) %>% distinct()

system.ci.var <- system.ci
system.ci.var$Model <- "CI.w"

# Data frame for the sd of the input data
# same as non-weighted data frame, since they are both
# based on the same input data (mu values)
class.ci.var <- class.ci
class.ci.var$Model <- "CI.w"

# Data frame for the sd of the input data
# same as non-weighted data frame, since they are both
# based on the same input data (mu values)
IUCN.ci.var <- IUCN.ci
IUCN.ci.var$Model <- "CI.w"

slopes$SE <- scales::rescale(slopes$slope_se, to = c(0, 1))
system.se <- slopes %>% group_by(system) %>%
  mutate(sd = sd(SE)) %>% dplyr::select(system, sd) %>% distinct()
system.se$Model <- "se"

slope.class$SE <- scales::rescale(slope.class$slope_se, to = c(0, 1))
class.se <- slope.class %>% group_by(Class) %>%
  mutate(sd = sd(SE)) %>% dplyr::select(Class, sd) %>% distinct()
class.se$Model <- "se"

slopes.IUCN$SE <- scales::rescale(slopes.IUCN$slope_se, to = c(0, 1))
IUCN.se <- slopes.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(SE)) %>% dplyr::select(Red.List.status, sd) %>% distinct()
IUCN.se$Model <- "se"

# SD
# System
LPI.long <- LPI.long %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))

LPI.long$sdpop1 <- scales::rescale(LPI.long$sdpop, to = c(0, 1))
LPIraw <- LPI.long %>% dplyr::select(id, species, sdpop, sdpop1, Country.list, 
                              system, biome, Class) %>%
  distinct()

system.sd <- LPI.long %>% group_by(system) %>%
  mutate(sd = sd(scalepop)) %>% dplyr::select(system, sd) %>%
  distinct()
system.sd$Model <- "sd"

# Biome
LPI.long.biome <- filter(LPI.long, biome != "Unknown" & biome != "Mangroves" & 
                           biome != "Oceanic islands")

LPI.long.biome$sdpop1 <- scales::rescale(LPI.long.biome$sdpop, to = c(0, 1))
LPIraw.biome <- LPI.long.biome %>% dplyr::select(id, species, sdpop, sdpop1, Country.list, 
                                          system, biome, Class) %>%
  distinct()


LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical floodplain rivers and wetland complexes" = 
                               "Tropical wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical coastal rivers" = 
                               "Tropical wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical coniferous forests" = 
                               "Tropical and subtropical forests")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical moist broadleaf forests" = 
                               "Tropical and subtropical forests")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical dry broadleaf forests" = 
                               "Tropical and subtropical forests")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical upland rivers" = 
                               "Tropical wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Tropical and subtropical grasslands savannas and shrublands" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Flooded grasslands and savannas" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate coastal rivers" = 
                               "Temperate wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate floodplain rivers and wetlands" = 
                               "Temperate wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate broadleaf and mixed forests" = 
                               "Temperate forests")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate upland rivers" = 
                               "Temperate wetlands and rivers")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate coniferous forests" = 
                               "Temperate forests")
LPIraw.biome$biome <- recode(LPIraw.biome$biome, 
                             "Temperate upwelling" = 
                               "Temperate wetlands and rivers")

# Taxa
LPI.long.class <- filter(LPI.long, Class== "Aves" | Class == "Mammalia" | 
                           Class == "Reptilia" | Class == "Amphibia" |
                           Class == "Actinopterygii" | Class == "Elasmobranchii")
LPI.long.class <- LPI.long.class %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))

LPI.long.class$sdpop1 <- scales::rescale(LPI.long.class$sdpop, to = c(0, 1))

LPIraw.class <- LPI.long.class %>% dplyr::select(id, species, sdpop, sdpop1, Country.list, 
                                          system, biome, Class) %>%
  distinct()

class.sd <- LPI.long.class %>% group_by(Class) %>%
  mutate(sd = sd(scalepop)) %>% dplyr::select(Class, sd) %>%
  distinct()
class.sd$Model <- "sd"

#class.sd <- LPIraw.class %>% group_by(Class) %>%
#  mutate(sd = sd(sdpop1)) %>% dplyr::select(Class, sd) %>%
#  distinct()
#class.sd$Model <- "sd"

# IUCN
LPIraw.IUCN <- inner_join(LPI.long, IUCN, by = "species")
# Exclude DD (Data deficient species) and EW (Extinct in the wild)
LPIraw.IUCN <- filter(LPIraw.IUCN, Red.List.status != "DD" & 
                        Red.List.status != "EW")

# Reorder levels to match increasing risk
LPIraw.IUCN$Red.List.status <- factor(LPIraw.IUCN$Red.List.status, 
                                      levels = c("LC", "NT", "VU", "EN", "CR"),
                                      labels = c("Least concern", "Near threatened", 
                                                 "Vulnerable", "Endangered",
                                                 "Critically endangered"))

LPIraw.IUCN <- LPIraw.IUCN %>% drop_na(Red.List.status)

LPIraw.IUCN <- LPIraw.IUCN %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))

LPIraw.IUCN$sdpop1 <- scales::rescale(LPIraw.IUCN$sdpop, to = c(0, 1))
LPIraw.IUCN <- LPIraw.IUCN %>% dplyr::select(id, species, sdpop, sdpop1, Country.list, 
                                      system, biome, Class, Red.List.status) %>%
  distinct()

IUCN.sd <- LPIraw.IUCN %>% drop_na(Red.List.status)
IUCN.sd <- IUCN.sd %>% group_by(Red.List.status) %>%
  mutate(sd = sd(sdpop1)) %>% dplyr::select(Red.List.status, sd) %>%
  distinct()

IUCN.sd$Model <- "sd"

# Global hab spec
LPIraw.habspec <- inner_join(LPI.long, habspec, by = "species")

LPIraw.habspec$sdpop1 <- scales::rescale(LPIraw.habspec$sdpop, to = c(0, 1))

# Global mean pop size
LPI.mean.pop2 <- LPI.mean.pop
id_species <- mus.popsize %>% dplyr::select(species, id)
LPI.mean.pop2 <- left_join(LPI.mean.pop2, id_species, by = "id")
LPI.mean.pop2 <- LPI.mean.pop2 %>% dplyr::select(id, sdpop, meanpop, species) %>%
  distinct() %>% drop_na(species)
LPI.mean.pop2$sdpop1 <- scales::rescale(LPI.mean.pop2$sdpop, to = c(0, 1))
LPI.mean.pop2$log.meanpop <- log(LPI.mean.pop2$meanpop)
LPI.mean.pop2 <- LPI.mean.pop2 %>% filter(is.finite(log.meanpop))

# Global bird range
ranges_simple <- b_m_ranges %>% dplyr::select(species, range) %>% distinct()
LPIraw.bird.mammal.range <- left_join(LPI.long, ranges_simple, by = "species")
LPIraw.bird.mammal.range$sdpop1 <- scales::rescale(LPIraw.bird.mammal.range$sdpop, to = c(0, 1))
LPIraw.bird.mammal.range <- LPIraw.bird.mammal.range %>% dplyr::select(id, sdpop1, range, species) %>%
  distinct()
LPIraw.bird.mammal.range <- filter(LPIraw.bird.mammal.range, range != "EXTINCT")
LPIraw.bird.mammal.range$range <- as.numeric(as.character(LPIraw.bird.mammal.range$range))
LPIraw.bird.mammal.range <- LPIraw.bird.mammal.range %>% 
  mutate(log.range = log(range)) %>% drop_na(log.range)

# Data frame for the sd of the input data
mus$sigma.2s <- scales::rescale(mus$sigma.2, to = c(0, 1))

system.ci.sigma <- mus %>% group_by(system) %>%
  mutate(sd = sd(sigma.2s)) %>% dplyr::select(system, sd) %>% distinct()

system.ci.sigma$Model <- "sigma"

# Data frame for the sd of the input data
mu.class$sigma.2s <- scales::rescale(mu.class$sigma.2, to = c(0, 1))
taxa.ci.sigma <- mu.class %>% group_by(Class) %>%
  mutate(sd = sd(sigma.2s)) %>% dplyr::select(Class, sd) %>% distinct()

taxa.ci.sigma$Model <- "sigma"

# Data frame for the sd of the input data
mus.biome$sigma.2s <- scales::rescale(mus.biome$sigma.2, to = c(0, 1))
biome.ci.sigma <- mus.biome %>% group_by(biome) %>%
  mutate(sd = sd(sigma.2s)) %>% dplyr::select(biome, sd) %>% distinct()

biome.ci.sigma$Model <- "sigma"

mus.IUCN$sigma.2s <- scales::rescale(mus.IUCN$sigma.2, to = c(0, 1))
IUCN.ci.sigma <- mus.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(sigma.2s)) %>% dplyr::select(Red.List.status, sd) %>% distinct()

IUCN.ci.sigma$Model <- "sigma"

# UK scale ----
uk.mus <- mus %>% filter(grepl("United Kingdom", Country.list))
uk.mus.habspec <- inner_join(uk.mus, habspec, by = "species") %>% 
  drop_na(hab_spec)
uk.mus.popsize <- inner_join(uk.mus, mus.popsize, by = "id") %>%
  dplyr::select(-mu.y)
colnames(uk.mus.popsize)[27] <- "mu"
uk.mus.popsize <- left_join(uk.mus.popsize, id_species, by = "id")

# Get geographic range data
ranges <- ranges[, -1]
colnames(ranges) <- c("species", "km2_range")
uk.mus.range <- inner_join(uk.mus, ranges, by = "species") %>%
  mutate(log.range = log(km2_range)) %>%
  filter(is.finite(log.range))

uk.mus.habspec.prof <- inner_join(uk.mus, habspec_prof, by = "species") %>%
  drop_na(hab_spec)

# IUCN status
uk.mus.IUCN <- mus.IUCN %>% filter(grepl("United Kingdom", Country.list))

# Differences across taxa in the UK
uk.mus.class <- filter(uk.mus, Class == "Aves" | Class == "Mammalia" | 
                         Class == "Reptilia" | Class == "Amphibia" |
                         Class == "Actinopterygii" | Class == "Elasmobranchii")

# ** slope
uk.slopes <- slopes %>% filter(grepl("United Kingdom", Country.list))
uk.slopes.ranges <- inner_join(uk.slopes, ranges, by = "species")
uk.slopes.ranges <- uk.slopes.ranges %>% 
  mutate(log.range = log(km2_range)) %>%
  filter(is.finite(log.range))

# Population size
uk.slopes.popsize <- inner_join(uk.slopes, slopes.popsize, by = "id")%>%
  dplyr::select(-slope.y)
colnames(uk.slopes.popsize)[14] <- "slope"
uk.slopes.popsize <- left_join(uk.slopes.popsize, id_species, by = "id")

# Habitat specificity
uk.slopes.habspec <- inner_join(uk.slopes, habspec, by = "species") %>% 
  drop_na(hab_spec)
uk.slopes.habspec.prof <- inner_join(uk.slopes, habspec_prof, by = "species") %>% 
  drop_na(hab_spec)

# IUCN status
uk.slopes.IUCN <- slopes.IUCN %>% filter(grepl("United Kingdom", Country.list))

# Differences across taxa in the UK
uk.slopes.class <- filter(uk.slopes, Class == "Aves" | Class == "Mammalia" | 
                            Class == "Reptilia" | Class == "Amphibia" |
                            Class == "Actinopterygii" | Class == "Elasmobranchii")

# ** CI 
# Geographic range
uk.mus.range$CI <- (uk.mus.range$uCI - uk.mus.range$lCI)/2

# Population size
uk.mus.popsize$CI <- (uk.mus.popsize$uCI.x - uk.mus.popsize$lCI.x)/2

# Habitat specificity
uk.mus.habspec$CI <- (uk.mus.habspec$uCI - uk.mus.habspec$lCI)/2

# IUCN status
uk.mus.IUCN$CI <- (uk.mus.IUCN$uCI - uk.mus.IUCN$lCI)/2

# Differences across systems in the UK
uk.mus$CI <- (uk.mus$uCI - uk.mus$lCI)/2

# ** SD

UK.LPIraw <- LPIraw %>%  filter(grepl("United Kingdom", Country.list))
ranges_uk <- uk.mus.range %>% dplyr::select(id, km2_range)

UK.LPIraw.range <- left_join(UK.LPIraw, ranges_uk, by = "id")
UK.LPIraw.range$sdpop1 <- scales::rescale(UK.LPIraw.range$sdpop, to = c(0, 1))

uk_popsize_simple <- uk.mus.popsize %>% dplyr::select(id, log.meanpop)
UK.LPIraw.popsize <- left_join(UK.LPIraw, uk_popsize_simple, by = "id")
# Remove the NAs - those were for the populations that didn't have "true"
# population abundances, e.g. they were measured using an index
UK.LPIraw.popsize <- UK.LPIraw.popsize %>% drop_na(log.meanpop)
UK.LPIraw.popsize$sdpop1 <- scales::rescale(UK.LPIraw.popsize$sdpop, to = c(0, 1))

UK.LPIraw.habspec <- left_join(UK.LPIraw, habspec, by = "species")
UK.LPIraw.habspec <- UK.LPIraw.habspec %>% drop_na(hab_spec)
UK.LPIraw.habspec$sdpop1 <- scales::rescale(UK.LPIraw.habspec$sdpop, to = c(0, 1))

UK.LPIraw.IUCN <- left_join(UK.LPIraw, uk.mus.IUCN, by = "species")
UK.LPIraw.IUCN <- UK.LPIraw.IUCN %>% drop_na(Red.List.status)
UK.LPIraw.IUCN$sdpop1 <- scales::rescale(UK.LPIraw.IUCN$sdpop, to = c(0, 1))

UK.LPIraw$sdpop1 <- scales::rescale(UK.LPIraw$sdpop, to = c(0, 1))

# Geographic range
UK.LPIraw.range$log.range <- log(UK.LPIraw.range$km2_range)
UK.LPIraw.range <- UK.LPIraw.range %>% drop_na(log.range)

# Figures ----
# Load all models
load("~/RarityHub/data/output/global_models_18thFeb.RData")
load("~/RarityHub/data/output/models_uk18thFeb.RData")
list2env(model.list.global, globalenv())
list2env(model.list.uk, globalenv())

# Load effect sizes
uk_outputs18thFeb <- read.csv("data/output/uk_outputs18thFeb.csv")
global_outputs18thFeb <- read.csv("data/output/global_outputs18thFeb.csv")

# ** Rarity figure ----

# Global birds and mammals
mcmc_preds_range_global <- as.data.frame(predict(rarity.range.birds.mammals,
                                                 type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_range_global$range <- b_m_ranges$range

(range.plot.birds.mammals <- ggplot() +
    geom_point(data = b_m_ranges, aes(x = log(range), y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_range_global, aes(x = log(range), ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_range_global, aes(x = log(range), y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                       limits = c(-0.235, 0.235),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20),
                       limits = c(0, 20)) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', 'Log geographic range ' ~ (km^2))), 
         y = bquote(atop(italic(mu), ' '))))

# Global bird and mammal range and sigma
mcmc_preds_range_global_sigma <- as.data.frame(predict(rarity.range.birds.mammals.sigma,
                                                       type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_range_global_sigma$range <- b_m_ranges$range

(range.plot.birds.mammals.sigma <- ggplot() +
    geom_point(data = b_m_ranges, aes(x = log(range), y = sigma.2),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_range_global_sigma, aes(x = log(range), ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_range_global_sigma, aes(x = log(range), y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20),
                       limits = c(0, 20)) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', 'Log geographic range ' ~ (km^2))), 
         y = bquote(atop(italic(sigma^2), ' '))))

# Pop size global
mcmc_preds_popsize_global <- as.data.frame(predict(rarity.popsize.global,
                                                   type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_popsize_global$log.meanpop <- mus.popsize$log.meanpop

(popsize.plot.global <- ggplot() +
    geom_point(data = mus.popsize, aes(x = log.meanpop, y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_popsize_global, aes(x = log.meanpop,
                                                      ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_popsize_global, aes(x = log.meanpop, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                       limits = c(-0.235, 0.235),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20, 30),
                       limits = c(-10, 30)) +
    theme_LPI() + 
    labs(x = "\nLog population size", 
         y = bquote(atop(italic(mu), ' '))))

# Pop size global tau
mcmc_preds_popsize_global_sigma <- as.data.frame(predict(rarity.popsize.global.sigma,
                                                         type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_popsize_global_sigma$log.meanpop <- mus.popsize$log.meanpop

(popsize.plot.global.sigma <- ggplot() +
    geom_point(data = mus.popsize, aes(x = log.meanpop, y = sigma.2),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_popsize_global_sigma, aes(x = log.meanpop,
                                                            ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_popsize_global_sigma, aes(x = log.meanpop, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20, 30),
                       limits = c(-10, 30)) +
    theme_LPI() + 
    labs(x = "\nLog population size", 
         y = bquote(atop(italic(sigma), ' '))))

# Habitat specifity versus population change global
mcmc_preds_habspec_global <- as.data.frame(predict(rarity.hab.global,
                                                   type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_habspec_global$hab_spec <- mus.habspec$hab_spec

(habspec.plot.global <- ggplot() +
    geom_point(data = mus.habspec, aes(x = hab_spec, y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_habspec_global, aes(x = hab_spec,
                                                      ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_habspec_global, aes(x = hab_spec, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                       limits = c(-0.235, 0.235),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48),
                       limits = c(0, 48)) +
    theme_LPI() + 
    labs(x = "\nNumber of habitats", 
         y = bquote(atop(italic(mu), ' '))))

# Habitat specifity versus sigma global
mcmc_preds_habspec_global_sigma <- as.data.frame(predict(rarity.hab.global.sigma,
                                                         type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_habspec_global_sigma$hab_spec <- mus.habspec$hab_spec

(habspec.plot.global.sigma <- ggplot() +
    geom_point(data = mus.habspec, aes(x = hab_spec, y = sigma.2),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_habspec_global_sigma, aes(x = hab_spec,
                                                            ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_habspec_global_sigma, aes(x = hab_spec, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48),
                       limits = c(0, 48)) +
    theme_LPI() + 
    labs(x = "\nNumber of habitats", 
         y = bquote(atop(italic(sigma^2), ' '))))

# Adding effect size plots
all_outputs_5th_april <- read.csv("data/output/all_outputs_5th_April.csv")

rarity.effects.uk <- all_outputs_5th_april %>%
  filter(Model.name %in% c("Geographic range (UK) - mu",
                           "Geographic range (UK) - weighted",                        
                           "Geographic range (UK) - slope",
                           "Mean population size (UK) - mu",                          
                           "Mean population size (UK) - weighted",                    
                           "Mean population size (UK) - slope",
                           "Habitat specificity (UK) - mu",                           
                           "Habitat specificity (UK) - weighted",                     
                           "Habitat specificity (UK) - slope"),
         Variable != "(Intercept)" & Variable != "units")

rarity.effects.uk$Variable <- recode(rarity.effects.uk$Variable, 
                                     "log(km2_range)" = "log.range")

rarity.effects.uk$Variable <- recode(rarity.effects.uk$Variable, 
                                     "log(meanpop)" = "log.meanpop")

rarity.effects.uk <- rarity.effects.uk %>% separate(Model.name,
                                                    c("Name", "Model"), " - ")

rarity.effects.uk$Model <- recode(rarity.effects.uk$Model, 
                                  "mu" = "State-space")
rarity.effects.uk$Model <- recode(rarity.effects.uk$Model, 
                                  "weighted" = "Weighted")
rarity.effects.uk$Model <- recode(rarity.effects.uk$Model, 
                                  "slope" = "Linear")

# Checking sample sizes
length(unique(uk.mus.range$id))
length(unique(uk.mus.popsize$id))
length(unique(uk.mus.habspec$id))

rarity.effects.uk$Variable <- factor(rarity.effects.uk$Variable, 
                                     levels = c("log.range",
                                                "log.meanpop", 
                                                "hab_spec"),
                                     labels=c("Geographic range (412)",
                                              "Mean population size (175)", 
                                              "Habitat specificity (457)"))

# Standardising the effect sizes using the sd of the raw data for each model
rarity.effects.uk$sd <- NA
rarity.effects.uk[1, 11] <- sd(uk.mus.range$mu)
rarity.effects.uk[2, 11] <- sd(uk.mus.range$mu)
rarity.effects.uk[3, 11] <- sd(uk.slopes.ranges$slope)
rarity.effects.uk[4, 11] <- sd(uk.mus.popsize$mu)
rarity.effects.uk[5, 11] <- sd(uk.mus.popsize$mu)
rarity.effects.uk[6, 11] <- sd(uk.slopes.popsize$slope)
rarity.effects.uk[7, 11] <- sd(uk.mus.habspec$mu)
rarity.effects.uk[8, 11] <- sd(uk.mus.habspec$mu)
rarity.effects.uk[9, 11] <- sd(uk.slopes.habspec$slope)
rarity.effects.uk$std.mean <- rarity.effects.uk$Posterior.mean/rarity.effects.uk$sd
rarity.effects.uk$std.lower <- rarity.effects.uk$Lower.95..CI/rarity.effects.uk$sd
rarity.effects.uk$std.upper <- rarity.effects.uk$Upper.95..CI/rarity.effects.uk$sd

(eff_plot_uk1 <- ggplot(rarity.effects.uk, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, 
                        shape = Model, colour = Model), 
                    position = position_dodge(0.7)) + 
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    scale_y_continuous(breaks = c(0.06, 0.03, 0, -0.03, -0.06), 
                       limits = c(-0.07, 0.07),
                       labels = c("0.06", "0.03", "0", "-0.03", "-0.06")) +
    theme_LPI2() +
    theme(legend.position = c(0.8, 0.2),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Coefficient\n", title = "(g) Rarity effects (UK)\n") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    guides(shape = F, colour = F))

rarity.effects.global <- all_outputs_5th_april %>%
  filter(Model.name %in% c("Geographic range (birds) - mu",                           
                           "Geographic range (birds) - weighted",                     
                           "Geographic range (birds) - slope",
                           "Mean population size (global) - mu",                      
                           "Mean population size (global) - weighted",                
                           "Mean population size (global) - slope",
                           "Habitat specificity (global) - mu",                       
                           "Habitat specificity (global) - weighted",                 
                           "Habitat specificity (global) - slope"),
         Variable != "(Intercept)" & Variable != "units")

rarity.effects.global$Variable <- recode(rarity.effects.global$Variable, 
                                         "log(km2_range)" = "log.range")

rarity.effects.global$Variable <- recode(rarity.effects.global$Variable, 
                                         "log(range)" = "log.range")

rarity.effects.global$Variable <- recode(rarity.effects.global$Variable, 
                                         "log(meanpop)" = "log.meanpop")

rarity.effects.global <- rarity.effects.global %>% separate(Model.name,
                                                            c("Name", "Model"), " - ")

rarity.effects.global$Model <- recode(rarity.effects.global$Model, 
                                      "mu" = "State-space")
rarity.effects.global$Model <- recode(rarity.effects.global$Model, 
                                      "weighted" = "Weighted")
rarity.effects.global$Model <- recode(rarity.effects.global$Model, 
                                      "slope" = "Linear")

# Checking sample sizes
length(unique(bird_mu_ranges$id))
length(unique(mus.popsize$id))
length(unique(mus.habspec$id))

rarity.effects.global$Variable <- factor(rarity.effects.global$Variable, 
                                         levels = c("log.range",
                                                    "log.meanpop", 
                                                    "hab_spec"),
                                         labels=c("Geographic range (5381)",
                                                  "Mean population size (4310)", 
                                                  "Habitat specificity (7901)"))

# Standardising the effect sizes using the sd of the raw data for each model
rarity.effects.global$sd <- NA
rarity.effects.global[1, 11] <- sd(bird_mu_ranges$mu)
rarity.effects.global[2, 11] <- sd(bird_mu_ranges$mu)
rarity.effects.global[3, 11] <- sd(bird_slopes_ranges$slope)
rarity.effects.global[4, 11] <- sd(mus.popsize$mu)
rarity.effects.global[5, 11] <- sd(mus.popsize$mu)
rarity.effects.global[6, 11] <- sd(slopes.popsize$slope)
rarity.effects.global[7, 11] <- sd(mus.habspec$mu)
rarity.effects.global[8, 11] <- sd(mus.habspec$mu)
rarity.effects.global[9, 11] <- sd(slopes.habspec$slope)
rarity.effects.global$std.mean <- rarity.effects.global$Posterior.mean/rarity.effects.global$sd
rarity.effects.global$std.lower <- rarity.effects.global$Lower.95..CI/rarity.effects.global$sd
rarity.effects.global$std.upper <- rarity.effects.global$Upper.95..CI/rarity.effects.global$sd

(eff_plot_global <- ggplot(rarity.effects.global, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, 
                        shape = Model, colour = Model), 
                    position = position_dodge(0.7), size = 0.7) + 
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    scale_y_continuous(breaks = c(0.06, 0.03, 0, -0.03, -0.06), 
                       limits = c(-0.07, 0.07),
                       labels = c("0.06", "0.03", "0", "-0.03", "-0.06")) +
    theme_LPI() +
    theme(legend.position = c(0.85, 0.85),
          axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nMetric") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))


rarity.mu.panel <- grid.arrange(range.plot.birds, popsize.plot.global,
                                habspec.plot.global, eff_plot_global,
                                range.plot.birds.sigma,
                                popsize.plot.global.sigma,
                                habspec.plot.global.sigma, eff_plot_global_ci,
                                ncol = 4)

ggsave(rarity.mu.panel, filename = "figures/Fig3_rarity_effects.pdf", device = "pdf",
       dpi = 300, width = 20, height = 10)

# IUCN effects UK scale
IUCN.effects.uk <- all_outputs_5th_april %>%
  filter(Model.name %in% c("Red List status (UK) - mu",                           
                           "Red List status (UK) - weighted",                     
                           "Red List status (UK) - slope"),
         Variable != "units")

IUCN.effects.uk$Variable <- gsub("Red.List.status", "", 
                                 IUCN.effects.uk$Variable, fixed = TRUE)

IUCN.effects.uk <- IUCN.effects.uk %>% separate(Model.name,
                                                c("Name", "Model"), " - ")

IUCN.effects.uk$Model <- recode(IUCN.effects.uk$Model, 
                                "mu" = "State-space")
IUCN.effects.uk$Model <- recode(IUCN.effects.uk$Model, 
                                "weighted" = "Weighted")
IUCN.effects.uk$Model <- recode(IUCN.effects.uk$Model, 
                                "slope" = "Linear")

# Checking sample sizes
uk.mus.IUCN %>% group_by(Red.List.status) %>% tally

IUCN.effects.uk$Variable <- factor(IUCN.effects.uk$Variable, 
                                   levels = c("Least concern",
                                              "Near threatened", 
                                              "Vulnerable",
                                              "Endangered",
                                              "Critically endangered"),
                                   labels=c("Least concern (405)",
                                            "Near threatened (32)", 
                                            "Vulnerable (26)",
                                            "Endangered (1)",
                                            "Critically endangered (3)"))

# Standardising the effect sizes using the sd of the raw data for each IUCN category
uk.IUCN.sd <- uk.mus.IUCN %>% group_by(Red.List.status) %>%
  summarise(sd = sd(mu))
uk.IUCN.slope.sd <- uk.slopes.IUCN %>% group_by(Red.List.status) %>%
  summarise(sd = sd(slope))

IUCN.sd <- rbind(uk.IUCN.sd, uk.IUCN.sd, uk.IUCN.slope.sd)

IUCN.effects.uk <- cbind(IUCN.effects.uk, IUCN.sd)

IUCN.effects.uk$std.mean <- IUCN.effects.uk$Posterior.mean/IUCN.effects.uk$sd
IUCN.effects.uk$std.lower <- IUCN.effects.uk$Lower.95..CI/IUCN.effects.uk$sd
IUCN.effects.uk$std.upper <- IUCN.effects.uk$Upper.95..CI/IUCN.effects.uk$sd

# Because the sample size is very low for EN and CR, those can't be standardised
# For now putting the raw effects sizes back in
IUCN.effects.uk[4, 13] <- IUCN.effects.uk[4, 5]
IUCN.effects.uk[4, 14] <- IUCN.effects.uk[4, 6]
IUCN.effects.uk[4, 15] <- IUCN.effects.uk[4, 7]
IUCN.effects.uk[9, 13] <- IUCN.effects.uk[9, 5]
IUCN.effects.uk[9, 14] <- IUCN.effects.uk[9, 6]
IUCN.effects.uk[9, 15] <- IUCN.effects.uk[9, 7]
IUCN.effects.uk[14, 13] <- IUCN.effects.uk[14, 5]
IUCN.effects.uk[14, 14] <- IUCN.effects.uk[14, 6]
IUCN.effects.uk[14, 15] <- IUCN.effects.uk[14, 7]

IUCN.effects.uk[5, 13] <- IUCN.effects.uk[5, 5]
IUCN.effects.uk[5, 14] <- IUCN.effects.uk[5, 6]
IUCN.effects.uk[5, 15] <- IUCN.effects.uk[5, 7]
IUCN.effects.uk[10, 13] <- IUCN.effects.uk[10, 5]
IUCN.effects.uk[10, 14] <- IUCN.effects.uk[10, 6]
IUCN.effects.uk[10, 15] <- IUCN.effects.uk[10, 7]
IUCN.effects.uk[15, 13] <- IUCN.effects.uk[15, 5]
IUCN.effects.uk[15, 14] <- IUCN.effects.uk[15, 6]
IUCN.effects.uk[15, 15] <- IUCN.effects.uk[15, 7]

#IUCN.effects.uk <- IUCN.effects.uk %>% filter(Red.List.status != "Endangered")

(eff_plot_uk2 <- ggplot(IUCN.effects.uk, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Variable,
                        shape = Model), 
                    position = position_dodge(0.7)) + 
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ffa15b", "#ff000e")) +
    scale_y_continuous(breaks = c(0.8, 0.4, 0, -0.4, -0.8), 
                       limits = c(-0.84, 0.84),
                       labels = c("0.8", "0.4", "0", "-0.4", "-0.8")) +
    theme_LPI2() +
    theme(legend.position = c(0.3, 0.2),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Coefficient\n", title = "(i) Red List status effects (UK)\n") +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    guides(colour = F, shape = F))

# Global IUCN effects
IUCN.effects.global <- all_outputs_5th_april %>%
  filter(Model.name %in% c("Red List status (global) - mu",                           
                           "Red List status (global) - weighted",                     
                           "Red List status (global) - slope"),
         Variable != "units")

IUCN.effects.global$Variable <- gsub("Red.List.status", "", 
                                     IUCN.effects.global$Variable, fixed = TRUE)

IUCN.effects.global <- IUCN.effects.global %>% separate(Model.name,
                                                        c("Name", "Model"), " - ")

IUCN.effects.global$Model <- recode(IUCN.effects.global$Model, 
                                    "mu" = "State-space")
IUCN.effects.global$Model <- recode(IUCN.effects.global$Model, 
                                    "weighted" = "Weighted")
IUCN.effects.global$Model <- recode(IUCN.effects.global$Model, 
                                    "slope" = "Linear")

# Checking sample sizes
mus.IUCN %>% group_by(Red.List.status) %>% tally

IUCN.effects.global$Variable <- factor(IUCN.effects.global$Variable, 
                                       levels = c("Least concern",
                                                  "Near threatened", 
                                                  "Vulnerable",
                                                  "Endangered",
                                                  "Critically endangered"),
                                       labels=c("Least concern (6523)",
                                                "Near threatened (505)", 
                                                "Vulnerable (582)",
                                                "Endangered (304)",
                                                "Critically endangered (152)"))

# Standardising the effect sizes using the sd of the raw data for each IUCN category
global.IUCN.sd <- mus.IUCN %>% group_by(Red.List.status) %>%
  summarise(sd = sd(mu))
global.IUCN.slope.sd <- slopes.IUCN %>% group_by(Red.List.status) %>%
  summarise(sd = sd(slope))

IUCN.sd.global <- rbind(global.IUCN.sd, global.IUCN.sd, global.IUCN.slope.sd)

IUCN.effects.global <- cbind(IUCN.effects.global, IUCN.sd.global)

IUCN.effects.global$std.mean <- IUCN.effects.global$Posterior.mean/IUCN.effects.global$sd
IUCN.effects.global$std.lower <- IUCN.effects.global$Lower.95..CI/IUCN.effects.global$sd
IUCN.effects.global$std.upper <- IUCN.effects.global$Upper.95..CI/IUCN.effects.global$sd

(eff_plot_IUCN <- ggplot(IUCN.effects.global, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Variable,
                        shape = Model),
                    position = position_dodge(0.7), size = 0.7) + 
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ffa15b", "#ff000e")) +
    scale_y_continuous(breaks = c(0.8, 0.4, 0, -0.4, -0.8), 
                       limits = c(-0.84, 0.84),
                       labels = c("0.8", "0.4", "0", "-0.4", "-0.8")) +
    theme_LPI3() +
    theme(legend.position = c(0.25, 0.85),
          axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nRed List status") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") + 
    guides(colour = F))

ggsave(eff_plot_IUCN, filename = "figures/Fig4_IUCN_global_effects.pdf", device = "pdf",
       dpi = 300, width = 5, height = 5)

rarity.panel <- grid.arrange(range.plot, range.plot.birds,
                             pop.size.plot, pop.size.plot.global,
                             hab.spec.plot, hab.spec.plot.global,
                             eff_plot_uk1, eff_plot_global, 
                             eff_plot_uk2, eff_plot_IUCN, ncol = 2,
                             heights = c(0.184, 0.184, 0.184, 0.224, 0.224))

ggsave(rarity.panel, filename = "figures/updated/Figure1_updated.png", 
       width = 11, height = 28)

# ** IUCN figure ----
# Map

mus.IUCN$Red.List.status <- factor(mus.IUCN$Red.List.status,
                                   levels = c("Least concern",
                                              "Near threatened",
                                              "Vulnerable",
                                              "Endangered",
                                              "Critically endangered"),
                                   labels = c("Least concern",
                                              "Near threatened",
                                              "Vulnerable",
                                              "Endangered",
                                              "Critically endangered"))

mus.IUCN <- arrange(mus.IUCN, Red.List.status)

(world.map.IUCN <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(data = mus.IUCN, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                    colour = Red.List.status, size = duration),
               alpha = 0.6) +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ff8a32", "#ff000e")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Duration")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

ggsave(world.map.IUCN, filename = "figures/SI_IUCN_map.pdf", device = "pdf",
       dpi = 300, height = 8, width = 12)

(duration.IUCN <- ggplot(mus.IUCN, aes(duration, fill = Red.List.status, colour = Red.List.status)) +
    geom_line(stat = "density", aes(y = ..count..), size = 2, alpha = 0.7) +
    theme_LPI() +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ff8a32", "#ff000e")) +
    scale_fill_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                 "#ff8a32", "#ff000e")) +
    geom_vline(xintercept = mean(mus[mus.IUCN$Red.List.status == "Least concern",]$duration),
               size = 1, colour = "#4bb33c") +
    geom_vline(xintercept = mean(mus[mus.IUCN$Red.List.status == "Near threatened",]$duration),
               size = 1, colour = "#d5db00") +
    geom_vline(xintercept = mean(mus[mus.IUCN$Red.List.status == "Vulnerable",]$duration),
               size = 1, colour = "#ffef00") +
    geom_vline(xintercept = mean(mus[mus.IUCN$Red.List.status == "Endangered",]$duration),
               size = 1, colour = "#ff8a32") +
    geom_vline(xintercept = mean(mus[mus.IUCN$Red.List.status == "Critically endangered",]$duration),
               size = 1, colour = "#ff000e") +
    scale_y_continuous(limits = c(0, 300)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45)) +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")) +
    labs(x = "\nYears", y = "Number of time series\n") +
    geom_vline(xintercept = 5, colour = "grey30", size = 1, linetype = "dashed") +
    guides(colour = F, fill = F))

ggsave(duration.IUCN, filename = "figures/SI_IUCN_duration.pdf", device = "pdf",
       dpi = 300, height = 5, width = 5)


# Distributions
(rain_IUCN <- 
    ggplot(data = mus.IUCN, 
           aes(x = Red.List.status, y = mu, fill = Red.List.status)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mu, color = Red.List.status), 
               position = position_jitter(width = .15), size = .5, alpha = 0.3) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(mu), ' ')), x = "\nRed List status") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ff8a32", "#ff000e")) +
    scale_fill_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                 "#ff8a32", "#ff000e")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(-0.24, 0.24),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_IUCN, filename = "figures/Fig4_IUCN_mu.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)

mus.IUCN2 <- filter(mus.IUCN, sigma.2 < 1)

(rain_IUCN_sigma <- 
    ggplot(data = mus.IUCN2, 
           aes(x = Red.List.status, y = sigma.2, fill = Red.List.status)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = sigma.2, color = Red.List.status), 
               position = position_jitter(width = .15), size = .5, alpha = 0.3) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(sigma^2), ' ')), x = "\nRed List status") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ff8a32", "#ff000e")) +
    scale_fill_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                 "#ff8a32", "#ff000e")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.8),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_IUCN_sigma, filename = "figures/Fig4_IUCN_sigma.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)


# IUCN threats
(threat_plot <- ggplot(top_threats, aes(x = reorder(title, n), y = n, 
                                        colour = reorder(title, n), fill = reorder(title, n))) +
    geom_bar(stat = "identity") + 
    scale_colour_viridis_d(option = "magma", direction = -1,
                           begin = 0.3, end = 0.9) +
    scale_fill_viridis_d(option = "magma", direction = -1, 
                         begin = 0.3, end = 0.9) +
    coord_flip() +
    theme_LPI() +
    guides(colour = FALSE, fill = FALSE) +
    labs(x = "Threat\n", y = "\nNumber of species"))

ggsave(threat_plot, filename = "figures/threats.pdf", device = "pdf", dpi = 300,
       height = 4.5, width = 18)

(threat_mu_plot <- ggplot(top_threats_mus, aes(x = mu, y = title, fill = title)) +
    geom_density_ridges() +
    theme_LPI() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey10") +
    labs(x = bquote(atop(' ', italic(mu))), y = "Threat\n") +
    scale_fill_viridis_d(option = "magma", direction = -1, 
                         begin = 0.3, end = 0.8) +
    guides(fill = FALSE))

ggsave(threat_mu_plot, filename = "figures/threats_mu.pdf", device = "pdf", dpi = 300,
       height = 5, width = 13)

# Calculate predictions
mcmc_preds_n_threats <- as.data.frame(predict(threat.n.m,
                                              type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_n_threats$n <- threats_sum_species_mus$n


(threat_n_mu <- ggplot() +
    geom_point(data = threats_sum_species_mus, aes(x = n, y = mu), colour = "grey70") +
    geom_ribbon(data = mcmc_preds_n_threats, aes(x = n, ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_n_threats, aes(x = n, y = fit), 
              colour = "black", size = 1) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    labs(x = "\nNumber of threats", y = bquote(atop(italic(mu), ' '))) +
    scale_y_continuous(limits = c(-0.24, 0.24),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    theme_LPI())

ggsave(threat_n_mu, filename = "figures/Fig4threat_n_mu.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)



# IUCN map
(world.IUCN.map <- ggplot(mus.IUCN, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                        colour = Red.List.status)) +
    borders("worldHires",  colour = "gray75", fill = "gray75", size = 0.3) +
    theme_map() +
    geom_point(alpha = 0.6, size = 2) +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ff8a32", "#ff000e")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "italic"),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size = 20, vjust = 1, hjust = 0)) +
    labs(colour = "Red\n", 
         title = "(a) Distribution of records"))


# IUCN map
# Reordering levels so that the green doesn't plot over  everything
mus.IUCN$Red.List.status <- factor(mus.IUCN$Red.List.status, 
                                   levels = c("Critically endangered",
                                              "Endangered",
                                              "Near threatened", 
                                              "Vulnerable",
                                              "Least concern"),
                                   labels=c("Critically endangered",
                                            "Endangered",
                                            "Near threatened", 
                                            "Vulnerable",
                                            "Least concern"))

mus.IUCN <- arrange(mus.IUCN, desc(Red.List.status))

(world.IUCN.map <- ggplot(mus.IUCN, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                        colour = Red.List.status)) +
    borders("worldHires",  colour = "gray75", fill = "gray75", size = 0.3) +
    theme_map() +
    geom_point(alpha = 0.6, size = 2) +
    scale_colour_manual(values = c("#ff000e", "#ff8a32", "#d5db00",
                                   "#ffef00", "#4bb33c")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16, face = "italic"),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size = 20, vjust = 1, hjust = 0)) +
    labs(colour = "Red List status  "))

IUCN.effects.panel <- grid.arrange(eff_plot_uk2, 
                                   eff_plot_IUCN, ncol = 2)
ggsave(IUCN.effects.panel, filename = "figures/updated/Figure2_effects2.png",
       height = 7.1, width = 11.1)

IUCN.panel <- grid.arrange(world.IUCN.map, IUCN.effects.panel,
                           nrow = 2, heights = c(0.47, 0.53))

ggsave(world.IUCN.map, filename = "figures/updated/SI_IUCN_map.png",
       height = 6, width = 10)

ggsave(IUCN.panel, filename = "figures/updated/Figure2_IUCN.png", 
       width = 11.1, height = 13.1)

# ** Lat/biome figure ----
# Global map of variation in population change
# rearranging values, so that they are all visible on map
mus$abs.mu<- abs(mus$mu)
mus <- arrange(mus, abs.mu)

world <- map_data("world")
world <- world[world$region != "Antarctica",]

map.coords <- mus %>% dplyr::select(Decimal.Longitude, Decimal.Latitude, 
                                    system, End) %>% distinct()

(world.map <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(data = mus, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                               colour = mu, size = duration),
               alpha = 0.6) +
    scale_colour_viridis() +
    #scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Duration")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

ggsave(world.map, filename = "figures/mu_map.pdf", device = "pdf",
       dpi = 600, height = 8, width = 12)

mus <- arrange(mus, sigma.2)

(world.map.sigma <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(data = mus, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                               colour = sigma.2, size = duration),
               alpha = 0.6) +
    scale_colour_viridis() +
    #scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Duration")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

ggsave(world.map.sigma, filename = "figures/SI_sigma_map.pdf", device = "pdf",
       dpi = 600, height = 8, width = 12)

slopes <- arrange(slopes, abs(slope))

(world.map.slopes <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(data = slopes, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                  colour = slope, size = lengthyear),
               alpha = 0.6) +
    scale_colour_viridis() +
    #scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Duration")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

ggsave(world.map.slopes, filename = "figures/SI_slopes_map.pdf", device = "pdf",
       dpi = 600, height = 8, width = 12)

# Panel
SI_map_panel <- grid.arrange(world.map, world.map.slopes, world.map.sigma, nrow = 3)
ggsave(SI_map_panel, filename = "figures/SI_maps.pdf", device = "pdf",
       dpi = 300, height = 16, width = 12)


# Latitudinal pop change
(latitude.plot <- ggplot(mus, aes(x = mu, y = Decimal.Latitude)) +
    geom_point(aes(colour = mu), alpha = 0.7) +
    #geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_viridis() +
    scale_x_continuous(breaks = c(0.2, 0, -0.2), 
                       labels = c("0.2", "0", "-0.2")) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', italic(mu))), y = "Latitude\n") +
    scale_y_continuous(limits = c(-80, 80), 
                       breaks = c(-80, -40, 0, 40, 80)) +
    guides(colour = F))

ggsave(latitude.plot, filename = "figures/Fig1_mu.pdf", device = "pdf",
       dpi = 300, height = 8, width = 4)

mus <- arrange(mus, sigma.2)
mus.sigma <- mus %>% filter(sigma.2 < 1)
# Latitudinal fluctuations sigma
(latitude.plot.sigma <- ggplot(mus.sigma, aes(x = sigma.2, y = Decimal.Latitude)) +
    geom_point(aes(colour = sigma.2), alpha = 0.7) +
    # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    #geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_viridis() +
    scale_x_reverse(limits = c(0.65, 0),
                    breaks = c(0, 0.3, 0.6), 
                    labels = c("0", "0.3", "0.6")) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', italic(sigma))), y = "Latitude\n") +
    scale_y_continuous(limits = c(-80, 80), position = "right", 
                       breaks = c(-80, -40, 0, 40, 80)) + 
    guides(colour = F))

ggsave(latitude.plot.sigma, filename = "figures/Fig1_sigma.pdf", device = "pdf",
       dpi = 300, height = 8, width = 4)

# Biome effects
biome.effects.global <- all_outputs_5th_april %>%
  filter(Model.name %in% c("Biome (global) - mu",                           
                           "Biome (global) - weighted",                     
                           "Biome (global) - slope"),
         Variable != "units")

biome.effects.global$Variable <- gsub("biome", "", 
                                      biome.effects.global$Variable, fixed = TRUE)

biome.effects.global <- biome.effects.global %>% separate(Model.name,
                                                          c("Name", "Model"), " - ")

biome.effects.global$Model <- recode(biome.effects.global$Model, 
                                     "mu" = "State-space")
biome.effects.global$Model <- recode(biome.effects.global$Model, 
                                     "weighted" = "Weighted")
biome.effects.global$Model <- recode(biome.effects.global$Model, 
                                     "slope" = "Linear")

# Checking sample sizes
mus.biome %>% group_by(biome) %>% tally

biome.effects.global$Variable <- factor(biome.effects.global$Variable, 
                                        levels = c("Boreal forests/taiga",
                                                   "Deserts and xeric shrublands",
                                                   "Trop. and subtrop. grasslands savannas and shrublands",
                                                   "Large lakes",
                                                   "Mediterranean forests woodlands and scrub",
                                                   "Montane freshwaters",
                                                   "Montane grasslands and shrublands",
                                                   "Polar freshwaters",
                                                   "Polar seas",
                                                   "Temperate forests",
                                                   "Temperate wetlands and rivers",
                                                   "Temperate grasslands savannas and shrublands",
                                                   "Tropical wetlands and rivers",
                                                   "Tropical and subtropical forests",
                                                   "Tropical coral",
                                                   "Tundra",
                                                   "Xeric freshwaters and endorheic basins"),
                                        labels=c("Boreal forests (1289)",
                                                 "Deserts (46)",
                                                 "Tropical savannas (216)",
                                                 "Large lakes (225)",
                                                 "Mediterranean forests (225)",
                                                 "Montane freshwaters (23)",
                                                 "Montane grasslands (36)",
                                                 "Polar freshwaters (362)",
                                                 "Polar seas (26)",
                                                 "Temperate forests (1785)",
                                                 "Temperate wetlands (1712)",
                                                 "Temperate grasslands (321)",
                                                 "Tropical wetlands (184)",
                                                 "Tropical forests (178)",
                                                 "Tropical coral (102)",
                                                 "Tundra (201)",
                                                 "Xeric freshwaters (52)"))

# Standardising the effect sizes using the sd of the raw data for each biome category
global.biome.sd <- mus.biome %>% group_by(biome) %>%
  summarise(sd = sd(mu))

global.biome.slope.sd <- slopes.biome %>% group_by(biome) %>%
  summarise(sd = sd(slope))

biome.sd.global <- rbind(global.biome.sd, global.biome.sd, global.biome.slope.sd)

biome.effects.global <- cbind(biome.effects.global, biome.sd.global)

biome.effects.global$std.mean <- biome.effects.global$Posterior.mean/biome.effects.global$sd
biome.effects.global$std.lower <- biome.effects.global$Lower.95..CI/biome.effects.global$sd
biome.effects.global$std.upper <- biome.effects.global$Upper.95..CI/biome.effects.global$sd

(eff_plot_biome <- ggplot(biome.effects.global, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Model,
                        shape = Model),
                    position = position_dodge(0.7), size = 1) + 
    scale_colour_manual(values = c("grey70", "grey40", "black")) +
    scale_y_continuous(breaks = c(1, 0.5, 0, -0.5, -1), 
                       limits = c(-1.18, 1.18),
                       labels = c("1.0", "0.5", "0", "-0.5", "-1.0")) +
    theme_LPI2() +
    theme(legend.position = c(0.06, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Standardised effect size\n") +
    geom_hline(yintercept = 0, linetype = "dashed"))

map.lat <- grid.arrange(world.mu.map, latitude.plot, widths = c(0.66, 0.34))

map.biome <- grid.arrange(map.lat, eff_plot_biome, nrow = 2)

ggsave(eff_plot_biome, filename = "figures/Figure1_biome_effects.pdf", device = "pdf", dpi = 300,
       width = 18, height = 6)

# ** System/taxa figure ----
# System map
(world.map.system <- ggplot() +
   geom_map(map = world, data = world,
            aes(long, lat, map_id = region), 
            color = "gray80", fill = "gray80", size = 0.3) +
   coord_proj("+proj=wintri") +
   theme_map() +
   geom_point(data = mus, aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                              colour = system, size = duration),
              alpha = 0.6) +
   scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
   scale_y_continuous(limits = c(-80, 80)) +
   guides(size = guide_legend(title = "Duration")) +
   theme(legend.position = "bottom",
         legend.title = element_text(size = 16),
         legend.text = element_text(size = 10),
         legend.justification = "top"))

(rain_system <- 
    ggplot(data = mus, 
           aes(x = system, y = mu, fill = system)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mu, color = system), 
               position = position_jitter(width = .15), size = .5, alpha = 0.1) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(mu), ' ')), x = "\nSystem") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_fill_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(-0.24, 0.24),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_system, filename = "figures/Fig1_system_mu.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)

(rain_system_sigma <- 
    ggplot(data = mus, 
           aes(x = system, y = sigma.2, fill = system)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = sigma.2, color = system), 
               position = position_jitter(width = .15), size = .5, alpha = 0.05) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(sigma^2), ' ')), x = "\nSystem") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_fill_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.8),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_system_sigma, filename = "figures/Fig1_system_mu_sigma.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)

mu.class$Class <- factor(mu.class$Class,
                         levels = c("Elasmobranchii", "Actinopterygii",
                                    "Amphibia", "Aves", "Mammalia", "Reptilia"),
                         labels = c("Elasmobranchii", "Actinopterygii",
                                    "Amphibia", "Aves", "Mammalia", "Reptilia"))

# Taxa map
(world.map.taxa <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray80", fill = "gray80", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(data = mu.class[mu.class$Class == "Aves",], aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                                               size = duration),
               alpha = 0.6, colour = "#d8b70a") +
    geom_point(data = mu.class[mu.class$Class != "Aves",], aes(x = Decimal.Longitude, y = Decimal.Latitude, 
                                                               colour = Class, size = duration),
               alpha = 0.6) +
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b", "#972d15", "#dc863b")) +
    scale_y_continuous(limits = c(-80, 80)) +
    guides(size = guide_legend(title = "Duration")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 10),
          legend.justification = "top"))

map.panel <- grid.arrange(world.map.system, world.map.taxa, nrow = 2)
ggsave(map.panel, filename = "figures/Figure_SI_maps2.pdf", device = "pdf",
       height = 16, width = 12, dpi = 300)

(rain_taxa <- 
    ggplot(data = mu.class, 
           aes(x = Class, y = mu, fill = Class)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = mu, color = Class), 
               position = position_jitter(width = .15), size = .5, alpha = 0.2) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(mu), ' ')), x = "\nTaxa") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    scale_fill_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(-0.24, 0.24),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_taxa, filename = "figures/Fig2_taxa_mu.pdf", device = "pdf", dpi = 600,
       height = 5, width = 5)

(rain_taxa_sigma <- 
    ggplot(data = mu.class, 
           aes(x = Class, y = sigma.2, fill = Class)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = sigma.2, color = Class), 
               position = position_jitter(width = .15), size = .5, alpha = 0.05) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = bquote(atop(italic(sigma^2), ' ')), x = "\nTaxa") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    scale_fill_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.8),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_taxa_sigma, filename = "figures/Fig2_taxa_mu_sigma.pdf", device = "pdf", dpi = 600,
       height = 5, width = 5)

(duration.taxa <- ggplot(mu.class, aes(duration, fill = Class, colour = Class)) +
    geom_line(stat = "density", aes(y = ..count..), size = 2) +
    theme_LPI() +
    scale_fill_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b", "#d8b70a", "#972d15", "#dc863b")) +
    scale_y_continuous(limits = c(0, 450)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45)) +
    geom_vline(xintercept = mean(mus[mus$Class == "Actinopterygii",]$duration), size = 1, colour = "#262F46") +
    geom_vline(xintercept = mean(mus[mus$Class == "Elasmobranchii",]$duration), size = 1, colour = "#96CDCD") +
    geom_vline(xintercept = mean(mus[mus$Class == "Amphibia",]$duration), size = 1, colour = "#02401b") +
    geom_vline(xintercept = mean(mus[mus$Class == "Aves",]$duration), size = 1, colour = "#d8b70a") +
    geom_vline(xintercept = mean(mus[mus$Class == "Mammalia",]$duration), size = 1, colour = "#972d15") +
    geom_vline(xintercept = mean(mus[mus$Class == "Reptilia",]$duration), size = 1, colour = "#dc863b") +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")) +
    labs(x = "\nYears", y = "Number of time series\n", 
         title = expression(paste(bold("c   \n"), "Monitoring duration across taxa\n"))) +
    geom_vline(xintercept = 5, colour = "grey30", size = 1, linetype = "dashed") +
    guides(colour = F, fill = F))

(duration.system <- ggplot(mus, aes(duration, fill = system, colour = system)) +
    #geom_histogram(binwidth = 1, alpha = 0.2, position = "identity") +
    geom_line(stat = "density", aes(y = ..count..), size = 2) +
    theme_LPI() +
    scale_fill_manual(values = c("#a2a475", "#273046", "#0b775e")) +
    scale_colour_manual(values = c("#a2a475", "#273046", "#0b775e")) +
    scale_y_continuous(limits = c(0, 450)) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45)) +
    geom_vline(xintercept = mean(mus[mus$system == "Freshwater",]$duration), size = 1, colour = "#0b775e") +
    geom_vline(xintercept = mean(mus[mus$system == "Marine",]$duration), size = 1, colour = "#273046") +
    geom_vline(xintercept = mean(mus[mus$system == "Terrestrial",]$duration), size = 1, colour = "#a2a475") +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          plot.title = element_text(hjust = 0)) +
    labs(x = "\nYears", y = "Number of time series\n", 
         title = expression(paste(bold("b   \n"), "Monitoring duration across realms\n"))) +
    geom_vline(xintercept = 5, colour = "grey30", size = 1, linetype = "dashed") +
    guides(colour = F, fill = F))

(duration.all <- ggplot(mus, aes(duration)) +
    #geom_histogram(binwidth = 1, alpha = 0.2, position = "identity") +
    geom_line(stat = "density", aes(y = ..count..), size = 2) +
    theme_LPI() +
    scale_fill_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_y_continuous(limits = c(0, 450)) +
    geom_vline(xintercept = mean(mus$duration), size = 1) +
    scale_x_continuous(breaks = c(5, 15, 25, 35, 45)) +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          plot.title = element_text(hjust = 0.99)) +
    labs(x = "\nYears", y = "Number of time series\n", 
         title = expression(paste(bold("a   \n"), "Monitoring duration across all populations\n"))) +
    geom_vline(xintercept = 5, colour = "grey30", size = 1, linetype = "dashed") +
    guides(colour = F, fill = F))

duration.panel <- grid.arrange(duration.all, duration.system, duration.taxa,
                               ncol = 3)
ggsave(duration.panel, filename = "figures/SI_duration.pdf", device = "pdf",
       dpi = 300, width = 15, height = 5)

# Effect sizes
taxa.effects.global <- global_outputs18thFeb %>%
  filter(Model.name %in% c("Taxa - mu",                           
                           "Taxa - weighted",                     
                           "Taxa - slope"),
         Variable != "units")

taxa.effects.global$Variable <- gsub("Class", "", 
                                     taxa.effects.global$Variable, fixed = TRUE)

taxa.effects.global <- taxa.effects.global %>% separate(Model.name,
                                                        c("Name", "Model"), " - ")

taxa.effects.global$Model <- recode(taxa.effects.global$Model, 
                                    "mu" = "State-space")
taxa.effects.global$Model <- recode(taxa.effects.global$Model, 
                                    "weighted" = "Weighted")
taxa.effects.global$Model <- recode(taxa.effects.global$Model, 
                                    "slope" = "Linear")

# Checking sample sizes
mu.class %>% group_by(Class) %>% tally

taxa.effects.global$Variable <- factor(taxa.effects.global$Variable, 
                                       levels = c("Elasmobranchii",
                                                  "Actinopterygii",
                                                  "Amphibia", 
                                                  "Aves",
                                                  "Mammalia",
                                                  "Reptilia"),
                                       labels=c("Elasmobranchii (127)",
                                                "Actinopterygii (1626)",
                                                "Amphibia (193)", 
                                                "Aves (5854)",
                                                "Mammalia (1158)",
                                                "Reptilia (322)"))

# Standardising the effect sizes using the sd of the raw data for each taxa category
global.taxa.sd <- mu.class %>% group_by(Class) %>%
  summarise(sd = sd(mu))
global.taxa.slope.sd <- slope.class %>% group_by(Class) %>%
  summarise(sd = sd(slope))

taxa.sd.global <- rbind(global.taxa.sd, global.taxa.sd, global.taxa.slope.sd)

taxa.effects.global <- cbind(taxa.effects.global, taxa.sd.global)

taxa.effects.global$std.mean <- taxa.effects.global$Posterior.mean/taxa.effects.global$sd
taxa.effects.global$std.lower <- taxa.effects.global$Lower.95..CI/taxa.effects.global$sd
taxa.effects.global$std.upper <- taxa.effects.global$Upper.95..CI/taxa.effects.global$sd


(eff_plot_taxa <- ggplot(taxa.effects.global, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, 
                        shape = Model, colour = Variable),
                    position = position_dodge(0.7), size = 0.7) + 
    scale_y_continuous(breaks = c(0.6, 0.3, 0, -0.3, -0.6), 
                       limits = c(-0.69, 0.69),
                       labels = c("0.6", "0.3", "0", "-0.3", "-0.6")) +
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b",
                                   "#d8b70a", "#972d15", "#dc863b")) +
    theme_LPI3() +
    theme(legend.position = c(0.8, 0.27)) +
    theme(axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nTaxa") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    guides(colour = F))

ggsave(eff_plot_taxa, filename = "figures/Fig2_taxa_effects.pdf", device = "pdf",
       dpi = 600, width = 5, height = 5)

# System effect sizes
system.effects.global <- all_outputs_5th_april %>%
  filter(Model.name %in% c("System (global) - mu",                           
                           "System (global) - weighted",                     
                           "System (global) - slope"),
         Variable != "units")

system.effects.global$Variable <- gsub("system", "", 
                                       system.effects.global$Variable, fixed = TRUE)

system.effects.global <- system.effects.global %>% separate(Model.name,
                                                            c("Name", "Model"), " - ")

system.effects.global$Model <- recode(system.effects.global$Model, 
                                      "mu" = "State-space")
system.effects.global$Model <- recode(system.effects.global$Model, 
                                      "weighted" = "Weighted")
system.effects.global$Model <- recode(system.effects.global$Model, 
                                      "slope" = "Linear")

# Checking sample sizes
mus %>% group_by(system) %>% tally

system.effects.global$Variable <- factor(system.effects.global$Variable, 
                                         levels = c("Freshwater", 
                                                    "Marine",
                                                    "Terrestrial"),
                                         labels=c("Freshwater (2570)",
                                                  "Marine (2418)", 
                                                  "Terrestrial (4298)"))

# Standardising the effect sizes using the sd of the raw data for each system category
global.system.sd <- mus %>% group_by(system) %>%
  summarise(sd = sd(mu))
global.system.slope.sd <- slopes %>% group_by(system) %>%
  summarise(sd = sd(slope))

system.sd.global <- rbind(global.system.sd, global.system.sd, global.system.slope.sd)

system.effects.global <- cbind(system.effects.global, system.sd.global)

system.effects.global$std.mean <- system.effects.global$Posterior.mean/system.effects.global$sd
system.effects.global$std.lower <- system.effects.global$Lower.95..CI/system.effects.global$sd
system.effects.global$std.upper <- system.effects.global$Upper.95..CI/system.effects.global$sd

(eff_plot_system <- ggplot(system.effects.global, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, 
                        shape = Model, colour = Variable),
                    position = position_dodge(0.7), size = 0.7) + 
    scale_y_continuous(breaks = c(0.4, 0.2, 0, -0.2, -0.4), 
                       limits = c(-0.4, 0.4),
                       labels = c("0.4", "0.2", "0", "-0.2", "-0.4")) +
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    theme_LPI3() +
    theme(legend.position = c(0.3, 0.2)) +
    theme(axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nRealm") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    guides(colour = F))

ggsave(eff_plot_system, filename = "figures/Fig1_system_effects.pdf", device = "pdf",
       dpi = 300, width = 5, height = 5)

global.effects <- grid.arrange(eff_plot_system, eff_plot_taxa, ncol = 2)
duration.panel <- grid.arrange(duration.system, duration.taxa, ncol = 2)
duration.effects.panel <- grid.arrange(duration.panel, global.effects, 
                                       heights = c(0.45, 0.55))

ggsave(duration.effects.panel, filename = "figures/updated/Figure4_effects.png",
       height = 11, width = 10)

# ** Phylogeny ----
birdtree <- read.nexus("output.nex")

bird_mu_ranges2$species <- as.factor(bird_mu_ranges2$species)

bird_mus <- bird_mu_ranges2 %>% dplyr::group_by(species) %>%
  dplyr::summarise(mean.mu = mean(mu), mean.sigma = mean(sigma.2))
bird_mus <- bird_mus %>% filter(mean.sigma < 1)

# Format bird pop change so that it lines up with species names in tree
bird_mus$species <- gsub(" ", "_", bird_mus$species, fixed = TRUE)
bird_mus$node <- NA
tipnode <- seq_along(birdtree[[1]]$tip.label)
names(tipnode) <- birdtree[[1]]$tip.label
bird_mus$node <- tipnode[bird_mus$species] ## convert the tip label to tip node number
i <- is.na(bird_mus$node)
bird_mus$node[i] = bird_mus$species[i] ## fill in internal node number

(tree <- ggtree(birdtree[[1]], layout = "circular") +
    geom_tiplab(size=0))

matrix <- as.matrix(bird_mus[, "mean.mu"])
matrix.names <- bird_mus$species
rownames(matrix) <- matrix.names

(heat.tree1 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(-0.25, 0.25))+
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree1, filename = "figures/Fig2_bird_mu_tree.pdf", device = "pdf",
       dpi = 300, height = 10, width = 10)

matrix <- as.matrix(bird_mus[, "mean.sigma"])
matrix.names <- bird_mus$species
rownames(matrix) <- matrix.names


(heat.tree2 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(0, 0.25))  +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree2, filename = "figures/Fig2_bird_sigma_tree.pdf", device = "pdf",
       dpi = 300, height = 10, width = 10)


(phylo_density <- ggplot(bird_mus) + 
    geom_density(aes(x = mean.mu), fill = "#d8b70a", colour = "#d8b70a",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    labs(x = "\n\nMean population trend", 
         y = "Density\n",
         title = "(e) Population trends\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density, filename = "figures/updated/bird_density.png", 
       height = 6, width = 5)

bird_mus2 <- bird_mu_ranges2
bird_mus2$CI <- (bird_mus2$uCI - bird_mus2$lCI)/2
bird_mus2$CI.scaled <- scales::rescale(bird_mus2$CI, to = c(0,1))
bird_mus2 <- bird_mus2 %>% group_by(species) %>%
  summarise(mean.CI.scaled = mean(CI.scaled))

bird_mus2$species <- gsub(" ", "_", bird_mus2$species, fixed = TRUE)

bird_mus2$node <- NA
tipnode <- seq_along(birdtree[[1]]$tip.label)
names(tipnode) <- birdtree[[1]]$tip.label
bird_mus2$node <- tipnode[bird_mus2$species] ## convert the tip label to tip node number
i <- is.na(bird_mus2$node)
bird_mus2$node[i] = bird_mus2$species[i] ## fill in internal node number

(tree <- ggtree(birdtree[[1]], layout = "circular") +
    geom_tiplab(size = 0))

matrix <- as.matrix(bird_mus2[, "mean.CI.scaled"])
matrix.names <- bird_mus2$species
rownames(matrix) <- matrix.names

(heat.tree2 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis() +
    labs(title = "(e) Population fluctuations") +
    theme(plot.title = element_text(size = 20, vjust = 1, hjust = 0)))

ggsave(heat.tree2, filename = "figures/updated/bird_tree4.png", height = 5, width = 5)


(phylo_density_fl <- ggplot(bird_mus2) + 
    geom_density(aes(x = mean.CI.scaled), fill = "#d8b70a", colour = "#d8b70a",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6),
                       labels = c("0", "0.2", "0.4", "0.6")) +
    labs(x = "\n\nMean population fluctuations", 
         y = "Density\n",
         title = "(k) Population fluctuations\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density_fl, filename = "figures/updated/bird_density_fl.png", 
       height = 5, width = 6)

# Amphibian tree
amphtree <- read.nexus("output_amp.nex")

amp_mus <- mus %>% filter(Class == "Amphibia") %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean.mu = mean(mu), mean.sigma = mean(sigma.2))

amp_mus$species <- gsub(" ", "_", amp_mus$species, fixed = TRUE)

amp_mus$node <- NA
tipnode <- seq_along(amphtree[[1]]$tip.label)
names(tipnode) <- amphtree[[1]]$tip.label
amp_mus$node <- tipnode[amp_mus$species] ## convert the tip label to tip node number
i <- is.na(amp_mus$node)
amp_mus$node[i] = amp_mus$species[i] ## fill in internal node number

(tree <- ggtree(amphtree[[1]], layout = "circular") +
    geom_tiplab(size=0))

matrix <- as.matrix(amp_mus[, "mean.mu"])
matrix.names <- amp_mus$species
rownames(matrix) <- matrix.names

(heat.tree3 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(-0.25, 0.25))+
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree3, filename = "figures/Fig2_amp_mu_tree.pdf", device = "pdf",
       dpi = 300, height = 10, width = 10)

matrix <- as.matrix(amp_mus[, "mean.sigma"])
matrix.names <- amp_mus$species
rownames(matrix) <- matrix.names

(heat.tree4 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(0, 0.25)) +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree4, filename = "figures/Fig2_amp_sigma_tree.pdf", device = "pdf",
       dpi = 300, height = 10, width = 10)


(phylo_density2 <- ggplot(amp_mus) + 
    geom_density(aes(x = mean.mu), fill = "#02401b", colour = "#02401b",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(limits = c(-0.2, 0.2),
                       breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    labs(x = "\n\nMean population trend", 
         y = "Density\n",
         title = "(d) Population trends\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density2, filename = "figures/updated/amp_density.png", 
       height = 6, width = 5)

amp_mus2 <- mus %>% filter(Class == "Amphibia")
amp_mus2$CI <- (amp_mus2$uCI - amp_mus2$lCI)/2
amp_mus2$CI.scaled <- scales::rescale(amp_mus2$CI, to = c(0,1))
amp_mus2 <- amp_mus2 %>% group_by(species) %>%
  summarise(mean.CI.scaled = mean(CI.scaled))

amp_mus2$species <- gsub(" ", "_", amp_mus2$species, fixed = TRUE)

amp_mus2$node <- NA
tipnode <- seq_along(amphtree[[1]]$tip.label)
names(tipnode) <- amphtree[[1]]$tip.label
amp_mus2$node <- tipnode[amp_mus2$species] ## convert the tip label to tip node number
i <- is.na(amp_mus2$node)
amp_mus2$node[i] = amp_mus2$species[i] ## fill in internal node number

(tree <- ggtree(amphtree[[1]], layout = "circular") +
    geom_tiplab(size=0))

matrix <- as.matrix(amp_mus2[, "mean.CI.scaled"])
matrix.names <- amp_mus2$species
rownames(matrix) <- matrix.names

(heat.tree4 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis() +
    labs(title = "(d) Population fluctuations") +
    theme(plot.title = element_text(size = 20, vjust = 1, hjust = 0)))

(phylo_density_amp_fl <- ggplot(amp_mus2) + 
    geom_density(aes(x = mean.CI.scaled), fill = "#02401b", colour = "#02401b",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6),
                       labels = c("0", "0.2", "0.4", "0.6")) +
    labs(x = "\n\nMean population fluctuations", 
         y = "Density\n",
         title = "(j) Population fluctuations\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density_amp_fl, filename = "figures/updated/amp_density2.png", 
       height = 6, width = 5)

# Reptiles
repttree <- read.nexus("output_rept.nex")

rept_mus <- mus %>% filter(Class == "Reptilia") %>% group_by(species) %>%
  dplyr::summarise(mean.mu = mean(mu), mean.sigma = mean(sigma.2))

rept_mus$species <- gsub(" ", "_", rept_mus$species, fixed = TRUE)

rept_mus$node <- NA
tipnode <- seq_along(repttree[[1]]$tip.label)
names(tipnode) <- repttree[[1]]$tip.label
rept_mus$node <- tipnode[rept_mus$species] ## convert the tip label to tip node number
i <- is.na(rept_mus$node)
rept_mus$node[i] = rept_mus$species[i] ## fill in internal node number

(tree <- ggtree(repttree[[1]], layout = "circular") +
    geom_tiplab(size=0))

matrix <- as.matrix(rept_mus[, "mean.mu"])
matrix.names <- rept_mus$species
rownames(matrix) <- matrix.names

(heat.tree5 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(-0.25, 0.25))+
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree5, filename = "figures/Fig2_rept_mu_tree.pdf", device = "pdf", dpi = 300,
       height = 10, width = 10)

matrix <- as.matrix(rept_mus[, "mean.sigma"])
matrix.names <- rept_mus$species
rownames(matrix) <- matrix.names

(heat.tree6 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis(limits = c(0, 0.25))+
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1.0, 'cm'),
          legend.text = element_text(margin = margin(t = 10))))

ggsave(heat.tree6, filename = "figures/Fig2_rept_sigma_tree.pdf", device = "pdf", dpi = 300,
       height = 10, width = 10)


(phylo_density3 <- ggplot(rept_mus) + 
    geom_density(aes(x = mean.mu), fill = "#dc863b", colour = "#dc863b",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    labs(x = "\n\nMean population trend", 
         y = "Density\n",
         title = "(f) Population trends\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density3, filename = "figures/updated/rept_density.png", 
       height = 6, width = 5)

rept_mus2 <- mus %>% filter(Class == "Reptilia")
rept_mus2$CI <- (rept_mus2$uCI - rept_mus2$lCI)/2
rept_mus2$CI.scaled <- scales::rescale(rept_mus2$CI, to = c(0,1))
rept_mus2 <- rept_mus2 %>% group_by(species) %>%
  summarise(mean.CI.scaled = mean(CI.scaled))

rept_mus2$species <- gsub(" ", "_", rept_mus2$species, fixed = TRUE)

rept_mus2$node <- NA
tipnode <- seq_along(repttree[[1]]$tip.label)
names(tipnode) <- repttree[[1]]$tip.label
rept_mus2$node <- tipnode[rept_mus2$species] ## convert the tip label to tip node number
i <- is.na(rept_mus2$node)
rept_mus2$node[i] = rept_mus2$species[i] ## fill in internal node number

(tree <- ggtree(repttree[[1]], layout = "circular") +
    geom_tiplab(size = 0))

matrix <- as.matrix(rept_mus2[, "mean.CI.scaled"])
matrix.names <- rept_mus2$species
rownames(matrix) <- matrix.names

(heat.tree6 <- gheatmap(tree, matrix, offset = 1, width = 0.1,
                        colnames = F) + 
    scale_fill_viridis() +
    labs(title = "(f) Population fluctuations") +
    theme(plot.title = element_text(size = 20, vjust = 1, hjust = 0)))

(phylo_density_rept_fl <- ggplot(rept_mus2) + 
    geom_density(aes(x = mean.CI.scaled), fill = "#dc863b", colour = "#dc863b",
                 alpha = 0.6) +
    theme_LPI() +
    scale_y_continuous(limits = c(0, 14),
                       breaks = c(0, 4 , 8, 12),
                       labels = c("0", "4", "8", "12")) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    labs(x = "\n\nMean population fluctuations", 
         y = "Density\n",
         title = "(l) Population fluctuations\n") +
    geom_vline(xintercept = 0, linetype = "dashed"))

ggsave(phylo_density_rept_fl, filename = "figures/updated/rept_density2.png", 
       height = 6, width = 5)

# Panel
(phylo_panel <- grid.arrange(heat.tree3, heat.tree1, heat.tree5, ncol = 3))

ggsave(phylo_panel, filename = "figures/updated/phylo_panel5.png",
       width = 15, height = 5.5)

(phylo_panel2 <- grid.arrange(heat.tree4, heat.tree2, heat.tree6, ncol = 3))

ggsave(phylo_panel2, filename = "figures/updated/phylo_panel6.png",
       width = 15, height = 5.5)

density_panel1 <- grid.arrange(phylo_density2, phylo_density, 
                               phylo_density3, ncol = 3)

ggsave(density_panel1, filename = "figures/updated/phylo_density1.png",
       height = 5.5, width = 15)

density_panel2 <- grid.arrange(phylo_density_amp_fl, phylo_density_fl, 
                               phylo_density_rept_fl, ncol = 3)

ggsave(density_panel2, filename = "figures/updated/phylo_density2.png",
       height = 5.5, width = 15)

# Mammals
mamtree <- read.nexus("mammals.txt")

# Find which species to include
m_species <- mus %>% filter(Class == "Mammalia") %>%
  dplyr::select(species) %>% distinct()

m_species <- m_species %>% filter(!species %in% c("Erethizon dorsata",
                                                  "Spermophilus franklinii",
                                                  "Pagophilus groenlandicus",
                                                  "Spermophilus tridecemlineatus",
                                                  "Taurotragus oryx",
                                                  "Dicrostonyx vinogradovi",
                                                  "Lemmus portenkoi",
                                                  "Damaliscus korrigum",
                                                  "Saiga borealis",
                                                  "Spermophilus columbianus",
                                                  "Monachus schauinslandi",
                                                  "Otaria flavescens",
                                                  "Spermophilus parryii",
                                                  "Hydrochoeris hydrochaeris",
                                                  "Cynogale bennettii",
                                                  "Lasiurus cinereus"))

m_species$species <- gsub(" ", "_", m_species$species, fixed = TRUE)
species <- m_species$species

#write.csv(m_species, file = "data/input/mam_names.csv")

#pruned.tree <- drop.tip(amphtree, setdiff(amphtree$tip.label, species))

#tree_species <- as.data.frame(amphtree$tip.label)

ggtree(mamtree) +
  geom_treescale() + 
  theme_tree2()

m_mus <- mus %>% filter(Class == "Mammalia") %>%
  filter(!species %in% c("Erethizon dorsata",
                         "Spermophilus franklinii",
                         "Pagophilus groenlandicus",
                         "Spermophilus tridecemlineatus",
                         "Taurotragus oryx",
                         "Dicrostonyx vinogradovi",
                         "Lemmus portenkoi",
                         "Damaliscus korrigum",
                         "Saiga borealis",
                         "Spermophilus columbianus",
                         "Monachus schauinslandi",
                         "Otaria flavescens",
                         "Spermophilus parryii",
                         "Hydrochoeris hydrochaeris",
                         "Cynogale bennettii",
                         "Lasiurus cinereus")) %>%
  group_by(species) %>%
  summarise(mean.mu = mean(mu))

m_mus$species <- gsub(" ", "_", m_mus$species, fixed = TRUE)

m_mus$node <- NA
tipnode <- seq_along(mamtree$tip.label)
names(tipnode) <- mamtree$tip.label
m_mus$node <- tipnode[m_mus$species] ## convert the tip label to tip node number
i <- is.na(m_mus$node)
m_mus$node[i] = m_mus$species[i] ## fill in internal node number

(tree <- ggtree(mamtree, aes(color = mean.mu), layout = "circular") %<+% m_mus +
    geom_tiplab(size=0) + 
    scale_color_viridis() + 
    theme(legend.position = "bottom"))

ggsave(tree, filename = "figures/updated/amp_tree.png",
       height = 6, width = 6)

matrix <- as.matrix(m_mus[, "mean.mu"])
matrix.names <- m_mus$species
rownames(matrix) <- matrix.names

(heat.tree <- gheatmap(tree, matrix, offset = 1, width = 0.1) + scale_fill_viridis())

ggsave(heat.tree, filename = "figures/updated/mam_heat.png", height = 6, width = 6)


# ** Fluctuations ----

# Figures with the confidence interval
# Geographic range
(range.plot.ci <- ggplot(uk.mus, aes(x = log(km2_range), y = CI)) +
   geom_point(alpha = 0.3, colour = "black") +
   scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75), 
                      labels = c("0", "0.25", "0.50", "0.75")) +
   theme_LPI() + 
   labs(x = bquote('\nLog geographic range ' ~ (km^2)), y = "Confidence interval\n", 
        title = "(a) Geographic range\n"))

# Local population size (logged)
(pop.size.plot.ci <- ggplot(uk.mus, aes(x = log(meanpop), y = CI)) +
    geom_point(alpha = 0.3, colour = "black") +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75), 
                       labels = c("0", "0.25", "0.50", "0.75")) +
    theme_LPI() +
    labs(x = "\nLog population size", y = "Confidence interval\n", 
         title = "(b) Mean population size\n"))

# Habitat specifity versus population change
(hab.spec.plot.ci <- ggplot(uk.mus, aes(x = hab_spec, y = CI)) +
    geom_point(alpha = 0.3, colour = "black") +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75), 
                       labels = c("0", "0.25", "0.50", "0.75")) +
    theme_LPI() +
    labs(x = "\nNumber of habitats", y = "Confidence interval\n", 
         title = "(c) Habitat specificity\n"))

# IUCN status
(IUCN.plot.UK.ci <- ggplot(uk.mus3, aes(x = IUCN, y = CI, colour = IUCN)) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.3)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75), 
                       labels = c("0", "0.25", "0.50", "0.75")) +
    theme_LPI2() + 
    theme(axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00", "#ffa15b", "#ff000e")) +
    labs(title = "(e) Red List status\n", y = "Population change\n") +
    guides(colour = F))

# Effect sizes from the different models
uk.effects <- read.csv("data/CI_effects.csv")

uk.effects1 <- uk.effects[c(1:3), ]

uk.effects1$Metric <- factor(uk.effects1$Metric, 
                             levels = c("Geographic range", "Mean population size",
                                        "Habitat specificity"),
                             labels=c("Geographic range", "Mean population size",
                                      "Habitat specificity"))

(eff_plot_uk_ci <- ggplot(uk.effects1, aes(x = Metric, y = mean)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), 
                    position = position_dodge(0.7)) + 
    scale_colour_manual(values = c("grey80", "grey60", "grey40", "black")) +
    scale_y_continuous(breaks = c(-0.01, 0, 0.01), limits = c(-0.012, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    theme_LPI2() +
    theme(legend.position = c(0.8, 0.2),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Coefficient\n", title = "(a) Rarity effects\n") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    guides(shape = F, colour = F))

uk.effects2 <- uk.effects[c(4:8), ]

# Reorder levels to match increasing risk
uk.effects2$Metric <- factor(uk.effects2$Metric, 
                             levels = c("Least concern", "Near threatened",
                                        "Vulnerable", "Endangered",
                                        "Critically endangered"),
                             labels=c("Least concern", "Near threatened", "Vulnerable",
                                      "Endangered", "Critically endangered"))


(eff_plot_uk_ci2 <- ggplot(uk.effects2, aes(x = Metric, y = mean)) +
    geom_pointrange(aes(ymin = lower, ymax = upper, colour = Metric), 
                    position = position_dodge(0.7)) + 
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ffa15b", "#ff000e")) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6), 
                       labels = c("0", "0.2", "0.4", "0.6")) +
    theme_LPI2() +
    theme(legend.position = c(0.3, 0.2),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Coefficient\n", title = "(b) Red List status effects\n") +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    guides(colour = F))

ci.panel.UK <- grid.arrange(eff_plot_uk_ci, eff_plot_uk_ci2, ncol = 2)
ggsave(ci.panel.UK, filename = "figures/Figure5_ci_effects.png", width = 10, height = 6)


# Table of models with mean rate sd
dataList <- list(IUCN.m.se, IUCN.m.sd, IUCN.m.ci, IUCN.m.ci.var, IUCN.m.ci.sigma)

# Create a list of input model names
dataListNames <- list("se",
                      "sd",
                      "CI",
                      "CI.w",
                      "sigma")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.sd <- as.data.frame(do.call(rbind, readyList))

mcmc.outputs.sd <- mcmc.outputs.sd %>% dplyr::select(modelName, variable:effect)

colnames(mcmc.outputs.sd) <- c("Model name", "Variable", "Posterior mean",
                               "Lower 95% CI", "Upper 95% CI", "Effective sample size",
                               "pMCMC", "Effect")



IUCN.effects <- mcmc.outputs.sd %>% filter(`Model name` == "sd" |
                                             `Model name` == "CI" |
                                             `Model name` == "se" |
                                             `Model name` == "CI.w" |
                                             `Model name` == "sigma")

IUCN.effects <- IUCN.effects %>% filter(Variable != "units")

IUCN.sdevs <- rbind(IUCN.ci, IUCN.ci.var, IUCN.se, IUCN.sd, IUCN.ci.sigma)
colnames(IUCN.effects)[1] <- "Model" 
IUCN.effects$Model <- gsub("Red List status fluctuations - ", "", IUCN.effects$Model, fixed = TRUE)
IUCN.effects$Variable <- gsub("Red.List.status", "", IUCN.effects$Variable, fixed = TRUE)
IUCN.effects$Variable <- gsub(".x", "", IUCN.effects$Variable, fixed = TRUE)
IUCN.effects$merge.id <- paste0(IUCN.effects$Model, IUCN.effects$Variable)
IUCN.sdevs$merge.id <- paste0(IUCN.sdevs$Model, IUCN.sdevs$Red.List.status)
IUCN.effects.merged <- merge(IUCN.effects, IUCN.sdevs, by = "merge.id")

IUCN.effects.merged$std.mean <- IUCN.effects.merged$`Posterior mean`/IUCN.effects.merged$sd
IUCN.effects.merged$std.lower <- IUCN.effects.merged$`Lower 95% CI`/IUCN.effects.merged$sd
IUCN.effects.merged$std.upper <- IUCN.effects.merged$`Upper 95% CI`/IUCN.effects.merged$sd

# Rearrange the IUCN categories 
IUCN.effects.merged$Variable <- factor(IUCN.effects.merged$Variable, 
                                       levels = c("Least concern",
                                                  "Near threatened",
                                                  "Vulnerable", "Endangered", 
                                                  "Critically endangered"),
                                       labels=c("Least concern (6523)",
                                                "Near threatened (505)", 
                                                "Vulnerable (582)",
                                                "Endangered (304)",
                                                "Critically endangered (152)"))

IUCN.effects.merged$Model.x <- factor(IUCN.effects.merged$Model.x, 
                                      levels = c("se",
                                                 "CI", 
                                                 "CI.w",
                                                 "sigma",
                                                 "sd"),
                                      labels = c("SE",
                                                 "CI",
                                                 "CI Weighted", 
                                                 "Sigma",
                                                 "SD"))

(eff_plot_IUCN_ci <- ggplot(IUCN.effects.merged, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, 
                        colour = Variable, shape = Model.x), 
                    position = position_dodge(0.7), size = 0.7) + 
    scale_colour_manual(values = c("#4bb33c", "#d5db00", "#ffef00",
                                   "#ffa15b", "#ff000e")) +
    scale_y_continuous(limits = c(0, 2),
                       breaks = c(0, 0.5, 1.0, 1.5, 2.0),
                       labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
    theme_LPI3() +
    theme(legend.position = c(0.25, 0.85),
          axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nRed List status") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") + 
    guides(colour = F))

ggsave(eff_plot_IUCN_ci, filename = "figures/Fig4_IUCN_effects_fl.pdf", device = "pdf",
       dpi = 300, width = 4.5, height = 4.5)

# Testing fluctuations for GLOBAL system models
# Table of models with mean rate sd
dataList <- list(system.m.se, system.m.sd, system.m.ci, system.m.ci.var, system.m.sigma)

# Create a list of input model names
dataListNames <- list("se",
                      "sd",
                      "CI",
                      "CI.w",
                      "sigma")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.sd <- as.data.frame(do.call(rbind, readyList))

mcmc.outputs.sd <- mcmc.outputs.sd %>% dplyr::select(modelName, variable:effect)

colnames(mcmc.outputs.sd) <- c("Model name", "Variable", "Posterior mean",
                               "Lower 95% CI", "Upper 95% CI", "Effective sample size",
                               "pMCMC", "Effect")

system.effects <- mcmc.outputs.sd %>% filter(`Model name` == "sd" |
                                               `Model name` == "CI" |
                                               `Model name` == "se" |
                                               `Model name` == "CI.w" |
                                               `Model name` == "sigma")

system.effects <- system.effects %>% filter(Variable != "units")

system.sdevs <- rbind(system.ci, system.ci.var, system.se, system.sd, system.ci.sigma)
colnames(system.effects)[1] <- "Model" 
system.effects$Variable <- gsub("system", "", system.effects$Variable, fixed = TRUE)
system.effects$merge.id <- paste0(system.effects$Model, system.effects$Variable)
system.sdevs$merge.id <- paste0(system.sdevs$Model, system.sdevs$system)
system.effects.merged <- merge(system.effects, system.sdevs, by = "merge.id")

system.effects.merged$std.mean <- system.effects.merged$`Posterior mean`/system.effects.merged$sd
system.effects.merged$std.lower <- system.effects.merged$`Lower 95% CI`/system.effects.merged$sd
system.effects.merged$std.upper <- system.effects.merged$`Upper 95% CI`/system.effects.merged$sd

system.effects.merged$Model.x <- factor(system.effects.merged$Model.x, 
                                        levels = c("se",
                                                   "CI", 
                                                   "CI.w",
                                                   "sigma",
                                                   "sd"),
                                        labels = c("SE",
                                                   "CI",
                                                   "CI Weighted", 
                                                   "Sigma",
                                                   "SD"))

system.effects.merged$Variable <- factor(system.effects.merged$Variable, 
                                         levels = c("Freshwater", 
                                                    "Marine",
                                                    "Terrestrial"),
                                         labels=c("Freshwater (2570)",
                                                  "Marine (2418)", 
                                                  "Terrestrial (4298)"))


(eff_plot_system_ci <- ggplot(system.effects.merged, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Variable,
                        shape = Model.x), size = 0.7, 
                    position = position_dodge(0.7)) + 
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_y_continuous(limits = c(-0.5, 1.5),
                       breaks = c(-0.5, 0, 0.5, 1.0, 1.5),
                       labels = c("-0.5", "0", "0.5", "1.0", "1.5")) +   
    theme_LPI3() +
    theme(legend.position = c(0.3, 0.3),
          axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nRealm") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    guides(colour = F))

ggsave(eff_plot_system_ci, filename = "figures/Fig1_system_effects_fl.pdf", device = "pdf",
       dpi = 300, width = 5, height = 5)


# Testing fluctuations for GLOBAL taxa models
# Table of models with mean rate sd
dataList <- list(class.m.se, class.m.sd, class.m.ci, class.m.ci.var, class.m.sigma)

# Create a list of input model names
dataListNames <- list("se",
                      "sd",
                      "CI",
                      "CI.w",
                      "sigma")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.sd <- as.data.frame(do.call(rbind, readyList))

mcmc.outputs.sd <- mcmc.outputs.sd %>% dplyr::select(modelName, variable:effect)

colnames(mcmc.outputs.sd) <- c("Model name", "Variable", "Posterior mean",
                               "Lower 95% CI", "Upper 95% CI", "Effective sample size",
                               "pMCMC", "Effect")

taxa.effects <- mcmc.outputs.sd %>% filter(`Model name` == "sd" |
                                             `Model name` == "CI" |
                                             `Model name` == "se" |
                                             `Model name` == "CI.w" |
                                             `Model name` == "sigma")

taxa.effects <- taxa.effects %>% filter(Variable != "units")

taxa.sdevs <- rbind(class.ci, class.ci.var, class.se, class.sd, taxa.ci.sigma)
colnames(taxa.effects)[1] <- "Model" 
taxa.effects$Variable <- gsub("Class", "", taxa.effects$Variable, fixed = TRUE)
taxa.effects$merge.id <- paste0(taxa.effects$Model, taxa.effects$Variable)
taxa.sdevs$merge.id <- paste0(taxa.sdevs$Model, taxa.sdevs$Class)
taxa.effects.merged <- merge(taxa.effects, taxa.sdevs, by = "merge.id")

taxa.effects.merged$std.mean <- taxa.effects.merged$`Posterior mean`/taxa.effects.merged$sd
taxa.effects.merged$std.lower <- taxa.effects.merged$`Lower 95% CI`/taxa.effects.merged$sd
taxa.effects.merged$std.upper <- taxa.effects.merged$`Upper 95% CI`/taxa.effects.merged$sd

taxa.effects.merged$Model.x <- factor(taxa.effects.merged$Model.x, 
                                      levels = c("se",
                                                 "CI", 
                                                 "CI.w",
                                                 "sigma",
                                                 "sd"),
                                      labels=c("SE",
                                               "CI",
                                               "CI Weighted", 
                                               "Sigma",
                                               "SD"))

taxa.effects.merged$Variable <- factor(taxa.effects.merged$Variable, 
                                       levels = c("Elasmobranchii",
                                                  "Actinopterygii",
                                                  "Amphibia", 
                                                  "Aves",
                                                  "Mammalia",
                                                  "Reptilia"),
                                       labels=c("Elasmobranchii (127)",
                                                "Actinopterygii (1626)",
                                                "Amphibia (193)", 
                                                "Aves (5854)",
                                                "Mammalia (1158)",
                                                "Reptilia (322)"))

(eff_plot_taxa_ci <- ggplot(taxa.effects.merged, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Variable,
                        shape = Model.x), 
                    position = position_dodge(0.7), size = 0.7) + 
    scale_colour_manual(values = c("#96CDCD", "#262F46", "#02401b",
                                   "#d8b70a", "#972d15", "#dc863b")) +
    scale_y_continuous(limits = c(-0.3, 2.02),
                       breaks = c(0, 0.5, 1.0, 1.5, 2.0),
                       labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
    theme_LPI3() +
    theme(legend.position = c(0.3, 0.2),
          axis.line.x = element_line(),
          axis.line.y = element_line()) +
    labs(y = "Standardised effect size\n", x = "\nTaxa") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    guides(colour = F))

ggsave(eff_plot_taxa_ci, filename = "figures/Fig2_taxa_effects_fl.pdf", device = "pdf",
       dpi = 600, width = 5, height = 5)

fluctuations <- all_outputs_5th_april[c(344:358, 368:382, 392:406),]

rarity.effects.global.fl <- fluctuations %>%
  filter(Variable != "(Intercept)" & Variable != "units")

rarity.effects.global.fl$Variable <- recode(rarity.effects.global.fl$Variable, 
                                            "log(km2_range)" = "log.range")

rarity.effects.global.fl$Variable <- recode(rarity.effects.global.fl$Variable, 
                                            "log(range)" = "log.range")

rarity.effects.global.fl$Variable <- recode(rarity.effects.global.fl$Variable, 
                                            "log(meanpop)" = "log.meanpop")

rarity.effects.global.fl <- rarity.effects.global.fl %>% separate(Model.name,
                                                                  c("Name", "Model"), "fluctuations ")

# Checking sample sizes
length(unique(bird_mu_ranges$id))
length(unique(mus.popsize$id))
length(unique(mus.habspec$id))

rarity.effects.global.fl$Variable <- factor(rarity.effects.global.fl$Variable, 
                                            levels = c("log.range",
                                                       "log.meanpop", 
                                                       "hab_spec"),
                                            labels=c("Geographic range (5381)",
                                                     "Mean population size (4310)", 
                                                     "Habitat specificity (7901)"))

# Standardising the effect sizes using the sd of the raw data for each model
rarity.effects.global.fl$sd <- NA
rarity.effects.global.fl[1, 11] <- sd(bird_mu_ranges$CI)
rarity.effects.global.fl[2, 11] <- sd(bird_mu_ranges$CI)
rarity.effects.global.fl[3, 11] <- sd(bird_mu_ranges$sigma.2)
rarity.effects.global.fl[4, 11] <- sd(bird_slopes_ranges$SE)
rarity.effects.global.fl[5, 11] <- sd(LPIraw.bird.range$sdpop1)
rarity.effects.global.fl[6, 11] <- sd(mus.popsize$CI)
rarity.effects.global.fl[7, 11] <- sd(mus.popsize$CI)
rarity.effects.global.fl[8, 11] <- sd(mus.popsize$sigma.2)
rarity.effects.global.fl[9, 11] <- sd(slopes.popsize$SE)
rarity.effects.global.fl[10, 11] <- sd(LPI.mean.pop2$sdpop1, na.rm = T)
rarity.effects.global.fl[11, 11] <- sd(mus.habspec$CI)
rarity.effects.global.fl[12, 11] <- sd(mus.habspec$CI)
rarity.effects.global.fl[13, 11] <- sd(mus.habspec$sigma.2)
rarity.effects.global.fl[14, 11] <- sd(slopes.habspec$SE)
rarity.effects.global.fl[15, 11] <- sd(LPIraw.habspec$sdpop1)

rarity.effects.global.fl$std.mean <- rarity.effects.global.fl$Posterior.mean/rarity.effects.global.fl$sd
rarity.effects.global.fl$std.lower <- rarity.effects.global.fl$Lower.95..CI/rarity.effects.global.fl$sd
rarity.effects.global.fl$std.upper <- rarity.effects.global.fl$Upper.95..CI/rarity.effects.global.fl$sd

rarity.effects.global.fl$Model <- factor(rarity.effects.global.fl$Model,
                                         levels = c("SE", "CI",
                                                    "CI weighted", "sigma",
                                                    "SD"),
                                         labels = c("SE", "CI",
                                                    "CI weighted", "Sigma",
                                                    "SD"))

(eff_plot_global_ci <- ggplot(rarity.effects.global.fl, aes(x = Variable, y = std.mean,
                                                            colour = Model, shape = Model)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper), 
                    position = position_dodge(0.7), size = 0.7) + 
    scale_colour_manual(values = c("grey80", "grey65", "grey50", "grey35", "black")) +
    scale_y_continuous(breaks = c(-0.08, -0.04, 0, 0.04, 0.08), 
                       limits = c(-0.085, 0.085),
                       labels = c("-0.08", "-0.04", "0", "0.04", "0.08")) +
    theme_LPI() +
    theme(legend.position = c(0.85, 0.85),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Standardised effect size\n", x = "\nMetric") +
    geom_hline(yintercept = 0, linetype = "dashed"))

ci.panel.sd1 <- grid.arrange(eff_plot_global_ci, eff_plot_IUCN_ci, ncol = 2)
ggsave(ci.panel.sd1, filename = "figures/updated/CI_collage1.jpg", height = 6.5, 
       width = 10)

ci.panel.sd2 <- grid.arrange(eff_plot_system_ci, eff_plot_taxa_ci,  ncol = 2)
ggsave(ci.panel.sd2, filename = "figures/updated/CI_collage3.jpg", height = 6, 
       width = 10)

# ** Correlation plots ----
rarity.correlation <- inner_join(uk.mus.popsize, uk.mus.range,
                                 by = "id") %>% inner_join(., uk.mus.habspec,
                                                           by = "id") %>%
  dplyr::select(log.range, log.meanpop, hab_spec)

colnames(rarity.correlation) <- c("Log geographic range", "Log mean population size", "Habitat specificity")


rarity.corr <- cor(rarity.correlation)
uk.corr <- corrplot(rarity.corr, method = "color",
                    number.cex = 1.7, type = "upper",
                    col = brewer.pal(n = 8, name = "PuOr"),
                    tl.col = "black", addCoef.col = "black",
                    diag = FALSE)

rarity.correlation.global <- inner_join(b_m_ranges, mus.popsize,
                                        by = "id") %>% inner_join(., mus.habspec,
                                                                  by = "id") %>%
  dplyr::select(range, log.meanpop, hab_spec)

rarity.correlation.global$range <- log(rarity.correlation.global$range)
rarity.correlation.global <- rarity.correlation.global %>% drop_na(range)
colnames(rarity.correlation.global) <- c("Log geographic range", "Log mean population size", "Habitat specificity")

rarity.corr.global <- cor(rarity.correlation.global)
global.corr <- corrplot(rarity.corr.global, method = "circle", type = "upper")

global.corr <- corrplot(rarity.corr.global, method = "color",
                        number.cex = 1.7, type = "upper",
                        col = brewer.pal(n = 8, name = "PuOr"),
                        tl.col = "black", addCoef.col = "black",
                        diag = FALSE)


# SI figures ----

# ** Sampling units ----
# Pop change split by most popular methods of collecting data
methods <- mus %>% dplyr::select(id, mu, Units) %>%
  group_by(Units) %>% tally %>% arrange(desc(n))

# Renaming some of the units categories that are essentially the same
# e.g. individuals and number of individuals

mus$Units2 <- mus$Units
mus$Units2 <- recode(mus$Units2, 
                     "Number of individuals" = 
                       "Individuals")
mus$Units2 <- recode(mus$Units2, 
                     "Number of nests" = 
                       "Nests")
mus$Units2 <- recode(mus$Units2, 
                     "Number of pairs" = 
                       "Pairs")
mus$Units2 <- recode(mus$Units2, 
                     "Number of breeding pairs" = 
                       "Pairs")
mus$Units2 <- recode(mus$Units2, 
                     "Breeding pairs" = 
                       "Pairs")
mus$Units2 <- recode(mus$Units2, 
                     "Annual index" = 
                       "Index")
mus$Units2 <- recode(mus$Units2, 
                     "Annual Index" = 
                       "Index")
mus$Units2 <- recode(mus$Units2, 
                     "Annual index of population change" = 
                       "Index")
mus$Units2 <- recode(mus$Units2, 
                     "Annual index of abundance" = 
                       "Index")
mus$Units2 <- recode(mus$Units2, 
                     "Population index" = 
                       "Index")

methods2 <- mus %>% dplyr::select(id, mu, Units2) %>%
  group_by(Units2) %>% tally %>% arrange(desc(n))

mus.methods <- mus %>% filter(Units2 == "Index" |
                                Units2 == "Individuals" |
                                Units2 == "Pairs" |
                                Units2 == "Nests" |
                                Units2 == "Population estimate")

mus.methods$Units2 <- factor(mus.methods$Units2, 
                             levels = c("Index",
                                        "Individuals",
                                        "Pairs", "Nests", 
                                        "Population estimate"),
                             labels=c("Index (2827)",
                                      "Individuals (1167)",
                                      "Pairs (505)", "Nests (224)", 
                                      "Population estimate (97)"))


(mus.methods.plot <- ggplot(mus.methods, aes(x = Units2, y = mu)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8,
                     fill = "grey30", ) +
    geom_point(aes(y = mu), 
               position = position_jitter(width = .15),
               size = 2, alpha = 0.05) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    theme_LPI2() + 
    theme(axis.title.x = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")) +
    labs(y = bquote(atop(italic(mu), ' '))) +
    guides(colour = F))

ggsave(mus.methods.plot, filename = "figures/SI_methods.pdf", width = 9, height = 9,
       dpi = 300, device ="pdf")

# ** UK maps ----
uk.map <- filter(uk.mus, Decimal.Latitude < 60 & Decimal.Longitude > -10)

(uk.mu.map <- ggplot(uk.map, aes(x = Decimal.Longitude, y = Decimal.Latitude, colour = mu)) +
    borders("worldHires", region = "UK:Great Britain",  colour = "#bfbfbf", 
            fill = "#bfbfbf", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_viridis(breaks = c(0.14, 0, -0.14),
                         labels = c("0.14", "0", "-0.14")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "italic"),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size=20, vjust=1, hjust=0)) +
    labs(colour = bquote(italic(mu))))

uk.map2 <- filter(uk.slopes, Decimal.Latitude < 60 & Decimal.Longitude > -10)

(uk.slope.map <- ggplot(uk.map2, aes(x = Decimal.Longitude, y = Decimal.Latitude, colour = slope)) +
    borders("worldHires", region = "UK:Great Britain",  colour = "#bfbfbf", 
            fill = "#bfbfbf", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_viridis(breaks = c(0.30, 0, -0.30),
                         labels = c("0.30", "0", "-0.30")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size=20, vjust=1, hjust=0)) +
    labs(colour = "Slope"))

(uk.sigma.map <- ggplot(uk.map, aes(x = Decimal.Longitude, y = Decimal.Latitude, colour = sigma.2)) +
    borders("worldHires", region = "UK:Great Britain",  colour = "#bfbfbf", 
            fill = "#bfbfbf", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_viridis(breaks = c(0.1, 0.4, 0.7),
                         labels = c("0.1", "0.4", "0.7")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 16, face = "italic"),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size=20, vjust=1, hjust=0)) +
    labs(colour = bquote(italic(sigma^2))))

(uk.se.map <- ggplot(uk.map2, aes(x = Decimal.Longitude, y = Decimal.Latitude, colour = slope_se)) +
    borders("worldHires", region = "UK:Great Britain",  colour = "#bfbfbf", 
            fill = "#bfbfbf", size = 0.3) +
    coord_proj("+proj=wintri") +
    theme_map() +
    geom_point(alpha = 0.8, size = 2) +
    scale_colour_viridis(breaks = c(0.05, 0.15, 0.25),
                         labels = c("0.05", "0.15", "0.25")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.justification = "top",
          plot.title = element_text(size=20, vjust=1, hjust=0)) +
    labs(colour = "Slope SE"))

uk.pop.change.panel <- grid.arrange(uk.mu.map, uk.sigma.map,
                                    uk.slope.map, uk.se.map, ncol =2)
ggsave(uk.pop.change.panel, filename = "figures/SI_uk_maps.png", height = 10, width = 10)

# ** Hab spec UK profiling method ----
mcmc_preds_habspec_uk2 <- as.data.frame(predict(rarity.hab.prof.uk,
                                                type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_habspec_uk2$hab_spec <- uk.mus.habspec.prof$hab_spec

(hab.spec.plot2 <- ggplot() +
    geom_point(data = uk.mus.habspec.prof, aes(x = hab_spec, y = mu),
               size = 4, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(data = mcmc_preds_habspec_uk2, aes(x = hab_spec,
                                                   ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_habspec_uk2, aes(x = hab_spec, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0.20, 0.10, 0, -0.10, -0.20), 
                       limits = c(-0.25, 0.25),
                       labels = c("0.20", "0.10", "0", "-0.10", "-0.20")) +
    theme_LPI() +
    labs(x = "\nNumber of habitats", 
         y = bquote(atop(italic(mu), ' '))))

ggsave(hab.spec.plot2, filename = "figures/habspec_prof.pdf", device = "pdf", dpi = 300,
       height = 9, width = 9)

# ** Duration ----
(mu.duration <- ggplot(mus, aes(x = duration, y = mu)) +
   geom_pointrange(aes(ymin = lCI, ymax = uCI), alpha = 0.1, size = 0.3,
                   colour = "black") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   scale_y_continuous(limits = c(-1, 1),
                      breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                      labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
   theme_LPI() +
   labs(y = bquote(atop(italic(mu), ' ')), x = "\nYears"))

(slope.duration <- ggplot(slopes, aes(x = lengthyear, y = slope)) +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.1,
                    colour = "black", size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(limits = c(-1, 1),
                       breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
                       labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    theme_LPI() +
    labs(y = "Slope\n", x = "\nYears"))

(sigma.duration <- ggplot(mus[mus$sigma.2 < 1,], aes(x = duration, y = sigma.2)) +
    # removed one outlier point
    geom_point(alpha = 0.1, colour = "black") +
    #geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                       limits = c(0, 1),
                       labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
    theme_LPI() +
    labs(y = bquote(atop(italic(sigma^2), ' ')), x = "\nYears"))

duration.panel <- grid.arrange(mu.duration, slope.duration, sigma.duration,
                                ncol = 3)
ggsave(duration.panel, filename = "figures/SI_duration_panel.pdf", device ="pdf",
       height = 5, width = 15, dpi = 300)

# ** Fluctuations biome ----
# Biome effects
biome.effects.global.fl <- global_outputs18thFeb %>%
  filter(Model.name %in% c("Biome - fluctuations sigma",
                           "Biome - fluctuations CI",                       
                           "Biome - fluctuations CI weighted",
                           "Biome - fluctuations SE",                      
                           "Biome - fluctuations SD") &
         Variable != "units")

biome.effects.global.fl$Variable <- gsub("biome", "", 
                                      biome.effects.global.fl$Variable, fixed = TRUE)

biome.effects.global.fl <- biome.effects.global.fl %>% separate(Model.name,
                                                          c("Name", "Model"), " - ")

biome.effects.global.fl$Model <- recode(biome.effects.global.fl$Model, 
                                     "mu" = "State-space")
biome.effects.global.fl$Model <- recode(biome.effects.global.fl$Model, 
                                     "weighted" = "Weighted")
biome.effects.global.fl$Model <- recode(biome.effects.global.fl$Model, 
                                     "slope" = "Linear")

# Checking sample sizes
mus.biome %>% group_by(biome) %>% tally

biome.effects.global.fl$Variable <- factor(biome.effects.global.fl$Variable, 
                                        levels = c("Boreal forests/taiga",
                                                   "Deserts and xeric shrublands",
                                                   "Trop. and subtrop. grasslands savannas and shrublands",
                                                   "Large lakes",
                                                   "Mediterranean forests woodlands and scrub",
                                                   "Montane freshwaters",
                                                   "Montane grasslands and shrublands",
                                                   "Polar freshwaters",
                                                   "Polar seas",
                                                   "Temperate forests",
                                                   "Temperate wetlands and rivers",
                                                   "Temperate grasslands savannas and shrublands",
                                                   "Tropical wetlands and rivers",
                                                   "Tropical and subtropical forests",
                                                   "Tropical coral",
                                                   "Tundra",
                                                   "Xeric freshwaters and endorheic basins"),
                                        labels = c("Boreal forests (1289)",
                                                 "Deserts (46)",
                                                 "Tropical savannas (216)",
                                                 "Large lakes (225)",
                                                 "Mediterranean forests (225)",
                                                 "Montane freshwaters (23)",
                                                 "Montane grasslands (36)",
                                                 "Polar freshwaters (362)",
                                                 "Polar seas (26)",
                                                 "Temperate forests (1785)",
                                                 "Temperate wetlands (1712)",
                                                 "Temperate grasslands (321)",
                                                 "Tropical wetlands (184)",
                                                 "Tropical forests (178)",
                                                 "Tropical coral (102)",
                                                 "Tundra (201)",
                                                 "Xeric freshwaters (52)"))

# Standardising the effect sizes using the sd of the raw data for each biome category
mus.biome$sigma.2s <- scales::rescale(mus.biome$sigma.2, to = c(0, 1))
global.biome.sigma.sd <- mus.biome %>% group_by(biome) %>%
  summarise(sd = sd(sigma.2s))

mus.biome$CI <- (mus.biome$uCI - mus.biome$lCI)/2
mus.biome$CIs <- scales::rescale(mus.biome$CI, to = c(0, 1))
global.biome.ci.sd <- mus.biome %>% group_by(biome) %>%
  summarise(sd = sd(CIs))

slopes.biome$SEs <- scales::rescale(slopes.biome$slope_se, to = c(0, 1))
global.biome.se.sd <- slopes.biome %>% group_by(biome) %>%
  summarise(sd = sd(SEs))

LPI.long.biome <- filter(LPI.long, biome != "Unknown" & biome != "Mangroves" & 
                           biome != "Oceanic islands")

LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical floodplain rivers and wetland complexes" = 
                               "Tropical wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical coastal rivers" = 
                               "Tropical wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical coniferous forests" = 
                               "Tropical and subtropical forests")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical moist broadleaf forests" = 
                               "Tropical and subtropical forests")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical dry broadleaf forests" = 
                               "Tropical and subtropical forests")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical upland rivers" = 
                               "Tropical wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Tropical and subtropical grasslands savannas and shrublands" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Flooded grasslands and savannas" = 
                               "Trop. and subtrop. grasslands savannas and shrublands")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate coastal rivers" = 
                               "Temperate wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate floodplain rivers and wetlands" = 
                               "Temperate wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate broadleaf and mixed forests" = 
                               "Temperate forests")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate upland rivers" = 
                               "Temperate wetlands and rivers")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate coniferous forests" = 
                               "Temperate forests")
LPI.long.biome$biome <- recode(LPI.long.biome$biome, 
                             "Temperate upwelling" = 
                               "Temperate wetlands and rivers")

global.biome.sd.sd <- LPI.long.biome %>% group_by(biome) %>%
  summarise(sd = sd(scalepop))

biome.sd.global <- rbind(global.biome.sigma.sd, global.biome.ci.sd, 
                         global.biome.ci.sd, global.biome.se.sd,
                         global.biome.sd.sd)

biome.effects.global.fl <- cbind(biome.effects.global.fl, biome.sd.global)

biome.effects.global.fl$std.mean <- biome.effects.global.fl$Posterior.mean/biome.effects.global.fl$sd
biome.effects.global.fl$std.lower <- biome.effects.global.fl$Lower.95..CI/biome.effects.global.fl$sd
biome.effects.global.fl$std.upper <- biome.effects.global.fl$Upper.95..CI/biome.effects.global.fl$sd

(eff_plot_biome_fl <- ggplot(biome.effects.global.fl, aes(x = Variable, y = std.mean)) +
    geom_pointrange(aes(ymin = std.lower, ymax = std.upper, colour = Model,
                        shape = Model),
                    position = position_dodge(0.7), size = 1) + 
    scale_colour_manual(values = c("grey80", "grey60", "grey40", "grey20", "black")) +
    #scale_y_continuous(breaks = c(1, 0.5, 0, -0.5, -1), 
    #                   limits = c(-1.18, 1.18),
    #                   labels = c("1.0", "0.5", "0", "-0.5", "-1.0")) +
    theme_LPI2() +
    theme(legend.position = c(0.1, 0.8),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")) +
    labs(y = "Standardised effect size\n") +
    geom_hline(yintercept = 0, linetype = "dashed"))

ggsave(eff_plot_biome_fl, filename = "figures/SI_biome_effects_fl.pdf", device = "pdf", dpi = 300,
       width = 18, height = 6)

# ** Threats and fluctuations ----
(threat_sigma_plot <- ggplot(top_threats_mus, aes(x = sigma.2, y = title, fill = title)) +
   geom_density_ridges() +
   theme_LPI() +
   scale_x_continuous(limits = c(0, 0.6)) +
   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey10") +
   labs(x = bquote(atop(' ', italic(sigma^2))), y = "Threat\n") +
   scale_fill_viridis_d(option = "magma", direction = -1, 
                        begin = 0.3, end = 0.8) +
   guides(fill = FALSE))

ggsave(threat_sigma_plot, filename = "figures/SI_threats_sigma.pdf", device = "pdf", dpi = 300,
       height = 5, width = 9)

# Calculate predictions
mcmc_preds_n_threats_sigma <- as.data.frame(predict(threat.n.m.sigma,
                                              type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_n_threats_sigma$n <- threats_sum_species_mus$n

(threat_n_sigma <- ggplot() +
    geom_point(data = threats_sum_species_mus, aes(x = n, y = sigma.2), colour = "grey70") +
    geom_ribbon(data = mcmc_preds_n_threats_sigma, aes(x = n, ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_n_threats_sigma, aes(x = n, y = fit), 
              colour = "black", size = 1) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    labs(x = "\nNumber of threats", y = bquote(atop(italic(mu), ' '))) +
    scale_y_continuous(limits = c(0, 0.6)) +
    theme_LPI())

ggsave(threat_n_sigma, filename = "figures/SI_threat_n_sigma.pdf", device = "pdf", dpi = 300,
       height = 5, width = 5)

# ** Number of observations and duration ----
(records_duration <- ggplot(mus, aes(x = duration, y = End)) +
  geom_point(alpha = 0.2, colour = "grey30", size = 4) +
   theme_LPI() +
   geom_smooth(method = "lm", colour = "aquamarine4", size = 1.5) +
   labs(x = "\nDuration (years)", y = "Survey points (n)\n"))

ggsave(records_duration, filename = "figures/SI_duration_records.pdf",
       device = "pdf", dpi = 300, height = 9, width = 9)

# ** Rarity taxa interaction ----
rarity_int <- global_outputs18thFeb %>%
  filter(str_detect(Model.name, pattern = "interaction") &
           str_detect(Variable, pattern = ":"))

rarity_int <- rarity_int[-1, ]

rarity_int$Variable <- gsub(":", "*", rarity_int$Variable, fixed = TRUE)
rarity_int$Variable <- gsub("Class", "", rarity_int$Variable, fixed = TRUE)
rarity_int$Variable <- gsub("log.meanpop", "Mean population size", rarity_int$Variable, fixed = TRUE)
rarity_int$Variable <- gsub("hab_spec", "Habitat specificity", rarity_int$Variable, fixed = TRUE)

rarity_int$test <- "A"
rarity_int[rarity_int$Variable == "Mean population size*Aves",]$test <- "B"
rarity_int[rarity_int$Variable == "Mean population size*Mammalia",]$test <- "B"

(rarity_int_plot <- ggplot(rarity_int, aes(x = fct_rev(Variable), y = Posterior.mean)) +
    geom_pointrange(aes(ymin = `Lower.95..CI`, ymax = `Upper.95..CI`,
                        fill = test),
                    size = 1.5, shape = 21) +
    scale_fill_manual(values = c("grey30", "aquamarine4")) +
    coord_flip() +
    theme_LPI() +
    labs(x = NULL, y = "\nPosterior mean") +
    guides(fill = FALSE) +
    scale_y_continuous(limits = c(-0.013, 0.013),
                       breaks = c(-0.01, 0, 0.01),
                       labels = c("-0.01", "0", "0.01")) +
    geom_hline(yintercept = 0, linetype = "dashed"))

ggsave(rarity_int_plot, filename = "figures/SI_rarity_int.pdf",
       device = "pdf", dpi = 300, height = 5, width = 9)

# ** Phylogeny ----
bird_phylo <- read.csv("data/output/bird_VCV_trend_raw.csv")
bird_phylo$metric <- factor(bird_phylo$metric,
                            levels = c("phylo", "species", "units"),
                            labels = c("Phylogeny", "Species", "Sigma"))

# Phylogeny and species effects
(bird_plot_trend <- ggplot(bird_phylo, aes(x = metric, y = value)) +
    geom_violin(colour = "#d8b70a", fill = "#d8b70a", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003), 
                       limits = c(0, 0.003),
                       labels = c("0", "0.001", "0.002", "0.003")) +
    theme_LPI() +
    theme(legend.position = c(0.8, 0.2),
          axis.title.x = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    labs(y = "Posterior mean\n") +
    geom_hline(yintercept = 0, linetype = "dashed"))

bird_phylo_fl <- read.csv("data/output/bird_VCV_fl_raw.csv")
bird_phylo_fl$metric <- factor(bird_phylo_fl$metric,
                            levels = c("phylo", "species", "units"),
                            labels = c("Phylogeny", "Species", "Sigma"))

(bird_plot_fl <- ggplot(bird_phylo_fl, aes(x = metric, y = value)) +
    geom_violin(colour = "#d8b70a", fill = "#d8b70a", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.001, 0.002, 0.003), 
                       limits = c(0, 0.003),
                       labels = c("0", "0.001", "0.002", "0.003")) +
    theme_LPI() +
    labs(y = "Posterior mean\n", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))

amp_phylo <- read.csv("data/output/amp_VCV_trend_raw.csv")
amp_phylo$metric <- factor(amp_phylo$metric,
                               levels = c("phylo", "species", "units"),
                               labels = c("Phylogeny", "Species", "Sigma"))

(amp_plot_trend <- ggplot(amp_phylo, aes(x = metric, y = value)) +
    geom_violin(colour = "#02401b", fill = "#02401b", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03), 
                       limits = c(0, 0.03),
                       labels = c("0", " 0.01", " 0.02", " 0.03")) +
    theme_LPI() +
    labs(y = "Posterior mean\n", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))

amp_phylo_fl <- read.csv("data/output/amp_VCV_fl_raw.csv")
amp_phylo_fl$metric <- factor(amp_phylo_fl$metric,
                           levels = c("phylo", "species", "units"),
                           labels = c("Phylogeny", "Species", "Sigma"))

(amp_plot_fl <- ggplot(amp_phylo_fl, aes(x = metric, y = value)) +
    geom_violin(colour = "#02401b", fill = "#02401b", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.03, 0.06, 0.09), 
                       limits = c(0, 0.09),
                       labels = c("0", " 0.03", " 0.06", " 0.09")) +
    theme_LPI() +
    labs(y = "Posterior mean\n", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))


rept_phylo <- read.csv("data/output/rept_VCV_trend_raw.csv")
rept_phylo$metric <- factor(rept_phylo$metric,
                           levels = c("phylo", "species", "units"),
                           labels = c("Phylogeny", "Species", "Sigma"))

(rept_plot_trend <- ggplot(rept_phylo, aes(x = metric, y = value)) +
    geom_violin(colour = "#dc863b", fill = "#dc863b", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.02, 0.04, 0.06), 
                       limits = c(0, 0.06),
                       labels = c("0", " 0.02", " 0.04", " 0.06")) +
    theme_LPI() +
    labs(y = "Posterior mean\n", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))

rept_phylo_fl <- read.csv("data/output/rept_VCV_fl_raw.csv")
rept_phylo_fl$metric <- factor(rept_phylo_fl$metric,
                            levels = c("phylo", "species", "units"),
                            labels = c("Phylogeny", "Species", "Sigma"))

(rept_plot_fl <- ggplot(rept_phylo_fl, aes(x = metric, y = mean)) +
    geom_violin(colour = "#dc863b", fill = "#dc863b", alpha = 0.6) + 
    scale_y_continuous(breaks = c(0, 0.04, 0.08, 0.12), 
                       limits = c(0, 0.12),
                       labels = c("0", " 0.04", " 0.08", " 0.12")) +
    theme_LPI() +
    labs(y = "Posterior mean\n", x = NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30"))

phylo_SI_panel <- grid.arrange(bird_plot_trend, bird_plot_fl,
                               amp_plot_trend, amp_plot_fl,
                               rept_plot_trend, rept_plot_fl)

ggsave(phylo_SI_panel, filename = "figures/SI_phylo.pdf", device = "pdf",
       dpi = 300, height = 15, width = 12)

# ** Left/Right truncation ----
global_mus_left_trunc <- read.csv("data/output/global_mus_left_trunc.csv")
global_mus_right_trunc <- read.csv("data/output/global_mus_right_trunc.csv")

colnames(mus)
colnames(global_mus_left_trunc)
mus <- mus %>% dplyr::select(1:34)

mus$Model <- "Full time-series"
global_mus_left_trunc$Model <- "Left-truncated"
global_mus_right_trunc$Model <- "Right-truncated"

all_mus <- rbind(mus, global_mus_left_trunc, global_mus_right_trunc)

(truncation <- ggplot(all_mus, aes(x = mu, colour = Model)) +
    geom_line(stat = "density", aes(linetype = Model), size = 2) +
    scale_colour_viridis_d(begin = 0, end = 0.9) +
    theme_LPI() +
    scale_x_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    labs(y = "Density\n", x = bquote(atop(' ', italic(mu)))))

ggsave(truncation, filename = "figures/SI_truncation.pdf",
       device = "pdf", dpi = 300, height = 7, width = 9)

# ** UK rarity
# Geographic range
mcmc_preds_range_uk <- as.data.frame(predict(rarity.range.uk,
                                             type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_range_uk$range <- uk.mus.range$km2_range

(range.plot <- ggplot() +
    geom_point(data = uk.mus.range, aes(x = log(km2_range), y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_range_uk, aes(x = log(range), ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_range_uk, aes(x = log(range), y = fit), 
              colour = "black", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(breaks = c(0.2, 0.1, 0, -0.1, -0.2), 
                       limits = c(-0.25, 0.25),
                       labels = c("0.2", "0.1", "0", "-0.1", "-0.2")) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', 'Log geographic range ' ~ (km^2))), 
         y = bquote(atop(mu), '')))

# Local population size (logged) versus population change
mcmc_preds_popsize_uk <- as.data.frame(predict(rarity.popsize.uk,
                                               type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_popsize_uk$log.meanpop <- uk.mus.popsize$log.meanpop

(pop.size.plot <- ggplot() +
    geom_point(data = uk.mus.popsize, aes(x = log.meanpop, y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_popsize_uk, aes(x = log.meanpop,
                                                  ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_popsize_uk, aes(x = log.meanpop, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                       limits = c(-0.235, 0.235),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    theme_LPI() + 
    labs(x = "\nLog population size", 
         y = bquote(atop(italic(mu), ' '))))

# Habitat specifity versus population change
mcmc_preds_habspec_uk <- as.data.frame(predict(rarity.hab.uk,
                                               type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_habspec_uk$hab_spec <- uk.mus.habspec$hab_spec

(hab.spec.plot <- ggplot() +
    geom_point(data = uk.mus.habspec, aes(x = hab_spec, y = mu),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_habspec_uk, aes(x = hab_spec,
                                                  ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_habspec_uk, aes(x = hab_spec, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2), 
                       limits = c(-0.235, 0.235),
                       labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    theme_LPI() + 
    labs(x = "\nNumber of habitats", 
         y = bquote(atop(italic(mu), ' '))))

# For fluctuations
# Geographic range
mcmc_preds_range_uk_fl <- as.data.frame(predict(rarity.range.uk.sigma,
                                             type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_range_uk_fl$range <- uk.mus.range$km2_range

(range.plot.fl <- ggplot() +
    geom_point(data = uk.mus.range, aes(x = log(km2_range), y = sigma.2),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_range_uk, aes(x = log(range), ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_range_uk, aes(x = log(range), y = fit), 
              colour = "black", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    theme_LPI() + 
    labs(x = bquote(atop(' ', 'Log geographic range ' ~ (km^2))), 
         y = bquote(atop(sigma^2), '')))

# Local population size (logged) versus population change
mcmc_preds_popsize_uk_fl <- as.data.frame(predict(rarity.popsize.uk.sigma,
                                               type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_popsize_uk_fl$log.meanpop <- uk.mus.popsize$log.meanpop

(pop.size.plot.fl <- ggplot() +
    geom_point(data = uk.mus.popsize, aes(x = log.meanpop, y = sigma.2.x),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_popsize_uk_fl, aes(x = log.meanpop,
                                                  ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_popsize_uk_fl, aes(x = log.meanpop, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    theme_LPI() + 
    labs(x = "\nLog population size", 
         y = bquote(atop(italic(sigma^2), ' '))))

# Habitat specifity versus population change
mcmc_preds_habspec_uk_fl <- as.data.frame(predict(rarity.hab.uk.sigma,
                                               type = "terms", interval = "confidence"))

# Add the range values to the predictions data frame
mcmc_preds_habspec_uk_fl$hab_spec <- uk.mus.habspec$hab_spec

(hab.spec.plot.fl <- ggplot() +
    geom_point(data = uk.mus.habspec, aes(x = hab_spec, y = sigma.2),
               size = 1, colour = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = mcmc_preds_habspec_uk_fl, aes(x = hab_spec,
                                                  ymin = lwr, ymax = upr), 
                fill = "black", alpha = 0.5) +
    geom_line(data = mcmc_preds_habspec_uk_fl, aes(x = hab_spec, y = fit), 
              colour = "black", size = 1) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), 
                       limits = c(-0.03, 0.8),
                       labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
    theme_LPI() + 
    labs(x = "\nNumber of habitats", 
         y = bquote(atop(italic(sigma^2), ' '))))

uk.rarity.panel <- grid.arrange(range.plot, pop.size.plot, hab.spec.plot, 
                                range.plot.fl, pop.size.plot.fl, hab.spec.plot.fl,
                                ncol = 3)
ggsave(uk.rarity.panel, filename = "figures/SI_UK_rarity.pdf",
       device = "pdf", dpi = 300, height = 10, width = 15)


# Declining populations

mus.neg <- filter(mus, uCI < 0)
mus.pos <- filter(mus, lCI > 0)

length(mus.neg$mu) / length(mus$mu) # 15% are declining

length(mus.pos$mu) / length(mus$mu) # 18% are increasing

100-15-18
# 67% no net change in abundance