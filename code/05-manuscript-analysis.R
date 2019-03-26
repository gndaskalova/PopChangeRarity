# How does vertebrate population change differ across biomes and taxa?
# How does rarity influence population change in UK vertebrates?

# Gergana Daskalova (gndaskalova@gmail.com)

# 17th Feb 2019

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
library(proj4)
library(ggalt)
library(RColorBrewer)
library(ggridges)
library(ape)

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

# The code below doesn't work on a Mac because of the text encoding
# The code was run on a Windows machine and the output saved
# LPI.long.pop2 <- LPI.long.pop %>% 
#   dplyr::select(Units, id, pop, Country.list) %>% 
#   filter(str_detect(tolower(Units), pattern = ("count"))) %>%
#   filter(!str_detect(tolower(Units), pattern = ("%")))
# 
# LPI.long.pop3 <- LPI.long.pop %>% 
#   dplyr::select(Units, id, pop, Country.list) %>% 
#   filter(str_detect(tolower(Units), pattern = ("individual"))) %>%
#   filter(!str_detect(tolower(Units), pattern = ("%")))
# 
# LPI.long.pop4 <- LPI.long.pop %>% 
#   dplyr::select(Units, id, pop, Country.list) %>% 
#   filter(str_detect(tolower(Units), pattern = ("number"))) %>%
#   filter(!str_detect(tolower(Units), pattern = ("%")))
# 
# LPI.long.pop5 <- LPI.long.pop %>% 
#   dplyr::select(Units, id, pop, Country.list) %>% 
#   filter(str_detect(tolower(Units), pattern = ("nest"))) %>%
#   filter(!str_detect(tolower(Units), pattern = ("%")))
# 
# LPI.long.pop6 <- LPI.long.pop %>% 
#   dplyr::select(Units, id, pop, Country.list) %>% 
#   filter(str_detect(tolower(Units), pattern = ("pair"))) %>%
#   filter(!str_detect(tolower(Units), pattern = ("%")))
# 
# LPI.long.pop.all <- rbind(LPI.long.pop2, LPI.long.pop3,
#                           LPI.long.pop3, LPI.long.pop4,
#                           LPI.long.pop5, LPI.long.pop6) %>%
#   distinct()
# 
# LPI.mean.pop <- LPI.long.pop.all %>%
#   group_by(Country.list, id) %>%   # group rows so that each group is one population
#   summarise(meanpop = mean(pop),
#             sdpop = sd(pop),
#             se = std.error(pop)) %>%  # Create column for mean population
#   ungroup() %>% dplyr::select(id, meanpop, sdpop) %>% distinct() %>%
#   drop_na(meanpop)

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

# Global models ----
# In models with categorical variables, the intercept was set to 0
# so that we can compare across all levels of the categorical variable
# otherwise e.g. Aves becomes the intercept and we we comparing
# mammals with birds, not whether each category has experienced net pop change

# ** mu ----
# Latitude model
latitude.m <- MCMCglmm(mu ~ Decimal.Latitude, data = mus, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(latitude.m)
# Check model convergence
# plot(latitude.m$VCV)
# autocorr(latitude.m$VCV)

# Duration model
duration.m <- MCMCglmm(mu ~ duration - 1, data = mus, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(duration.m)
# plot(duration.m$VCV)
# autocorr(duration.m$VCV)

# System model
system.m <- MCMCglmm(mu ~ system - 1, data = mus, random = ~species,
                     nitt = 200000, burnin = 20000)
summary(system.m)
# Check model convergence
# plot(system.m$VCV)
# autocorr(system.m$VCV)

# Test for multimodal distribution
mus.terr <- filter(mus, system == "Terrestrial")
dip.test(mus.terr$mu)
mus.fresh <- filter(mus, system == "Freshwater")
dip.test(mus.terr$mu)
mus.mar <- filter(mus, system == "Marine")
dip.test(mus.mar$mu)
# They are all at least bimodal

# Biome model
biome.m <- MCMCglmm(mu ~ biome - 1, data = mus.biome, random = ~species,
                    nitt = 200000, burnin = 20000)
summary(biome.m)
# Check model convergence
# plot(biome.m$VCV)
# autocorr(biome.m$VCV)

# Taxa model
class.m <- MCMCglmm(mu ~ Class - 1, data = mu.class, random = ~species,
                    nitt = 200000, burnin = 20000)
summary(class.m)
# Check model convergence
# plot(class.m$VCV)
# autocorr(class.m$VCV)

# IUCN model
IUCN.m <- MCMCglmm(mu ~ Red.List.status - 1, data = mus.IUCN, random = ~species,
                   nitt = 200000, burnin = 20000)
summary(IUCN.m)
# Check model convergence
# plot(IUCN.m$VCV)
# autocorr(IUCN.m$VCV)

# Model for type of threats
threat.type.m <- MCMCglmm(mu ~ title - 1, data = top_threats_mus, random = ~species,
                          nitt = 200000, burnin = 20000)
summary(threat.type.m)
# Check model convergence
# plot(threat.type.m$VCV)
# autocorr(threat.type.m$VCV)

# Number of threats and pop trend
threats_lpi2 <- filter(threats_lpi, !title %in% c("Unspecified species",
                                                  "Type Unknown/Unrecorded",
                                                  "Scale Unknown/Unrecorded",
                                                  "Motivation Unknown/Unrecorded"))

threats_sum_species <- threats_lpi2 %>% group_by(species) %>%
  tally() %>% arrange(desc(n))

threats_sum_species_mus <- inner_join(species_mus, threats_sum_species,
                                      by = "species")

# Model for number of threats
threat.n.m <- MCMCglmm(mu ~ n, data = threats_sum_species_mus, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(threat.n.m)
# Check model convergence
# plot(threat.n.m$VCV)
# autocorr(threat.n.m$VCV)

# Global hab spec model
rarity.hab.global <- MCMCglmm(mu ~ hab_spec, data = mus.habspec, random = ~species,
                       nitt = 200000, burnin = 20000) 
summary(rarity.hab.global)
# plot(rarity.hab.global$VCV)
# autocorr(rarity.hab.global$VCV)

# Global mean pop size model
rarity.popsize.global <- MCMCglmm(mu ~ log.meanpop, data = mus.popsize, random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global)
# plot(rarity.popsize.global$VCV)
# autocorr(rarity.popsize.global$VCV)

# Birds and mammals
b_m_ranges$log.range <- log(b_m_ranges$range)
b_m_ranges <- b_m_ranges %>% drop_na(log.range)

rarity.range.birds.mammals <- MCMCglmm(mu ~ log(range), data = b_m_ranges, 
                               nitt = 200000, burnin = 20000) 
summary(rarity.range.birds.mammals)
# plot(rarity.range.birds.mammals$VCV)
# autocorr(rarity.range.birds.mammals$VCV)

# Taxa interaction
# Global hab spec model
rarity.hab.global.taxa <- MCMCglmm(mu ~ hab_spec*Class, data = mus.habspec, random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.taxa)
# plot(rarity.hab.global.taxa$VCV)
# autocorr(rarity.hab.global.taxa$VCV)

# Global mean pop size model
# Add taxa back to the data frame

taxa <- mus %>% dplyr::select(id, Class) %>% distinct()
mus.popsize <- left_join(mus.popsize, taxa, by = "id")

rarity.popsize.global.taxa <- MCMCglmm(mu ~ log.meanpop*Class, data = mus.popsize, random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.taxa)
# plot(rarity.popsize.global.taxa$VCV)
# autocorr(rarity.popsize.global.taxa$VCV)

# Interaction for mammals and birds with geographic range
rarity.range.birds.mammals.int <- MCMCglmm(mu ~ log(range)*Class, data = b_m_ranges,  random = ~species,
                                       nitt = 200000, burnin = 20000) 
summary(rarity.range.birds.mammals.int)
# plot(rarity.range.birds.mammals.int$VCV)
# autocorr(rarity.range.birds.mammals.int$VCV)

# Methods model
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
                                      "Pairs (505)", "Nests(224)", 
                                      "Population estimate (97)"))

units.m <- MCMCglmm(mu ~ Units2, data = mus.methods, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(units.m)
# Check model convergence
# plot(units.m$VCV)
# autocorr(units.m$VCV)

# ** mu by tau ----
# weighted by observation error
# System model
weight1 <- mus$tau.2
system.m.var <- MCMCglmm(mu ~ system - 1, mev = weight1, data = mus, random = ~species,
                         nitt = 200000, burnin = 20000)
summary(system.m.var)
# Check model convergence
# plot(system.m.var$VCV)
# autocorr(system.m.var$VCV)

# Biome model
weight2 <- mus.biome$tau.2
biome.m.var <- MCMCglmm(mu ~ biome - 1, mev = weight2,  data = mus.biome,  random = ~species,
                        nitt = 200000, burnin = 20000)
summary(biome.m.var)
# Check model convergence
# plot(biome.m.var$VCV)
# autocorr(biome.m.var$VCV)

# Taxa model
weight3 <- mu.class$tau.2
class.m.var <- MCMCglmm(mu ~ Class - 1, mev = weight3, data = mu.class,  random = ~species,
                        nitt = 200000, burnin = 20000)
summary(class.m.var)
# Check model convergence
# plot(class.m.var$VCV)
# autocorr(class.m.var$VCV)

# IUCN model
weight4 <- mus.IUCN$tau.2
IUCN.m.var <- MCMCglmm(mu ~ Red.List.status - 1, mev = weight4, data = mus.IUCN,  random = ~species,
                       nitt = 200000, burnin = 20000)
summary(IUCN.m.var)
# Check model convergence
# plot(IUCN.m.var$VCV)
# autocorr(IUCN.m.var$VCV)

# Global hab spec model
weight5 <- mus.habspec$tau.2
rarity.hab.global.var <- MCMCglmm(mu ~ hab_spec, mev = weight5, data = mus.habspec,  random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.var)
# plot(rarity.hab.global.var$VCV)
# autocorr(rarity.hab.global.var$VCV)

# Global mean pop size model
weight6 <- mus.popsize$tau.2
rarity.popsize.global.var <- MCMCglmm(mu ~ log.meanpop, mev =  weight6, data = mus.popsize,  random = ~species,
                                      nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.var)
# plot(rarity.popsize.global.var$VCV)
# autocorr(rarity.popsize.global.var$VCV)

# Global bird and mammal range model
weight7 <- b_m_ranges$tau.2
rarity.range.birds.mammals.var <- MCMCglmm(mu ~ log(range), mev = weight7, data = b_m_ranges,  random = ~species,
                               nitt = 200000, burnin = 20000) 
summary(rarity.range.birds)
# plot(rarity.range.birds.var$VCV)
# autocorr(rarity.range.birds.var$VCV)

# ** slope ----
# System model
system.m2 <- MCMCglmm(slope ~ system - 1, data = slopes,  random = ~species,
                      nitt = 200000, burnin = 20000)
summary(system.m2)
# Check model convergence
# plot(system.m2$VCV)
# autocorr(system.m2$VCV)

# Biome model
biome.m2 <- MCMCglmm(slope ~ biome - 1, data = slopes.biome,  random = ~species,
                     nitt = 200000, burnin = 20000)
summary(biome.m2)
# Check model convergence
# plot(biome2.m$VCV)
# autocorr(biome.m2$VCV)

# Taxa model
class.m2 <- MCMCglmm(slope ~ Class - 1, data = slope.class,  random = ~species,
                     nitt = 200000, burnin = 20000)
summary(class.m2)
# Check model convergence
# plot(class.m2$VCV)
# autocorr(class.m2$VCV)

# IUCN model
IUCN.m2 <- MCMCglmm(slope ~ Red.List.status - 1, data = slopes.IUCN,  random = ~species,
                    nitt = 200000, burnin = 20000)
summary(IUCN.m2)
# Check model convergence
# plot(IUCN.m2$VCV)
# autocorr(IUCN.m2$VCV)

# Global hab spec model
rarity.hab.global2 <- MCMCglmm(slope ~ hab_spec, data = slopes.habspec,  random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.hab.global2)
# plot(rarity.hab.global2$VCV)
# autocorr(rarity.hab.global2$VCV)

# Global mean pop size model
slopes.popsize <- na.omit(slopes.popsize)
slopes.popsize$log.meanpop <- log(slopes.popsize$meanpop)
slopes.popsize <- slopes.popsize %>% filter(is.finite(log.meanpop))

rarity.popsize.global2 <- MCMCglmm(slope ~ log.meanpop, data = slopes.popsize, random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global2)
# plot(rarity.popsize.global2$VCV)
# autocorr(rarity.popsize.global2$VCV)

rarity.range.birds.mammals2 <- MCMCglmm(slope ~ log.range, data = b_m_slopes_ranges, random = ~species,
                                   nitt = 200000, burnin = 20000) 
summary(rarity.range.birds2)
# plot(rarity.range.birds2$VCV)
# autocorr(rarity.range.birds2$VCV)

# ** CI ----
# confidence intervals around mu

# System model
mus$CI <- (mus$uCI - mus$lCI)/2
mus$CI <- rescale(mus$CI, to = c(0, 1))
system.m.ci <- MCMCglmm(CI ~ system - 1, data = mus, random = ~species,
                         nitt = 200000, burnin = 20000)
summary(system.m.ci)
# Check model convergence
# plot(system.m.ci$VCV)
# autocorr(system.m.ci$VCV)

system.ci <- mus %>% group_by(system) %>%
  mutate(sd = sd(CI)) %>% select(system, sd) %>% distinct()

system.ci$Model <- "CI"

mus$system <- factor(mus$system, 
                         levels = c("Terrestrial", "Marine", "Freshwater"),
                         labels = c("Terrestrial", "Marine", "Freshwater"))

duration.system.m <- MCMCglmm(duration ~ system, data = mus, random = ~species,
                                nitt = 200000, burnin = 20000)

summary(duration.system.m)


# Biome model
mus.biome$CI <- (mus.biome$uCI - mus.biome$lCI)/2
mus.biome$CI <- rescale(mus.biome$CI, to = c(0, 1))

biome.m.ci <- MCMCglmm(CI ~ biome - 1,  data = mus.biome, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(biome.m.ci)
# Check model convergence
# plot(biome.m.ci$VCV)
# autocorr(biome.m.ci$VCV)

# Taxa model
mu.class$CI <- (mu.class$uCI - mu.class$lCI)/2
mu.class$CI <- scales::rescale(mu.class$CI, to = c(0, 1))

class.m.ci <- MCMCglmm(CI ~ Class - 1, data = mu.class, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(class.m.ci)
# Check model convergence
# plot(class.m.ci$VCV)
# autocorr(class.m.ci$VCV)

class.ci <- mu.class %>% group_by(Class) %>%
  mutate(sd = sd(CI)) %>% select(Class, sd) %>% distinct()

class.ci$Model <- "CI"

# Check duration across taxa
duration.class.m <- MCMCglmm(duration ~ Class, data = mu.class, random = ~species,
                       nitt = 200000, burnin = 20000)

summary(duration.class.m)

# IUCN model
mus.IUCN$CI <- (mus.IUCN$uCI - mus.IUCN$lCI)/2
mus.IUCN$CI <- rescale(mus.IUCN$CI, to = c(0, 1))

IUCN.m.ci <- MCMCglmm(CI ~ Red.List.status - 1, data = mus.IUCN, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(IUCN.m.ci)
# Check model convergence
# plot(IUCN.m.ci$VCV)
# autocorr(IUCN.m.ci$VCV)

IUCN.ci <- mus.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(CI)) %>% select(Red.List.status, sd) %>% distinct()

IUCN.ci$Model <- "CI"

duration.IUCN.m <- MCMCglmm(duration ~ Red.List.status, data = mus.IUCN, random = ~species,
                                 nitt = 200000, burnin = 20000)
summary(duration.IUCN.m)


# Global hab spec model
mus.habspec$CI <- (mus.habspec$uCI - mus.habspec$lCI)/2
mus.habspec$CI <- rescale(mus.habspec$CI, to = c(0, 1))

rarity.hab.global.ci <- MCMCglmm(CI ~ hab_spec, data = mus.habspec, random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.ci)
# plot(rarity.hab.global.ci$VCV)
# autocorr(rarity.hab.global.ci$VCV)

# Global mean pop size model
mus.popsize$CI <- (mus.popsize$uCI - mus.popsize$lCI)/2
mus.popsize$CI <- rescale(mus.popsize$CI, to = c(0, 1))

rarity.popsize.global.ci <- MCMCglmm(CI ~ log.meanpop, data = mus.popsize, random = ~species,
                                      nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.ci)
# plot(rarity.popsize.global.ci$VCV)
# autocorr(rarity.popsize.global.ci$VCV)

# Global bird and mammal range model
b_m_ranges$CI <- (b_m_ranges$uCI - b_m_ranges$lCI)/2
b_m_ranges$CI <- rescale(b_m_ranges$CI, to = c(0, 1))

bird.mammal.range.m.ci <- MCMCglmm(CI ~ log(range), data = b_m_ranges, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(bird.mammal.range.m.ci)
# Check model convergence
# plot(bird.mammal.range.m.ci$VCV)
# autocorr(bird.mammal.range.m.ci$VCV)

# ** CI by tau ----
# System model
system.m.ci.var <- MCMCglmm(CI ~ system - 1, data = mus, mev = weight1, random = ~species,
                            nitt = 200000, burnin = 20000)
summary(system.m.ci.var)
# Check model convergence
# plot(system.m.ci.var$VCV)
# autocorr(system.m.ci.var$VCV)

# Data frame for the sd of the input data
system.ci <- mus %>% group_by(system) %>%
  mutate(sd = sd(CI)) %>% select(system, sd) %>% distinct()

system.ci.var <- system.ci
system.ci.var$Model <- "CI.w"

# Biome model
biome.m.ci.var <- MCMCglmm(CI ~ biome - 1,  data = mus.biome,  mev = weight2, random = ~species,
                           nitt = 200000, burnin = 20000)
summary(biome.m.ci.var)
# Check model convergence
# plot(biome.m.ci.var$VCV)
# autocorr(biome.m.ci.var$VCV)

# Taxa model
class.m.ci.var <- MCMCglmm(CI ~ Class - 1, data = mu.class, mev = weight3, random = ~species,
                           nitt = 200000, burnin = 20000)
summary(class.m.ci.var)
# Check model convergence
# plot(class.m.ci.var$VCV)
# autocorr(class.m.ci.var$VCV)

# Data frame for the sd of the input data
# same as non-weighted data frame, since they are both
# based on the same input data (mu values)
class.ci.var <- class.ci
class.ci.var$Model <- "CI.w"

# IUCN model
IUCN.m.ci.var <- MCMCglmm(CI ~ Red.List.status - 1, mev = weight4, data = mus.IUCN, random = ~species,
                          nitt = 200000, burnin = 20000)
summary(IUCN.m.ci.var)
# Check model convergence
# plot(IUCN.m.ci.var$VCV)
# autocorr(IUCN.m.ci.var$VCV)

# Data frame for the sd of the input data
# same as non-weighted data frame, since they are both
# based on the same input data (mu values)
IUCN.ci.var <- IUCN.ci
IUCN.ci.var$Model <- "CI.w"

# Global hab spec model
rarity.hab.global.ci.var <- MCMCglmm(CI ~ hab_spec, mev = weight5, data = mus.habspec, random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.ci.var)
# plot(rarity.hab.global.ci.var$VCV)
# autocorr(rarity.hab.global.ci.var$VCV)

# Global mean pop size model
rarity.popsize.global.ci.var <- MCMCglmm(CI ~ log.meanpop, mev =  weight6, data = mus.popsize, random = ~species,
                                      nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.ci.var)
# plot(rarity.popsize.global.ci.var$VCV)
# autocorr(rarity.popsize.global.ci.var$VCV)

# Bird range model
bird.mammal.range.m.ci.var <- MCMCglmm(CI ~ log(range), mev = weight7, data = b_m_ranges, random = ~species,
                            nitt = 200000, burnin = 20000)
summary(bird.range.m.ci.var)
# Check model convergence
# plot(bird.range.m.ci.var$VCV)
# autocorr(bird.range.m.ci.var$VCV)

# ** SE ----
# Using the se from the slope estimates
# System model
slopes$SE <- rescale(slopes$slope_se, to = c(0, 1))
system.m.se <- MCMCglmm(SE ~ system - 1, data = slopes, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(system.m.se)
# Check model convergence
# plot(system.m.se$VCV)
# autocorr(system.m.se$VCV)

system.se <- slopes %>% group_by(system) %>%
  mutate(sd = sd(SE)) %>% select(system, sd) %>% distinct()
system.se$Model <- "se"

# Biome model
slopes.biome$SE <- rescale(slopes.biome$slope_se, to = c(0, 1))

biome.m.se <- MCMCglmm(SE ~ biome - 1,  data = slopes.biome, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(biome.m.se)
# Check model convergence
# plot(biome.m.se$VCV)
# autocorr(biome.m.se$VCV)

# Taxa model
slope.class$SE <- rescale(slope.class$slope_se, to = c(0, 1))

class.m.se <- MCMCglmm(SE ~ Class - 1, data = slope.class, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(class.m.se)
# Check model convergence
# plot(class.m.se$VCV)
# autocorr(class.m.se$VCV)

class.se <- slope.class %>% group_by(Class) %>%
  mutate(sd = sd(SE)) %>% select(Class, sd) %>% distinct()
class.se$Model <- "se"

# IUCN model
slopes.IUCN$SE <- rescale(slopes.IUCN$slope_se, to = c(0, 1))

IUCN.m.se <- MCMCglmm(SE ~ Red.List.status - 1, data = slopes.IUCN, random = ~species,
                      nitt = 200000, burnin = 20000)
summary(IUCN.m.se)
# Check model convergence
# plot(IUCN.m.se$VCV)
# autocorr(IUCN.m.se$VCV)
IUCN.se <- slopes.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(SE)) %>% select(Red.List.status, sd) %>% distinct()
IUCN.se$Model <- "se"

# Global hab spec model
slopes.habspec$SE <- rescale(slopes.habspec$slope_se, to = c(0, 1))

rarity.hab.global.se <- MCMCglmm(SE ~ hab_spec, data = slopes.habspec, random = ~species,
                                     nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.se)
# plot(rarity.hab.global.se$VCV)
# autocorr(rarity.hab.global.se$VCV)

# Global mean pop size model
slopes.popsize$SE <- rescale(slopes.popsize$slope_se, to = c(0, 1))

rarity.popsize.global.se <- MCMCglmm(SE ~ log.meanpop, data = slopes.popsize, random = ~species,
                                         nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.se)
# plot(rarity.popsize.global.se$VCV)
# autocorr(rarity.popsize.global.se$VCV)

# Global bird range model
b_m_slopes_ranges$SE <- rescale(b_m_slopes_ranges$slope_se, to = c(0, 1))

bird.mammal.range.m.se <- MCMCglmm(SE ~ log(range), data = b_m_slopes_ranges, random = ~species,
                            nitt = 200000, burnin = 20000)
summary(bird.mammal.range.m.se)
# Check model convergence
# plot(bird.mammal.range.m.se$VCV)
# autocorr(bird.mammal.range.m.se$VCV)

# ** SD ----
# Using the sd of the scaled raw pop data

# System model
LPI.long <- LPI.long %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))

LPI.long$sdpop1 <- scales::rescale(LPI.long$sdpop, to = c(0, 1))
LPIraw <- LPI.long %>% select(id, species, sdpop, sdpop1, Country.list, 
                              system, biome, Class) %>%
  distinct()

system.m.sd <- MCMCglmm(sdpop1 ~ system - 1, data = LPIraw, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(system.m.sd)
# Check model convergence
# plot(system.m.sd$VCV)
# autocorr(system.m.sd$VCV)

system.sd <- LPI.long %>% group_by(system) %>%
  mutate(sd = sd(scalepop)) %>% dplyr::select(system, sd) %>%
  distinct()
system.sd$Model <- "sd"

# Biome model
LPI.long.biome <- filter(LPI.long, biome != "Unknown" & biome != "Mangroves" & 
                             biome != "Oceanic islands")

LPI.long.biome <- LPI.long.biome %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))
LPI.long.biome$sdpop1 <- rescale(LPI.long.biome$sdpop, to = c(0, 1))
LPIraw.biome <- LPI.long.biome %>% select(id, species, sdpop, sdpop1, Country.list, 
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

biome.m.sd <- MCMCglmm(sdpop1 ~ biome - 1,  data = LPIraw.biome, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(biome.m.sd)
# Check model convergence
# plot(biome.m.sd$VCV)
# autocorr(biome.m.sd$VCV)

# Taxa model
LPI.long.class <- filter(LPI.long, Class== "Aves" | Class == "Mammalia" | 
                             Class == "Reptilia" | Class == "Amphibia" |
                             Class == "Actinopterygii" | Class == "Elasmobranchii")
LPI.long.class <- LPI.long.class %>% group_by(id) %>%
  mutate(sdpop = sd(scalepop))

LPI.long.class$sdpop1 <- rescale(LPI.long.class$sdpop, to = c(0, 1))

LPIraw.class <- LPI.long.class %>% select(id, species, sdpop, sdpop1, Country.list, 
                                          system, biome, Class) %>%
  distinct()


class.m.sd <- MCMCglmm(sdpop1 ~ Class - 1, data = LPIraw.class, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(class.m.sd)
# Check model convergence
# plot(class.m.sd$VCV)
# autocorr(class.m.sd$VCV)

class.sd <- LPI.long.class %>% group_by(Class) %>%
  mutate(sd = sd(scalepop)) %>% dplyr::select(Class, sd) %>%
  distinct()
class.sd$Model <- "sd"

# IUCN model
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

# Calculate the sd of the raw data for each category
IUCN.sd <- inner_join(LPI.long, IUCN, by = "species")
# Exclude DD (Data deficient species) and EW (Extinct in the wild)
IUCN.sd <- filter(IUCN.sd, Red.List.status != "DD" & 
                        Red.List.status != "EW")

# Reorder levels to match increasing risk
IUCN.sd$Red.List.status <- factor(IUCN.sd$Red.List.status, 
                                      levels = c("LC", "NT", "VU", "EN", "CR"),
                                      labels = c("Least concern", "Near threatened", 
                                                 "Vulnerable", "Endangered",
                                                 "Critically endangered"))

IUCN.sd <- IUCN.sd %>% drop_na(Red.List.status)
IUCN.sd <- IUCN.sd %>% group_by(Red.List.status) %>%
  mutate(sd = sd(scalepop)) %>% dplyr::select(Red.List.status, sd) %>%
  distinct()

IUCN.sd$Model <- "sd"

# Filtering the data frame for the model
LPIraw.IUCN$sdpop1 <- rescale(LPIraw.IUCN$sdpop, to = c(0, 1))
LPIraw.IUCN <- LPIraw.IUCN %>% select(id, species, sdpop, sdpop1, Country.list, 
                                          system, biome, Class, Red.List.status) %>%
  distinct()

IUCN.m.sd <- MCMCglmm(sdpop1 ~ Red.List.status - 1, data = LPIraw.IUCN, random = ~species,
                      nitt = 200000, burnin = 20000)
summary(IUCN.m.sd)
# Check model convergence
# plot(IUCN.m.sd$VCV)
# autocorr(IUCN.m.sd$VCV)

# Global hab spec model
LPIraw.habspec <- inner_join(LPI.long, habspec, by = "species")

LPIraw.habspec$sdpop1 <- rescale(LPIraw.habspec$sdpop, to = c(0, 1))
LPIraw.habspec <- LPIraw.habspec %>% select(id, species, sdpop, sdpop1, Country.list, 
                                      system, biome, Class, hab_spec) %>%
  distinct()

rarity.hab.global.sd <- MCMCglmm(sdpop1 ~ hab_spec, data = LPIraw.habspec, random = ~species,
                                 nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.sd)
# plot(rarity.hab.global.sd$VCV)
# autocorr(rarity.hab.global.sd$VCV)

# Global mean pop size model
LPI.mean.pop2 <- LPI.mean.pop
id_species <- mus.popsize %>% dplyr::select(species, id)
LPI.mean.pop2 <- left_join(LPI.mean.pop2, id_species, by = "id")
LPI.mean.pop2 <- LPI.mean.pop2 %>% dplyr::select(id, sdpop, meanpop, species) %>%
  distinct() %>% drop_na(species)
LPI.mean.pop2$sdpop1 <- rescale(LPI.mean.pop2$sdpop, to = c(0, 1))
LPI.mean.pop2$log.meanpop <- log(LPI.mean.pop2$meanpop)
LPI.mean.pop2 <- LPI.mean.pop2 %>% filter(is.finite(log.meanpop))

prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, 
                                  alpha.mu = 0, alpha.v = 10000)))

rarity.popsize.global.sd <- MCMCglmm(sdpop1 ~ log.meanpop, data = LPI.mean.pop2, random = ~species,
                                     nitt = 200000, burnin = 20000, prior = prior1) 
summary(rarity.popsize.global.sd)
# plot(rarity.popsize.global.sd$VCV)
# autocorr(rarity.popsize.global.sd$VCV)

# Global bird range model
ranges_simple <- b_m_ranges %>% dplyr::select(species, range) %>% distinct()
LPIraw.bird.mammal.range <- left_join(LPI.long, ranges_simple, by = "species")
LPIraw.bird.mammal.range$sdpop1 <- rescale(LPIraw.bird.mammal.range$sdpop, to = c(0, 1))
LPIraw.bird.mammal.range <- LPIraw.bird.mammal.range %>% dplyr::select(id, sdpop1, range, species) %>%
  distinct()
LPIraw.bird.mammal.range <- filter(LPIraw.bird.mammal.range, range != "EXTINCT")
LPIraw.bird.mammal.range$range <- as.numeric(as.character(LPIraw.bird.mammal.range$range))
LPIraw.bird.mammal.range <- LPIraw.bird.mammal.range %>% 
  mutate(log.range = log(range)) %>% drop_na(log.range)

rarity.bird.mammal.range.sd <- MCMCglmm(sdpop1 ~ log(range), data = LPIraw.bird.mammal.range, random = ~species,
                                     nitt = 200000, burnin = 20000) 
summary(rarity.bird.mammal.range.sd)
# plot(rarity.bird.mammal.range.sd$VCV)
# autocorr(rarity.bird.mammal.range.sd$VCV)

# ** sigma ----
# Duration model
mus$sigma.2 <- rescale(mus$sigma.2, to = c(0, 1))
duration.m.sigma <- MCMCglmm(sigma.2 ~ duration - 1, data = mus, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(duration.m)
# plot(duration.m$VCV)
# autocorr(duration.m$VCV)

# System model
system.m.sigma <- MCMCglmm(sigma.2 ~ system - 1, data = mus, random = ~species,
                            nitt = 200000, burnin = 20000)
summary(system.m.sigma)
# Check model convergence
# plot(system.m.sigma$VCV)
# autocorr(system.m.sigma$VCV)

# Data frame for the sd of the input data
system.ci.sigma <- mus %>% group_by(system) %>%
  mutate(sd = sd(sigma.2)) %>% select(system, sd) %>% distinct()

system.ci.sigma$Model <- "sigma"

# Taxa model
mu.class$sigma.2 <- rescale(mu.class$sigma.2, to = c(0, 1))
class.m.sigma <- MCMCglmm(sigma.2 ~ Class - 1, data = mu.class,  random = ~species,
                           nitt = 200000, burnin = 20000)
summary(class.m.sigma)
# Check model convergence
# plot(class.m.sigma$VCV)
# autocorr(class.m.sigma$VCV)

# Data frame for the sd of the input data
taxa.ci.sigma <- mu.class %>% group_by(Class) %>%
  mutate(sd = sd(sigma.2)) %>% select(Class, sd) %>% distinct()

taxa.ci.sigma$Model <- "sigma"

# Biome model
mus.biome$sigma.2 <- rescale(mus.biome$sigma.2, to = c(0, 1))
biome.m.sigma <- MCMCglmm(sigma.2 ~ biome - 1, data = mus.biome,  random = ~species,
                    nitt = 200000, burnin = 20000)
summary(biome.m.sigma)
# Check model convergence
# plot(biome.m.sigma$VCV)
# autocorr(biome.m.sigma$VCV)

# Data frame for the sd of the input data
biome.ci.sigma <- mus.biome %>% group_by(biome) %>%
  mutate(sd = sd(sigma.2)) %>% select(biome, sd) %>% distinct()

biome.ci.sigma$Model <- "sigma"

# Global hab spec model
mus.habspec$sigma.2 <- rescale(mus.habspec$sigma.2, to = c(0, 1))
rarity.hab.global.sigma <- MCMCglmm(sigma.2 ~ hab_spec, data = mus.habspec,  random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.hab.global.sigma)
# plot(rarity.hab.global$VCV)
# autocorr(rarity.hab.global$VCV)

# Global mean pop size model
mus.popsize$sigma.2 <- rescale(mus.popsize$sigma.2, to = c(0, 1))
rarity.popsize.global.sigma <- MCMCglmm(sigma.2 ~ log.meanpop, data = mus.popsize,  random = ~species,
                                  nitt = 200000, burnin = 20000) 
summary(rarity.popsize.global.sigma)
# plot(rarity.popsize.global$VCV)
# autocorr(rarity.popsize.global$VCV)

# Global bird and mammal range
b_m_ranges$sigma.2 <- rescale(b_m_ranges$sigma.2, to = c(0, 1))

rarity.range.birds.mammals.sigma <- MCMCglmm(sigma.2 ~ log(range), data = b_m_ranges,  random = ~species,
                               nitt = 200000, burnin = 20000) 
summary(rarity.range.birds.mammals.sigma)
# plot(rarity.range.birds.mammals.sigma$VCV)
# autocorr(rarity.range.birds.mammals.sigma$VCV)

# IUCN model
mus.IUCN$sigma.2 <- rescale(mus.IUCN$sigma.2, to = c(0, 1))
IUCN.m.ci.sigma <- MCMCglmm(sigma.2 ~ Red.List.status - 1, data = mus.IUCN,  random = ~species,
                      nitt = 200000, burnin = 20000)
summary(IUCN.m.ci.sigma)
# Check model convergence
# plot(IUCN.m.ci.sigma$VCV)
# autocorr(IUCN.m.ci.sigma$VCV)

IUCN.ci.sigma <- mus.IUCN %>% group_by(Red.List.status) %>%
  mutate(sd = sd(sigma.2)) %>% select(Red.List.status, sd) %>% distinct()

IUCN.ci.sigma$Model <- "sigma"

# Model for type of threats
threat.type.m.sigma <- MCMCglmm(sigma.2 ~ title - 1, data = top_threats_mus, random = ~species,
                          nitt = 200000, burnin = 20000)
summary(threat.type.m.sigma)
# Check model convergence
# plot(threat.type.m.sigma$VCV)
# autocorr(threat.type.m.sigma$VCV)

# Model for number of threats
threat.n.m.sigma <- MCMCglmm(sigma.2 ~ n, data = threats_sum_species_mus, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(threat.n.m.sigma)
# Check model convergence
# plot(threat.type.m.sigma$VCV)
# autocorr(threat.type.m.sigma$VCV)

# ** Summary output table ----
# Save all model outputs
model.list.global <- list(latitude.m, system.m, system.m.var, system.m2, system.m.sigma,
                  system.m.ci, system.m.ci.var,
                  system.m.se, system.m.sd,
                  biome.m, biome.m.var, biome.m2, biome.m.sigma,
                  biome.m.ci, biome.m.ci.var,
                  biome.m.se, biome.m.sd,
                  class.m, class.m.var, class.m2, class.m.sigma,
                  class.m.ci, class.m.ci.var,
                  class.m.se, class.m.sd,
                  duration.m, duration.m.sigma, duration.system.m,
                  duration.class.m, duration.IUCN.m, units.m,
                  rarity.range.birds.mammals, rarity.range.birds.mammals.int,
                  rarity.range.birds.mammals.var,
                  rarity.range.birds.mammals2, 
                  rarity.range.birds.mammals.sigma,
                  bird.mammal.range.m.ci, bird.mammal.range.m.ci.var,
                  bird.mammal.range.m.se, rarity.bird.mammal.range.sd,
                  rarity.popsize.global, rarity.popsize.global.taxa,
                  rarity.popsize.global.var, rarity.popsize.global2,
                  rarity.popsize.global.sigma,
                  rarity.popsize.global.ci, rarity.popsize.global.ci.var,
                  rarity.popsize.global.se, rarity.popsize.global.sd,
                  rarity.hab.global, rarity.hab.global.taxa, 
                  rarity.hab.global.var, rarity.hab.global2,
                  rarity.hab.global.sigma,
                  rarity.hab.global.ci, rarity.hab.global.ci.var,
                  rarity.hab.global.se, rarity.hab.global.sd,
                  IUCN.m, IUCN.m.var, IUCN.m2, IUCN.m.ci.sigma,
                  IUCN.m.ci, IUCN.m.ci.var,
                  IUCN.m.se, IUCN.m.sd,
                  threat.type.m, threat.n.m,
                  threat.type.m.sigma, threat.n.m.sigma)

names(model.list.global) <- c("latitude.m", "system.m", "system.m.var",
                               "system.m2", "system.m.sigma",
                               "system.m.ci", "system.m.ci.var",
                               "system.m.se", "system.m.sd",
                               "biome.m", "biome.m.var", "biome.m2", "biome.m.sigma",
                               "biome.m.ci", "biome.m.ci.var",
                               "biome.m.se", "biome.m.sd",
                               "class.m", "class.m.var", "class.m2", "class.m.sigma",
                               "class.m.ci", "class.m.ci.var",
                               "class.m.se", "class.m.sd",
                               "duration.m", "duration.m.sigma", "duration.system.m",
                               "duration.class.m", "duration.IUCN.m", "units.m",
                               "rarity.range.birds.mammals", "rarity.range.birds.mammals.int",
                               "rarity.range.birds.mammals.var",
                               "rarity.range.birds.mammals2", 
                               "rarity.range.birds.mammals.sigma",
                               "bird.mammal.range.m.ci", "bird.mammal.range.m.ci.var",
                               "bird.mammal.range.m.se", "rarity.bird.mammal.range.sd",
                               "rarity.popsize.global", "rarity.popsize.global.taxa",
                               "rarity.popsize.global.var", "rarity.popsize.global2",
                               "rarity.popsize.global.sigma",
                               "rarity.popsize.global.ci", "rarity.popsize.global.ci.var",
                               "rarity.popsize.global.se", "rarity.popsize.global.sd",
                               "rarity.hab.global", "rarity.hab.global.taxa", "rarity.hab.global.var",
                               "rarity.hab.global2",
                               "rarity.hab.global.sigma",
                               "rarity.hab.global.ci", "rarity.hab.global.ci.var",
                               "rarity.hab.global.se", "rarity.hab.global.sd",
                               "IUCN.m", "IUCN.m.var", "IUCN.m2", "IUCN.m.ci.sigma",
                               "IUCN.m.ci", "IUCN.m.ci.var",
                               "IUCN.m.se", "IUCN.m.sd",
                               "threat.type.m", "threat.n.m",
                               "threat.type.m.sigma", "threat.n.m.sigma")

# save(model.list.global, file = "data/output/global_models_18thFeb.RData")

# Create a list of input model names (global scale)
dataListNames <- list("Latitude - mu", 
                      "Realm - mu", 
                       "Realm - weighted", 
                       "Realm - slope", 
                      "Realm - fluctuations sigma", 
                       "Realm - fluctuations CI",                       
                       "Realm - fluctuations CI weighted", 
                       "Realm - fluctuations SE", 
                       "Realm - fluctuations SD", 
                       "Biome - mu", 
                       "Biome - slope", 
                       "Biome - weighted", 
                      "Biome - fluctuations sigma",
                       "Biome - fluctuations CI",
                       "Biome - fluctuations CI weighted",
                       "Biome - fluctuations SE",
                       "Biome - fluctuations SD",
                       "Taxa - mu", 
                       "Taxa - weighted", 
                       "Taxa - slope", 
                      "Taxa - fluctuations sigma", 
                       "Taxa - fluctuations CI", 
                       "Taxa - fluctuations CI weighted",
                       "Taxa - fluctuations SE", 
                       "Taxa - fluctuations SD", 
                      "Duration - mu",
                      "Duration - sigma",
                      "Duration - system",
                      "Duration - taxa",
                      "Duration - IUCN conservation status",
                      "Sampling units - mu",
                       "Geographic range (birds/mammals) - mu",
                      "Geographic range (birds/mammals) - mu*taxa interaction",
                       "Geographic range (birds/mammals) - weighted",
                       "Geographic range (birds/mammals) - slope",
                      "Geographic range (birds/mammals) - fluctuations sigma",
                       "Geographic range (birds/mammals) - fluctuations CI", 
                       "Geographic range (birds/mammals) - fluctuations CI weighted",
                       "Geographic range (birds/mammals)  - fluctuations SE",
                       "Geographic range (birds/mammals) - fluctuations SD",
                       "Mean population size - mu", 
                      "Mean population size - mu*taxa interaction",
                       "Mean population size - weighted", 
                       "Mean population size - slope", 
                      "Mean population size - fluctuations sigma",
                       "Mean population size - fluctuations CI", 
                       "Mean population size - fluctuations CI weighted",
                       "Mean population size - fluctuations SE",
                       "Mean population size - fluctuations SD",
                       "Habitat specificity - mu", 
                      "Habitat specificity - mu*taxa interaction",
                       "Habitat specificity - weighted", 
                       "Habitat specificity - slope", 
                      "Habitat specificity - fluctuations sigma",
                       "Habitat specificity - fluctuations CI", 
                       "Habitat specificity - fluctuations CI weighted",
                       "Habitat specificity - fluctuations SE", 
                       "Habitat specificity - fluctuations SD",
                       "IUCN conservation status - mu", 
                       "IUCN conservation status - weighted", 
                       "IUCN conservation status - slope", 
                      "IUCN conservation status - fluctuations sigma",
                       "IUCN conservation status - fluctuations CI", 
                       "IUCN conservation status - fluctuations CI weighted", 
                       "IUCN conservation status - fluctuations SE", 
                       "IUCN conservation status - fluctuations SD",
                      "IUCN threat type - mu",
                      "IUCN threat number - mu",
                      "IUCN threat type - fluctuations sigma",
                      "IUCN threat number - fluctuations sigma")
 
# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(model.list.global, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)
  
# Turn the list of data.frames into one big data.frame
mcmc.outputs.all.global <- as.data.frame(do.call(rbind, readyList))
 
mcmc.outputs.all.global <- mcmc.outputs.all.global %>% dplyr::select(modelName, variable:effect)
  
colnames(mcmc.outputs.all.global) <- c("Model name", "Variable", "Posterior mean",
                                   "Lower 95% CI", "Upper 95% CI", "Effective sample size",
                                   "pMCMC", "Effect")
  
write.csv(mcmc.outputs.all.global, file = "data/output/global_outputs18thFeb.csv")
  
mcmc.outputs.all.global$`Model name`[duplicated(mcmc.outputs.all.global$`Model name`)] <- " "
  
mcmc.outputs.all.global$`Variable` <- gsub("system", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("Class", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("Red.List.status", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("biome", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("title", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("units", "sigma", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("hab_spec", "Habitat specificity", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("log.meanpop", "log(Mean population size)", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("log(range)", "log(Geographic range)", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("&", "/", mcmc.outputs.all.global$`Variable`, fixed = TRUE)
mcmc.outputs.all.global$`Variable` <- gsub("Units2", "", mcmc.outputs.all.global$`Variable`, fixed = TRUE)

mcmc.outputs.all.global$`Effective sample size` <- round(mcmc.outputs.all.global$`Effective sample size`)
rownames(mcmc.outputs.all.global) <- c()
mcmc.outputs.all.global[429, 2] <- "Number of threats"
mcmc.outputs.all.global[443, 2] <- "Number of threats"

stargazer(mcmc.outputs.all.global, type = "html", rownames = FALSE,
          summary = FALSE, digits = 3)
 
 
# UK models ----
# ** mu ----
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

# Geographic range
rarity.range.uk <- MCMCglmm(mu ~ log(km2_range), data = uk.mus.range, random = ~species,
                         nitt = 200000, burnin = 20000) 
summary(rarity.range.uk)
# plot(rarity.range.uk$VCV)
# autocorr(rarity.range.uk$VCV)

# Population size
rarity.popsize.uk <- MCMCglmm(mu ~ log.meanpop, data = uk.mus.popsize, random = ~species,
                           nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk)
# plot(rarity.popsize.uk$VCV)
# autocorr(rarity.popsize.uk$VCV)

# Habitat specificity
rarity.hab.uk <- MCMCglmm(mu ~ hab_spec, data = uk.mus.habspec, random = ~species,
                       nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk)
# plot(rarity.hab.uk$VCV)
# autocorr(rarity.hab.uk$VCV)

# Habitat specificity - profiling method
rarity.hab.prof.uk <- MCMCglmm(mu ~ hab_spec, data = uk.mus.habspec.prof, random = ~species,
                            nitt = 200000, burnin = 20000) 
summary(rarity.hab.prof.uk)
# plot(rarity.hab.prof.uk$VCV)
# autocorr(rarity.hab.prof.uk$VCV)

# IUCN status
uk.mus.IUCN <- mus.IUCN %>% filter(grepl("United Kingdom", Country.list))

UK.IUCN <- MCMCglmm(mu ~ Red.List.status - 1, data = uk.mus.IUCN, random = ~species,
                    nitt = 200000, burnin = 20000) 
summary(UK.IUCN)
# plot(UK.IUCN$VCV)
# autocorr(UK.IUCN$VCV)

# Differences across systems in the UK
UK.system <- MCMCglmm(mu ~ system - 1, data = uk.mus, random = ~species,
                      nitt = 200000, burnin = 20000) 
summary(UK.system)
# plot(UK.system$VCV)
# autocorr(UK.system$VCV)

# Differences across taxa in the UK
uk.mus.class <- filter(uk.mus, Class == "Aves" | Class == "Mammalia" | 
                        Class == "Reptilia" | Class == "Amphibia" |
                        Class == "Actinopterygii" | Class == "Elasmobranchii")

UK.taxa <- MCMCglmm(mu ~ Class - 1, data = uk.mus.class, random = ~species,
                    nitt = 200000, burnin = 20000)
summary(UK.taxa)
# plot(UK.taxa$VCV)
# plot(UK.taxa$Sol)
# autocorr(UK.taxa$VCV)

# ** mu by tau ----
# weighted by observation error
# Geographic range
uk.weight1 <- uk.mus.range$tau.2

rarity.range.uk.var <- MCMCglmm(mu ~ log.range, mev = uk.weight1, data = uk.mus.range, random = ~species,
                             nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.var)
# plot(rarity.range.uk.var$VCV)
# autocorr(rarity.range.uk.var$VCV)

# Population size
uk.weight.pop <- uk.mus.popsize$tau.2

rarity.popsize.uk.var <- MCMCglmm(mu ~ log.meanpop, mev = uk.weight.pop, random = ~species,
                               data = uk.mus.popsize, nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.var)
# plot(rarity.popsize.uk.var$VCV)
# autocorr(rarity.popsize.uk.var$VCV)

# Habitat specificity
uk.weight2 <- uk.mus.habspec$tau.2

rarity.hab.uk.var <- MCMCglmm(mu ~ hab_spec, mev = uk.weight2, 
                           data = uk.mus.habspec, random = ~species,
                           nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.var)
# plot(rarity.hab.var.uk$VCV)
# autocorr(rarity.hab.var.uk$VCV)

# IUCN status
uk.weight3 <- uk.mus.IUCN$tau.2

UK.IUCN.var <- MCMCglmm(mu ~ Red.List.status - 1, mev = uk.weight3, 
                        data = uk.mus.IUCN, random = ~species,
                        nitt = 200000, burnin = 20000) 
summary(UK.IUCN.var)
# plot(UK.IUCN.var$VCV)
# autocorr(UK.IUCN.var$VCV)

# Differences across systems in the UK
uk.weight.system <- uk.mus$tau.2
UK.system.var <- MCMCglmm(mu ~ system - 1, mev = uk.weight.system,
                          data = uk.mus, random = ~species,
                          nitt = 200000, burnin = 20000) 
summary(UK.system.var)
# plot(UK.system.var$VCV)
# autocorr(UK.system.var$VCV)

# Differences across taxa in the UK
uk.weight4 <- uk.mus.class$tau.2

UK.taxa.var <- MCMCglmm(mu ~ Class - 1, mev = uk.weight4, 
                        data = uk.mus.class, random = ~species,
                        nitt = 200000, burnin = 20000)
summary(UK.taxa.var)
# plot(UK.taxa.var$VCV)
# autocorr(UK.taxa.var$VCV)

# ** slope ----
uk.slopes <- slopes %>% filter(grepl("United Kingdom", Country.list))
uk.slopes.ranges <- inner_join(uk.slopes, ranges, by = "species")
uk.slopes.ranges <- uk.slopes.ranges %>% 
  mutate(log.range = log(km2_range)) %>%
  filter(is.finite(log.range))

# Geographic range
rarity.range.uk2 <- MCMCglmm(slope ~ log.range, 
                          data = uk.slopes.ranges, random = ~species,
                          nitt = 200000, burnin = 20000) 
summary(rarity.range.uk2)
# plot(rarity.range.uk2$VCV)
# autocorr(rarity.range.uk2$VCV)

# Population size
uk.slopes.popsize <- inner_join(uk.slopes, slopes.popsize, by = "id")%>%
  dplyr::select(-slope.y)
colnames(uk.slopes.popsize)[14] <- "slope"
uk.slopes.popsize <- left_join(uk.slopes.popsize, id_species, by = "id")

rarity.popsize.uk2 <- MCMCglmm(slope ~ log(meanpop), data = uk.slopes.popsize, 
                            random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk2)
# plot(rarity.popsize.uk2$VCV)
# autocorr(rarity.popsize.uk2$VCV)

# Habitat specificity
uk.slopes.habspec <- inner_join(uk.slopes, habspec, by = "species") %>% 
  drop_na(hab_spec)
uk.slopes.habspec.prof <- inner_join(uk.slopes, habspec_prof, by = "species") %>% 
  drop_na(hab_spec)

rarity.hab.uk2 <- MCMCglmm(slope ~ hab_spec, data = uk.slopes.habspec, 
                           random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk2)
# plot(rarity.hab.uk2$VCV)
# autocorr(rarity.hab.uk2$VCV)

# Habitat specificity - profiling method
rarity.hab.uk.prof2 <- MCMCglmm(slope ~ hab_spec, data = uk.slopes.habspec.prof, 
                                random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.prof2)
# plot(rarity.hab.uk.prof2$VCV)
# autocorr(rarity.hab.uk.prof2$VCV)

# IUCN status
uk.slopes.IUCN <- slopes.IUCN %>% filter(grepl("United Kingdom", Country.list))

UK.IUCN2 <- MCMCglmm(slope ~ Red.List.status - 1, data = uk.slopes.IUCN, 
                     random = ~species, nitt = 200000, burnin = 20000)
summary(UK.IUCN2)
# plot(UK.IUCN2$VCV)
# autocorr(UK.IUCN2$VCV)

# Differences across systems in the UK
UK.system2 <- MCMCglmm(slope ~ system - 1, data = uk.slopes, 
                       random = ~species, nitt = 200000, burnin = 20000) 
summary(UK.system2)
# plot(UK.system2$VCV)
# autocorr(UK.system2$VCV)

# Differences across taxa in the UK
uk.slopes.class <- filter(uk.slopes, Class == "Aves" | Class == "Mammalia" | 
                            Class == "Reptilia" | Class == "Amphibia" |
                            Class == "Actinopterygii" | Class == "Elasmobranchii")

UK.class2 <- MCMCglmm(slope ~ Class - 1, data = uk.slopes.class, 
                      random = ~species, nitt = 200000, burnin = 20000)
summary(UK.class2)
# plot(UK.class2$VCV)
# autocorr(UK.class2$VCV)

# ** CI ----
# Variance models (confidence intervals around mu)
# Geographic range
uk.mus.range$CI <- (uk.mus.range$uCI - uk.mus.range$lCI)/2
uk.mus.range$CI <- scales::rescale(uk.mus.range$CI, to = c(0, 1))

rarity.range.uk.ci <- MCMCglmm(CI ~ log.range, data = uk.mus.range, random = ~species,
                            nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.ci)
# plot(rarity.range.uk.ci$VCV)
# autocorr(rarity.range.uk.ci$VCV)

# Population size
uk.mus.popsize$CI <- (uk.mus.popsize$uCI.x - uk.mus.popsize$lCI.x)/2
uk.mus.popsize$CI <- scales::rescale(uk.mus.popsize$CI, to = c(0, 1))

rarity.popsize.uk.ci <- MCMCglmm(CI ~ log.meanpop, data = uk.mus.popsize, random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.ci)
# plot(rarity.popsize.uk.ci$VCV)
# autocorr(rarity.popsize.uk.ci$VCV)

# Habitat specificity
uk.mus.habspec$CI <- (uk.mus.habspec$uCI - uk.mus.habspec$lCI)/2
uk.mus.habspec$CI <- scales::rescale(uk.mus.habspec$CI, to = c(0, 1))

rarity.hab.uk.ci <- MCMCglmm(CI ~ hab_spec, data = uk.mus.habspec, random = ~species,
                          nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.ci)
# plot(rarity.hab.uk.ci$VCV)
# autocorr(rarity.hab.uk.ci$VCV)

# IUCN status
uk.mus.IUCN$CI <- (uk.mus.IUCN$uCI - uk.mus.IUCN$lCI)/2
uk.mus.IUCN$CI <- scales::rescale(uk.mus.IUCN$CI, to = c(0, 1))

UK.IUCN.ci <- MCMCglmm(CI ~ Red.List.status - 1, data = uk.mus.IUCN, random = ~species,
                       nitt = 200000, burnin = 20000) 
summary(UK.IUCN.ci)
# plot(UK.IUCN$VCV)
# autocorr(UK.IUCN$VCV)

# Differences across systems in the UK
uk.mus$CI <- (uk.mus$uCI - uk.mus$lCI)/2
uk.mus$CI <- scales::rescale(uk.mus$CI, to = c(0, 1))

UK.system.ci <- MCMCglmm(CI ~ system - 1, data = uk.mus, random = ~species,
                         nitt = 200000, burnin = 20000) 
summary(UK.system.ci)
# plot(UK.system.ci$VCV)
# autocorr(UK.system.ci$VCV)

# Differences across taxa in the UK
UK.taxa.ci <- MCMCglmm(CI ~ Class - 1, data = uk.mus.class, random = ~species,
                       nitt = 200000, burnin = 20000)
summary(UK.taxa.ci)
# plot(UK.taxa.ci$VCV)
# plot(UK.taxa.ci$Sol)
# autocorr(UK.taxa.ci$VCV)

# ** CI by tau ----
# Weighted variance models 
# Geographic range
rarity.range.uk.ci.var <- MCMCglmm(CI ~ log.range, data = uk.mus.range, random = ~species,
                                mev = uk.weight1, nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.ci.var)
# plot(rarity.range.uk.ci.var$VCV)
# autocorr(rarity.range.uk.ci.var$VCV)

# Population size
rarity.popsize.uk.ci.var <- MCMCglmm(CI ~ log.meanpop, data = uk.mus.popsize, random = ~species,
                                  mev = uk.weight.pop, nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.ci.var)
# plot(rarity.popsize.uk.ci.var$VCV)
# autocorr(rarity.popsize.uk.ci.var$VCV)

# Habitat specificity
rarity.hab.uk.ci.var <- MCMCglmm(CI ~ hab_spec, data = uk.mus.habspec, mev = uk.weight2,
                                 random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.ci.var)
# plot(rarity.hab.uk.ci$VCV)
# autocorr(rarity.hab.uk.ci$VCV)

# IUCN status
UK.IUCN.ci.var <- MCMCglmm(CI ~ Red.List.status - 1, data = uk.mus.IUCN, random = ~species,
                           mev = uk.weight3, nitt = 200000, burnin = 20000) 
summary(UK.IUCN.ci.var)
# plot(UK.IUCN.ci.var$VCV)
# autocorr(UK.IUCN.ci.var$VCV)

# Differences across systems in the UK
UK.system.ci.var <- MCMCglmm(CI ~ system - 1, data = uk.mus, mev = uk.weight.system, 
                             random = ~species, nitt = 200000, burnin = 20000) 
summary(UK.system.ci.var)
# plot(UK.system.ci$VCV)
# autocorr(UK.system.ci$VCV)

# Differences across taxa in the UK
UK.taxa.ci.var <- MCMCglmm(CI ~ Class - 1, data = uk.mus.class, mev = uk.weight4,
                           random = ~species, nitt = 200000, burnin = 20000)
summary(UK.taxa.ci.var)
# plot(UK.taxa.ci.var$VCV)
# plot(UK.taxa.ci.var$Sol)
# autocorr(UK.taxa.ci.var$VCV)


# ** SE ----
# Using the se from the linear models

# Geographic range
uk.slopes.ranges$SE <- scales::rescale(uk.slopes.ranges$slope_se, to = c(0, 1))

rarity.range.uk.se <- MCMCglmm(SE ~ log.range, data = uk.slopes.ranges, random = ~species,
                            nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.se)
# plot(rarity.range.uk.se$VCV)
# autocorr(rarity.range.uk.se$VCV)

# Population size
uk.slopes.popsize$SE.UK <- scales::rescale(uk.slopes.popsize$slope_se.y, to = c(0, 1))

rarity.popsize.uk.se <- MCMCglmm(SE.UK ~ log(meanpop), data = uk.slopes.popsize, random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.se)
# plot(rarity.popsize.uk.se$VCV)
# autocorr(rarity.popsize.uk.se$VCV)

# Habitat specificity
uk.slopes.habspec$SE <- scales::rescale(uk.slopes.habspec$slope_se, to = c(0, 1))

rarity.hab.uk.se <- MCMCglmm(SE ~ hab_spec, data = uk.slopes.habspec, random = ~species,
                          nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.se)
# plot(rarity.hab.uk.se$VCV)
# autocorr(rarity.hab.uk.se$VCV)

# IUCN status
uk.slopes.IUCN$SE <- scales::rescale(uk.slopes.IUCN$slope_se, to = c(0, 1))

UK.IUCN.se <- MCMCglmm(SE ~ Red.List.status - 1, data = uk.slopes.IUCN, random = ~species,
                       nitt = 200000, burnin = 20000) 
summary(UK.IUCN.se)
# plot(UK.IUCN.se$VCV)
# autocorr(UK.IUCN.se$VCV)

# ** SD ----

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

rarity.range.uk.sd <- MCMCglmm(sdpop1 ~ log.range, data = UK.LPIraw.range, 
                               random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.sd)
# plot(rarity.range.rate.uk.sd$VCV)
# autocorr(rarity.range.rate.uk.sd$VCV)

# Population size
rarity.popsize.uk.sd <- MCMCglmm(sdpop1 ~ log.meanpop, data = UK.LPIraw.popsize, 
                                 random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.sd)
# plot(rarity.popsize.rate.uk.sd$VCV)
# autocorr(rarity.popsize.rate.uk.sd$VCV)

# Habitat specificity
rarity.hab.uk.sd <- MCMCglmm(sdpop1 ~ hab_spec, data = UK.LPIraw.habspec, 
                             random = ~species, nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.sd)
# plot(rarity.hab.uk.sd$VCV)
# autocorr(rarity.hab.uk.sd$VCV)

# IUCN status
UK.IUCN.sd <- MCMCglmm(sdpop1 ~ Red.List.status - 1, data = UK.LPIraw.IUCN, 
                       random = ~species, nitt = 200000, burnin = 20000) 
summary(UK.IUCN.sd)
# plot(UK.IUCN.sd$VCV)
# autocorr(UK.IUCN.sd$VCV)

# ** sigma ----
# Geographic range
uk.mus.range$sigma.2 <- scales::rescale(uk.mus.range$sigma.2, to = c(0, 1))

rarity.range.uk.sigma <- MCMCglmm(sigma.2 ~ log(km2_range), data = uk.mus.range, random = ~species,
                            nitt = 200000, burnin = 20000) 
summary(rarity.range.uk.sigma)
# plot(rarity.range.uk.sigma$VCV)
# autocorr(rarity.range.uk.sigma$VCV)

# Population size
uk.mus.popsize$sigma.2 <- scales::rescale(uk.mus.popsize$sigma.2.x, to = c(0, 1))

rarity.popsize.uk.sigma <- MCMCglmm(sigma.2 ~ log.meanpop, data = uk.mus.popsize, random = ~species,
                              nitt = 200000, burnin = 20000) 
summary(rarity.popsize.uk.sigma)
# plot(rarity.popsize.uk.sigma$VCV)
# autocorr(rarity.popsize.uk.sigma$VCV)

# Habitat specificity
uk.mus.habspec$sigma.2 <- scales::rescale(uk.mus.habspec$sigma.2, to = c(0, 1))

rarity.hab.uk.sigma <- MCMCglmm(sigma.2 ~ hab_spec, data = uk.mus.habspec, random = ~species,
                          nitt = 200000, burnin = 20000) 
summary(rarity.hab.uk.sigma)
# plot(rarity.hab.uk.sigma$VCV)
# autocorr(rarity.hab.uk.sigma$VCV)

# Habitat specificity - profiling method
uk.mus.habspec.prof$sigma.2 <- scales::rescale(uk.mus.habspec.profe$sigma.2, to = c(0, 1))

rarity.hab.prof.uk.sigma <- MCMCglmm(sigma.2 ~ hab_spec, data = uk.mus.habspec.prof, random = ~species,
                               nitt = 200000, burnin = 20000) 
summary(rarity.hab.prof.uk.sigma)
# plot(rarity.hab.prof.uk.sigma$VCV)
# autocorr(rarity.hab.prof.uk.sigma$VCV)

# IUCN status
uk.mus.IUCN$sigma.2 <- scales::rescale(uk.mus.IUCN$sigma.2, to = c(0, 1))

UK.IUCN.sigma <- MCMCglmm(sigma.2 ~ Red.List.status - 1, data = uk.mus.IUCN, random = ~species,
                    nitt = 200000, burnin = 20000) 
summary(UK.IUCN.sigma)
# plot(UK.IUCN.sigma$VCV)
# autocorr(UK.IUCN.sigma$VCV)

# Differences across systems in the UK
uk.mus$sigma.2 <- scales::rescale(uk.mus$sigma.2, to = c(0, 1))

UK.system.sigma <- MCMCglmm(sigma.2 ~ system - 1, data = uk.mus, random = ~species,
                      nitt = 200000, burnin = 20000) 
summary(UK.system.sigma)
# plot(UK.system.sigma$VCV)
# autocorr(UK.system.sigma$VCV)

# Differences across taxa in the UK
uk.mus.class$sigma.2 <- scales::rescale(uk.mus.class$sigma.2, to = c(0, 1))

UK.taxa.sigma <- MCMCglmm(sigma.2 ~ Class - 1, data = uk.mus.class, random = ~species,
                    nitt = 200000, burnin = 20000)
summary(UK.taxa.sigma)
# plot(UK.taxa.sigma$VCV)
# plot(UK.taxa.sigma$Sol)
# autocorr(UK.taxa.sigma$VCV)

# ** Summary output table ----
# Save model output objects
# Save all model outputs
model.list.uk <- list(UK.system, UK.system.var, UK.system2,
                      UK.system.sigma,
                      UK.system.ci, UK.system.ci.var,
                      UK.taxa, UK.taxa.var, UK.class2,
                      UK.taxa.sigma,
                      UK.taxa.ci, UK.taxa.ci.var,
                      rarity.range.uk, rarity.range.uk.var, rarity.range.uk2,
                      rarity.range.uk.sigma,
                      rarity.range.uk.ci, rarity.range.uk.ci.var, 
                      rarity.range.uk.se, rarity.range.uk.sd,
                      rarity.popsize.uk, rarity.popsize.uk.var, rarity.popsize.uk2,
                      rarity.popsize.uk.sigma,
                      rarity.popsize.uk.ci, rarity.popsize.uk.ci.var, 
                      rarity.popsize.uk.se, rarity.popsize.uk.sd,
                      rarity.hab.uk, rarity.hab.prof.uk, rarity.hab.uk.var, 
                      rarity.hab.uk2, rarity.hab.uk.prof2,
                      rarity.hab.uk.sigma, rarity.hab.prof.uk.sigma,
                      rarity.hab.uk.ci, rarity.hab.uk.ci.var,
                      rarity.hab.uk.se, rarity.hab.uk.sd,
                      UK.IUCN, UK.IUCN.var, UK.IUCN2,
                      UK.IUCN.sigma,
                      UK.IUCN.ci, UK.IUCN.ci.var,
                      UK.IUCN.se, UK.IUCN.sd)
 
names(model.list.uk) <- c("UK.system", "UK.system.var", "UK.system2",
                          "UK.system.sigma",
                          "UK.system.ci", "UK.system.ci.var",
                          "UK.taxa", "UK.taxa.var", "UK.class2",
                          "UK.taxa.sigma",
                          "UK.taxa.ci", "UK.taxa.ci.var",
                          "rarity.range.uk", "rarity.range.uk.var", "rarity.range.uk2",
                        "rarity.range.uk.sigma",
                        "rarity.range.uk.ci", "rarity.range.uk.ci.var", 
                        "rarity.range.uk.se", "rarity.range.uk.sd",
                        "rarity.popsize.uk", "rarity.popsize.uk.var", "rarity.popsize.uk2",
                        "rarity.popsize.uk.sigma",
                        "rarity.popsize.uk.ci", "rarity.popsize.uk.ci.var", 
                        "rarity.popsize.uk.se", "rarity.popsize.uk.sd",
                        "rarity.hab.uk", "rarity.hab.prof.uk", "rarity.hab.uk.var", 
                        "rarity.hab.uk2", "rarity.hab.uk.prof2",
                        "rarity.hab.uk.sigma", "rarity.hab.prof.uk.sigma",
                        "rarity.hab.uk.ci", "rarity.hab.uk.ci.var",
                        "rarity.hab.uk.se", "rarity.hab.uk.sd",
                        "UK.IUCN", "UK.IUCN.var", "UK.IUCN2",
                        "UK.IUCN.sigma",
                        "UK.IUCN.ci", "UK.IUCN.ci.var",
                        "UK.IUCN.se", "UK.IUCN.sd")

# save(model.list.uk, file = "data/output/models_uk18thFeb.RData")

# Create a list of input model names (global scale)
dataListNames.uk <- list("Realm - mu", 
                      "Realm - weighted", 
                      "Realm - slope", 
                      "Realm - fluctuations sigma", 
                      "Realm - fluctuations CI",                       
                      "Realm - fluctuations CI weighted", 
                      "Taxa - mu", 
                      "Taxa - weighted", 
                      "Taxa - slope", 
                      "Taxa - fluctuations sigma", 
                      "Taxa - fluctuations CI", 
                      "Taxa - fluctuations CI weighted",
                      "Geographic range (all) - mu",
                      "Geographic range (all) - weighted",
                      "Geographic range (all) - slope",
                      "Geographic range (all) - fluctuations sigma",
                      "Geographic range (all) - fluctuations CI", 
                      "Geographic range (all) - fluctuations CI weighted",
                      "Geographic range (all)  - fluctuations SE",
                      "Geographic range (all) - fluctuations SD",
                      "Mean population size - mu", 
                      "Mean population size - weighted", 
                      "Mean population size - slope", 
                      "Mean population size - fluctuations sigma",
                      "Mean population size - fluctuations CI", 
                      "Mean population size - fluctuations CI weighted",
                      "Mean population size - fluctuations SE",
                      "Mean population size - fluctuations SD",
                      "Habitat specificity - mu", 
                      "Habitat specificity (profiling) - mu",
                      "Habitat specificity - weighted", 
                      "Habitat specificity - slope",  
                      "Habitat specificity (profiling) - slope", 
                      "Habitat specificity - fluctuations sigma",
                      "Habitat specificity (profiling) - fluctuations sigma",
                      "Habitat specificity - fluctuations CI", 
                      "Habitat specificity - fluctuations CI weighted",
                      "Habitat specificity - fluctuations SE", 
                      "Habitat specificity - fluctuations SD",
                      "IUCN conservation status - mu", 
                      "IUCN conservation status - weighted", 
                      "IUCN conservation status - slope", 
                      "IUCN conservation status - fluctuations sigma",
                      "IUCN conservation status - fluctuations CI", 
                      "IUCN conservation status - fluctuations CI weighted", 
                      "IUCN conservation status - fluctuations SE", 
                      "IUCN conservation status - fluctuations SD")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList.uk <- mapply(cbind, lapply(model.list.uk, clean.MCMC), "modelName" = dataListNames.uk, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.all.uk <- as.data.frame(do.call(rbind, readyList.uk))

mcmc.outputs.all.uk <- mcmc.outputs.all.uk %>% dplyr::select(modelName, variable:effect)

colnames(mcmc.outputs.all.uk) <- c("Model name", "Variable", "Posterior mean",
                                       "Lower 95% CI", "Upper 95% CI", "Effective sample size",
                                       "pMCMC", "Effect")

write.csv(mcmc.outputs.all.uk, file = "data/output/uk_outputs18thFeb.csv")

mcmc.outputs.all.uk$`Model name`[duplicated(mcmc.outputs.all.uk$`Model name`)] <- " "

mcmc.outputs.all.uk$`Variable` <- gsub("system", "", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("Class", "", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("Red.List.status", "", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("biome", "", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("title", "", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("units", "sigma", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("hab_spec", "Habitat specificity", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("log.meanpop", "log(Mean population size)", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("log(range)", "log(Geographic range)", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)
mcmc.outputs.all.uk$`Variable` <- gsub("&", "/", mcmc.outputs.all.uk$`Variable`, fixed = TRUE)

mcmc.outputs.all.uk$`Effective sample size` <- round(mcmc.outputs.all.uk$`Effective sample size`)
rownames(mcmc.outputs.all.uk) <- c()

stargazer(mcmc.outputs.all.uk, type = "html", rownames = FALSE,
          summary = FALSE, digits = 3)

# Phylogeny models ----
# Note that here we include phylogeny models using one tree
# The models were reran 10 times using different trees
# The manuscript presents the results of the mean across the ten
# different models to account for phylogenetic uncertainty

# ** Birds ----
birds_mu <- mus %>% filter(Class == "Aves")
bird_mu_ranges <- inner_join(birds_mu, bird_ranges, by = "species")
bird_mu_ranges <- bird_mu_ranges %>% filter(range != "EXTINCT") %>%
  drop_na(range)
bird_mu_ranges$range <- as.numeric(as.character(bird_mu_ranges$range))
bird_mu_ranges <- bird_mu_ranges %>% 
  drop_na(range)

# Adding bird phylogeny
birdtree <- read.nexus("output.nex")

# Then convert to the type of tree MCMCglmm can use with this code using ape package:
inv.phylo <- inverseA(birdtree[[8]], nodes = "TIPS", scale = TRUE)

phylo.names <- inv.phylo$node.names
phylo.names <- gsub("_", " ", phylo.names, fixed = TRUE)
names <- as.data.frame(phylo.names)
colnames(names) <- "species"
mu_names <- as.data.frame(unique(bird_mu_ranges$species))
colnames(mu_names) <- "species"
test <- anti_join(mu_names, names, by = "species")

# Matching the tree labels to the species in the pop change database
bird_mu_ranges2 <- bird_mu_ranges %>% filter(!species %in% c("Thalasseus sandvicensis",
                                                             "Ardea alba",
                                                             "Dryobates minor",
                                                             "Antigone canadensis",
                                                             "Vermivora cyanoptera",
                                                             "Mareca strepera",
                                                             "Spatula clypeata",
                                                             "Falcipennis canadensis",
                                                             "Hydroprogne caspia",
                                                             "Tringa semipalmata",
                                                             "Larus smithsonianus",
                                                             "Hylatomus pileatus",
                                                             "Antrostomus vociferus",
                                                             "Spatula discors",
                                                             "Circus hudsonius",
                                                             "Spatula cyanoptera",
                                                             "Mareca americana",
                                                             "Mareca penelope",
                                                             "Leuconotopicus villosus",
                                                             "Dryobates pubescens",
                                                             "Gallinula galeata",
                                                             "Selasphorus calliope",
                                                             "Parkesia motacilla",
                                                             "Leuconotopicus albolarvatus",
                                                             "Leuconotopicus borealis",
                                                             "Melozone aberti",
                                                             "Parkesia noveboracensis",
                                                             "Dryobates scalaris",
                                                             "Dryobates nuttallii",
                                                             "Gelochelidon nilotica",
                                                             "Antrostomus carolinensis",
                                                             "Gallinago delicata",
                                                             "Melanitta deglandi",
                                                             "Thalasseus maximus",
                                                             "Melozone crissalis",
                                                             "Rhynchophanes mccownii",
                                                             "Sternula antillarum",
                                                             "Lyrurus tetrix",
                                                             "Sternula albifrons",
                                                             "Microcarbo pygmaeus",
                                                             "Anser caerulescens",
                                                             "Peucaea aestivalis",
                                                             "Peucaea cassinii",
                                                             "Spatula querquedula",
                                                             "Tringa brevipes",
                                                             "Antigone vipio",
                                                             "Microcarbo melanoleucos",
                                                             "Aerodramus maximus",
                                                             "Calidris pugnax",
                                                             "Leiopicus medius",
                                                             "Leptoptilos crumenifer",
                                                             "Alopecoenas xanthonurus",
                                                             "Thalassarche melanophris",
                                                             "Anser canagicus",
                                                             "Thinornis cucullatus",
                                                             "Onychoprion lunatus",
                                                             "Aquila fasciata",
                                                             "Calidris subruficollis",
                                                             "Ardea plumifera",
                                                             "Hydrobates furcatus",
                                                             "Thalasseus bergii",
                                                             "Leucogeranus leucogeranus",
                                                             "Ardenna grisea",
                                                             "Threskiornis moluccus",
                                                             "Bubo scandiacus",
                                                             "Onychoprion fuscatus",
                                                             "Zapornia parva",
                                                             "Saundersilarus saundersi",
                                                             "Anser rossii",
                                                             "Eolophus roseicapilla",
                                                             "Sternula superciliaris",
                                                             "Calidris virgata",
                                                             "Melanitta americana",
                                                             "Calidris pygmaea",
                                                             "Caprimulgus longipennis"))



#View one of the 100 phylogenetic trees obtained from Jetz.
plot(birdtree[[1]], cex = 0.1)
plot(birdtree[[75]], cex = 0.1)

#To test phylogenetic signal

##set priors
a <- 1000
prior_pa <- list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), 
                                nu=1, alpha.mu=0, 
                                alpha.V=diag(1)*a), 
                        G1=list(V=diag(1), nu=1, 
                                alpha.mu=0, 
                                alpha.V=diag(1)*a)))

#Model used to calculate phylogenetic signal
#phylo: to detect whether traits are similar enough that they follow Brownian Motion model
#species: variance explained by species-specific effects. Does each species have a different slope?

bird_mu_ranges2$phylo <- gsub("[[:blank:]]", "_", bird_mu_ranges2$species)
bird_mu_ranges2$phylo <- as.factor(bird_mu_ranges2$phylo)

rarity.range.birds.phyl <- MCMCglmm(mu ~ 1, random = ~phylo + species,
                                    data = bird_mu_ranges2, 
                                    prior = prior_pa,
                                    ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                                    nitt = 120000, burnin = 20000, thin = 10) 

summary(rarity.range.birds.phyl)
save(rarity.range.birds.phyl, file = "phylo_model.RData")
plot(rarity.range.birds.phyl$VCV)

# Fluctuations and phylogeny for birds
rept_sigmas$sigma.2s <- scales::rescale(rept_sigmas$sigma.2, to = c(0, 1))

bird.phyl.sigma1 <- MCMCglmm(sigma.2s ~ 1, random = ~phylo + species,
                             data = bird_mu_ranges2, 
                             prior = prior_pa,
                             ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                             nitt = 120000, burnin = 20000, thin = 10) 

summary(bird.phyl.sigma1)
plot(bird.phyl.sigma1$VCV)
save(bird.phyl.sigma1, file = "data/output/r_phylo_sigma1.RData")

# ** Reptiles ----
# Global amphibians
rept_sigmas <- mus %>% filter(Class == "Reptilia")

# Adding amphibian phylogeny
reptree <- read.nexus("data/output_rept.nex")

if(any(is.ultrametric(reptree)) == FALSE) {
  reptree <- lapply(reptree, chronoMPL)
  class(reptree) <- "multiPhylo"
}

# Then convert to the type of tree MCMCglmm can use with this code using ape package:
inv.phylo <- inverseA(reptree[[1]], nodes = "TIPS", scale = TRUE)

phylo.names <- inv.phylo$node.names
phylo.names <- gsub("_", " ", phylo.names, fixed = TRUE)
names <- as.data.frame(phylo.names)
colnames(names) <- "species"
mu_names <- as.data.frame(unique(rept_sigmas$species))
colnames(mu_names) <- "species"
test <- anti_join(mu_names, names, by = "species")

rept_sigmas <- rept_sigmas %>% filter(!species %in% c("Natator depressus",
                                                      "Eretmochelys imbricata",
                                                      "Pseudemydura umbrina",
                                                      "Dermochelys coriacea",
                                                      "Caretta caretta",
                                                      "Lepidochelys olivacea",
                                                      "Chelonia mydas",
                                                      "Lepidochelys kempii",
                                                      "Alligator mississippiensis",
                                                      "Crocodylus porosus",
                                                      "Glyptemys insculpta",
                                                      "Chelydra serpentina",
                                                      "Terrapene ornata",
                                                      "Crocodylus acutus",
                                                      "Crocodylus novaeguineae",
                                                      "Alligator sinensis",
                                                      "Testudo hermanni",
                                                      "Testudo graeca",
                                                      "Graptemys geographica",
                                                      "Malaclemys terrapin",
                                                      "Emydoidea blandingii",
                                                      "Batagur baska",
                                                      "Crocodylus niloticus",
                                                      "Gopherus agassizii",
                                                      "Terrapene carolina",
                                                      "Chrysemys picta",
                                                      "Clemmys guttata",
                                                      "Pseudemys peninsularis",
                                                      "Pseudemys nelsoni",
                                                      "Sternotherus minor",
                                                      "Sternotherus odoratus"))

#To test phylogenetic signal

##set priors
a <- 1000
prior_pa <- list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), 
                                nu=1, alpha.mu=0, 
                                alpha.V=diag(1)*a), 
                        G1=list(V=diag(1), nu=1, 
                                alpha.mu=0, 
                                alpha.V=diag(1)*a)))

#Model used to calculate phylogenetic signal
#phylo: to detect whether traits are similar enough that they follow Brownian Motion model
#species: variance explained by species-specific effects. Does each species have a different slope?

rept_sigmas$phylo <- gsub("[[:blank:]]", "_", rept_sigmas$species)
rept_sigmas$phylo <- as.factor(rept_sigmas$phylo)

# Population trends and phylogeny for reptiles
rept.phyl.mu1 <- MCMCglmm(mu ~ 1, random = ~phylo + species,
                          data = rept_sigmas, 
                          prior = prior_pa,
                          ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                          nitt = 120000, burnin = 20000, thin = 10) 

summary(rept.phyl.mu1)
plot(rept.phyl.mu1$VCV)
save(rept.phyl.mu1, file = "data/output/r_phylo_mu1.RData")

# Fluctuations and phylogeny for reptiles
rept_sigmas$sigma.2s <- scales::rescale(rept_sigmas$sigma.2, to = c(0, 1))

rept.phyl.sigma1 <- MCMCglmm(sigma.2s ~ 1, random = ~phylo + species,
                             data = rept_sigmas, 
                             prior = prior_pa,
                             ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                             nitt = 120000, burnin = 20000, thin = 10) 

summary(rept.phyl.sigma1)
plot(rept.phyl.sigma1$VCV)
save(rept.phyl.sigma1, file = "data/output/r_phylo_sigma1.RData")

# ** Amphibians ----
# Global amphibians
amp_mus <- mus %>% filter(Class == "Amphibia")

# Adding amphibian phylogeny
amptree <- read.nexus("data/output_amp.nex")

if(any(is.ultrametric(amptree)) == FALSE) {
  amptree <- lapply(amptree, chronoMPL)
  class(amptree) <- "multiPhylo"
}

# Then convert to the type of tree MCMCglmm can use with this code using ape package:
inv.phylo <- inverseA(amptree[[1]], nodes = "TIPS", scale = TRUE)

phylo.names <- inv.phylo$node.names
phylo.names <- gsub("_", " ", phylo.names, fixed = TRUE)
names <- as.data.frame(phylo.names)
colnames(names) <- "species"
mu_names <- as.data.frame(unique(amp_mus$species))
colnames(mu_names) <- "species"
test <- anti_join(mu_names, names, by = "species")

amp_mus <- amp_mus %>% filter(!species %in% c("Lithobates sylvaticus",
                                              "Lithobates clamitans",
                                              "Lithobates sphenocephalus",
                                              "Lithobates catesbeianus",
                                              "Lithobates septentrionalis",
                                              "Lithobates pipiens",
                                              "Lithobates palustris",
                                              "Lithobates yavapaiensis",
                                              "Lithobates berlandieri"))

#To test phylogenetic signal

##set priors
a <- 1000
prior_pa <- list(R=list(V=diag(1), nu=0.002), 
                 G=list(G1=list(V=diag(1), 
                                nu=1, alpha.mu=0, 
                                alpha.V=diag(1)*a), 
                        G1=list(V=diag(1), nu=1, 
                                alpha.mu=0, 
                                alpha.V=diag(1)*a)))

#Model used to calculate phylogenetic signal
#phylo: to detect whether traits are similar enough that they follow Brownian Motion model
#species: variance explained by species-specific effects. Does each species have a different slope?

amp_mus$phylo <- gsub("[[:blank:]]", "_", amp_mus$species)
amp_mus$phylo <- as.factor(amp_mus$phylo)

# Trends and amphibian phylogeny
amp.phyl1 <- MCMCglmm(mu ~ 1, random = ~phylo + species,
                      data = amp_mus, 
                      prior = prior_pa,
                      ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                      nitt = 120000, burnin = 20000, thin = 10) 

summary(amp.phyl1)
plot(amp.phyl1$VCV)
save(amp.phyl1, file = "data/output/a_phylo_model1.RData")

# Fluctuations and phylogeny for amphibians
amp_mus$sigma.2s <- scales::rescale(amp_mus$sigma.2, to = c(0, 1))

amp.phyl.sigma1 <- MCMCglmm(sigma.2s ~ 1, random = ~phylo + species,
                            data = amp_mus, 
                            prior = prior_pa,
                            ginverse = list(phylo = inv.phylo$Ainv), pr = TRUE,
                            nitt = 120000, burnin = 20000, thin = 10) 

summary(amp.phyl.sigma1)
plot(amp.phyl.sigma1)

save(amp.phyl.sigma1, file = "data/output/a_phylo_sigma1.RData")

# UK species names table ----
species.names <- uk.mus %>% dplyr::select(species, id) %>% 
  group_by(species) %>% tally
write.csv(species.names, file = "data/species_names.csv")