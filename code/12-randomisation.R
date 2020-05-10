# Randomisation sensitivity analysis

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

df1 %>% 
  group_by(class) %>% 
  mutate(value = value[sample(row_number())]) 


# Packages ----
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(broom)
library(viridis)

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

LPI.long.sorted <- arrange(LPI.long, id)

test3 <- LPI.long.sorted %>% 
  group_by(id) %>% 
  mutate(scalepop = scalepop[sample(row_number())]) 


# ** Models ----
# Run linear models of abundance trends over time for each population and extract model coefficients
LPI.models.random <- test3 %>%
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

hist(LPI.models.random$slope)

real_slopes <- read.csv("data/input/global_slopes.csv")
hist(real_slopes$slope)

LPI.models.random$category <- "Randomised"
real_slopes$category <- "Real data"

slope_comparison <- bind_rows(LPI.models.random, real_slopes)

(comparison <- ggplot(slope_comparison, aes(x = slope, fill = category)) +
  geom_density(alpha = 0.6) +
    theme_LPI() +
    labs(x = "\nPopulation trend (slope)", y = "Density\n") +
    scale_fill_viridis_d(option = "magma", end = 0.7))

ggsave(comparison, filename = "figures/randomisation.pdf", height = 10, width = 10)
