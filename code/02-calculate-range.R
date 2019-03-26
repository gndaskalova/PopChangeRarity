# Calculating sensible (excluding obvious outliers) species ranges from GBIF occurrence data
# John Godlee (johngodlee@gmail.com)
# 2018_03_26

# Prep ----
# Packages
library(data.table)
library(dplyr)
library(CoordinateCleaner)
library(geosphere)
library(ggplot2)
library(gridExtra)
library(parallel)
library(tidyr)
library(scrubr)

# Set working directory to the location of the source file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in GBIF data ----
# gbif_1 <- fread(file = 'GD_occurrence.csv', 
#                 skip = 1, sep = "\t", 
#                 na.strings = c("", "NA", NA, "null", "\\N"))
# save(gbif_1, file = "gbif_1.RData")

load("D:/Gergana_Daskalova/gbif_1.RData") # Note this file is very large and it cannot
# be stored on GitHub

# The file represents all GBIF occurrences for the species which have monitored populations
# in the UK

# Clean data ----
# Clean it
gbif_colnames <- gbif_1 %>%
  select(V1, V4, V5, V6, V7, V8, V9, V10, V13, V16, V17, V22) %>%
  rename(gbif_id = V1,
         kingdom = V4,
         phylum = V5,
         order = V6,
         class = V7,
         family = V8,
         genus = V9,
         species_binomial = V10,
         species_binomial_auth = V13,
         decimallatitude = V16, 
         decimallongitude = V17,
         year = V22) %>%
  filter(!is.na(decimallongitude),
         !is.na(decimallatitude))  # Remove rows with NA in lat lon, breaks CoordinateCleaner


# Identify outliers - batch processed because otherwise crashes supershrub :(

gbif_list <- split(gbif_colnames, factor(gbif_colnames$species_binomial))

gbif_list_spl1 <- gbif_list[1:50]
gbif_spl1 <- do.call("rbind", gbif_list_spl1)

gbif_list_spl2 <- gbif_list[51:100]
gbif_spl2 <- do.call("rbind", gbif_list_spl2)

gbif_list_spl3 <- gbif_list[101:150]
gbif_spl3 <- do.call("rbind", gbif_list_spl3)

gbif_list_spl4 <- gbif_list[151:200]
gbif_spl4 <- do.call("rbind", gbif_list_spl4)

gbif_list_spl5 <- gbif_list[201:231]
gbif_spl5 <- do.call("rbind", gbif_list_spl5)


gbif_clean_spl1 <- gbif_spl1 %>%
  coord_impossible(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_incomplete(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_unlikely(lat = "decimallatitude", lon = "decimallongitude") %>%
  mutate(flag = CleanCoordinates(., 
                   value = "flags", 
                   lon = "decimallongitude",
                   lat = "decimallatitude",
                   species = "species_binomial",
                   outliers = FALSE,
                   zeros = TRUE,
                   GBIF = TRUE,
                   institutions = TRUE,
                   capitals = FALSE,
                   urban = FALSE,
                   duplicates = FALSE,
                   seas = FALSE)) %>%
  filter(flag == "FALSE") %>%
  filter_("decimallatitude!=decimallongitude") %>%  # Where lat==long
  filter(round(decimallongitude) != decimallongitude,  # Where there are no decimal places
         round(decimallatitude) != decimallatitude) %>%
  group_by(species_binomial) %>%
  filter(n() > 20)
  

gbif_clean_spl2 <- gbif_spl2 %>%
  coord_impossible(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_incomplete(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_unlikely(lat = "decimallatitude", lon = "decimallongitude") %>%
  mutate(flag = CleanCoordinates(., 
                                 value = "flags", 
                                 lon = "decimallongitude",
                                 lat = "decimallatitude",
                                 species = "species_binomial",
                                 outliers = FALSE,
                                 zeros = TRUE,
                                 GBIF = TRUE,
                                 institutions = TRUE,
                                 capitals = FALSE,
                                 urban = FALSE,
                                 duplicates = FALSE,
                                 seas = FALSE)) %>%
  filter(flag == "FALSE") %>%
  filter_("decimallatitude!=decimallongitude") %>%  # Where lat==long
  filter(round(decimallongitude) != decimallongitude,  # Where there are no decimal places
         round(decimallatitude) != decimallatitude) %>%
  group_by(species_binomial) %>%
  filter(n() > 20)

gbif_clean_spl3 <- gbif_spl3 %>%
  coord_impossible(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_incomplete(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_unlikely(lat = "decimallatitude", lon = "decimallongitude") %>%
  mutate(flag = CleanCoordinates(., 
                                 value = "flags", 
                                 lon = "decimallongitude",
                                 lat = "decimallatitude",
                                 species = "species_binomial",
                                 outliers = FALSE,
                                 zeros = TRUE,
                                 GBIF = TRUE,
                                 institutions = TRUE,
                                 capitals = FALSE,
                                 urban = FALSE,
                                 duplicates = FALSE,
                                 seas = FALSE)) %>%
  filter(flag == "FALSE") %>%
  filter_("decimallatitude!=decimallongitude") %>%  # Where lat==long
  filter(round(decimallongitude) != decimallongitude,  # Where there are no decimal places
         round(decimallatitude) != decimallatitude) %>%
  group_by(species_binomial) %>%
  filter(n() > 20)

gbif_clean_spl4 <- gbif_spl4 %>%
  coord_impossible(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_incomplete(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_unlikely(lat = "decimallatitude", lon = "decimallongitude") %>%
  mutate(flag = CleanCoordinates(., 
                                 value = "flags", 
                                 lon = "decimallongitude",
                                 lat = "decimallatitude",
                                 species = "species_binomial",
                                 outliers = FALSE,
                                 zeros = TRUE,
                                 GBIF = TRUE,
                                 institutions = TRUE,
                                 capitals = FALSE,
                                 urban = FALSE,
                                 duplicates = FALSE,
                                 seas = FALSE)) %>%
  filter(flag == "FALSE") %>%
  filter_("decimallatitude!=decimallongitude") %>%  # Where lat==long
  filter(round(decimallongitude) != decimallongitude,  # Where there are no decimal places
         round(decimallatitude) != decimallatitude) %>%
  group_by(species_binomial) %>%
  filter(n() > 20)

gbif_clean_spl5 <- gbif_spl5 %>%
  coord_impossible(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_incomplete(lat = "decimallatitude", lon = "decimallongitude") %>%
  coord_unlikely(lat = "decimallatitude", lon = "decimallongitude") %>%
  mutate(flag = CleanCoordinates(., 
                                 value = "flags", 
                                 lon = "decimallongitude",
                                 lat = "decimallatitude",
                                 species = "species_binomial",
                                 outliers = FALSE,
                                 zeros = TRUE,
                                 GBIF = TRUE,
                                 institutions = TRUE,
                                 capitals = FALSE,
                                 urban = FALSE,
                                 duplicates = FALSE,
                                 seas = FALSE)) %>%
  filter(flag == "FALSE") %>%
  filter_("decimallatitude!=decimallongitude") %>%  # Where lat==long
  filter(round(decimallongitude) != decimallongitude,  # Where there are no decimal places
         round(decimallatitude) != decimallatitude) %>%
  group_by(species_binomial) %>%
  filter(n() > 20)


gbif_clean <- rbind(gbif_clean_spl1, gbif_clean_spl2, gbif_clean_spl3, gbif_clean_spl4, gbif_clean_spl5)

# Rename object for range analysis
GBIF <- gbif_clean

# save(GBIF, file = "GBIF.Rdata")

# Ensure lat and long are encoded as numbers
GBIF$decimallatitude <- as.numeric(GBIF$decimallatitude)
GBIF$decimallongitude <- as.numeric(GBIF$decimallongitude)

# Remove top 90% and lower 10 quantile of lat and long for each species ----
GBIF_quant <- GBIF %>% group_by(species_binomial) %>%
  filter(decimallatitude < quantile(GBIF$decimallatitude, 0.995) &
           decimallatitude > quantile(GBIF$decimallatitude, 0.005)) %>%
  filter(decimallongitude < quantile(GBIF$decimallongitude, 0.995) &
           decimallongitude > quantile(GBIF$decimallongitude, 0.005))

# save(GBIF_quant, file = "GBIF_quant.Rdata")

# Convex hull range extent ----
# This method performs better with species whose ranges span the international date line or the poles
# The previous method vastly over-estimated these species ranges by wrapping them round
# Create list of data per species
GBIF <- GBIF_quant
GBIF_list <- split(GBIF, GBIF$species_binomial)
#save(GBIF_list, file = "GBIF_list.Rdata")

# Identify rows which contain occurrences bounding the convex hull of each species 
GBIF_chull_list <- lapply(GBIF_list, function(x) chull(x$decimallongitude, x$decimallatitude))

# Subset GBIF_list to only include those rows identified as bounding the convex hull
GBIF_chull_subset <- Map(function(x, y) x[y, ], GBIF_list, GBIF_chull_list)

# Calculate the polygon area of each convex hull
chull_km2_range <- 0.000001*unlist(lapply(GBIF_chull_subset, function(x) areaPolygon(x[,11:10])))

# Make summary data frame
GBIF_ranges <- data.frame(names(chull_km2_range))
GBIF_ranges$chull_km2_range <- chull_km2_range
colnames(GBIF_ranges) <- c("species", "range")

# Comparison with ranges before quantile filtering ----
colnames(gbif_ranges_before_quant)[2] <- "species"
range.comparison <- inner_join(GBIF_ranges, 
                               gbif_ranges_before_quant, 
                               by = "species")

(range.comp <- ggplot(range.comparison, aes(x = range, 
                                           y = chull_km2_range)) +
  geom_point() +
  geom_smooth(method = "lm") +
  coord_equal() + 
  theme_bw())

# Save range estimates ----
write.csv(GBIF_ranges, file = "data/input/gbif_ranges_clean.csv")