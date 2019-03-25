# Using `rredlist` to extract habitat specificity data from IUCN Red List
# John Godlee (johngodlee@gmail.com)
# 2018_03_14

# The habitat data is contained in:
# Full Account -> Classification -> Habitat, on the IUCN website
# This seems like a better system to classify the species by, 
# it has distinct categories which are consistent between species, but it may have worse coverage for some species than the "Habitat and Ecology" text we used last time
# http://www.iucnredlist.org/technical-documents/classification-schemes/habitats-classification-scheme-ver3

# Packages ----
library(rredlist)
library(dplyr)

# Load data ----

# IUCN Redlist API key
IUCN_REDLIST_KEY <- "5095ce9b03e5bae012cde5cc112dd0baf09092fa4a723d042142a8c21ea3ad00"

# Get species list
mus_names <- read.csv("data/global_mus_scaled.csv") %>% 
  mutate(species = paste(Genus, Species, sep = ' ')) %>%
  dplyr::select(species)


# Get habitat data for each species ----

# Convert iucn_sp_names to list, 1 name per list entry
iucn_sp_list <- as.list(as.character(mus_names$species))

# Download species habitats to list
habitat_sp <- lapply(iucn_sp_list, function(i) {rl_habitats(name = i, key = IUCN_REDLIST_KEY)})

save(habitat_sp, file = "global_sp.RData")

# Transform list to tidy data frame with list of habitats per species ----

# Add column to each `result` data frame of species name
habitat_name <- lapply(habitat_sp, function(i) {
  i$result$species = rep(i$name, times = length(i$result$code)); 
  return(i)
})

# Collapse list to data frame
habitat_df <- lapply(habitat_name, "[[", 2) %>%
  bind_rows()

# Get unique habitats for each species
habitat_df_unique <- habitat_df %>%
  group_by(species, habitat) %>%
  tally() %>%
  dplyr::select(species, habitat)

# Create summary data frame of number of habitats per species
habitat_df_count <- habitat_df_unique %>%
  group_by(species) %>%
  tally()


# Write list of habitats data frame to csv ----
write.csv(habitat_df_unique, "data/iucn_sp_habitats_global.csv")
write.csv(habitat_df_count, "data/iucn_sp_habitats_count_global.csv")