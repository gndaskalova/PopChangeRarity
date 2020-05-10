# looking into possible bands of mus of the same value

library(tidyverse)

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
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
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

# Load data ----
LPI <- read.csv("data/input/LPIdata_Feb2016.csv")  # raw LPI data
mus <- read.csv("data/input/global_mus_scaled.csv")  # overall population change using mu
slopes <- read.csv("data/input/global_slopes.csv")  # overall pop change using slopes

# Add duration to the data frame from the state-space models
durations <- slopes %>% dplyr::select(id, lengthyear)
colnames(durations) <- c("id", "duration")
mus <- left_join(mus, durations, by = "id")

# Extract the populations with the same lat and long
mus$coord <- paste(mus$Decimal.Latitude, mus$Decimal.Longitude, sep = " ")
length(mus$coord)  # 9286
length(unique(mus$coord))  # 3033
3033/9286
# 6253 geographic locations have more than one population

# Look if the mu values are exactly the same at some of the bands
mus_020 <- mus %>% filter(mu > 0.199999999999 & mu < 0.2000001)

# extract raw data for those 13 populations

# Transform from wide to long format
LPI.long <- gather(data = LPI, key = "year", value = "pop", 26:70)
# Get rid of the X in front of years
LPI.long$year <- parse_number(LPI.long$year)

# Create new column with genus and species together
LPI.long$species <- paste(LPI.long$Genus, LPI.long$Species)

# Remove duplicate rows
LPI.long <- LPI.long %>% distinct() %>% filter(is.finite(pop))

# Get just the populatations with count data
LPI.long.pop <- LPI.long
LPI.long <- LPI.long %>%
  drop_na(pop) %>%
  group_by(id) %>%   # group rows so that each group is one population
  mutate(scalepop = scales::rescale(pop, to = c(-1, 1))) %>%
  filter(length(unique(year)) > 5) %>%
  drop_na(scalepop)

pop_ids <- mus_020 %>% dplyr::select(id)
pop_ids <- left_join(pop_ids, LPI.long, by = "id")
pop_ids$species_id <- paste(pop_ids$species, pop_ids$id, sep = "_")

(pop_trends <- ggplot(pop_ids, aes(x = year, y = pop)) +
  geom_point() +
  facet_wrap(~species_id, scales = "free") +
    theme_LPI() +
    labs(x = "Year", y = "Population abundance"))

ggsave(pop_trends, filename = "figures/pop_trends_020.png", height = 12, width = 15)

mus_0025 <- mus %>% filter(mu > 0.02499999 & mu < 0.02500000001)

pop_ids2 <- mus_0025 %>% dplyr::select(id)
pop_ids2 <- left_join(pop_ids2, LPI.long, by = "id")
pop_ids2$species_id <- paste(pop_ids2$species, pop_ids2$id, sep = "_")

# Get references for studies 468, 10193, 17803
refs <- LPI %>% filter(id %in% c("468", "10193", "17803"))
refs <- distinct(refs)
refs <- refs %>% dplyr::select(id, Data.source.citation)
write.csv(refs, file = "data/output/refs_tableS8.csv")

# Get refs for studies with low variance
refs2 <- slopes %>% filter(slope_se < 0.001) %>% dplyr::select(id)
refs2 <- left_join(refs2, LPI, by = "id") %>% dplyr::select(id, Data.source.citation) %>% distinct()

write.csv(refs2, file = "data/output/refs_tableS7.csv")

(pop_trends2 <- ggplot(pop_ids2, aes(x = year, y = pop)) +
    geom_point() +
    facet_wrap(~species_id, scales = "free") +
    theme_LPI() +
    labs(x = "Year", y = "Population abundance"))

ggsave(pop_trends2, filename = "figures/pop_trends_0025.png", height = 7, width = 12)

# Remaking a violin figure using slopes
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

(rain_system2 <- 
    ggplot(data = slopes, 
           aes(x = system, y = slope, fill = system)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y = slope, color = system), 
               position = position_jitter(width = .15), size = .5, alpha = 0.1) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8) +
    labs(y = "Slope\n", x = "\nSystem") +
    guides(fill = FALSE, color = FALSE) +
    scale_colour_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    scale_fill_manual(values = c("#0b775e", "#273046", "#a2a475")) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    theme_bw() +
  #  scale_y_continuous(limits = c(-0.24, 0.24),
   #                    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    #                   labels = c("-0.2", "-0.1", "0", "0.1", "0.2")) +
    raincloud_theme +
    theme_LPI3())

ggsave(rain_system2, filename = "figures/slopes_system.png", height = 5, width = 5)
