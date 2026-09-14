# script for filtering G4 and G6 protozoa data at the genus level
# Workflow:
#   1. for duplicated samples (two buffers) - taking the sample with higher total reads
#   2. removing samples with total reads < 500 for both G4 and G6
#   3. filtering for potentially parasitic protozoans
#   4. aggregating G4 and G6
#   5. adding the un-infected hosts

# output:
# protozoa_long_genus.csv - edgelist of individual hosts and parasites

#############################################################################

# Load necessary libraries
library(tidyverse)

rm(list=ls())


#############################################################################
# analysis
#############################################################################

# raw protozoa data
protozoa_data_raw <- read_csv("data/data_raw/parasites/NONHUMAN_PARASITE_ALL_VILLAGE_G4-6.csv")

# small mammals data
data_mammals <- read_csv("data/data_raw/small_mammals/Terrestrial_Mammals.csv") %>% 
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>% 
  dplyr::select(host_ID, field_identification, village, habitat_type, season) %>%
  dplyr::rename(host_species = field_identification, grid = habitat_type) %>% 
  mutate(season = factor(season, levels = c("1","2","3"))) %>% 
  mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth","agriculture","flooded_rice","agroforest","village")))


# organizing the table
protozoa_data <- protozoa_data_raw %>% 
  filter(str_starts(Sample_Name, "TMR-")) %>% 
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name))) %>% 
  group_by(host_ID) %>% # for duplicated samples, taking the one with the higher number of reads
  slice_max(order_by = reads_G4, n = 1, with_ties = FALSE) %>%
  ungroup()

# list of host with valid data (having enough reads)
full_host_ids <- protozoa_data %>% 
  filter(reads_G4>500 | reads_G6>500) %>% 
  select(host_ID)

# long format
protozoa_data_long <- protozoa_data %>% 
  pivot_longer(
    cols = matches("_(G4|G6)$") & !matches("^reads_", ignore.case = TRUE),
    names_to = c("final_id", "type"),
    names_pattern = "^(.*)_(G4|G6)$",
    values_to = "reads"
  ) %>% 
  filter(reads>0) %>% 
  select(-Sample_Name)

# only the parasitic protozoa (at the genus level)
protozoa_sp <- c("Hypotrichomonas","Entamoeba","Hexamastix","Tritrichomonas","Tetratrichomonas","Pentatrichomonas","Eimeria","Blastocystis","Cryptosporidium","Balantidium")

# filtering
protozoa_data_long_filtered <- protozoa_data_long %>% 
  filter(!(reads_G4 < 500 & type == "G4")) %>% # filter samples with less than 500 total reads
  filter(!(reads_G6 < 500 & type == "G6")) %>% 
  mutate(final_id = case_when(
    final_id == "Pentatrichomonas_hominis" ~ "Pentatrichomonas",
    final_id == "Neobalantidium_coli"      ~ "Balantidium",
    TRUE ~ final_id
  )) %>% 
  filter(final_id %in% protozoa_sp) %>% 
  group_by(host_ID, final_id) %>% 
  summarise(link=1) %>% 
  ungroup() 

# un-infected hosts
group_zero <- anti_join(full_host_ids, protozoa_data_long_filtered, by = "host_ID") %>% 
  mutate(link = 0)

# Add the uninfected hosts
protozoa_long_genus <- bind_rows(protozoa_data_long_filtered, group_zero) %>% 
  mutate(type = "Protozoa") %>% 
  left_join(data_mammals, by="host_ID")

write_csv(protozoa_long_genus, "data/data_processed/protozoa_long_genus.csv")
