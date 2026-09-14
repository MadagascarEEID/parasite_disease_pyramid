# script for combining, filtering, and curating parasite data

library(tidyverse)

rm(list=ls())


#############################################################################
# parasite data
#############################################################################
# reading the data and filtering to host species and village

# defining the species of interest
sm_species <- c("Rattus rattus", "Microgale brevicaudata")

# protozoa
data_prt <- read_csv("data/data_processed/protozoa_long_genus.csv") %>% 
  filter(host_species %in% sm_species, village == "Mandena") %>% 
  #mutate(otu_ID = final_id) %>% 
  select(host_ID, final_id, type, link, host_species, village, grid)

# nematode
data_nc <- read_csv("data/data_processed/NC_long_otu_lulu.csv") %>% 
  filter(host_species %in% sm_species, village == "Mandena") %>% 
  select(host_ID, final_id, type, link, host_species, village, grid) %>% 
  mutate(type = "Nematode")

# microbiome
data_microbiome <- read_csv("data/data_processed/data_microbiome_rra0.001_th5000_2.csv")


#############################################################################
# filtering
#############################################################################

# filtering hosts with full information (protozoa + nematode + microbiome)

# matching g4-6, nc, and microbiome host IDs
host_match <- base::intersect(base::intersect(data_nc$host_ID, data_prt$host_ID), data_microbiome$host_ID)

# filtering protozoa data
data_prt_filtered <- data_prt %>% 
  filter(host_ID %in% host_match)

# filtering nematode data
data_nc_filtered <- data_nc %>% 
  filter(host_ID %in% host_match)

# final parasite table
data_parasite <- data_prt_filtered %>% 
  bind_rows(data_nc_filtered)

data_parasite <- data_parasite %>%
  mutate(final_id = if_else(final_id %in% c("Trichostrongylidae_1", "Trichostrongylidae_2"),
                            "Trichostrongylidae", final_id)) %>%   # merging 2 Trichostrongylidae OTUs
  mutate(final_id = if_else(final_id =="Heligmosomoidea_1",
                            "Heligmosomoidea", final_id)) %>% 
  distinct(host_ID,final_id,host_species,village,grid,link,type)

###########
# filtering parasites with lower than minimum degree
parasite_min <- 10

parasite_degree <- data_parasite %>%
  group_by(final_id) %>%
  summarise(n = sum(link))

table(parasite_degree$n)

data_parasite_filtered <- data_parasite %>%
  left_join(parasite_degree, by=c("final_id")) %>% 
  filter(n >= parasite_min)


###########
# Keep one row per host with host-level metadata
host_info <- data_parasite %>%
  distinct(host_ID, host_species, village, grid)

# Keep only real parasite links for the wide parasite matrix
parasite_links <- data_parasite_filtered %>%
  filter(!is.na(final_id)) %>%
  select(host_ID, final_id, link) %>%
  distinct() %>% 
  group_by(host_ID, final_id) %>%
  summarise(link = max(link), .groups = "drop")

# Wide format: one row per host, one column per parasite final_id
data_parasite_wide <- host_info %>%
  left_join(parasite_links, by = "host_ID") %>%
  pivot_wider(
    names_from = final_id,
    values_from = link,
    values_fill = 0,
    values_fn = max
  ) %>%
  select(-any_of("NA"))

# parasite metadata
parasite_info <- data_parasite %>% 
  distinct(final_id, type)

# Long format, including uninfected hosts
data_parasite_long <- data_parasite_wide %>%
  pivot_longer(
    cols = -c(host_ID, host_species, village, grid),
    names_to = "final_id",
    values_to = "link"
  ) %>%
  left_join(parasite_info, by = "final_id")


###########
# saving the data
write_csv(data_parasite_wide, "data/data_processed/parasites_wide_10_lulu.csv")
write_csv(data_parasite_long, "data/data_processed/parasites_long_10_lulu.csv")

