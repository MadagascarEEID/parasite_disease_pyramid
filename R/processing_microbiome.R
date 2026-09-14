## ============================================================================
##  Gut microbiome (16S) metabarcoding: filter and normalize ASVs into a
##  long-format host x ASV relative-abundance table for downstream analysis.
##
##  OVERVIEW
##  --------
##  Step 0  Load host metadata and the per-sample ASV table, keep one sample per
##          host.
##  Step 1  Keep only the target host species.
##  Step 2  Remove non-bacterial, chloroplast, and mitochondrial ASVs (taxonomy).
##  Step 3  Within each sample, zero ASVs below a 0.1% relative-abundance floor,
##          then drop ASVs that do not reach the floor in a minimum number of
##          samples.
##  Step 4  Remove samples (hosts) with low total reads (< 5000).
##  Step 5  Recompute per-sample relative abundances and write the long-format
##          output.
##
##  INPUT FILES
##  -----------
##  data/data_raw/small_mammals/Terrestrial_Mammals.csv        host metadata
##  data/data_raw/microbiome/merged_full_sample_table.csv      per-sample ASV counts
##  data/data_raw/microbiome/ASVs_taxonomy_new.tsv             ASV taxonomy
##
##  OUTPUT FILE
##  -----------
##  data/data_processed/data_microbiome_rra0.001_th5000_2.csv
##          Long format: one row per host x ASV with relative abundance > 0.
##          Columns: host metadata (host_ID, host_species, village, grid, season),
##          total_reads, asv_ID, and reads (within-sample relative abundance).
##
## ============================================================================

library(tidyverse)
library(magrittr)   

rm(list = ls())

## helper: drop ASV columns that are zero in every sample (leaves metadata intact)
drop_empty_asv <- function(df) {
  empty_asv <- df %>%
    dplyr::select(starts_with("ASV_")) %>%
    dplyr::select(where(~ all(.x == 0))) %>%
    names()
  dplyr::select(df, -all_of(empty_asv))
}


## ============================================================================
## STEP 0 — configuration and data loading
## ============================================================================

## target host species
sm_species <- c("Rattus rattus", "Microgale brevicaudata")

## host metadata
data_mammals <- read_csv("data/data_raw/small_mammals/Terrestrial_Mammals.csv")

## per-sample ASV counts. keep true samples only
data_asv <- read_csv("data/data_raw/microbiome/merged_full_sample_table.csv")
data_asv %<>% filter(sample_type == "SAMPLE")

## Numeric host IDs are extracted from mixed-format identifiers. Restrict the
## mammal table to TMR animals *before* extracting the number, so an animal from
## another series whose number happens to coincide cannot be matched by mistake.
data_mammals %<>%
  filter(grepl("^TMR", animal_id)) %>%
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id)))

## one sample per host (the one with the most reads)
data_asv %<>% mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name))) %>%
  group_by(host_ID) %>%
  slice_max(order_by = unfiltered_reads, n = 1, with_ties = FALSE) %>%
  ungroup()

## host metadata for the small mammals that have microbiome data
data_sm <- semi_join(data_mammals, data_asv, by = "host_ID") %>%
  dplyr::select(host_ID, field_identification, village, habitat_type, season) %>%
  dplyr::rename(host_species = field_identification, grid = habitat_type) %>%
  mutate(season = factor(season, levels = c("1","2","3"))) %>%
  mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth",
                                        "agriculture","flooded_rice","agroforest","village")))

## require one metadata record per host. otherwise the join below would duplicate
## a host's ASV row once per duplicate record (e.g. recaptures)
if (any(duplicated(data_sm$host_ID))) {
  stop("Duplicated host_ID in data_sm (>1 mammal record per host) -- the ASV join would ",
       "duplicate rows. Resolve recaptures/duplicates before joining.")
}

## join host metadata to ASV counts
data_asv_f <- data_asv %>%
  dplyr::select(host_ID, unfiltered_reads, contains("ASV"))

dat <- left_join(data_sm, data_asv_f, by = "host_ID")


## ============================================================================
## STEP 1 — keep target host species
## ============================================================================
dat1 <- dat %>%
  filter(host_species %in% sm_species) %>%
  drop_empty_asv()


## ============================================================================
## STEP 2 — remove non-bacterial, chloroplast and mitochondrial ASVs
## ============================================================================
tax <- read_delim("data/data_raw/microbiome/ASVs_taxonomy_new.tsv") %>%
  dplyr::rename(asv_ID = ASV)

## ASVs to remove: not Bacteria, chloroplast, mitochondria, or unassigned kingdom
tax_exclude <- tax %>%
  filter(asv_ID %in% colnames(dat1)) %>%
  filter(Kingdom != "Bacteria" | Order == "Chloroplast" | Family == "Mitochondria" | is.na(Kingdom))

dat2 <- dat1 %>%
  select(-all_of(tax_exclude$asv_ID))


## ============================================================================
## STEP 3 — within-sample relative-abundance floor + minimum prevalence
## ============================================================================

asv_rel_reads_th <- 0.001   # within-sample detection floor (0.1%)

## keep an ASV only if it reaches the floor in at least this many samples.
##   1 = keep if at or above the floor in any single sample
##   2 = drop single-sample detections
min_prevalence <- 2

## Relative abundance is computed against the BACTERIAL (post-taxonomy) per-sample
## total, consistent with the final relative abundances below -- not against the
## raw unfiltered reads, which still included the reads removed in Step 2.
bact_total <- rowSums(dplyr::select(dat2, starts_with("ASV_")))

## Some samples have zero bacterial reads. the division below is guarded so that
## 0/0 -> 0 rather than NaN. These samples stay all-zero and are removed by the
## depth filter (Step 4).
n_empty <- sum(bact_total == 0)
if (n_empty > 0) message(n_empty, " sample(s) with zero bacterial reads (dropped later by the depth filter).")

## apply the floor: counts -> relative abundance -> zero below floor -> back to counts
dat3 <- dat2 %>%
  mutate(across(starts_with("ASV_"), ~ ifelse(bact_total > 0, .x / bact_total, 0))) %>%
  mutate(across(starts_with("ASV_"), ~ ifelse(.x < asv_rel_reads_th, 0, .x))) %>%
  mutate(across(starts_with("ASV_"), ~ .x * bact_total))

## drop ASVs that do not reach the floor in at least `min_prevalence` samples
asv_prevalence <- dat3 %>%
  dplyr::select(starts_with("ASV_")) %>%
  summarise(across(everything(), ~ sum(.x > 0)))
keep_asv <- names(asv_prevalence)[as.numeric(asv_prevalence[1, ]) >= min_prevalence]
dat3 <- dat3 %>% dplyr::select(!starts_with("ASV_"), all_of(keep_asv))

## per-sample total reads after filtering. drop the now-outdated raw read count
host_total_reads <- dat3 %>%
  select(starts_with("ASV_")) %>%
  rowSums()

dat3_updated <- dat3 %>%
  mutate(total_reads = host_total_reads) %>%
  select(-unfiltered_reads)


## ============================================================================
## STEP 4 — remove low-depth samples (< 5000 total reads)
## ============================================================================

## diagnostic: distribution of per-sample total reads, with the threshold marked
dat3_updated %>%
  ggplot(aes(x = total_reads)) +
  geom_histogram(binwidth = 1000) +
  facet_wrap(~host_species) +
  geom_vline(xintercept = 5000) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, color = 'black'),
        title = element_text(size = 20), strip.text.x = element_text(size = 12)) +
  labs(x = "Total Reads", y = "Count")

total_reads_th <- 5000
dat4 <- dat3_updated %>%
  filter(total_reads >= total_reads_th) %>%
  drop_empty_asv()

message("hosts removed by depth filter: ",
        length(unique(dat3_updated$host_ID)) - length(unique(dat4$host_ID)))


## ============================================================================
## STEP 5 — recompute relative abundances and write long-format output
## ============================================================================

## per-sample totals on the depth-filtered data, then convert counts -> relative abundance
host_total_reads <- dat4 %>%
  select(starts_with("ASV_")) %>%
  rowSums()

dat4_final <- dat4 %>%
  mutate(total_reads = host_total_reads) %>%
  mutate(across(starts_with("ASV_"), ~ .x / total_reads))

## long format: one row per host x ASV, keeping present ASVs only
dat4_final_long <- dat4_final %>%
  pivot_longer(cols = starts_with("ASV_"), names_to = "asv_ID", values_to = "reads") %>%
  filter(reads > 0)

write_csv(dat4_final_long, "data/data_processed/data_microbiome_rra0.001_th5000_2.csv")
