## ============================================================================
##  Nematode ITS2 (NC1/NC2) metabarcoding: ASVs -> 97% OTUs -> LULU curation
##  -> host x OTU table for downstream analysis.
##
##  OVERVIEW
##  --------
##  Step 0  Load host metadata and the per-sample ASV table. Keep one sample per
##          host. in low-depth samples (<= 500 reads) ASV abundances are set to zero.
##  Step 1  Cluster ASVs into operational taxonomic units (OTUs) at 97% sequence
##          similarity with DECIPHER::Clusterize.
##  Step 2  Curate OTUs with LULU, which merges an OTU into a more abundant,
##          co-occurring OTU when they are highly similar in sequence, to reduce
##          over-splitting of closely related ITS2 variants.
##  Step 3  Assign each OTU the taxonomy of its representative ASV.
##  Step 4  Remove non-parasite taxa and build the final long-format host x OTU
##          table, retaining uninfected hosts as absences (link = 0).
##
##  INPUT FILES
##  -----------
##  data/data_raw/small_mammals/Terrestrial_Mammals.csv        host metadata
##  data/data_raw/parasites/ALL_NC_1pct_RRA_NIH_filtered_samples.csv
##                                                             per-sample ASV
##                                                             relative abundances
##                                                             and total read counts
##  data/data_raw/parasites/ASV_NC_all_village_1pct.fa         ASV sequences
##  data/data_raw/parasites/NCASVs_taxonomy80.tsv              ASV taxonomy
##
##  OUTPUT FILES
##  ------------
##  data/data_processed/NC_long_otu_lulu.csv        MAIN OUTPUT: one row per
##          host x OTU, plus one row per uninfected host. link = 1 (present) /
##          0 (absent); columns include otu_ID, genus, final_id and host metadata.
##
##  data/data_processed/lulu_out/
##      NC_OTU97_cluster_map.csv        ASV -> OTU assignment (all 97% OTUs)
##      NC_OTU97_table_uncurated.csv    OTU x host table before LULU
##      NC_OTU_curated_table.csv        OTU x host table after LULU
##      NC_curated_OTU_taxonomy.csv     taxonomy of each curated OTU
##      NC_lulu_merge_report.csv        OTUs merged by LULU and their taxonomy
##      lulu_result.rds                 full LULU result object
##
## ============================================================================

library(Biostrings) 
library(DECIPHER)
library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(lulu)
library(pwalign)

rm(list = ls())


## ============================================================================
## STEP 0 — configuration and data loading
## ============================================================================

fa_path <- "data/data_processed/NC_filtered.fa"          # ASV sequences kept for clustering
out_dir <- "data/data_processed/lulu_out"                # intermediate OTU/LULU outputs
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## host species of interest (all sampled small-mammal species)
sm_species <- c("Microgale brevicaudata","Setifer setosus","Rattus rattus","Mus musculus","Suncus etruscus","Microgale longicaudata","Suncus murinus",
                "Eliurus webbi","Eliurus minor","Tenrec ecaudatus","Nesomys audeberti","Eliurus cf. tanala","Nesomys cf. audeberti","Microgale prolixacaudata",
                "Eliurus cf. webbi","Microgale cf. prolixacaudata","Nesomys cf. rufus","Microgale parvula","Microgale cf. talazaci","Nesogale talazaci",
                "Eliurus cf. grandidieri","Nesogale cf. talazaci")

## host metadata (species, village, habitat, season), keyed by numeric host_ID
data_mammals <- read_csv("data/data_raw/small_mammals/Terrestrial_Mammals.csv") %>%
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>%
  dplyr::select(host_ID, field_identification, village, habitat_type, season) %>%
  dplyr::rename(host_species = field_identification, grid = habitat_type) %>%
  mutate(season = factor(season, levels = c("1","2","3"))) %>%
  mutate(grid = factor(grid, levels = c("semi-intact_forest","secondary_forest","brushy_regrowth",
                                        "agriculture","flooded_rice","agroforest","village")))

## per-sample ASV table: keep target host species and take one sample per host
## (the one with the most reads). In samples with 500 or fewer total reads, all
## ASV abundances are set to zero, so these samples contribute no nematode
## detections. ASVs left empty across all samples are then dropped.
NC_wide_asv <- read_csv("data/data_raw/parasites/ALL_NC_1pct_RRA_NIH_filtered_samples.csv") %>%
  filter(Species %in% sm_species) %>%
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name))) %>%
  group_by(host_ID) %>%
  slice_max(order_by = reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(c(host_ID, reads, starts_with("ASV"))) %>%
  filter(!is.na(reads)) %>%
  mutate(across(-1, ~ ifelse(reads <= 500, 0, .))) %>%
  select(-which(colSums(., na.rm = TRUE) == 0))

## ASV columns to cluster (all columns except host_ID and reads)
nc_asv_to_cluster <- colnames(NC_wide_asv[-c(1:2)])


## ============================================================================
## STEP 1 — cluster ASVs into 97% OTUs (DECIPHER)
## ============================================================================

## write the ASV sequences that passed filtering, then read them back in
seq_fa  <- ape::read.FASTA(file = "data/data_raw/parasites/ASV_NC_all_village_1pct.fa")
seq_fa2 <- seq_fa[names(seq_fa) %in% nc_asv_to_cluster]
ape::write.FASTA(seq_fa2, file = fa_path)

seqs <- readDNAStringSet(fa_path)
seqs <- OrientNucleotides(seqs, processors = NULL)       # put all sequences on the same strand

set.seed(123)                                            # reproducible clustering
clusters <- Clusterize(
  seqs,
  cutoff                   = 0.03,                        # 0.03 distance = 97% similarity
  method                   = "overlap",
  penalizeGapLetterMatches = NA,
  invertCenters            = TRUE,                        # mark cluster representatives (negative id)
  processors               = NULL)

## tidy cluster table: one row per ASV, with its OTU id and a flag for the
## representative (centroid) ASV of each OTU
NC_OTUs_97 <- clusters %>%
  rownames_to_column("asv_ID") %>%
  dplyr::rename(otu_ID = cluster) %>%
  mutate(type   = "NC",
         center = otu_ID < 0,
         otu_ID = paste0("nc_", abs(otu_ID)))
message("OTUs at 97%: ", dplyr::n_distinct(NC_OTUs_97$otu_ID),
        "  (from ", nrow(NC_OTUs_97), " ASVs)")


## ============================================================================
## STEP 2 — LULU post-clustering curation
## ============================================================================

## 2a. representative sequence of each OTU, named by OTU id
centre   <- NC_OTUs_97 %>% filter(center)
rep_seqs <- seqs[centre$asv_ID] 
names(rep_seqs) <- centre$otu_ID

## 2b. sum ASV relative abundances within each OTU to get an OTU x host table
M   <- as.matrix(NC_wide_asv[, setdiff(colnames(NC_wide_asv), c("host_ID", "reads"))])
rownames(M) <- NC_wide_asv$host_ID
otu_of_asv  <- setNames(NC_OTUs_97$otu_ID, NC_OTUs_97$asv_ID)[colnames(M)]
stopifnot(!anyNA(otu_of_asv))
otutab <- as.data.frame(rowsum(t(M), group = otu_of_asv))[names(rep_seqs), ]
write.csv(cbind(otu_ID = rownames(otutab), otutab),
          file.path(out_dir, "NC_OTU97_table_uncurated.csv"), row.names = FALSE)

## 2c. build the pairwise-similarity match list LULU needs.
## For every pair of OTU representatives above a k-mer similarity prefilter, we
## compute the aligned percent identity. An overlap alignment and a >= 90%
## coverage requirement accommodate the wide range of ITS2 lengths (~270-550 bp).
match_list_R <- function(rep_seqs, min_id = 84, min_cov = 0.9, kmer = 5,
                         prefilter = 0.60, gapOpening = 5, gapExtension = 2) {
  ids <- names(rep_seqs); L <- width(rep_seqs)
  ## k-mer cosine similarity to shortlist candidate pairs (avoids all-vs-all alignment)
  K   <- oligonucleotideFrequency(rep_seqs, width = kmer)
  Kn  <- K / sqrt(rowSums(K^2))
  cs  <- tcrossprod(Kn); diag(cs) <- 0
  cand <- which(cs >= prefilter, arr.ind = TRUE)
  cand <- cand[cand[, 1] < cand[, 2], , drop = FALSE]
  ## align each shortlisted pair and compute percent identity + coverage
  sm  <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  pa  <- pwalign::pairwiseAlignment(rep_seqs[cand[, 1]], rep_seqs[cand[, 2]],
                    type = "overlap", substitutionMatrix = sm,
                    gapOpening = gapOpening, gapExtension = gapExtension)
  idv  <- pwalign::pid(pa, type = "PID1")
  covg <- pmin(1, nchar(pa) / pmin(L[cand[, 1]], L[cand[, 2]]))
  keep <- idv >= min_id & covg >= min_cov
  ## LULU expects both directions of each pair
  ml <- data.frame(OTUid = ids[cand[keep, 1]], hit = ids[cand[keep, 2]], match = idv[keep])
  rbind(ml, data.frame(OTUid = ml$hit, hit = ml$OTUid, match = ml$match))
}
matchlist <- match_list_R(rep_seqs)

## 2d. run LULU. An OTU is merged into a more abundant "parent" OTU when they
## share >= 84% identity, co-occur in >= 90% of the samples where the OTU is
## found, and have a compatible abundance ratio.
res <- lulu(otutab, matchlist,
            minimum_ratio_type           = "avg",
            minimum_ratio                = 0.5,
            minimum_match                = 84,
            minimum_relative_cooccurence = 0.9)
message("LULU: kept ", res$curated_count, " OTUs; merged ", res$discarded_count)
saveRDS(res, file.path(out_dir, "lulu_result.rds"))


## ============================================================================
## STEP 3 — OTU taxonomy and LULU merge report
## ============================================================================

rank_cols <- c("Phylum","Class","Order","Family","Genus","Species")
tax <- read_tsv("data/data_raw/parasites/NCASVs_taxonomy80.tsv") %>%
  select(ASV, all_of(rank_cols)) %>% dplyr::rename(asv_ID = ASV)

## each OTU takes the taxonomy of its representative ASV
otu_tax <- centre %>% left_join(tax, by = "asv_ID") %>% mutate(rep_ASV = asv_ID)

curated     <- res$curated_table %>% rownames_to_column("otu_ID")
curated_tax <- otu_tax %>% filter(otu_ID %in% res$curated_otus)
write.csv(NC_OTUs_97,  file.path(out_dir, "NC_OTU97_cluster_map.csv"),   row.names = FALSE)
write.csv(curated,     file.path(out_dir, "NC_OTU_curated_table.csv"),   row.names = FALSE)
write.csv(curated_tax, file.path(out_dir, "NC_curated_OTU_taxonomy.csv"), row.names = FALSE)

## report of merged OTUs, with the taxonomy of both the merged OTU and its parent
merge_report <- res$otu_map %>%
  rownames_to_column("otu_ID") %>% filter(curated == "merged") %>%
  transmute(daughter = otu_ID, parent = parent_id, spread, total) %>%
  left_join(otu_tax, by = c("daughter" = "otu_ID")) %>%
  rename_with(~paste0("daughter_", .), all_of(c("rep_ASV", rank_cols))) %>%
  left_join(otu_tax, by = c("parent" = "otu_ID")) %>%
  rename_with(~paste0("parent_", .), all_of(c("rep_ASV", rank_cols))) %>%
  mutate(same_genus   = daughter_Genus   == parent_Genus,
         same_species = daughter_Species == parent_Species)
write.csv(merge_report, file.path(out_dir, "NC_lulu_merge_report.csv"), row.names = FALSE)


## ============================================================================
## STEP 4 — final long-format host x OTU table
## ============================================================================

## map each original 97% OTU to the curated OTU it belongs to after LULU
otu_to_parent <- res$otu_map %>% rownames_to_column("otu_ID") %>%
  dplyr::select(otu_ID, merged_otu_ID = parent_id)

## label each curated OTU:
##   genus       - genus of the representative ASV (Viannaia relabelled to
##                 Heligmosomoidea to correct a reference-database misannotation)
##   base_taxon  - lowest confidently assigned rank, floored at genus
##   final_taxon - base_taxon, numbered when several OTUs share the same label
parent_total <- tibble(merged_otu_ID = rownames(res$curated_table),
                       otu_total     = rowSums(res$curated_table))
parent_taxon <- otu_tax %>% filter(otu_ID %in% res$curated_otus) %>%
  transmute(merged_otu_ID = otu_ID,
            genus      = if_else(Genus == "Viannaia", "Heligmosomoidea", Genus),
            base_taxon = coalesce(if_else(Genus == "Viannaia", "Heligmosomoidea", Genus),
                                  Family, Order, Class, Phylum, "Unassigned")) %>%
  mutate(has_genus = !is.na(genus)) %>%
  left_join(parent_total, by = "merged_otu_ID") %>%
  group_by(base_taxon) %>% arrange(desc(otu_total), .by_group = TRUE) %>%
  mutate(final_taxon = if (n() > 1) paste(base_taxon, row_number(), sep = "_") else base_taxon) %>%
  ungroup() %>% dplyr::select(merged_otu_ID, genus, base_taxon, has_genus, final_taxon)

## attach curated-OTU id and labels back to every ASV
NC_OTUs_97_new <- NC_OTUs_97 %>% dplyr::select(asv_ID, otu_ID) %>%
  left_join(otu_to_parent, by = "otu_ID") %>%
  left_join(parent_taxon,  by = "merged_otu_ID") %>%
  dplyr::select(asv_ID, otu_ID = merged_otu_ID, genus, base_taxon, has_genus, final_taxon)

## taxa to exclude: genera that are not parasites of small mammals
nc_exclude <- c("Caenorhabditis", "Oscheius", "Helicotylenchus",
                "Syngamus", "Ancylostoma", "Necator")

## host x curated-OTU relative abundances (present OTUs only), excluding non-parasites
core_otu <- NC_wide_asv %>%
  pivot_longer(starts_with("ASV"), names_to = "asv_ID", values_to = "relative_reads") %>%
  filter(relative_reads > 0) %>%
  left_join(NC_OTUs_97_new, by = "asv_ID") %>%
  filter(!(base_taxon %in% nc_exclude)) %>%
  group_by(host_ID, otu_ID, genus, base_taxon, has_genus, final_taxon) %>%
  summarise(relative_reads = sum(relative_reads), .groups = "drop")

## every retained host. those with no parasite OTU become absences (link = 0)
NC_host <- tibble(host_ID = NC_wide_asv$host_ID) %>%
  left_join(data_mammals, by = "host_ID") %>% mutate(type = "NC", link = 0)

## present host x OTU records (link = 1), combined with uninfected hosts (link = 0)
NC_long_otu <- core_otu %>%
  mutate(link = 1, type = "NC",
         final_taxon_otu_ID = paste(otu_ID, final_taxon, sep = "_")) %>%
  dplyr::rename(final_id = final_taxon) %>%
  left_join(data_mammals, by = "host_ID")
NC_long_otu2 <- bind_rows(NC_long_otu, anti_join(NC_host, NC_long_otu, by = "host_ID"))

write_csv(NC_long_otu2, "data/data_processed/NC_long_otu_lulu.csv")
