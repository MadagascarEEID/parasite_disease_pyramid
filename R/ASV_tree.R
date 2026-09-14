## ============================================================================
##  Build a phylogenetic tree of the microbiome ASVs for downstream
##  phylogenetic diversity metrics (e.g. Faith's PD, UniFrac).
##
##  OVERVIEW
##  --------
##  Step 0  Read the ASV sequences and the list of ASVs retained by the
##          microbiome filtering pipeline, and subset the sequences to them.
##  Step 1  Align the sequences (DECIPHER) and, optionally, mask poorly aligned
##          (highly gapped) columns.
##  Step 2  Build an approximate maximum-likelihood tree (FastTree, GTR+Gamma),
##          midpoint-root it, and save it.
##  Step 3  Plot the tree, coloured by phylum, as a visual check.
##
##  INPUT FILES
##  -----------
##  data/data_raw/microbiome/ASV_merged_full.fa                ASV sequences
##  data/data_processed/data_microbiome_rra0.001_th5000_2.csv  retained ASVs
##                                                             (from processing_microbiome.R)
##  data/data_raw/microbiome/ASVs_taxonomy_new.tsv             ASV taxonomy (for plotting)
##
##  OUTPUT FILES
##  ------------
##  data/data_processed/ASV_tree_ML_rooted.rds     MAIN OUTPUT: midpoint-rooted
##                                                 tree (R object) for downstream use
##  data/data_processed/ASV_tree_fasttree.nwk      FastTree output (Newick)
##  data/data_raw/microbiome/ASV_merged_filtered.fa  filtered ASV sequences
##  data/data_raw/microbiome/ASV_alignment.fasta     multiple-sequence alignment
##
## ============================================================================

library(tidyverse)
library(magrittr)
library(ape)
library(phangorn)
library(ggtree)
library(Biostrings)

rm(list = ls())


## ============================================================================
## STEP 0 — load sequences and subset to the retained ASVs
## ============================================================================

## full ASV sequences
seqs_all <- Biostrings::readDNAStringSet("data/data_raw/microbiome/ASV_merged_full.fa")

## ASVs retained by the microbiome filtering pipeline
data_asv_filtered <- read_csv("data/data_processed/data_microbiome_rra0.001_th5000_2.csv")
asv_ids <- unique(data_asv_filtered$asv_ID)

## every retained ASV must have a sequence, otherwise the tree would silently
## omit tips
missing_asv <- setdiff(asv_ids, names(seqs_all))
if (length(missing_asv) > 0) {
  stop(length(missing_asv), " filtered ASV(s) are absent from the fasta, e.g. ",
       paste(head(missing_asv), collapse = ", "))
}

## subset by name: filters and orders the sequences in one correspondence-safe step
seqs <- seqs_all[asv_ids]

## keep the filtered sequences as an artifact
Biostrings::writeXStringSet(seqs, filepath = "data/data_raw/microbiome/ASV_merged_filtered.fa", format = "fasta")

cat("ASVs in tree:", length(seqs), "\n")


## ============================================================================
## STEP 1 — align and (optionally) mask the alignment
## ============================================================================

aligned <- DECIPHER::AlignSeqs(seqs)
# DECIPHER::BrowseSeqs(aligned)   # uncomment to visually inspect the alignment

aligned_dna <- as(aligned, "DNAStringSet")

## Mask poorly aligned / highly gapped columns before tree building. such columns
## add noise to branch lengths, which drive phylogenetic diversity metrics. Set
## do_masking <- FALSE to build on the raw alignment.
do_masking <- TRUE
if (do_masking) {
  masked <- DECIPHER::MaskAlignment(aligned,
                                    maxFractionGaps = 0.4,   # drop columns > 40% gaps
                                    threshold       = 1,     # information-content threshold
                                    showPlot        = FALSE)
  aligned_dna <- as(masked, "DNAStringSet")                 # coercion removes masked columns
  cat("Alignment length: ", width(aligned)[1],
      " -> after masking: ", width(aligned_dna)[1], " columns\n", sep = "")
}

## write the alignment for FastTree
Biostrings::writeXStringSet(aligned_dna, filepath = "data/data_raw/microbiome/ASV_alignment.fasta", format = "fasta")


## ============================================================================
## STEP 2 — build, root, and save the tree
## ============================================================================

## approximate maximum-likelihood tree, GTR+Gamma model, nucleotides (-nt)
tree_cmd <- paste(
  "FastTree -nt -gtr -gamma",
  "data/data_raw/microbiome/ASV_alignment.fasta",
  "> data/data_processed/ASV_tree_fasttree.nwk"
)
status <- system(tree_cmd)

## fail loudly if FastTree did not run or did not produce an output file
if (status != 0 || !file.exists("data/data_processed/ASV_tree_fasttree.nwk")) {
  stop("FastTree failed (exit status ", status, "). ",
       "Check that FastTree is installed and on PATH (for large trees, FastTreeDbl/MP are recommended).")
}

## read back, midpoint-root, and save
tree        <- read.tree("data/data_processed/ASV_tree_fasttree.nwk")
tree_rooted <- phangorn::midpoint(tree)
saveRDS(tree_rooted, "data/data_processed/ASV_tree_ML_rooted.rds")


## ============================================================================
## STEP 3 — plot the tree (visual check)
## ============================================================================

tax <- read_delim("data/data_raw/microbiome/ASVs_taxonomy_new.tsv") %>%
  dplyr::rename(label = ASV) %>%
  select(-any_of("...1")) %>%          # drop the unnamed index column if present
  relocate(label, .before = 1)

p <- ggtree(tree_rooted, layout = "fan") %<+% tax +
  geom_tippoint(aes(color = Phylum), size = 1) +
  geom_tiplab(label = "", size = 0) +
  ggtitle("Midpoint-rooted ASV phylogenetic tree") +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

p
