#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_ADDSPECIES
# (modules/local/dada2_addspecies.nf) with conf/modules.config
# DADA2_ADDSPECIES ext.args defaults (n, tryRC; allowMultiple via
# params.dada_addspecies_allowmultiple).
# Runs once per ASV fasta chunk inside the merged dada2_taxonomy
# rule; output name mirrors upstream "${taxtable.baseName}.species.tsv".
#
# Usage: Rscript dada2_add_species.R <taxtable.rds> <database>
#     <taxlevels_csv> <seed> <allow_multiple> <try_rc> <out_dir>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
taxtable_file <- argv[1]
database <- argv[2]
taxlevels <- unlist(strsplit(argv[3], ","))
seed <- argv[4]
allow_multiple <- as.logical(argv[5])
try_rc <- as.logical(argv[6])
out_dir <- argv[7]

# add species to taxlevels if not already present
if ( !"Species" %in% taxlevels ) {
    taxlevels <- c(taxlevels, "Species")
}

set.seed(as.integer(seed)) # Initialize random number generator for reproducibility

taxtable <- readRDS(taxtable_file)
taxtable <- data.frame(taxtable)

# remove Species to prevent bias
taxa_nospecies <- taxtable[,!colnames(taxtable) %in% 'Species']

tx <- addSpecies(taxa_nospecies, database, n = 1e5, allowMultiple = allow_multiple, tryRC = try_rc, verbose=TRUE)

tmp <- data.frame(row.names(tx))
expected_order <- c("ASV_ID",taxlevels,"confidence")
taxa <- as.data.frame( subset(tx, select = expected_order) )
taxa$sequence <- tmp[,1]
row.names(taxa) <- row.names(tmp)

# correct column name for exact species match
colnames(taxa)[which(names(taxa) == "Species")] <- "Species_exact"

# combine with added Species
if ( "Species" %in% colnames(taxtable) ) {
    taxa_export <- data.frame(append(taxa, list(Species=taxtable$Species), after=match("Genus", names(taxa))))
} else {
    taxa_export <- taxa
}

taxtable_base <- sub("\\.rds$", "", basename(taxtable_file))
write.table(taxa_export, file = file.path(out_dir, paste0(taxtable_base, ".species.tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')
