#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_TAXONOMY
# (modules/local/dada2_taxonomy.nf) with conf/modules.config
# DADA2_TAXONOMY ext.args defaults (minBoot, tryRC, seed).
# Runs once per ASV fasta chunk inside the merged dada2_taxonomy
# rule; output names mirror upstream "${fasta.baseName}${outfile}".
#
# Usage: Rscript dada2_assign_taxonomy.R <fasta> <database>
#     <outfile_suffix> <taxlevels_csv> <seed> <threads> <min_boot>
#     <try_rc> <out_dir>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
fasta <- argv[1]
database <- argv[2]
outfile <- argv[3]
taxlevels <- unlist(strsplit(argv[4], ","))
seed <- argv[5]
threads <- as.integer(argv[6])
min_boot <- argv[7]
try_rc <- as.logical(argv[8])
out_dir <- argv[9]

set.seed(as.integer(seed)) # Initialize random number generator for reproducibility

seq <- getSequences(fasta, collapse = TRUE, silence = FALSE)
taxa <- assignTaxonomy(seq, database, taxLevels = taxlevels, minBoot = as.integer(min_boot), tryRC = try_rc, multithread = threads, verbose=TRUE, outputBootstraps = TRUE)

# (1) Make a data frame, add ASV_ID from seq
tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa$tax), row.names = names(seq))

# (2) Set confidence to the bootstrap for the most specific taxon
# extract columns with taxonomic values
tax <- tx[ , grepl( "tax." , names( tx ) ) ]
# find first occurrence of NA
res <- max.col(is.na(tax), ties = "first")
# correct if no NA is present in column to NA
if(any(res == 1)) is.na(res) <- (res == 1) & !is.na(tax[[1]])
# find index of last entry before NA, which is the bootstrap value
res <- res-1
# if NA choose last entry
res[is.na(res)] <- ncol(tax)
# extract bootstrap values
boot <- tx[ , grepl( "boot." , names( tx ) ) ]
boot$last_tax <- res
valid_boot <- apply(boot,1,function(x) x[x[length(x)]][1]/100 )
# replace missing bootstrap values (NA) with 0
valid_boot[is.na(valid_boot)] <- 0
# add bootstrap values to column confidence
tx$confidence <- valid_boot

# (3) Reorder columns before writing to file
expected_order <- c("ASV_ID",paste0("tax.",taxlevels),"confidence","sequence")
expected_order <- intersect(expected_order,colnames(tx))
taxa_export <- subset(tx, select = expected_order)
colnames(taxa_export) <- sub("tax.", "", colnames(taxa_export))
rownames(taxa_export) <- names(seq)

fasta_base <- sub("\\.fasta$", "", basename(fasta))
write.table(taxa_export, file = file.path(out_dir, paste0(fasta_base, outfile, ".tsv")), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

# Save a version with rownames for addSpecies
taxa_export <- cbind( ASV_ID = tx$ASV_ID, taxa$tax, confidence = tx$confidence)
saveRDS(taxa_export, file.path(out_dir, paste0(fasta_base, outfile, ".rds")))
