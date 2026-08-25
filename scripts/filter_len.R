#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module FILTER_LEN
# (modules/local/filter_len.nf) with the FILTER_LEN_ITSX ext.args
# (ext.min_len_asv = 50, max_len_asv default 1000000). Runs in a
# scratch workdir; all output names are fixed (ASV_seqs.len.fasta,
# ASV_table.len.tsv, stats.len.tsv, ASV_len_orig.tsv, ASV_len_filt.tsv).
#
# Usage: Rscript filter_len.R <fasta> <ASV_table.tsv>
suppressPackageStartupMessages(library(Biostrings))

argv <- commandArgs(trailingOnly = TRUE)
fasta <- argv[1]
table_file <- argv[2]

min_len_asv <- 50
max_len_asv <- 1000000

# read abundance file, first column is ASV_ID
table <- read.table(file = table_file, sep = "\t", comment.char = "", header = TRUE)
colnames(table)[1] <- "ASV_ID"

# read fasta file of ASV sequences
input_seq <- readDNAStringSet(fasta)
input_seq <- data.frame(ID = names(input_seq), sequence = paste(input_seq))

# filter
filtered_seq <- input_seq[nchar(input_seq$sequence) %in% min_len_asv:max_len_asv, ]
id_list <- filtered_seq[, "ID", drop = FALSE]
filtered_table <- merge(table, id_list, by.x = "ASV_ID", by.y = "ID", all.x = FALSE, all.y = TRUE)

# report
distribution_before <- table(nchar(input_seq$sequence))
distribution_before <- data.frame(Length = names(distribution_before), Counts = as.vector(distribution_before))
distribution_after <- table(nchar(filtered_seq$sequence))
distribution_after <- data.frame(Length = names(distribution_after), Counts = as.vector(distribution_after))

# write
write.table(filtered_table, file = "ASV_table.len.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE, na = "")
write.table(data.frame(s = sprintf(">%s\n%s", filtered_seq$ID, filtered_seq$sequence)), "ASV_seqs.len.fasta", col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")
write.table(distribution_before, file = "ASV_len_orig.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE, na = "")
write.table(distribution_after, file = "ASV_len_filt.tsv", row.names = FALSE, sep = "\t", col.names = TRUE, quote = FALSE, na = "")

# stats
stats <- as.data.frame(t(rbind(colSums(table[-1]), colSums(filtered_table[-1]))))
stats$ID <- rownames(stats)
colnames(stats) <- c("lenfilter_input", "lenfilter_output", "sample")
write.table(stats, file = "stats.len.tsv", row.names = FALSE, sep = "\t")
