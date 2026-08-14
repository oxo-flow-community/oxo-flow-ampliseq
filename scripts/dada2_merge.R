#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_MERGE
# (modules/local/dada2_merge.nf). Default single-run path:
# upstream globs *.stats.tsv / *.ASVtable.rds from its cwd; the
# port receives the two run-level files as positional argv.
#
# Usage: Rscript dada2_merge.R <stats.tsv> <ASVtable.rds> <out_dir>
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(digest))

argv <- commandArgs(trailingOnly = TRUE)
stats_file <- argv[1]
asvtab_file <- argv[2]
out_dir <- argv[3]

# combine stats files (single-run default path: one file)
if (!exists("stats")){ stats <- read.csv(stats_file, header=TRUE, sep="\t") }
write.table( stats, file = file.path(out_dir, "DADA2_stats.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

# combine dada-class objects (single-run default path: readRDS)
files <- asvtab_file
if ( length(files) == 1 ) {
    ASVtab = readRDS(files)
} else {
    ASVtab <- mergeSequenceTables(tables=files, repeats = "error", orderBy = "abundance", tryRC = FALSE)
}
saveRDS(ASVtab, file.path(out_dir, "DADA2_table.rds"))

df <- t(ASVtab)
colnames(df) <- gsub('_1.filt.fastq.gz', '', colnames(df))
colnames(df) <- gsub('.filt.fastq.gz', '', colnames(df))
df <- data.frame(sequence = rownames(df), df, check.names=FALSE)
# Create an md5 sum of the sequences as ASV_ID and rearrange columns
df$ASV_ID <- sapply(df$sequence, digest, algo='md5', serialize = FALSE)
df <- df[,c(ncol(df),3:ncol(df)-1,1)]

# order ASV table by ASV_ID (md5sum of sequence), previously sorted decreasing by total counts (ASVs of same counts were randomly sorted)
df <- df[order(df$ASV_ID),]

# file to publish
write.table(df, file = file.path(out_dir, "DADA2_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = '')

# Write fasta file with ASV sequences to file
write.table(data.frame(s = sprintf(">%s\n%s", df$ASV_ID, df$sequence)), file.path(out_dir, 'ASV_seqs.fasta'), col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

# Write ASV file with ASV abundances to file
df$sequence <- NULL
write.table(df, file = file.path(out_dir, "ASV_table.tsv"), sep="\t", row.names = FALSE, quote = FALSE, na = '')
