#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_MERGE
# (modules/local/dada2_merge.nf). Upstream globs *.stats.tsv /
# *.ASVtable.rds from its cwd and merges them (unique-rbind stats,
# mergeSequenceTables when several tables). The port receives the
# run-level files as argv because all rules share one workdir; a
# single-run call and a multi-run (merge_runs = true) call take the
# same code path.
#
# Usage: Rscript dada2_merge.R "<stats.tsv> [stats.tsv ...]"
#     "<ASVtable.rds> [ASVtable.rds ...]" <out_dir>
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(digest))

argv <- commandArgs(trailingOnly = TRUE)
# the shell passes "f1 f2 " (trailing space from `tr '\n' ' '`), which
# strsplit turns into a trailing empty element — drop it before use
stats_files <- strsplit(argv[1], " ")[[1]]
asvtab_files <- strsplit(argv[2], " ")[[1]]
stats_files <- stats_files[nzchar(stats_files)]
asvtab_files <- asvtab_files[nzchar(asvtab_files)]
out_dir <- argv[3]

# combine stats files (upstream: unique rbind over sorted *.stats.tsv)
for (data in sort(stats_files)) {
    temp <- read.csv(data, header = TRUE, sep = "\t")
    if (!exists("stats")) {
        stats <- temp
    } else {
        stats <- unique(rbind(stats, temp))
    }
}
write.table(stats, file = file.path(out_dir, "DADA2_stats.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

# combine dada-class objects
files <- asvtab_files
if (length(files) == 1) {
    ASVtab <- readRDS(files[1])
} else {
    ASVtab <- mergeSequenceTables(tables = files, repeats = "error", orderBy = "abundance", tryRC = FALSE)
}
saveRDS(ASVtab, file.path(out_dir, "DADA2_table.rds"))

df <- t(ASVtab)
colnames(df) <- gsub("_1.filt.fastq.gz", "", colnames(df))
colnames(df) <- gsub(".filt.fastq.gz", "", colnames(df))
df <- data.frame(sequence = rownames(df), df, check.names = FALSE)
# Create an md5 sum of the sequences as ASV_ID and rearrange columns
df$ASV_ID <- sapply(df$sequence, digest, algo = "md5", serialize = FALSE)
df <- df[, c(ncol(df), 3:ncol(df) - 1, 1)]

# order ASV table by ASV_ID (md5sum of sequence), previously sorted decreasing by total counts (ASVs of same counts were randomly sorted)
df <- df[order(df$ASV_ID), ]

# file to publish
write.table(df, file = file.path(out_dir, "DADA2_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# Write fasta file with ASV sequences to file
write.table(data.frame(s = sprintf(">%s\n%s", df$ASV_ID, df$sequence)), file.path(out_dir, "ASV_seqs.fasta"), col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")

# Write ASV file with ASV abundances to file
df$sequence <- NULL
write.table(df, file = file.path(out_dir, "ASV_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
