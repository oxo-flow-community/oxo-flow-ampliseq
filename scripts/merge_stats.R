#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module MERGE_STATS
# (modules/local/merge_stats.nf).
#
# Usage: Rscript merge_stats.R <cutadapt_summary.tsv>
#     <DADA2_stats.tsv> <overall_summary.tsv>
argv <- commandArgs(trailingOnly = TRUE)
x <- read.table(argv[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
y <- read.table(argv[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#merge
df <- merge(x, y, by = "sample", all = TRUE)

#write
write.table(df, file = argv[3], quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
