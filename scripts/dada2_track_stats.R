#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_STATS
# (modules/local/dada2_stats.nf). The port passes all input files
# as argv; the script picks the per-sample filter_stats files and the
# run-level rds files by suffix (upstream lists "./filter_and_trim_files"
# and reads rds by position).
#
# Usage: Rscript dada2_track_stats.R <prefix> <stats_out>
#     <*.filter_stats.tsv...> <_1.dada.rds> <_2.dada.rds>
#     <*.mergers.rds> <*.seqtab.rds>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
prefix <- argv[1]
stats_out <- argv[2]
files <- argv[-(1:2)]

# combine filter_and_trim files
filter_stats_files <- sort(grep("\\.filter_stats\\.tsv$", files, value = TRUE))
for (data in filter_stats_files){
    if (!exists("filter_and_trim")){ filter_and_trim <- read.csv(data, header=TRUE, sep="\t") }
    if (exists("filter_and_trim")){
        tempory <-read.csv(data, header=TRUE, sep="\t")
        filter_and_trim <-unique(rbind(filter_and_trim, tempory))
        rm(tempory)
    }
}
rownames(filter_and_trim) <- filter_and_trim$ID
filter_and_trim["ID"] <- NULL

# read data
dadaFs = readRDS(grep("_1\\.dada\\.rds$", files, value = TRUE))
dadaRs = readRDS(grep("_2\\.dada\\.rds$", files, value = TRUE))
mergers = readRDS(grep("\\.mergers\\.rds$", files, value = TRUE))
nochim = readRDS(grep("\\.seqtab\\.rds$", files, value = TRUE))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
get_acc <- function(x) sum(x$abundance[x$accept])
if ( nrow(filter_and_trim) == 1 ) {
    track <- cbind(filter_and_trim, getN(dadaFs), getN(dadaRs), getN(mergers), get_acc(mergers), rowSums(nochim))
} else {
    dadaFs_getN <- data.frame( sapply(dadaFs, getN) )
    dadaRs_getN <- data.frame( sapply(dadaRs, getN) )
    mergers_getN <- data.frame( sapply(mergers, getN) )
    mergers_acc <- data.frame( sapply(mergers, get_acc) )
    nochim_rowSums <- data.frame( rowSums(nochim) )
    track <- cbind(
        filter_and_trim[order(rownames(filter_and_trim)), ],
        dadaFs_getN[order(rownames(dadaFs_getN)), ],
        dadaRs_getN[order(rownames(dadaRs_getN)), ],
        mergers_getN[order(rownames(mergers_getN)), ],
        mergers_acc[order(rownames(mergers_acc)), ],
        nochim_rowSums[order(rownames(nochim_rowSums)), ] )
}
colnames(track) <- c("DADA2_input", "filtered", "denoisedF", "denoisedR", "denoisedPairs", "merged", "nonchim")
rownames(track) <- sub(pattern = "_1.fastq.gz$", replacement = "", rownames(track)) #this is when cutadapt is skipped!
track <- cbind(sample = sub(pattern = "(.*?)\\..*$", replacement = "\\1", rownames(track)), track)
write.table( track, file = stats_out, sep = "\t", row.names = FALSE, quote = FALSE, na = '')
