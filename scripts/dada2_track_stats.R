#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 3.1.0 module DADA2_STATS (2.18.0
# baseline; multi-sample track table synced to the 3.1.0 #1057 key-merge fix)
# (modules/local/dada2_stats/main.nf). The port passes all input files
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
# filter_and_trim rownames end in '_1.fastq.gz', DADA2/seqtab rownames in '_1.filt.fastq.gz'; reduce both to the sample id
# (upstream 3.1.0 #1057: the track table is built by merging on this key,
# not by positional cbind, which mislabeled rows when one sample name was
# a prefix of another)
normKey_ft <- function(x) sub(pattern = "(.*?)\\..*$", replacement = "\\1", sub(pattern = "_1.fastq.gz$", replacement = "", x))
normKey_nc <- function(x) sub(pattern = ".filt.fastq.gz$", replacement = "", sub(pattern = "_2.filt.fastq.gz$", replacement = "", sub(pattern = "_1.filt.fastq.gz$", replacement = "", x)))
if ( nrow(filter_and_trim) == 1 ) {
    track <- cbind(filter_and_trim, getN(dadaFs), getN(dadaRs), getN(mergers), get_acc(mergers), rowSums(nochim))
} else {
    # Merge on the sample key; a positional cbind mislabels counts when one sample name is a prefix of another
    ft_keys <- normKey_ft(rownames(filter_and_trim))
    samples <- normKey_nc(rownames(nochim))   # canonical sample order from the ASV table
    gN  <- function(x) as.numeric(unname(sapply(x, getN)))
    gNa <- function(x) as.numeric(unname(sapply(x, get_acc)))
    track <- data.frame(key = ft_keys, filter_and_trim, check.names = FALSE, stringsAsFactors = FALSE)
    track <- merge(track, data.frame(key = samples,
            denoisedF      = gN(dadaFs),
            denoisedR      = gN(dadaRs),
            denoisedPairs  = gN(mergers),
            merged         = gNa(mergers),
            nonchim        = as.numeric(unname(rowSums(nochim))),
            stringsAsFactors = FALSE),
        by = "key", all = TRUE, sort = FALSE)
    # keep the deterministic filter_and_trim row order
    track <- track[order(match(track$key, ft_keys)), ]
    rownames(track) <- track$key
    track$key <- NULL
}
colnames(track) <- c("DADA2_input", "filtered", "denoisedF", "denoisedR", "denoisedPairs", "merged", "nonchim")
rownames(track) <- sub(pattern = "_1.fastq.gz$", replacement = "", rownames(track)) #this is when cutadapt is skipped!
track <- cbind(sample = sub(pattern = "(.*?)\\..*$", replacement = "\\1", rownames(track)), track)
write.table( track, file = stats_out, sep = "\t", row.names = FALSE, quote = FALSE, na = '')
