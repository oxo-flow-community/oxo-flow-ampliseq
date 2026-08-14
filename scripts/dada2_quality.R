#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_QUALITY
# (modules/local/dada2_quality.nf).
# Upstream lists files via list.files(".") in an isolated work dir;
# the port passes the read files as argv (sorted inside R) because
# all oxo-flow rules share one workdir.
#
# Usage: Rscript dada2_quality.R <prefix> <qc_out_dir> <svg_out_dir>
#                                <args_out_dir> <tsv_out_dir>
#                                <args_string> <read.fastq.gz> [...]
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
qc_out_dir <- args[2]
svg_out_dir <- args[3]
args_out_dir <- args[4]
tsv_out_dir <- args[5]
args_str <- args[6]
readfiles <- sort(args[-(1:6)])

# make list of number of sequences
readfiles_length <- countLines(readfiles) / 4
sum_readfiles_length <- sum(readfiles_length)

#use only the first x files when read number gets above 2147483647, read numbers above that do not fit into an INT and crash the process!
if ( sum_readfiles_length > 2147483647 ) {
    max_files = length(which(cumsum(readfiles_length) <= 2147483647 ))
    write.table(max_files, file = file.path(qc_out_dir, paste0("WARNING Only ",max_files," of ",length(readfiles)," files and ",sum(readfiles_length[1:max_files])," of ",sum_readfiles_length," reads were used for ", prefix, " plotQualityProfile.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    readfiles <- readfiles[1:max_files]
} else {
    max_files <- length(readfiles)
    write.table(max_files, file = file.path(qc_out_dir, paste0(max_files," files were used for ", prefix, " plotQualityProfile.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
}

plot <- plotQualityProfile(readfiles, n = 5e+04, aggregate = TRUE)
data <- plot$data

# collect quality scores
quality_scores <- unique( plot$data$Score )  |> sort()
write.table(paste(quality_scores, sep='\n'), file = file.path(tsv_out_dir, paste0(prefix, "_QualityScores.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')

#aggregate data for each sequencing cycle
df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
cycles <- sort(unique(data$Cycle))
for (cycle in cycles) {
    subdata <- data[data[, "Cycle"] == cycle, ]
    score <- list()
    #convert to list to calculate median
    for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata$Score[j], subdata$Count[j])))}
    temp = data.frame(Cycle=cycle, Count=sum(subdata$Count), Median=median(score), stringsAsFactors=FALSE)
    df <- rbind(df, temp)
}

#write output
write.table( t(df), file = file.path(tsv_out_dir, paste0(prefix, "_qual_stats",".tsv")), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
pdf(file.path(qc_out_dir, paste0(prefix, "_qual_stats",".pdf")))
plot
dev.off()
svg(file.path(svg_out_dir, paste0(prefix, "_qual_stats",".svg")))
plot
dev.off()

write.table(paste0('plotQualityProfile\t', args_str, '\nmax_files\t', max_files), file = file.path(args_out_dir, paste0(prefix, "_plotQualityProfile.args.txt")), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
