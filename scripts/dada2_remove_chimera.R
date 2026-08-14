#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_RMCHIMERA
# (modules/local/dada2_rmchimera.nf) with conf/modules.config
# DADA2_RMCHIMERA ext.args defaults.
#
# Usage: Rscript dada2_remove_chimera.R <seqtab.rds> <out.rds>
#     <args_out> <threads> <samples_csv>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
seqtab_in <- argv[1]
rds_out <- argv[2]
args_out <- argv[3]
threads <- as.integer(argv[4])
samples <- unlist(strsplit(argv[5], ","))
no_samples <- length(samples)
first_sample <- samples[1]

seqtab = readRDS(seqtab_in)

args_str <- 'method="consensus", minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 2, minParentAbundance = 8, allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16'

#remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", minSampleFraction = 0.9, ignoreNNegatives = 1, minFoldParentOverAbundance = 2, minParentAbundance = 8, allowOneOff = FALSE, minOneOffParentDistance = 4, maxShift = 16, multithread=threads, verbose=TRUE)
if ( no_samples == 1 ) { rownames(seqtab.nochim) <- first_sample }
saveRDS(seqtab.nochim, rds_out)

write.table(paste0('removeBimeraDenovo\t', args_str), file = args_out, row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
