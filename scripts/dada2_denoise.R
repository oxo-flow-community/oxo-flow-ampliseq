#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_DENOISING
# (modules/local/dada2_denoising.nf) with conf/modules.config
# DADA2_DENOISING ext.args / ext.args2 defaults.
# The port wires the two strategy params (mergepairs_strategy ->
# justConcatenate, sample_inference -> pool) from config; the
# "consensus" merge branch is documented as not ported (README
# fidelity table) — only the default "merge" branch is transcribed.
#
# Usage: Rscript dada2_denoise.R <prefix> <rds_out_dir> <log_out_dir>
#     <args_out_dir> <errF.rds> <errR.rds> <seed> <quality_type>
#     <threads> <mergepairs_strategy> <sample_inference>
#     <filt.fastq.gz...>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
prefix <- argv[1]
rds_out_dir <- argv[2]
log_out_dir <- argv[3]
args_out_dir <- argv[4]
errF_path <- argv[5]
errR_path <- argv[6]
seed <- as.integer(argv[7])
quality_type <- argv[8]
threads <- as.integer(argv[9])
mergepairs_strategy <- argv[10]
sample_inference <- argv[11]
files <- argv[-(1:11)]

errF <- readRDS(errF_path)
errR <- readRDS(errR_path)

filtFs <- sort(grep("_1.filt.fastq.gz$", files, value = TRUE))
filtRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", filtFs, fixed = TRUE)

pool_arg <- if (sample_inference == "pseudo") "pseudo" else if (sample_inference == "pooled") TRUE else FALSE
pool_str <- if (sample_inference == "pseudo") "\"pseudo\"" else if (sample_inference == "pooled") "TRUE" else "FALSE"

args_str <- paste0("selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE,BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL,pool = ", pool_str)

# denoising
sink(file = file.path(log_out_dir, paste0(prefix, ".dada.log")))
if (quality_type == "Auto") {
    # Avoid using memory-inefficient derepFastq() if not necessary
    dadaFs <- dada(filtFs, err = errF, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE, BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL, pool = pool_arg, multithread = threads)
    dadaRs <- dada(filtRs, err = errR, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE, BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL, pool = pool_arg, multithread = threads)
} else {
    derepFs <- derepFastq(filtFs, qualityType = quality_type)
    dadaFs <- dada(derepFs, err = errF, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE, BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL, pool = pool_arg, multithread = threads)
    derepRs <- derepFastq(filtRs, qualityType = quality_type)
    dadaRs <- dada(derepRs, err = errR, selfConsist = FALSE, priors = character(0), DETECT_SINGLETONS = FALSE, GAPLESS = TRUE, GAP_PENALTY = -8, GREEDY = TRUE, KDIST_CUTOFF = 0.42, MATCH = 5, MAX_CLUST = 0, MAX_CONSIST = 10, MIN_ABUNDANCE = 1, MIN_FOLD = 1, MIN_HAMMING = 1, MISMATCH = -4, OMEGA_A = 1e-40, OMEGA_C = 1e-40, OMEGA_P = 1e-4, PSEUDO_ABUNDANCE = Inf, PSEUDO_PREVALENCE = 2, SSE = 2, USE_KMERS = TRUE, USE_QUALS = TRUE, VECTORIZED_ALIGNMENT = TRUE, BAND_SIZE = 16, HOMOPOLYMER_GAP_PENALTY = NULL, pool = pool_arg, multithread = threads)
}
saveRDS(dadaFs, file.path(rds_out_dir, paste0(prefix, "_1.dada.rds")))
saveRDS(dadaRs, file.path(rds_out_dir, paste0(prefix, "_2.dada.rds")))
sink(file = NULL)

# merge (default strategy "merge"; the "consensus" branch is not ported)
just_concat <- (mergepairs_strategy == "concatenate")
args2_str <- paste0("homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE,justConcatenate = ", if (just_concat) "TRUE" else "FALSE", ", match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, homo_gap = NULL, endsfree = TRUE, vec = FALSE, propagateCol = character(0), trimOverhang = FALSE, returnRejects = TRUE, justConcatenate = just_concat, match = 1, mismatch = -64, gap = -64, minOverlap = 12, maxMismatch = 0, verbose = TRUE)

saveRDS(mergers, file.path(rds_out_dir, paste0(prefix, ".mergers.rds")))

# make table
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(rds_out_dir, paste0(prefix, ".seqtab.rds")))

write.table(paste0('dada\t', args_str), file = file.path(args_out_dir, "dada.args.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
write.table(paste0('mergePairs\t', args2_str), file = file.path(args_out_dir, "mergePairs.args.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
