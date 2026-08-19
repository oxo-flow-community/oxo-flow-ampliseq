#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_ERR
# (modules/local/dada2_err.nf) with conf/modules.config DADA2_ERR
# ext.args defaults (loessErrfun, qualityType "Auto").
# Upstream lists files via list.files("."); the port passes files
# as argv (suffix-sorted inside R).
#
# Usage: Rscript dada2_learn_errors.R <prefix> <rds_out_dir>
#     <qc_out_dir> <svg_out_dir> <args_out_dir> <log_out_dir>
#     <seed> <threads> <quality_type> <filt.fastq.gz...>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
prefix <- argv[1]
rds_out_dir <- argv[2]
qc_out_dir <- argv[3]
svg_out_dir <- argv[4]
args_out_dir <- argv[5]
log_out_dir <- argv[6]
seed <- argv[7]
threads <- as.integer(argv[8])
quality_type <- argv[9]
files <- argv[-(1:9)]

set.seed(as.integer(seed)) # Initialize random number generator for reproducibility

args_str <- paste0("nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0,qualityType = \"", quality_type, "\",errorEstimationFunction = loessErrfun")

fnFs <- sort(grep("_1.filt.fastq.gz$", files, value = TRUE))
fnRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", fnFs, fixed = TRUE)

sink(file = file.path(log_out_dir, paste0(prefix, ".err.log")))
errF_raw <- learnErrors(fnFs, nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0, qualityType = quality_type, errorEstimationFunction = loessErrfun, multithread = threads, verbose = TRUE)
# the conda bioconductor-dada2 1.26 build's learnErrors returns a
# detailed LIST (err_out/err_in/trans), CRAN releases return the
# matrix directly. The RDS holds the matrix (the cross-rule contract);
# the local package's own consumers (plotErrors) take the raw return.
errF <- if (is.list(errF_raw) && !is.null(errF_raw$err_out)) errF_raw$err_out else errF_raw
saveRDS(errF, file.path(rds_out_dir, paste0(prefix, "_1.err.rds")))
errR_raw <- learnErrors(fnRs, nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0, qualityType = quality_type, errorEstimationFunction = loessErrfun, multithread = threads, verbose = TRUE)
errR <- if (is.list(errR_raw) && !is.null(errR_raw$err_out)) errR_raw$err_out else errR_raw
saveRDS(errR, file.path(rds_out_dir, paste0(prefix, "_2.err.rds")))
sink(file = NULL)

pdf(file.path(qc_out_dir, paste0(prefix, "_1.err.pdf")))
plotErrors(errF_raw, nominalQ = TRUE)
dev.off()
svg(file.path(svg_out_dir, paste0(prefix, "_1.err.svg")))
plotErrors(errF_raw, nominalQ = TRUE)
dev.off()

pdf(file.path(qc_out_dir, paste0(prefix, "_2.err.pdf")))
plotErrors(errR_raw, nominalQ = TRUE)
dev.off()
svg(file.path(svg_out_dir, paste0(prefix, "_2.err.svg")))
plotErrors(errR_raw, nominalQ = TRUE)
dev.off()

sink(file = file.path(qc_out_dir, paste0(prefix, "_1.err.convergence.txt")))
dada2:::checkConvergence(errF_raw)
sink(file = NULL)

sink(file = file.path(qc_out_dir, paste0(prefix, "_2.err.convergence.txt")))
dada2:::checkConvergence(errR_raw)
sink(file = NULL)

write.table(paste0('learnErrors\t', args_str), file = file.path(args_out_dir, "learnErrors.args.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
