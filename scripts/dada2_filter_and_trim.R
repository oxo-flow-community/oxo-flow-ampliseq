#!/usr/bin/env Rscript
# Transcription of nf-core/ampliseq 2.18.0 module DADA2_FILTNTRIM
# (modules/local/dada2_filtntrim.nf) with conf/modules.config
# DADA2_FILTNTRIM ext.args defaults.
# All paths/params come in as argv; the ID column uses basename()
# because upstream inputs are bare file names staged into an isolated
# work dir, while oxo-flow rules run in a shared workdir.
#
# Usage: Rscript dada2_filter_and_trim.R <fw_in> <fw_out> <rv_in> <rv_out>
#     <trunclenf> <trunclenr> <threads> <stats_out> <args_out>
#     <quality_type> <truncq> <min_len> <max_len> <max_ee>
suppressPackageStartupMessages(library(dada2))

argv <- commandArgs(trailingOnly = TRUE)
fw_in <- argv[1]
fw_out <- argv[2]
rv_in <- argv[3]
rv_out <- argv[4]
trunclenf <- as.integer(argv[5])
trunclenr <- as.integer(argv[6])
threads <- as.integer(argv[7])
stats_out <- argv[8]
args_out <- argv[9]
quality_type <- argv[10]
truncq <- as.integer(argv[11])
min_len <- as.integer(argv[12])
max_len <- argv[13]
max_ee <- as.integer(argv[14])
max_len_r <- if (max_len == "Inf") Inf else as.integer(max_len)

trunc_args <- paste0("truncLen = c(", trunclenf, ", ", trunclenr, ")")
args_str <- paste0(
    "maxN = 0, trimRight = 0, minQ = 0, rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, id.sep = \"\\s\", id.field = NULL, n = 1e+05, OMP = TRUE,qualityType = \"", quality_type, "\",truncQ = ", truncq, ",minLen = ", min_len, ",maxLen = ", max_len, ",trimLeft = 0,maxEE = c(", max_ee, ", ", max_ee, "),rm.phix = TRUE"
)

out <- filterAndTrim(fw_in, fw_out, rv_in, rv_out,
    truncLen = c(trunclenf, trunclenr),
    maxN = 0, trimRight = 0, minQ = 0, rm.lowcomplex = 0, orient.fwd = NULL, matchIDs = FALSE, id.sep = "\\s", id.field = NULL, n = 1e+05, OMP = TRUE,
    qualityType = quality_type,
    truncQ = truncq,
    minLen = min_len,
    maxLen = max_len_r,
    trimLeft = 0,
    maxEE = c(max_ee, max_ee),
    rm.phix = TRUE,
    compress = TRUE,
    multithread = threads,
    verbose = TRUE)
out <- cbind(out, ID = basename(row.names(out)))

# If no reads passed the filter, write an empty GZ file
if(out[2] == '0'){
    for(fp in c(fw_out, rv_out)){
        print(paste("Writing out an empty file:", fp))
        handle <- gzfile(fp, "w")
        write("", handle)
        close(handle)
    }
}

write.table( out, file = stats_out, sep = "\t", row.names = FALSE, quote = FALSE, na = '')
write.table(paste('filterAndTrim\t', trunc_args, ',', args_str, sep = ""), file = args_out, row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
