#!/usr/bin/env Rscript
# Verbatim from nf-core/ampliseq 2.18.0 bin/parse_dada2_taxonomy.r
# (used by the QIIME2_INTAX module), with one porting change: the
# output file can be given as an optional second argument (default
# "tax.tsv"), because all oxo-flow rules share one workdir.

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1 && length(args) != 2){
    stop("Usage: parse_dada2_taxonomy.r <ASV_tax_species.tsv> [<out.tsv>]")
}

tax_file <- args[1]
OUT <- if (length(args) == 2) args[2] else "tax.tsv"

# read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

# Join columns 2:ncol(.) - 1, the taxonomy ranks (sequence is the last)
r <- colnames(tax)[!colnames(tax) %in% c('ASV_ID', 'sequence')]
tax$taxonomy <- do.call(paste, c(tax[r], sep = ';'))

#write
print (paste("write",OUT))
write.table(tax[,c('ASV_ID', 'taxonomy')], file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
