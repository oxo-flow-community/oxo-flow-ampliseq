#!/usr/bin/env Rscript
# Verbatim port of SIDLE_FILTTAX's inline R (nf-core/ampliseq
# subworkflows/local/sidle/filttax.nf): filters the reconstructed
# taxonomy table to features present in the reconstructed counts table
# and writes both the filtered table and a merged (all ref rows) table.
#
# Usage: sidle_filttax.R <reconstruction_taxonomy.tsv> <reconstructed_feature-table.tsv> <out_filtered=reconstructed_taxonomy.tsv> <out_merged=reconstructed_merged.tsv>
# (arg order matches upstream SIDLE_FILTTAX(SIDLE_TAXRECON.out.tsv, SIDLE_TABLERECON.out.tsv):
#  table_tofilter = reconstruction taxonomy, table_ref = reconstructed counts table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: sidle_filttax.R <reconstruction_taxonomy.tsv> <reconstructed_feature-table.tsv> <out_filtered.tsv> <out_merged.tsv>")
}

# Reconstruction taxonomy tsv (header: "Feature ID", "Confidence", "Taxon";
# the "# Constructed from taxonomy file" comment line is dropped by R's
# default comment.char = "#").
df_tofilter <- read.table(args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(df_tofilter)[1] <- "ID"

# Reconstructed counts table (header: "#OTU ID" + sample names; the
# "# Constructed from biom file" comment line needs an explicit skip = 1
# and comment.char = "" because sample names may start with "#").
df_ref <- read.table(args[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     skip = 1, comment.char = "")
colnames(df_ref)[1] <- "ID"

df_merged <- merge(df_tofilter, df_ref, by="ID", all.x=FALSE, all.y=TRUE)
write.table(df_merged, file = args[4], row.names=FALSE, sep="\t")

df_filtered <- subset(df_tofilter, df_tofilter$ID %in% df_ref$ID)
write.table(df_filtered, file = args[3], row.names=FALSE, sep="\t")
