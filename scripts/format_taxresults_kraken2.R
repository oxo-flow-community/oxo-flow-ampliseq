#!/usr/bin/env Rscript
# Extracted verbatim from nf-core/ampliseq 2.18.0
# modules/local/format_taxresults_kraken2.nf (inline R script).
# Adaptation: argv = <prefix> <report> <classifiedreads> <taxlevels>,
# where taxlevels is the caller's "${LEVELS:-D,P,C,O,F,G,S}" expansion.
# The upstream versions.yml block is dropped (no oxo-flow equivalent).

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
report.file <- args[2]
creads.file <- args[3]
taxlevels <- args[4]

## report
report <- read.table(file = report.file, header = FALSE, sep = "\t")
colnames(report) <- c("percent","fragments_total","fragments","rank","taxid","taxon")

taxonomy_standard_labels <- unlist(strsplit(taxlevels,","))

# prepare empty result data frame
headers <- c("rank","taxid","taxon","taxonomy_rank","taxonomy_standard")
headers <- c(headers, taxonomy_standard_labels)
df = data.frame(matrix(nrow = 0, ncol = length(headers)))
colnames(df) = headers

# go through each line of report and aggregate information
taxonomy_rank <- c("root")
for (i in 1:nrow(report)) {
    tax <- report$taxon[i]
    taxnospaces <- gsub("^ *", "", tax)

    # (1) extract full taxonomy with taxa ranks
    nspaces = nchar(tax)- nchar(taxnospaces)
    indent = nspaces/2 +1 #+1 is to catch also entries without any indent, e.g. "root"
    # if this is a lower taxonomic rank, add it, otherwise reset and add
    if ( indent > length(taxonomy_rank) ) {
        taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
    } else {
        taxonomy_rank <- taxonomy_rank[1:indent-1]
        taxonomy_rank <- c(taxonomy_rank,paste0(report$rank[i],"__",taxnospaces))
    }
    taxonomy_rank_string <- paste(taxonomy_rank,collapse=";")

    # (2) filter taxonomy_rank to only contain entries with taxonomy_standard_labels+"__"
    taxonomy_standard = list()
    taxonomy_ranks = gsub( "__.*", "", taxonomy_rank )
    for (x in taxonomy_standard_labels) {
        taxonomy_ranks_match <- match(x, taxonomy_ranks)
        if ( !is.na(taxonomy_ranks_match) ) {
            taxonomy_clean <- gsub( ".*__", "", taxonomy_rank[taxonomy_ranks_match] )
            taxonomy_standard <- c(taxonomy_standard,taxonomy_clean)
        } else {
            taxonomy_standard <- c(taxonomy_standard,"")
        }
    }
    taxonomy_standard_string <- paste(taxonomy_standard,collapse=";")
    names(taxonomy_standard) <- taxonomy_standard_labels

    # (3) propagate all in results data frame
    results <- c(
        rank=report$rank[i],
        taxid=report$taxid[i],
        taxon=taxnospaces,
        taxonomy_rank=taxonomy_rank_string,
        taxonomy_standard=taxonomy_standard_string,
        taxonomy_standard)
    df <- rbind(df, results)
}
df$taxon_taxid <- paste0(df$taxon," (taxid ",df$taxid,")")
write.table(df, file = paste0(prefix, ".kraken2.keys.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

# merge with reads
creads <- read.table(file = creads.file, header = FALSE, sep = "\t")
colnames(creads) <- c("classified","ASV_ID","taxonomy","ASV_length","kmer_LCA")
if ( !all( creads$taxonomy %in% df$taxon_taxid ) ) { stop(paste(creads.file,"and",report.file,"dont share all IDs, exit"), call.=FALSE) }
merged <- merge(creads, df, by.x="taxonomy", by.y="taxon_taxid", all.x=TRUE, all.y=FALSE)
write.table(merged, file = paste0(prefix, ".kraken2.complete.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

# get downstream compatible table
merged_reduced <- subset(merged, select = c("ASV_ID",taxonomy_standard_labels,"taxonomy"))
colnames(merged_reduced) <- c("ASV_ID",taxonomy_standard_labels,"lowest_match")
write.table(merged_reduced, file = paste0(prefix, ".kraken2.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

# get QIIME2 downstream compatible table
qiime <- merged[c("ASV_ID", "taxonomy_standard")]
colnames(qiime) <- c("ASV_ID", "taxonomy")
write.table(qiime, file = paste0(prefix, ".kraken2.into-qiime2.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
