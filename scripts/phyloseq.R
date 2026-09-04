#!/usr/bin/env Rscript
# Port of nf-core/ampliseq 2.18.0 modules/local/phyloseq.nf inline R.
# The upstream script body is kept verbatim; the nf-core templated paths
# ($otu_tsv, $tax_tsv, $sam_tsv, $tree, $prefix) become positional argv so
# the oxo-flow port can invoke it as:
#   Rscript scripts/phyloseq.R <prefix> <tax_tsv> <otu_tsv> <sam_tsv> <tree>
# Pass a non-existent path as <tree> (e.g. "none.tree") when no tree is
# available — the file.exists() guard then skips it, exactly like upstream
# when no phylogeny is staged.
# versions.yml emission is skipped on purpose (engine versions.yml adoption
# is a separate, later change for this repo).

args <- commandArgs(trailingOnly = TRUE)
prefix   <- args[1]
tax_tsv  <- args[2]
otu_tsv  <- args[3]
sam_tsv  <- args[4]
tree     <- args[5]

suppressPackageStartupMessages(library(phyloseq))

otu_df  <- read.table(otu_tsv, sep="\t", header=TRUE, row.names=1)
tax_df  <- read.table(tax_tsv, sep="\t", header=TRUE, row.names=1)
otu_mat <- as.matrix(otu_df)
tax_mat <- as.matrix(tax_df)

OTU     <- otu_table(otu_mat, taxa_are_rows=TRUE)
TAX     <- tax_table(tax_mat)
phy_obj <- phyloseq(OTU, TAX)

if (file.exists(sam_tsv)) {
    sam_df  <- read.table(sam_tsv, sep="\t", header=TRUE, row.names=1)
    SAM     <- sample_data(sam_df)
    phy_obj <- merge_phyloseq(phy_obj, SAM)
}

if (file.exists(tree)) {
    TREE    <- read_tree(tree)
    phy_obj <- merge_phyloseq(phy_obj, TREE)
}

saveRDS(phy_obj, file = paste0(prefix, "_phyloseq.rds"))
