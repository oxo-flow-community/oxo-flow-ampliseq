#!/usr/bin/env Rscript
# Port of nf-core/ampliseq 2.18.0
# modules/local/treesummarizedexperiment.nf inline R. The upstream script
# body is kept verbatim; the nf-core templated paths ($otu_tsv, $tax_tsv,
# $sam_tsv, $tree, $prefix) become positional argv so the oxo-flow port can
# invoke it as:
#   Rscript scripts/treesummarizedexperiment.R <prefix> <tax_tsv> <otu_tsv> <sam_tsv> <tree>
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

suppressPackageStartupMessages(library(TreeSummarizedExperiment))

# Read otu table. It must be in a SimpleList as a matrix where rows
# represent taxa and columns samples.
otu_mat  <- read.table(otu_tsv, sep="\t", header=TRUE, row.names=1)
otu_mat <- as.matrix(otu_mat)
assays <- SimpleList(counts = otu_mat)
# Read taxonomy table. Correct format for it is DataFrame.
taxonomy_table  <- read.table(tax_tsv, sep="\t", header=TRUE, row.names=1)
taxonomy_table <- DataFrame(taxonomy_table)

# Match rownames between taxonomy table and abundance matrix.
taxonomy_table <- taxonomy_table[match(rownames(otu_mat), rownames(taxonomy_table)), ]

# Create TreeSE object.
tse <- TreeSummarizedExperiment(
    assays = assays,
    rowData = taxonomy_table
)

# If taxonomy table contains sequences, move them to referenceSeq slot
if (!is.null(rowData(tse)[["sequence"]])) {
    referenceSeq(tse) <- DNAStringSet( rowData(tse)[["sequence"]] )
    rowData(tse)[["sequence"]] <- NULL
}

# If provided, we add sample metadata as DataFrame object. rownames of
# sample metadata must match with colnames of abundance matrix.
if (file.exists(sam_tsv)) {
    sample_meta  <- read.table(sam_tsv, sep="\t", header=TRUE, row.names=1)
    sample_meta <- sample_meta[match(colnames(tse), rownames(sample_meta)), ]
    sample_meta  <- DataFrame(sample_meta)
    colData(tse) <- sample_meta
}

# If provided, we add phylogeny. The rownames in abundance matrix must match
# with node labels in phylogeny.
if (file.exists(tree)) {
    phylogeny <- ape::read.tree(tree)
    rowTree(tse) <- phylogeny
}

saveRDS(tse, file = paste0(prefix, "_TreeSummarizedExperiment.rds"))
