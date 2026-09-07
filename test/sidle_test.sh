#!/bin/bash
# SIDLE branch-flip test (dry-run only; the SIDLE env is heavy so no
# live run in CI). Mirrors the other branch blocks in run.sh: flip the
# SIDLE config keys onto a temp copy of main.oxoflow and check which
# rules the engine schedules.
set -e
cd "$(dirname "$0")/.."

if [ -z "$OXO" ]; then
    echo "OXO environment variable is not set"
    exit 1
fi

sed -e 's/^skip_qiime = .*/skip_qiime = false/' \
    -e 's/^skip_dada_taxonomy = .*/skip_dada_taxonomy = true/' \
    -e 's/^skip_qiime_downstream = .*/skip_qiime_downstream = false/' \
    -e 's/^skip_alpha_rarefaction = .*/skip_alpha_rarefaction = false/' \
    -e 's/^skip_diversity_indices = .*/skip_diversity_indices = false/' \
    -e 's/^run_qiime2 = .*/run_qiime2 = true/' \
    -e 's|^sidle_ref_tax_custom = .*|sidle_ref_tax_custom = "test/fixtures/sidle/db_taxonomy.txt"|' \
    -e 's|^sidle_ref_seq_custom = .*|sidle_ref_seq_custom = "test/fixtures/sidle/db_sequences.fasta"|' \
    -e 's|^sidle_ref_aln_custom = .*|sidle_ref_aln_custom = "test/fixtures/sidle/db_aligned.fasta"|' \
    -e 's|^sidle_ref_tree_custom = .*|sidle_ref_tree_custom = ""|' \
    -e 's|^sidle_regions = .*|sidle_regions = "V4,V1V3"|' \
    -e 's|^sidle_trim_lengths = .*|sidle_trim_lengths = "253,150"|' \
    -e 's|^sidle_fw_primers = .*|sidle_fw_primers = "GTGCCAGCMGCCGCGGTAA,GTGYCAGCMGCCGCGGTAA"|' \
    -e 's|^sidle_rv_primers = .*|sidle_rv_primers = "ATTAGAWACCCBNGTAGTCCGGACTG,GGACTACNVGGGTWTCTAAT"|' \
    main.oxoflow > .tmp.oxoflow
trap 'rm -f .tmp.oxoflow /tmp/oxo-dryrun-sidle-*.txt' EXIT

"$OXO" dry-run .tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-sidle-$$.txt 2>&1

# Positive: the SIDLE chain is scheduled, including one scatter
# instance per region for the scattered rules.
for r in sidle_db_prep sidle_indb sidle_indbaligned sidle_dbfilt sidle_in \
         sidle_trim_V4 sidle_trim_V1V3 sidle_dbextract_V4 sidle_dbextract_V1V3 \
         sidle_align_V4 sidle_align_V1V3 sidle_dbrecon sidle_tablerecon \
         sidle_taxrecon sidle_filttax qiime2_intax_sidle; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-sidle-$$.txt || { echo "SIDLE branch: expected ${r} scheduled"; exit 1; }
done

# The optional aligned-reference rules enable seq reconstruction; the
# tree branch stays off (no sidle_ref_tree_custom in the test), so the
# de-novo diversity tree still runs.
for r in sidle_seqrecon qiime2_diversity_tree; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-sidle-$$.txt || { echo "SIDLE branch: expected ${r} scheduled"; exit 1; }
done

# Negative: the SIDLE branch replaces the DADA2-table import and the
# alternative taxonomies; the SEPP tree branch stays off.
for r in qiime2_inasv qiime2_intax qiime2_intax_pplace qiime2_intax_sintax \
         qiime2_intax_kraken2 qiime2_intax_vsearch_lca qiime2_classify \
         sidle_treerecon qiime2_intree_sidle; do
    if grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-sidle-$$.txt; then echo "SIDLE branch: ${r} should NOT be scheduled"; exit 1; fi
done

trap - EXIT
rm -f .tmp.oxoflow

echo "SIDLE PASS"
