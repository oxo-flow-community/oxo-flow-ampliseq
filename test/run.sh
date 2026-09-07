#!/usr/bin/env bash
# Acceptance test for oxo-flow-ampliseq port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

# Branch-flip acceptance: alt-classification chains must schedule exactly one
# taxonomy branch's rules when its config keys are set (dry-run is plan-time,
# so no tool execution happens; scheduling lines go to stderr, hence 2>&1).

echo "==> branch-flip: PPLACE phylogenetic placement"
sed -e 's/^run_qiime2 = false/run_qiime2 = true/' \
    -e 's/^skip_phyloseq = true/skip_phyloseq = false/' \
    -e 's/^skip_tse = true/skip_tse = false/' \
    -e 's/^skip_report = true/skip_report = false/' \
    -e 's/^pplace_tree = ""/pplace_tree = "test\/fixtures\/pplace\/ref.tree"/' \
    -e 's/^pplace_aln = ""/pplace_aln = "test\/fixtures\/pplace\/ref.aln.fa"/' \
    -e 's/^pplace_taxonomy = ""/pplace_taxonomy = "test\/fixtures\/pplace\/ref.tax.tsv"/' \
    main.oxoflow > .tmp.oxoflow
trap 'rm -f .tmp.oxoflow /tmp/oxo-dryrun-*-$$.txt' EXIT
"$OXO" dry-run .tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-pplace-$$.txt 2>&1
for r in pplace_clustalo_align pplace_epang_split pplace_epang_place pplace_gappa_graft \
         pplace_gappa_assign pplace_gappa_heattree pplace_format_tax pplace_reformat_tax \
         qiime2_intax_pplace qiime2_intree_pplace phyloseq_pplace tse_pplace; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-pplace-$$.txt || { echo "PPLACE branch: expected ${r} scheduled"; exit 1; }
done
# pplace-pair wins over the default DADA2 intax/diversity-tree/classify paths
for r in qiime2_intax qiime2_diversity_tree qiime2_classify; do
    if grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-pplace-$$.txt; then echo "PPLACE branch: ${r} must not run"; exit 1; fi
done
trap - EXIT
rm -f .tmp.oxoflow

echo "==> branch-flip: Kraken2 classification"
sed -e 's/^run_qiime2 = false/run_qiime2 = true/' \
    -e 's/^skip_dada_taxonomy = false/skip_dada_taxonomy = true/' \
    -e 's/^skip_phyloseq = true/skip_phyloseq = false/' \
    -e 's/^skip_tse = true/skip_tse = false/' \
    -e 's/^skip_report = true/skip_report = false/' \
    -e 's/^kraken2_ref_tax_custom = ""/kraken2_ref_tax_custom = "test\/fixtures\/kraken2\/db"/' \
    main.oxoflow > .tmp.oxoflow
trap 'rm -f .tmp.oxoflow /tmp/oxo-dryrun-*-$$.txt' EXIT
"$OXO" dry-run .tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-kraken2-$$.txt 2>&1
for r in kraken2_db_prep kraken2_classify kraken2_format_taxresults qiime2_intax_kraken2 \
         qiime2_barplot phyloseq_kraken2 tse_kraken2; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-kraken2-$$.txt || { echo "Kraken2 branch: expected ${r} scheduled"; exit 1; }
done
grep -qE "^  [0-9]+\. qiime2_intax  \[run" /tmp/oxo-dryrun-kraken2-$$.txt && { echo "Kraken2 branch: qiime2_intax must not run"; exit 1; }
trap - EXIT
rm -f .tmp.oxoflow

echo "==> branch-flip: SINTAX classification"
sed -e 's/^run_qiime2 = false/run_qiime2 = true/' \
    -e 's/^skip_dada_taxonomy = false/skip_dada_taxonomy = true/' \
    -e 's/^skip_phyloseq = true/skip_phyloseq = false/' \
    -e 's/^skip_tse = true/skip_tse = false/' \
    -e 's/^skip_report = true/skip_report = false/' \
    -e 's/^sintax_ref_tax_custom = ""/sintax_ref_tax_custom = "test\/fixtures\/sintax\/ref.fa"/' \
    main.oxoflow > .tmp.oxoflow
trap 'rm -f .tmp.oxoflow /tmp/oxo-dryrun-*-$$.txt' EXIT
"$OXO" dry-run .tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-sintax-$$.txt 2>&1
for r in sintax_format_db sintax_classify sintax_format_tax qiime2_intax_sintax \
         phyloseq_sintax tse_sintax; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-sintax-$$.txt || { echo "SINTAX branch: expected ${r} scheduled"; exit 1; }
done
trap - EXIT
rm -f .tmp.oxoflow

echo "==> branch-flip: VSEARCH LCA classification"
sed -e 's/^run_qiime2 = false/run_qiime2 = true/' \
    -e 's/^skip_dada_taxonomy = false/skip_dada_taxonomy = true/' \
    -e 's/^skip_phyloseq = true/skip_phyloseq = false/' \
    -e 's/^skip_tse = true/skip_tse = false/' \
    -e 's/^skip_report = true/skip_report = false/' \
    -e 's/^vsearch_lca_ref_tax_custom = ""/vsearch_lca_ref_tax_custom = "test\/fixtures\/vsearch_lca\/ref.fa"/' \
    main.oxoflow > .tmp.oxoflow
trap 'rm -f .tmp.oxoflow /tmp/oxo-dryrun-*-$$.txt' EXIT
"$OXO" dry-run .tmp.oxoflow --samples first:1 > /tmp/oxo-dryrun-vsearch-$$.txt 2>&1
for r in vsearch_lca_format_db vsearch_lca_classify vsearch_lca_format_tax \
         qiime2_intax_vsearch_lca phyloseq_vsearch_lca tse_vsearch_lca; do
    grep -qE "^  [0-9]+\. ${r}  \[run" /tmp/oxo-dryrun-vsearch-$$.txt || { echo "VSEARCH branch: expected ${r} scheduled"; exit 1; }
done
trap - EXIT
rm -f .tmp.oxoflow

echo "==> branch-flip: SIDLE multi-region OTU picking"
OXO="$OXO" ./test/sidle_test.sh

echo "PASS"
