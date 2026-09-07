# oxo-flow-ampliseq — Amplicon sequencing (16S/ITS): DADA2 denoising, taxonomy assignment, QIIME2 diversity/ANCOM, PICRUSt, SBDI export, phyloseq/TSE objects and QC

[![CI](https://github.com/oxo-flow-community/oxo-flow-ampliseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-ampliseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`nf-core/ampliseq`](https://github.com/nf-core/ampliseq) @ `2.18.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Turn raw paired-end amplicon reads (16S or ITS) into a published, quality-checked
DADA2 analysis: FastQC quality control, cutadapt primer trimming with a trimming
summary, DADA2 denoising (quality profiles, automatic truncation lengths,
filterAndTrim, error model learning, denoising, paired-end merging, chimera
removal and read tracking), taxonomy assignment against the curated SBDI-GTDB
reference (assignTaxonomy + addSpecies) — or, alternatively, phylogenetic
placement (PPLACE: clustalo → epa-ng → gappa), Kraken2, SINTAX or VSEARCH LCA
taxonomy, each feeding the same QIIME2 taxonomy slot plus per-tool
phyloseq/TreeSummarizedExperiment objects —, an overall per-sample read-tracking
summary, and a MultiQC report. A QIIME2 taxa barplot over your sample metadata
is also ported but gated behind `run_qiime2 = true` (default `false` — the
QIIME2 container is ~20GB unpacked).
Optional branches (off by default) add the SBDI Sweden biodiversity submission
export, phyloseq / TreeSummarizedExperiment R objects, and an Rmd-based HTML
summary report that aggregates everything into one document.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.17.0**. Install the prebuilt release binary
(recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively, conda users may `conda install -c bioconda oxo-flow-cli` — note
the bioconda package may lag behind the release binary; other platform binaries
are available on the [releases page](https://github.com/Traitome/oxo-flow/releases).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-ampliseq.git
cd oxo-flow-ampliseq
```

### 3. Requirements

- **Input data** — raw paired-end FASTQ reads per sample
  (`<raw_dir>/<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`; config key
  `raw_dir`, default `test/fixtures/raw` — point it at your real data
  directory), a sample groups file
  listing the sample IDs (default `test/fixtures/groups.tsv`), and a sample
  metadata TSV used by the QIIME2 taxa barplot (config key `metadata_file`,
  default `test/fixtures/metadata.tsv`). The SBDI-GTDB taxonomy reference
  database is downloaded automatically by the workflow — you do not provide it.
- **Compute** — up to **10 CPUs / 20 GB** per rule: `dada2_denoising`
  (48 h time limit) and `dada2_taxonomy` (24 h time limit). Most rules need
  1–6 CPUs / 1–6 GB (FastQC, cutadapt, DADA2 quality/stats, QIIME2 imports).
- **Tool delivery** — containers with **pinned images**: every rule declares a
  pinned Docker image (Docker or Singularity at runtime), so no conda
  environment is needed. The host additionally needs `curl` for the
  taxonomy-database download rule and network access to figshare.

## Usage

```bash
# 1. install oxo-flow (see Installation)
# 2. prepare data: <raw_dir>/<sample>_R1.fastq.gz / _R2.fastq.gz
#    and set raw_dir in main.oxoflow [config] (default test/fixtures/raw/)
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 8
# 5. run a subset
oxo-flow run main.oxoflow -t multiqc --samples first:2
```

Configuration lives in the `[config]` section of `main.oxoflow`. The primer
sequences `FW_primer` / `RV_primer` default to empty strings (no adapter
trimming); cutadapt behavior is controlled by `cutadapt_min_overlap` and
`cutadapt_max_error_rate`. DADA2 filtering/denoising knobs (`trunc_qmin`,
`min_len`, `max_ee`, `sample_inference`, `mergepairs_strategy`, …) mirror the
upstream parameters, and `run_id` names the run-level outputs. `raw_dir`
is the directory the reads are read from (`<raw_dir>/<sample>_R1/_R2.fastq.gz`;
default `test/fixtures/raw` — point it at your real data directory, e.g.
`raw_dir = "raw"`). `metadata_file`
points at the sample metadata TSV for the barplot (default
`test/fixtures/metadata.tsv`). All `skip_*` flags map 1:1 to the upstream
`params.skip_*` — all default to `false`, i.e. the full default path runs;
`skip_fastqc` requires `skip_multiqc` too, and `skip_taxonomy` /
`skip_dada_taxonomy` additionally gate the QIIME2 taxonomy import and barplot.
**`run_qiime2`** (default `false`) gates all QIIME2 rules: they run in
the `quay.io/qiime2/amplicon` container (~20GB unpacked, ~25GB free disk
needed for the pull; no conda qiime2 exists on common mirrors). Upstream
runs qiime2 always — set `run_qiime2 = true` (with the disk) to enable.

Branches ported on top of the default path (all off by default, matching
upstream's `params` defaults):

- **ITS analysis** — `illumina_pe_its = true` adds the read-through
  cutadapt pass over the reverse-complemented primers (upstream
  CUTADAPT_READTHROUGH) and disables DADA2 truncation; `cut_its` (`none`
  | `full` | `its1` | `its2`) together with `its_extractor` (`itsx` |
  `itsxrust`) and `its_partial` extract ITS regions before the length
  filter (`filter_len_itsx`) and ITS-specific taxonomy assignment.
- **QIIME2 downstream analyses** (gated by `run_qiime2` plus
  `skip_qiime_downstream`, both default `false`): alpha/beta diversity
  (tree, rarefaction, core metrics, group significance, emperor plots,
  adonis — `skip_alpha_rarefaction`, `skip_diversity_indices`,
  `diversity_rarefaction_depth`, `qiime_adonis_formula`), abundance table
  exports (`skip_abundance_tables`, `tax_agglom_min`/`max`), ANCOM /
  ANCOM-BC / ANCOM-BC2 differential abundance (`ancom`, `ancombc`,
  `ancombc2`, `ancombc_formula`, `ancombc2_formula`,
  `ancombc_effect_size`, `ancombc_significance`), and a QIIME2 taxonomy
  classifier pipeline (`qiime_ref_taxonomy` +
  `qiime_ref_taxonomy_urls`, or a pre-trained `classifier`).
- **Multi-run merge** — `merge_runs = true` makes `dada2_merge` merge
  all run-level stats/ASV tables (upstream `mergeSequenceTables`
  branch, `repeats="error"`, `orderBy="abundance"`).
- **PICRUSt2** — `picrust = true` runs functional predictions from the
  DADA2 ASV table (see fidelity table for the source-branch deviation).
- **SBDI export** — `sbdiexport = true` writes the SBDI Sweden biodiversity
  submission tables (`results/SBDI/`: event, dna, emof, asv-table, plus a
  re-annotated taxonomy annotation table) from the DADA2 ASV and taxonomy
  tables. Upstream also defaults this to off.
- **phyloseq / TreeSummarizedExperiment objects** —
  `skip_phyloseq` / `skip_tse` (default `true` in the port; upstream builds
  them by default) write `results/phyloseq/dada2_phyloseq.rds` and
  `results/treesummarizedexperiment/dada2_TreeSummarizedExperiment.rds`
  from the ASV, taxonomy and metadata tables. Set them to `false` to enable.
- **Summary report** — `skip_report` (default `true` in the port; upstream
  renders it by default) renders `results/summary_report/summary_report.html`
  from `assets/report_template.Rmd`: QC plots, DADA2 filter/error stats,
  taxonomy reference info, the SBDI/phyloseq/TSE artifacts (when present) and
  the ITS-cut summary are passed as rmarkdown parameters. Set it to `false`
  to enable.
- **Alternative classification chains** — instead of the DADA2/SBDI-GTDB
  taxonomy (`skip_dada_taxonomy = true`), exactly one of these can take
  over (priority order mirrors upstream: pplace-pair > DADA2 > SINTAX >
  Kraken2 > native QIIME2 > VSEARCH LCA):
  - **PPLACE phylogenetic placement** — set the `pplace_tree` /
    `pplace_aln` / `pplace_taxonomy` triple (upstream `fasta` + `tree` +
    `tree_taxonomy`). Reads are aligned to the reference with clustalo,
    placed with epa-ng (jplace) and re-grafted with gappa; taxonomy comes
    from `gappa assign` over the placement taxonomy (`pplace_format_tax` /
    `pplace_reformat_tax` normalize the output). With `run_qiime2 = true`
    the placed taxonomy lands in the QIIME2 taxonomy slot
    (`qiime2_intax_pplace`) and the grafted tree replaces the de-novo
    diversity tree (`qiime2_intree_pplace`); the DADA2 intax /
    diversity-tree / classify paths are gated off. `skip_phyloseq` /
    `skip_tse` additionally produce `phyloseq_pplace.rds` /
    `tse_pplace.rds`.
  - **Kraken2** — point `kraken2_ref_tax_custom` at a prepared Kraken2
    database directory (pre-built, as upstream's custom-DB path);
    `kraken2_db_prep` stages it, `kraken2_classify` + 
    `kraken2_format_taxresults` produce the taxonomy table that
    `qiime2_intax_kraken2` imports into the QIIME2 taxonomy slot, and
    per-tool `phyloseq_kraken2` / `tse_kraken2` objects follow.
  - **SINTAX** — point `sintax_ref_tax_custom` at a reference FASTA;
    `sintax_format_db` reformats it (vbmoth format), `sintax_classify`
    runs vsearch `--sintax` and `sintax_format_tax` normalizes the output
    for `qiime2_intax_sintax` / `phyloseq_sintax` / `tse_sintax`.
  - **VSEARCH LCA** — point `vsearch_lca_ref_tax_custom` at a reference
    FASTA; `vsearch_lca_format_db` builds the reference DB,
    `vsearch_lca_classify` runs `--usearch_global --lca` and
    `vsearch_lca_format_tax` normalizes the output for
    `qiime2_intax_vsearch_lca` / `phyloseq_vsearch_lca` /
    `tse_vsearch_lca`. This is the only classification branch that can
    combine with the native QIIME2 classifier path (upstream priority:
    native wins for the QIIME2 slot; the vsearch objects are still
    produced).
  All branches respect `skip_taxonomy` (kills the whole taxonomy layer) and
  feed the same downstream consumers (barplot, summary report) as the
  DADA2 taxonomy.

## Source

Ported from **[nf-core/ampliseq](https://github.com/nf-core/ampliseq)**,
version `2.18.0` (commit `2723d4c298d48321594920d0324697e14d73ee94`, MIT
license, see `LICENSE.upstream`). Created 2026-08-15; this workflow may lag
behind upstream releases — check the commit above and the fidelity table below
for the exact ported state. Attribution details are in `NOTICE.md`.

## Fidelity

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| RENAME_RAW_DATA_FILES | `rename_raw_data_files` | nf-core/ubuntu 20.04 | identical command (soft links) |
| FASTQC | `fastqc` | fastqc 0.12.1 | identical command; upstream publishes only `*.html` — the port also declares the `*.zip` files because MultiQC consumes them |
| CUTADAPT_BASIC | `cutadapt` | cutadapt 5.2 | identical `ext.args` (`-O 3 -e 0.1 -g/-G --discard-untrimmed`) |
| CUTADAPT_SUMMARY | `cutadapt_summary` | python 3.8.3 | verbatim `bin/cutadapt_summary.py`, `paired_end` mode |
| CUTADAPT_SUMMARY_MERGE | `cutadapt_summary_merge` | — | copy action |
| DADA2_QUALITY1 / DADA2_QUALITY2 | `dada2_quality_fw`, `dada2_quality_rv`, `dada2_quality_fw_preprocessed`, `dada2_quality_rv_preprocessed` | r-base 4.0.3 / dada2 1.26.0 | upstream runs the same process twice per stage with different prefixes; the port splits each into one rule per prefix |
| TRUNCLEN | `trunclen_fw`, `trunclen_rv` | pandas 1.1.5 | verbatim `bin/trunclen.py` |
| DADA2_FILTNTRIM | `dada2_filtntrim` | dada2 1.26.0 | identical `filterAndTrim` args; args file renamed `{sample}.filterAndTrim.args.txt` (all oxo-flow rules share one workdir, upstream name would collide); `ID` column uses the file basename so read-tracking sample names match upstream |
| DADA2_ERR | `dada2_err` | dada2 1.26.0 | identical `learnErrors` args + `checkConvergence`/`plotErrors` outputs |
| DADA2_DENOISING | `dada2_denoising` | dada2 1.26.0 | identical `dada` (incl. `getDadaOpt` defaults) + `mergePairs` + `makeSequenceTable`; `params.sample_inference` wired to `pool =` (upstream only records it in the args file); retries=3 (upstream `error_retry`), 48h limit (`process_long`) |
| DADA2_RMCHIMERA | `dada2_rmchimera` | dada2 1.26.0 | identical `removeBimeraDenovo` args |
| DADA2_STATS | `dada2_stats` | dada2 1.26.0 | identical read-tracking table |
| DADA2_MERGE | `dada2_merge` | dada2 1.26.0 / digest 0.6.27 | both branches: single-run path, and `merge_runs = true` merges all run-level stats/ASV tables (unique rbind + `mergeSequenceTables` `repeats="error", orderBy="abundance", tryRC=FALSE`); the port passes the run-level files as argv (upstream globs `*.stats.tsv`/`*.ASVtable.rds` in its cwd) — see `scripts/dada2_merge.R` |
| ITSX_CUTASV | `itsx_cutasv` | itsx 1.1.3 | identical `ITSx` call (`--save_regions` from `cut_its`, `--partial` from `its_partial`); the config-dependent outfile is copied to the fixed `intermediates/itsx/ASV_ITS_seqs.fasta` so downstream rules have static inputs |
| ITSXRUST_CUTASV | `itsxrust_cutasv` | itsxrust 0.2.2 | identical `itsxrust extract` (HMM from `share/itsxrust/hmm/F.hmm`), same region outputs + `sed` header cleanup |
| FILTER_LEN_ITSX | `filter_len_itsx` | biostrings 2.66.0 (R 4.2 build; upstream 2.58.0 on R 4.0.3) | verbatim `bin/filter_len.R` with the ITSX `ext.args` (min 50 / max 1000000) |
| DADA2_TAXONOMY_WF (ITS) | `dada2_taxonomy_its` | dada2 1.26.0 | same chunked assignTaxonomy/addSpecies machinery as `dada2_taxonomy`, on the ITS-cut length-filtered fasta with the `.ASV_ITS_tax.<ref>` suffix, then `bin/add_full_sequence_to_taxfile.py` maps the taxonomy back onto the full ASV fasta — same published `ASV_tax.<ref>.tsv` / `ASV_tax_species.<ref>.tsv` paths as the default branch (the two rules are mutually exclusive via their gates) |
| QIIME2_INASV / QIIME2_INSEQ (ITS) | `qiime2_inasv_its`, `qiime2_inseq_its` | qiime2 2026.4 | same imports as the default variants but over the ITS-cut length-filtered table/seqs; they share the `intermediates/qiime2/table.qza` / `rep-seqs.qza` output paths with the default rules (mutually exclusive gates) |
| METADATA_ALL / METADATA_PAIRWISE | `qiime2_metadata_categories` | r-base 4.2 | verbatim `bin/metadata_all.r` / `bin/metadata_pairwise.r` in a conda rule (the qiime2 image has no Rscript); runs only when the QIIME2 downstream analyses need categories |
| QIIME2_TREE | `qiime2_diversity_tree` | qiime2 2026.4 | identical mafft → mask → fasttree → midpoint-root chain; rooted tree exported as `tree.nwk` |
| QIIME2_ALPHARAREFACTION | `qiime2_alphararefaction` | qiime2 2026.4 | identical: max-depth = min-read-count of `results/overall_summary.tsv` (via `bin/count_table_minmax_reads.py`, same file format as upstream MERGE_STATS_STD output) capped at 75000, steps 250 or maxdepth/20, 10 iterations |
| QIIME2_DIVERSITY_CORE | `qiime2_diversity_core` | qiime2 2026.4 | identical core-metrics (sampling depth = min-read-count with `diversity_rarefaction_depth` floor, `UNIFRAC_USE_GPU=N`); the qza outputs feed the alpha/beta/betaord/adonis rules, the 4 distance matrices are exported as tsv |
| QIIME2_DIVERSITY_ALPHA / BETA / BETAORD / ADONIS | `qiime2_diversity_alpha`, `qiime2_diversity_beta`, `qiime2_diversity_betaord`, `qiime2_diversity_adonis` | qiime2 2026.4 | identical commands (beta: `--p-pairwise` per distance × metadata category; adonis: `--p-n-jobs 1`, `--p-formula` per comma-separated `qiime_adonis_formula`); upstream fans channels out, the port loops over the category/formula lists in-shell — **data-dependent outputs**: the declared file sets assume the default fixture metadata ("group") / configured formulas |
| QIIME2_EXPORT_ABSOLUTE / RELASV / RELTAX | `qiime2_export_absolute`, `qiime2_export_relasv`, `qiime2_export_reltax` | qiime2 2026.4 | identical export + biom convert + collapse/relative-frequency loops over `tax_agglom_min`..`tax_agglom_max` (default 2..6, so the declared outputs cover that range) |
| QIIME2_ANCOM / ANCOMBC / ANCOMBC2 | `qiime2_ancom`, `qiime2_ancombc`, `qiime2_ancombc2` | qiime2 2026.4 | identical per-category filtering (`--p-where "${cat}<>''"`), ASV + per-level analyses, `<2`-taxa WARNING branch, ANCOMBC `--p-prv-cut 0.1 --p-lib-cut 500 --p-alpha 0.05 --p-conserve` + da-barplot thresholds, ANCOMBC2 `--p-p-adjust-method "holm" --p-prevalence-cutoff 0.1 --p-alpha 0.05` + `bin/ancombc_volcanoplot.r`; `ancombc_formula`/`ancombc2_formula` variants run on the unfiltered table like upstream. Upstream's error-ignore WARNING files become `<...>.WARNING.txt` next to the export dirs |
| QIIME2_PREPTAX (incl. EXTRACT + TRAIN) | `qiime2_preptax` | qiime2 2026.4 | identical: downloads the `qiime_ref_taxonomy_urls` qza pair, `bin/taxref_reformat_qiime_silva138.sh`, imports, `extract-reads` with `FW_primer`/`RV_primer`, `fit-classifier-naive-bayes` → `intermediates/qiime2/classifier.qza` |
| QIIME2_TAXONOMY (classify) | `qiime2_classify` | qiime2 2026.4 | identical `classify-sklearn --p-n-jobs` + tabulate + export to `results/qiime2/taxonomy/`; a user-supplied `classifier` is copied in-shell (skips training); in classifier mode the DADA2-taxonomy import (`qiime2_intax`) is gated off and the classifier taxonomy takes over the same `intermediates/qiime2/taxonomy.qza` path. The gate encodes the upstream dispatch priority: classify is skipped when the pplace-pair, SINTAX or Kraken2 branch is active (they rank above native QIIME2) but not by the VSEARCH LCA branch (ranks below — native wins the shared `taxonomy.qza` slot via the `qiime2_intax_vsearch_lca` gate while the vsearch per-tool objects are still produced) |
| PICRUST | `picrust` | picrust2 2.6.3 | identical `picrust2_pipeline.py -t epa-ng --remove_intermediate --in_traits EC,KO` + `add_descriptions.py` ×3 (EC/KO/METACYC); the upstream source-message file (filename == message text) is written as `picrust_message.txt`; resource hint process_high + process_medium_memory = 10 cpus / 50G |
| FASTA_NEWICK_EPANG_GAPPA as PPLACE_STANDARD (clustalo branch) | `pplace_clustalo_align`, `pplace_epang_split`, `pplace_epang_place`, `pplace_gappa_graft`, `pplace_gappa_heattree` | clustalo (mulled-v2-4cefc385…, upstream-declared Wave image), epa-ng 0.3.8, gappa 0.8.0 | upstream profile-aligns query to the reference alignment with clustalo (`--profile1`), splits with `epa-ng --split fullaln.fa refphylogeny.tree`, places with `epa-ng --ref-msa --tree --query` (jplace gz), grafts with `gappa examine graft` and renders the per-branch heat tree (`examine heat-tree` + the `grep '^ *At'` colours step) — all args identical; upstream's hmmer/mafft alignmethod branches and the `pplace_sheet` variant are not ported (the clustalo/standard path is the upstream default, `pplace_alnmethod = 'clustalo'` in nextflow.config) |
| FORMAT_PPLACETAX as PPLACEFORMATTAX_STANDARD | `pplace_format_tax`, `pplace_reformat_tax` | dada2 env (R) / pandas 1.1.5 | identical per-ASV max-LWR selection (single max → fewer taxopath entries → iterative suffix reduction) producing `*.per_query_unique.tsv` + `{pplace_name}.taxonomy.tsv` (upstream FORMAT_PPLACETAX, R script extracted verbatim); `pplace_reformat_tax` runs the PHYLOSEQ_INTAX helper `reformat_tax_for_phyloseq.py` (pandas 1.1.5 container matches upstream `phyloseq_intax.nf`) to build the per-tool phyloseq/TSE input |
| KRAKEN2_TAXONOMY_WF | `kraken2_db_prep`, `kraken2_classify`, `kraken2_format_taxresults` | kraken2 (Wave container `kraken2_coreutils_pigz`, upstream-declared), r-base 4.2.1 | upstream UNTARs a `.tar.gz`/`.tgz` DB or accepts a directory (non-tar/dir inputs error); the port's `kraken2_db_prep` stages a pre-built directory and writes a `db.ready` pointer. Classify args identical (`--use-names --confidence` from `kraken2_confidence`); FORMAT_TAXRESULTS_KRAKEN2 extracted verbatim (R) with the default `taxlevels = "D,P,C,O,F,G,S"`, emitting the `.keys/.complete/.tsv/.into-qiime2.tsv` table set |
| SINTAX_TAXONOMY_WF | `sintax_format_db`, `sintax_classify`, `sintax_format_tax` | biocontainers:v1.2.0_cv1 (upstream-declared `docker.io` ref), vsearch 2.31.0 | `sintax_format_db` replicates FORMAT_TAXONOMY_SINTAX's gzip magic-byte check (`head -c2 | od -An -t u1` == "31 139" → `cp -fL` else `gzip -c`) and the `ref_taxonomy_*.txt` "dbversion label: user_supplied" bookkeeping; classify args identical (`--gzip_decompress --sintax_cutoff 0.8 --randseed {seed} --tabbedout`, ext.args from `conf/modules.config:582`); the converter script (verbatim `bin/convert_sintax_output.py`) normalizes to the per-ASV table. Upstream's `ASV_tax_ITS_tax…` cut_its name switch does not exist in the port (single output path per rule) |
| VSEARCH_LCA_TAXONOMY_WF | `vsearch_lca_format_db`, `vsearch_lca_classify`, `vsearch_lca_format_tax` | biocontainers:v1.2.0_cv1 (upstream-declared `docker.io` ref), vsearch 2.31.0 | same FORMAT_TAXONOMY gzip check as the SINTAX rule; classify args identical (`--usearch_global --id {vsearch_lca_id} --gzip_decompress --top_hits_only --output_no_hits --maxaccepts/--maxrejects/--lca_cutoff/--query_cov from config --n_mismatch --notrunclates --lcaout`, ext.args from `conf/modules.config:631`); the converter script (verbatim `bin/convert_vsearch_lca_output.py`) expands the LCA lineages over the default `vsearch_lca_taxlevels` |
| QIIME2_INTAX (per-chain variants) | `qiime2_intax_pplace`, `qiime2_intax_sintax`, `qiime2_intax_kraken2`, `qiime2_intax_vsearch_lca` | qiime2 2026.1 | identical imports into the shared `intermediates/qiime2/taxonomy.qza` slot (each branch's gate is mutually exclusive per the upstream dispatch priority); pplace/sintax/vsearch use the same `bin/parse_dada2_taxonomy.r` path as `qiime2_intax`; the kraken2 variant replicates the upstream quirk where KRAKEN2_TAXONOMY_WF's qiime2_tsv keeps its header row and QIIME2_INTAX_KRAKEN2 runs an empty script (plain `cp`) — imported as HeaderlessTSV with the header consumed as a data row, exactly as upstream carries it |
| QIIME2_TREE (PPLACE replacement) | `qiime2_intree_pplace` | qiime2 2026.1 | on the PPLACE branch the grafted `intermediates/pplace/grafted.newick` (GAPPA_GRAFT output) is imported as the rooted tree and replaces the de-novo mafft→mask→fasttree `qiime2_diversity_tree` run, mirroring upstream's `ch_trees_for_qiime` switch |
| SIDLE_TAXONOMY_WF (custom-DB path) | `sidle_db_prep`, `sidle_indb`, `sidle_indbaligned`, `sidle_dbfilt`, `sidle_in`, `sidle_trim`, `sidle_dbextract`, `sidle_align`, `sidle_dbrecon`, `sidle_tablerecon`, `sidle_taxrecon`, `sidle_filttax`, `sidle_seqrecon`, `sidle_treerecon` | pipesidle 0.1.0-beta (qiime2 2021.4 + rescript + q2-sidle + q2-fragment-insertion/sepp) | multi-region OTU picking (`sidle_wf`): stage a custom full-length reference (tax+seq[, aln][, sepp tree qza] via `sidle_ref_*_custom`, fail-fast checks in `sidle_db_prep`), import (HeaderlessTSVTaxonomyFormat), rescript cull/taxa-filter, biom-convert + import the DADA2 ASV table/seqs, per-region trim-dada2-posthoc / extract-reads / prepare-extracted-region / align-regional-kmers (scattered over `sidle_regions`), reconstruct-database/-counts/-taxonomy/-fragment-rep-seqs + SEPP fragment insertion, prefilter + filter/merge taxonomy via `scripts/sidle_prefilter_tablerecon.sh` + `scripts/sidle_filttax.R` (verbatim). `sidle_min_counts` maps to the prefilter arg (upstream `task.ext.min_counts`, default 0). Only the `sidle_ref_tax_custom` custom-db path is ported — the FORMAT_TAXONOMY_SIDLE standard-db path, the nanopore branch, syncom controls and the PPLACE_SHEET/hmmer/mafft pplace variants are not ported |
| QIIME2_INTAX/INTREE (SIDLE replacement) | `qiime2_intax_sidle`, `qiime2_intree_sidle` | pipesidle 0.1.0-beta | the reconstructed taxonomy + table are imported into the SHARED `intermediates/qiime2/taxonomy.qza`/`table.qza` slots (the SIDLE table REPLACES the DADA2 ASV table downstream, as upstream — "any ASV postprocessing is not allowed"); the SEPP `results/sidle/reconstructed_tree.nwk` replaces the de-novo `qiime2_diversity_tree` on the full SIDLE branch (tree + aln + tax all set). When no SIDLE tree is configured the de-novo tree still runs on the DADA2 rep-seqs |
| — (not ported) | — | — | nanopore branch (`params.nanopore` — absent from the 2.18.0 codebase, docs only), syncom controls (`params.syncom` — absent from the 2.18.0 codebase), the PPLACE_SHEET variant (`params.pplace_sheet`, hmmer-based) and the hmmer/mafft `pplace_alnmethod` branches — not ported; the non-clustalo placement branches are non-default upstream |
| softwareVersionsToYAML + `versions.yml` collection (`pipeline_info/nf_core_ampliseq_software_mqc_versions.yml`, mixed into MultiQC inputs) | engine-native export: `oxo-flow report --versions-yml <file> main.oxoflow` | — | oxo-flow ≥ 0.17.0 exports an nf-core-style `versions.yml` derived statically from the workflow declarations: one entry per rule with the pinned container tag or conda env file + its sha256, plus a `references:` section fed by the workflow's `[[reference_db]]` blocks (SBDI-GTDB R11-RS232-1 here). Deviation: it is a standalone CI-diff artifact, not a per-process runtime capture — upstream records each tool's runtime version at execution time and mixes the collected file into MultiQC, while the export reflects the pinned versions in the definition (resolved runtime package versions depend on the execution environment). Per-rule `versions.yml` emission inside every command is deliberately not replicated (it would change every rule's command while the default plan stays byte-identical). |
| MERGE_STATS_STD | `merge_stats` | r-base 4.2 (envs/dada2.yaml pin; upstream declares the Wave image bioconductor-dada2_r-base_r-digest_tbb, no visible pin) | identical merge by `sample` |
| DB download (launcher) | `download_taxonomy_db` | curl | upstream downloads the reference DB in the Nextflow launcher (`file(url)`); the port makes it an explicit system-backend rule |
| FORMAT_TAXONOMY | `format_taxonomy` | nf-core/ubuntu 20.04 | verbatim `bin/taxref_reformat_sbdi-gtdb.sh`; runs in a scratch dir (the script globs `*`). Upstream declares the `biocontainers:v1.2.0_cv1` container but the port runs the script in `quay.io/nf-core/ubuntu:20.04` |
| DADA2_TAXONOMY + DADA2_ADDSPECIES + collectFile | `dada2_taxonomy` | dada2 1.26.0 | **merged**: upstream splits `ASV_seqs.fasta` into 10000-sequence chunks (`splitFasta by: 10000`) and runs assignTaxonomy + addSpecies per chunk, then concatenates chunk tables with header + sorted rows (`collectFile keepHeader, skip 1, sort`). The port replicates chunking with `awk` + per-chunk `Rscript` calls + `head`/`tail -n +2 | sort` concatenation — same chunk files, same args, same outputs. addSpecies resource hint (1 cpu/50G) becomes rule-level 10 cpus/20G, 24h limit |
| QIIME2_INASV | `qiime2_inasv` | qiime2 2026.4 | identical: biom convert + `tools import` `BIOMV210Format` |
| QIIME2_INSEQ | `qiime2_inseq` | qiime2 2026.4 | identical `FeatureData[Sequence]` import |
| QIIME2_INTAX | `qiime2_intax` | qiime2 2026.4 | verbatim `bin/parse_dada2_taxonomy.r` (porting change: output path is argv[2]) + `HeaderlessTSVTaxonomyFormat` import |
| QIIME2_BARPLOT | `qiime2_barplot` | qiime2 2026.4 | identical `taxa barplot` + `tools export` |
| MULTIQC | `multiqc` | multiqc 1.34 | identical command in a scratch dir (`multiqc` scans cwd `.`); verbatim `assets/multiqc_config.yml` |
| SBDIEXPORT | `sbdiexport` | r-base + SBDI export scripts (sbdiexport 1.2.1) | identical `sbdiexport()` call (paired mode, `FW_primer`/`RV_primer`, dada2 taxmethod); writes `results/SBDI/{event,dna,emof,asv-table}.tsv` |
| SBDIEXPORTREANNOTATE | `sbdiexportreannotate` | r-base + SBDI export scripts | identical re-annotation table (`annotation.tsv`); the barrnap-prediction arg is omitted (no barrnap branch in the port — the R script treats it as `NA`, same as upstream when no predictions exist) |
| PHYLOSEQ | `phyloseq` | phyloseq 1.50.0 (upstream biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0, env pin 1.50.0) | identical inline R (`make_phyloseq` path); prefix literal `dada2`, tree arg = nonexistent `none.tree` (no de-novo tree in the default branch — `file.exists` guard skips it, as upstream when no tree is staged). Per-tool variants (`phyloseq_pplace` / `phyloseq_kraken2` / `phyloseq_sintax` / `phyloseq_vsearch_lca`) consume the alternative-chain taxonomy tables through the same inline R (upstream PHYLOSEQ_INTAX pre-formats each chain's table with `reformat_tax_for_phyloseq.py` — the port folds that normalization into `pplace_reformat_tax` / `kraken2_format_taxresults` / the sintax/vsearch converter scripts, so the variants read the already-per-ASV tables) |
| TREESUMMARIZEDEXPERIMENT | `treesummarizedexperiment` | TreeSummarizedExperiment 2.10.0 (upstream biocontainers/bioconductor-treesummarizedexperiment:2.10.0--r43hdfd78af_0, env pin 2.10.0) | identical inline R; referenceSeq slot filled from the taxonomy `sequence` column; same `none.tree` convention on the default branch. Per-tool variants (`tse_pplace` / `tse_kraken2` / `tse_sintax` / `tse_vsearch_lca`) mirror the phyloseq ones; the alt-chain variants fall back to the PPLACE grafted tree (`intermediates/pplace/grafted.newick`) in-shell when present, else `none.tree` (upstream receives the grafted tree through the channel — the in-shell fallback is the port's static-input equivalent) |
| SUMMARY_REPORT | `summary_report` | r-base 4.2 + rmarkdown | identical `rmarkdown::render` of `assets/report_template.Rmd` with the upstream params-list contract (params_list_named, all string values single-quoted); SBDI/phyloseq/TSE/ITS sections are `[ -f ]`-conditional in-shell (their artifacts are not declared inputs — see deviations); `mqc_plot`/picrust sections omitted (see deviations); on the alt-classification chains the per-tool phyloseq/TSE RDS paths are passed as the same comma-joined `phyloseq=`/`tse=` params upstream builds from its mixed channel |

Other notes:

- `params.FW_primer`/`RV_primer` default to `null` upstream, which renders
  a literal `null` adapter into the cutadapt command; the port defaults
  them to empty strings (`-g ""`/`-G ""` = no 5'/3' adapter trimming), the
  behavior the upstream docs describe.
- All skip flags map 1:1 to `params.skip_*` (default false = full path).
  `skip_fastqc` requires `skip_multiqc` too — the MultiQC rule consumes the
  FastQC zips. `skip_taxonomy`/`skip_dada_taxonomy` additionally gate the
  QIIME2 taxonomy import and barplot, mirroring upstream's empty-taxonomy
  channel handling.
- `metadata_file` is a config key (default `test/fixtures/metadata.tsv`);
  upstream takes it from the samplesheet.
- The workflow declares an `[[reference_db]]` block for the SBDI-GTDB
  taxonomy (R11-RS232-1), which the engine's versions.yml export
  (`oxo-flow report --versions-yml <file> main.oxoflow`) surfaces in its
  `references:` section; upstream records the equivalent in its
  `workflow_manifest`/summary report instead.
- Deviations: the QIIME2 rules pin the container at
  `quay.io/qiime2/amplicon:2026.1` (upstream modules use `2026.4` — the
  version the port was built and live-tested against); `skip_qiime_downstream`
  scopes to the newly ported downstream rules (diversity, exports, ANCOM,
  classifier) and does not turn off the taxonomy import/barplot;
  data-dependent outputs (metadata categories, `tax_agglom_min`/`max`,
  adonis formulas, `<2`-taxa WARNING branches) declare the file set for the
  default fixture/parameters; PICRUSt always uses the DADA2 source
  (upstream switches to the QIIME2-filtered table when `run_qiime2` +
  abundance tables + a taxonomy are available — the port documents the
  DADA2 basis in `results/picrust/picrust_message.txt`); the rarefaction
  WARNING txt files (upstream `error_ignore` emits) are not declared as
  rule outputs; the phyloseq/TSE/summary-report gates default to `true`
  (i.e. off) while upstream runs them by default — flip
  `skip_phyloseq`/`skip_tse`/`skip_report` to `false` to match upstream;
  the summary report's optional sections (SBDI, phyloseq, TSE, ITS-cut) are
  `[ -f ]`-conditional in-shell rather than declared inputs, so a cold run
  with those gates enabled may race the report (re-run once the upstream
  rules finish — the engine's staleness check re-triggers the report);
  the summary report omits the upstream `mqc_plot` and picrust sections
  (`mqc_plot` has no ported counterpart wiring; the upstream picrust
  report section references a params list that is dead in 2.18.0);
  the report's `workflow_manifest_version` is the constant `'2.18.0'`
  and its taxonomy title is "Release R11-RS232-1" (upstream hardcodes the
  typo "R10"); no barrnap branch exists in the port, so
  `sbdiexportreannotate` omits the prediction-file arg (R treats it as
  `NA`, the same path upstream takes when no predictions exist).
- Alt-classification chains (PPLACE / Kraken2 / SINTAX / VSEARCH LCA),
  added 2026-09-07: `pplace_gappa_assign` hardcodes the optional-emits
  flags upstream enables conditionally (`--per-query-results --krona
  --sativa`) and `pplace_gappa_heattree` hardcodes the tree-format flags
  upstream passes via `ext.args` (`--write-nexus-tree
  --write-phyloxml-tree --write-svg-tree`); only the standard clustalo
  placement path is implemented (upstream default, see the not-ported
  row above); the Kraken2 DB is staged from a pre-built directory via a
  `db.ready` marker instead of upstream's UNTAR channel, and the Kraken2
  chain runs in the upstream-declared community.wave.seqera.io Wave
  containers while the SINTAX/VSEARCH format rules keep upstream's
  `docker.io/biocontainers:v1.2.0_cv1` declaration; the per-tool
  phyloseq/TSE variants read already-normalized tables (the
  PHYLOSEQ_INTAX `reformat_tax_for_phyloseq.py` step is folded into the
  chain's format rule) and the TSE variants fall back to the grafted
  tree in-shell rather than through a channel; `format_pplacetax.R`
  runs in a conda env (`envs/pplacetax.yaml`) instead of upstream's
  Wave container.
- SIDLE multi-region OTU picking, added 2026-09-08: the port implements
  the `sidle_ref_tax_custom` custom-database path only (the
  standard-DB FORMAT_TAXONOMY_SIDLE path is not ported); upstream
  requires `skip_report = true` when SIDLE runs (its report template
  has no SIDLE section) — the port's `summary_report` is
  `[ -f ]`-conditional and has no SIDLE section either. Deviations:
  per-region DADA2_SPLITREGIONS output is not replicated — SIDLE
  consumes the merged `results/dada2/ASV_table.tsv`/`ASV_seqs.fasta`
  (single-sample-inference path, the same files `qiime2_inasv` reads);
  `qiime2_inseq` still imports the DADA2 rep-seqs (upstream switches
  to the reconstructed fragments for phyloseq/TSE — the shared
  `rep-seqs.qza` slot keeps the DADA2 sequences); regions are
  processed in `config.sidle_regions` declaration order (upstream
  uses `toSortedList`); scattered outputs are named `{region}`-only
  (upstream appends the trim length, e.g. `V4_253`); and
  `sidle_min_counts` feeds the prefilter script directly (upstream
  `task.ext.min_counts`, default 0). All SIDLE rules run in the
  upstream-declared `docker://nf-core/pipesidle:0.1.0-beta` image
  (qiime2 2021.4 + q2-sidle 0.1.0-beta + sepp 4.3.10).

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` + `dry-run` (with `--samples first:1`) against
`main.oxoflow`; it must exit 0.

## License

Apache-2.0 (this workflow, see `LICENSE` and `NOTICE.md`); upstream
nf-core/ampliseq is MIT (`LICENSE.upstream`).
