# oxo-flow-ampliseq — Amplicon sequencing (16S/ITS): DADA2 denoising, taxonomy assignment, QIIME2 diversity/ANCOM, PICRUSt, SBDI export, phyloseq/TSE objects and QC

[![CI](https://github.com/oxo-flow-community/oxo-flow-ampliseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-ampliseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`nf-core/ampliseq`](https://github.com/nf-core/ampliseq) @ `2.18.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

Turn raw paired-end amplicon reads (16S or ITS) into a published, quality-checked
DADA2 analysis: FastQC quality control, cutadapt primer trimming with a trimming
summary, DADA2 denoising (quality profiles, automatic truncation lengths,
filterAndTrim, error model learning, denoising, paired-end merging, chimera
removal and read tracking), taxonomy assignment against the curated SBDI-GTDB
reference (assignTaxonomy + addSpecies), a QIIME2 taxa barplot over your sample
metadata, an overall per-sample read-tracking summary, and a MultiQC report.
Optional branches (off by default) add the SBDI Sweden biodiversity submission
export, phyloseq / TreeSummarizedExperiment R objects, and an Rmd-based HTML
summary report that aggregates everything into one document.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Install the prebuilt release binary
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
  (`raw/<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`), a sample groups file
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
# 2. prepare data: raw/<sample>_R1.fastq.gz / _R2.fastq.gz (see test/fixtures/raw/)
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
upstream parameters, and `run_id` names the run-level outputs. `metadata_file`
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
| QIIME2_TAXONOMY (classify) | `qiime2_classify` | qiime2 2026.4 | identical `classify-sklearn --p-n-jobs` + tabulate + export to `results/qiime2/taxonomy/`; a user-supplied `classifier` is copied in-shell (skips training); in classifier mode the DADA2-taxonomy import (`qiime2_intax`) is gated off and the classifier taxonomy takes over the same `intermediates/qiime2/taxonomy.qza` path |
| PICRUST | `picrust` | picrust2 2.6.3 | identical `picrust2_pipeline.py -t epa-ng --remove_intermediate --in_traits EC,KO` + `add_descriptions.py` ×3 (EC/KO/METACYC); the upstream source-message file (filename == message text) is written as `picrust_message.txt`; resource hint process_high + process_medium_memory = 10 cpus / 50G |
| — (not ported) | — | — | nanopore branch (`params.nanopore` — absent from the 2.18.0 codebase, docs only), syncom controls (`params.syncom` — absent from the 2.18.0 codebase), `versions.yml` per-module tool version files (the port pins versions in the env files / container tags instead) |
| MERGE_STATS_STD | `merge_stats` | r-base 4.0.3 | identical merge by `sample` |
| DB download (launcher) | `download_taxonomy_db` | curl | upstream downloads the reference DB in the Nextflow launcher (`file(url)`); the port makes it an explicit system-backend rule |
| FORMAT_TAXONOMY | `format_taxonomy` | biocontainers 1.2.0 | verbatim `bin/taxref_reformat_sbdi-gtdb.sh`; runs in a scratch dir (the script globs `*`) |
| DADA2_TAXONOMY + DADA2_ADDSPECIES + collectFile | `dada2_taxonomy` | dada2 1.26.0 | **merged**: upstream splits `ASV_seqs.fasta` into 10000-sequence chunks (`splitFasta by: 10000`) and runs assignTaxonomy + addSpecies per chunk, then concatenates chunk tables with header + sorted rows (`collectFile keepHeader, skip 1, sort`). The port replicates chunking with `awk` + per-chunk `Rscript` calls + `head`/`tail -n +2 | sort` concatenation — same chunk files, same args, same outputs. addSpecies resource hint (1 cpu/50G) becomes rule-level 10 cpus/20G, 24h limit |
| QIIME2_INASV | `qiime2_inasv` | qiime2 2026.4 | identical: biom convert + `tools import` `BIOMV210Format` |
| QIIME2_INSEQ | `qiime2_inseq` | qiime2 2026.4 | identical `FeatureData[Sequence]` import |
| QIIME2_INTAX | `qiime2_intax` | qiime2 2026.4 | verbatim `bin/parse_dada2_taxonomy.r` (porting change: output path is argv[2]) + `HeaderlessTSVTaxonomyFormat` import |
| QIIME2_BARPLOT | `qiime2_barplot` | qiime2 2026.4 | identical `taxa barplot` + `tools export` |
| MULTIQC | `multiqc` | multiqc 1.34 | identical command in a scratch dir (`multiqc` scans cwd `.`); verbatim `assets/multiqc_config.yml` |
| SBDIEXPORT | `sbdiexport` | r-base + SBDI export scripts (sbdiexport 1.2.1) | identical `sbdiexport()` call (paired mode, `FW_primer`/`RV_primer`, dada2 taxmethod); writes `results/SBDI/{event,dna,emof,asv-table}.tsv` |
| SBDIEXPORTREANNOTATE | `sbdiexportreannotate` | r-base + SBDI export scripts | identical re-annotation table (`annotation.tsv`); the barrnap-prediction arg is omitted (no barrnap branch in the port — the R script treats it as `NA`, same as upstream when no predictions exist) |
| PHYLOSEQ | `phyloseq` | phyloseq 1.52.1 | identical inline R (`make_phyloseq` path); prefix literal `dada2`, tree arg = nonexistent `none.tree` (no phylogeny branch in the port — `file.exists` guard skips it, as upstream when no tree is staged) |
| TREESUMMARIZEDEXPERIMENT | `treesummarizedexperiment` | TreeSummarizedExperiment 2.10.1 | identical inline R; referenceSeq slot filled from the taxonomy `sequence` column; same `none.tree` convention |
| SUMMARY_REPORT | `summary_report` | r-base 4.2 + rmarkdown | identical `rmarkdown::render` of `assets/report_template.Rmd` with the upstream params-list contract (params_list_named, all string values single-quoted); SBDI/phyloseq/TSE/ITS sections are `[ -f ]`-conditional in-shell (their artifacts are not declared inputs — see deviations); `mqc_plot`/picrust sections omitted (see deviations) |

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

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` + `dry-run` (with `--samples first:1`) against
`main.oxoflow`; it must exit 0.

## License

Apache-2.0 (this workflow, see `LICENSE` and `NOTICE.md`); upstream
nf-core/ampliseq is MIT (`LICENSE.upstream`).
