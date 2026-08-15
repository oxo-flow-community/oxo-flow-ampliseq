# oxo-flow-ampliseq — Amplicon sequencing (16S/ITS): DADA2 denoising, taxonomy assignment and QC

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
| DADA2_MERGE | `dada2_merge` | dada2 1.26.0 / digest 0.6.27 | single-run default path only; multi-run `mergeSequenceTables` branch not ported |
| MERGE_STATS_STD | `merge_stats` | r-base 4.0.3 | identical merge by `sample` |
| DB download (launcher) | `download_taxonomy_db` | curl | upstream downloads the reference DB in the Nextflow launcher (`file(url)`); the port makes it an explicit system-backend rule |
| FORMAT_TAXONOMY | `format_taxonomy` | biocontainers 1.2.0 | verbatim `bin/taxref_reformat_sbdi-gtdb.sh`; runs in a scratch dir (the script globs `*`) |
| DADA2_TAXONOMY + DADA2_ADDSPECIES + collectFile | `dada2_taxonomy` | dada2 1.26.0 | **merged**: upstream splits `ASV_seqs.fasta` into 10000-sequence chunks (`splitFasta by: 10000`) and runs assignTaxonomy + addSpecies per chunk, then concatenates chunk tables with header + sorted rows (`collectFile keepHeader, skip 1, sort`). The port replicates chunking with `awk` + per-chunk `Rscript` calls + `head`/`tail -n +2 | sort` concatenation — same chunk files, same args, same outputs. addSpecies resource hint (1 cpu/50G) becomes rule-level 10 cpus/20G, 24h limit |
| QIIME2_INASV | `qiime2_inasv` | qiime2 2026.4 | identical: biom convert + `tools import` `BIOMV210Format` |
| QIIME2_INSEQ | `qiime2_inseq` | qiime2 2026.4 | identical `FeatureData[Sequence]` import |
| QIIME2_INTAX | `qiime2_intax` | qiime2 2026.4 | verbatim `bin/parse_dada2_taxonomy.r` (porting change: output path is argv[2]) + `HeaderlessTSVTaxonomyFormat` import |
| QIIME2_BARPLOT | `qiime2_barplot` | qiime2 2026.4 | identical `taxa barplot` + `tools export` |
| MULTIQC | `multiqc` | multiqc 1.34 | identical command in a scratch dir (`multiqc` scans cwd `.`); verbatim `assets/multiqc_config.yml` |
| — (not ported) | — | — | ITS branch (`params.its`, default false), nanopore branch (`params.nanopore`, default false), syncom controls (`params.syncom`, default false), QIIME2 analyses beyond the barplot (diversity, ANCOM, classifier — `params.qiime2` default false), multi-run merge, `versions.yml` per-module tool version files, PICRUSt and other optional reports |

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

## Test

```bash
bash test/run.sh
```

Runs `validate` + `lint` + `dry-run` (with `--samples first:1`) against
`main.oxoflow`; it must exit 0.

## License

Apache-2.0 (this workflow, see `LICENSE` and `NOTICE.md`); upstream
nf-core/ampliseq is MIT (`LICENSE.upstream`).
