#!/usr/bin/env python3
# Verbatim from nf-core/ampliseq 2.18.0 bin/cutadapt_summary.py

# --- Import libraries, do initializations  ---#
import os
import re, sys
from sys import argv

usage = "Usage: cutadapt_summary.py <single_end/paired_end> cutadapt_log_*.txt"

# --- Check and read arguments ---#
if len(argv) < 3:
    exit(usage)
if argv[1] != "single_end" and argv[1] != "paired_end":
    exit(usage)

regexes = [
    r" -o (\S+) ",
    r"Total (?:read pairs|reads) processed:\s+([0-9,,]+)",
    r"Reverse-complemented:\s+([0-9,,]+)",
    r"(?:Pairs|Reads) written .+?:\s+([0-9,,]+)",
    r"(?:Pairs|Reads) written .+?:.*?\(([^)]+)",
]

columns = [
    "sample",
    "cutadapt_total_processed",
    "cutadapt_reverse_complemented",
    "cutadapt_passing_filters",
    "cutadapt_passing_filters_percent",
]

# --- Search each file using regex ---#
print("\t".join(columns))
for FILE in argv[2:]:
    with open(FILE) as x:
        results = []
        TEXT = x.read()
        for REGEX in regexes:
            match = re.search(REGEX, TEXT)
            if match:
                results.append(match.group(1))
            else:
                results.append("")

        # modify sample names (all before ".")
        # keep only the basename first: paths like "work/cutadapt/S1.basic_1.fastq.gz"
        # must yield the same key as merge_stats.R's sample names ("S1")
        results[0] = os.path.basename(results[0]).split(".", 1)[0]

        # output per file
        print("\t".join(results))
