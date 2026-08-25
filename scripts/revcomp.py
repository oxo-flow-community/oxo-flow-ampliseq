#!/usr/bin/env python3
"""Reverse-complement of a primer sequence (IUPAC codes).

Port of the upstream makeComplement function
(subworkflows/local/utils_nfcore_ampliseq_pipeline/main.nf:435)
used for the illumina_pe_its read-through cutadapt pass
(CUTADAPT_READTHROUGH: -a rv_primer_revcomp -A fw_primer_revcomp).

Usage: python3 revcomp.py <primer_sequence>
"""
import sys

COMP = {
    "A": "T", "T": "A", "U": "A", "G": "C", "C": "G",
    "Y": "R", "R": "Y", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
}

if len(sys.argv) != 2:
    sys.exit("Usage: revcomp.py <primer_sequence>")

print("".join(COMP.get(base, "X") for base in sys.argv[1].upper()[::-1]))
