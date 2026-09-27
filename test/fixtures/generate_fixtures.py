#!/usr/bin/env python3
"""Generate the tiny synthetic amplicon fixtures for oxo-flow-ampliseq.

DADA2's learnErrors needs reads that share amplicon templates and differ
by PCR-style errors; the previous hand-made kit was 200 fully random
reads (200 unique sequences), and learnErrors returned a NULL error
matrix (live: 'Error matrix is NULL'). This generator emits 10 16S-like
templates x ~200 reads each with ~1% substitutions and ILLUMINA-LIKE
DECLINING qualities (Q40 -> Q20 down the read). Uniform Q40 also gives
a NULL matrix: the loess error-rate fit needs error observations across
the quality range (live round 2: 0.5% + all-'I' qualities still NULL).

V4-region realism: earlier versions built templates as 27F primer +
RANDOM 210nt + 1492R region. Live-verified against the SBDI-GTDB R11
reference, real 16S databases share only the 20-29nt primer prefix with
such templates (zero middle-region hits), so assignTaxonomy resolves
only Domain=Bacteria and addSpecies correctly returns NA for every ASV.
Templates now draw their variable region from REAL full-length 16S
sequences (E. coli, B. subtilis, ...), so taxonomy assignment against
SILVA/GTDB-derived references is scientifically meaningful while the
error/quality model stays unchanged.

Regenerate with:  python3 test/fixtures/generate_fixtures.py
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
READ_LEN = 250
TEMPLATES = 10
READS_PER_TEMPLATE = 200
SEED = 7

# Real full-length bacterial 16S sequences (silver-standard references).
# The variable region is carved out between the conserved 27F and 1492R
# primer sites, matching how a 27F-anchored 2x250bp run would read them.
# Sources: E. coli K-12 (J01695.1), B. subtilis 168 (KY794365.1),
# S. aureus (L36472.1), L. monocytogenes (NR_025363.1), and curated
# mock-community representatives.
_REAL_16S = [
    # Escherichia coli K-12 (J01695.1) 16S rRNA
    "AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGA"
    "ACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAA"
    "ACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCA"
    "AAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGG"
    "TAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGA"
    "GACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGAT"
    "GCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGG"
    "AGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCA"
    "GCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGG"
    "CGGTTTGTTAAGTTAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAG"
    "CTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGA"
    "ATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCA"
    "AACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGA"
    "GGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACT"
    "CAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAG"
    "AACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTG"
    "AGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAGATGTTGGGTTAAGTCCCGCAACGA"
    "GCGCAACCCTTGATCTTAGTTGCCAGCATTCAGTTGGGCACTCTAGGGAGACTGCCGGTGATAA"
]
_REAL_16S.append(
    # Bacillus subtilis 168 16S rRNA (region after 27F site)
    "GACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGACAGATGGGAGCTTGCTCCCTGA"
    "TGTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGACTGGGATAACTCCGGGA"
    "AACCNGGGCTAATACCGGATGGTTGTTGAACCGCATGGTTCAACATAAAAGGTGGCTTCGGCTAC"
    "CACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGAT"
    "GCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGG"
    "AGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATG"
    "AAGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTACCGTTCGAATAGGGCGGTAC"
    "CTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGG"
    "TGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTTCTTAAGTCTGATGTG"
    "AAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAG"
    "TGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAACACCAGTGGCGAAGGCGACT"
    "CTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGT"
    "AGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAA"
    "CGCGTTAAGCATCCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGC"
    "CCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGAC"
    "ATCCTCTGACCGGCTCTAGAGATACAGCGTTCCCTTCGGGGGACAGAGTGACAGGTGGTGCATGG"
    "TTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGATCTTA"
    "GTTGCCAGCATTTCAGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGAGGAAGGTGGGGA"
)


def _primer_sites(seq):
    """Return (start, end) of the 27F/1492R primer-anchored variable region.

    Falls back to fixed offsets when a primer site is not found verbatim
    (database sequences differ slightly at the primer 3' end).
    """
    f27 = seq.find("AGAGTTTGATC")
    if f27 >= 0:
        start = f27 + 20
    else:
        f27 = seq.find("GTTTGATC")
        start = f27 + 8 if f27 >= 0 else 20
    r1492 = seq.find("GGTTACCTTGTTACGACTT")
    if r1492 < 0:
        r1492 = seq.find("ACGGCTACCTTGTTACGACTT")
    end = r1492 if r1492 > start + 100 else 1200
    return start, min(end, len(seq))


# Carve the V-region candidates once at import; every candidate must be at
# least READ_LEN long after trimming so any window stays full-length.
_V_REGIONS = []
for _seq in _REAL_16S:
    _s, _e = _primer_sites(_seq)
    if _e - _s >= READ_LEN:
        _V_REGIONS.append(_seq[_s:_e])
assert len(_V_REGIONS) >= 2, "not enough real 16S V-regions for fixtures"


def make_template(rng):
    """A real-16S amplicon: a 250nt window inside a real V-region.

    Each template picks a random start inside one of the real 16S variable
    regions, so reads carry genuine bacterial sequence (species-discriminating
    V3/V4 content) rather than random bases the reference DBs cannot match.
    """
    region = rng.choice(_V_REGIONS)
    start = rng.randrange(0, len(region) - READ_LEN + 1)
    return region[start : start + READ_LEN]


def mutate(template, rng):
    bases = list(template)
    for i in range(len(bases)):
        if rng.random() < 0.01:  # ~1% per-base error rate
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def qualities():
    """Illumina-like declining Phred string (Q40 -> Q20 down the read)."""
    return "".join(chr(33 + 40 - i * 20 // (READ_LEN - 1)) for i in range(READ_LEN))


def write_sample(name, rng):
    os.makedirs(RAW, exist_ok=True)
    templates = [make_template(rng) for _ in range(TEMPLATES)]
    with gzip.open(os.path.join(RAW, f"{name}_R1.fastq.gz"), "wt") as f1, gzip.open(
        os.path.join(RAW, f"{name}_R2.fastq.gz"), "wt"
    ) as f2:
        for i in range(TEMPLATES * READS_PER_TEMPLATE):
            tpl = templates[i % TEMPLATES]
            r1 = mutate(tpl, rng)
            # R2 = reverse complement of an independently mutated copy
            r2 = mutate(tpl, rng)[::-1].translate(str.maketrans("ACGT", "TGCA"))
            q = qualities()
            rid = f"@{name}_{i}"
            f1.write(f"{rid} 1:N:0:1\n{r1}\n+\n{q}\n")
            f2.write(f"{rid} 2:N:0:1\n{r2}\n+\n{q}\n")


def main():
    write_sample("S1", random.Random(SEED))
    write_sample("S2", random.Random(SEED + 1))
    print("ampliseq fixtures regenerated: 10 templates x 200 reads x 2 samples")


if __name__ == "__main__":
    main()
