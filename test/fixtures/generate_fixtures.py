#!/usr/bin/env python3
"""Generate the tiny synthetic amplicon fixtures for oxo-flow-ampliseq.

DADA2's learnErrors needs reads that share amplicon templates and differ
by PCR-style errors; the previous hand-made kit was 200 fully random
reads (200 unique sequences), and learnErrors returned a NULL error
matrix (live: 'Error matrix is NULL'). This generator emits 10 16S-like
templates x ~200 reads each with ~0.5% substitutions + Phred-40 quality
— the error model has real structure to learn.

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


def make_template(rng):
    """A 16S-like amplicon: conserved flanks around a variable region."""
    flank5 = "AGAGTTTGATCCTGGCTCAG"
    flank3 = "GGTTACCTTGTTACGACTT"
    v4 = "".join(rng.choice("ACGT") for _ in range(READ_LEN - len(flank5) - len(flank3)))
    return flank5 + v4 + flank3


def mutate(template, rng):
    bases = list(template)
    for i in range(len(bases)):
        if rng.random() < 0.005:  # ~0.5% per-base error rate
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


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
            q = "I" * READ_LEN
            rid = f"@{name}_{i}"
            f1.write(f"{rid} 1:N:0:1\n{r1}\n+\n{q}\n")
            f2.write(f"{rid} 2:N:0:1\n{r2}\n+\n{q}\n")


def main():
    write_sample("S1", random.Random(SEED))
    write_sample("S2", random.Random(SEED + 1))
    print("ampliseq fixtures regenerated: 10 templates x 200 reads x 2 samples")


if __name__ == "__main__":
    main()
