# Genome Assembly Quality Statistics: A Complete Guide from Scratch

---

## Part 1: Why Do We Need Assembly Statistics?

When you assemble a genome, you get a collection of DNA sequences called **contigs**
(contiguous sequences). In a perfect world, you'd get one contig per chromosome —
for *C. capitata* with 6 chromosomes, that would be 6 contigs. In reality, you get
hundreds or thousands, because the assembler can't resolve every repeat or
ambiguity in the genome.

So the first question is: **how good is my assembly?** This has three dimensions:

1. **Contiguity** — How long and continuous are the pieces? (N50, L50, etc.)
2. **Completeness** — Did I capture all the genes/sequence? (BUSCO)
3. **Accuracy** — Are the individual bases correct? (QV from YAK/Merqury)

Each measures something fundamentally different. A highly contiguous assembly
can still have base errors. A base-accurate assembly can still be missing chunks
of the genome. You need all three.

---

## Part 2: Contiguity Metrics

### 2.1 The Basic Metrics

**Total assembly size**: The sum of the lengths of all contigs. For *C. capitata*,
you'd expect ~600 Mbp. If you get 400 Mbp, a lot of the genome is missing.
If you get 900 Mbp, there might be spurious duplications.

**Number of contigs**: Fewer is generally better (ideal = number of chromosomes).
More contigs means more fragmentation.

**Largest contig**: The length of the single longest contig. Tells you the best-case
scenario of your assembly.

### 2.2 N50: The Most Important Assembly Statistic

N50 answers the question: **"If I pick a random base in the assembly, what is the
expected minimum length of the contig it falls on?"**

More precisely, N50 is defined as:

> Sort all contigs from longest to shortest. Starting from the longest,
> add up their lengths until the running sum reaches 50% of the total
> assembly size. The length of the contig that pushes you past the 50%
> mark is the N50.

#### Worked Example

Suppose your assembly has these contigs (sorted largest to smallest):

```
Contig A:  100 Mbp
Contig B:   80 Mbp
Contig C:   60 Mbp
Contig D:   40 Mbp
Contig E:   20 Mbp
─────────────────
Total:     300 Mbp
50% mark:  150 Mbp
```

Walk through from the top:
- After A: cumulative = 100 Mbp (not yet at 150)
- After B: cumulative = 180 Mbp (crossed 150!)

The contig that crossed the 50% threshold is **Contig B (80 Mbp)**.
Therefore, **N50 = 80 Mbp**.

#### What N50 Tells You

N50 captures the "typical" contig length, weighted by size. It's more
informative than the average or median because it's weighted — a few
huge contigs matter more than many tiny ones. A higher N50 means your
assembly is more contiguous.

#### A Common Misconception

N50 is NOT the median contig length. The median would be Contig C (60 Mbp)
in the example above. N50 weights by length, so it's pulled toward the
larger contigs.

### 2.3 L50: How Many Contigs to Cover Half the Assembly

L50 is the companion to N50. It answers: **"How many of the largest contigs
do you need to cover 50% of the assembly?"**

In the example above:
- You need Contigs A and B to reach 150 Mbp
- **L50 = 2**

Lower L50 is better — it means fewer, larger contigs make up the bulk
of your assembly.

### 2.4 N90 and L90

Same concept as N50/L50, but using the 90% threshold instead of 50%.

In our example:
- 90% of 300 Mbp = 270 Mbp
- A (100) + B (80) + C (60) + D (40) = 280 Mbp → crossed 270 at Contig D
- **N90 = 40 Mbp**, **L90 = 4**

N90/L90 tells you about the quality of the "tail" of your assembly.
For a chromosome-scale assembly, L90 should ideally equal the number
of chromosomes.

### 2.5 NG50: The Normalized Version

N50 has a subtle problem: it depends on the total assembly size.
If your assembler produces a 300 Mbp assembly of a 600 Mbp genome
(missing half the genome), N50 is calculated against 300 Mbp — making
it look artificially good.

**NG50** fixes this by using the **expected genome size** instead of the
assembly size as the denominator. For *C. capitata*, you'd use 600 Mbp.

This makes NG50 comparable across different assemblies of the same
genome, even if they vary in total size.

### 2.6 Summary Table

| Metric | What it measures | Better = |
|--------|-----------------|----------|
| Total size | How much genome is in the assembly | Close to expected genome size |
| # Contigs | Fragmentation | Lower |
| Largest contig | Best case contiguity | Higher |
| N50 | Weighted median contig length (50% threshold) | Higher |
| L50 | # contigs needed for 50% of assembly | Lower |
| N90 | Weighted contig length at 90% threshold | Higher |
| L90 | # contigs needed for 90% of assembly | Lower |
| NG50 | N50 normalized to expected genome size | Higher |

### 2.7 How We Computed Stats in the Script

The `awk` command in `03_assembly_stats.sh` does exactly the algorithm
described above. Here's a simplified, annotated version:

```awk
# Read a FASTA file and compute N50
awk '
/^>/ {
    # When we hit a header line (starts with >), save the previous contig length
    if (n > 0) {
        lengths[n] = len       # Store this contig's length
        total += len           # Add to running total
        if (len > max) max = len  # Track the longest
    }
    n++                        # Increment contig counter
    len = 0                    # Reset length for new contig
    next                       # Skip to next line
}
{
    len += length($0)          # Add this line's characters to current contig length
}
END {
    # Don't forget the last contig
    lengths[n] = len
    total += len
    if (len > max) max = len

    # Sort lengths descending (bubble sort — fine for small arrays)
    for (i = 1; i <= n; i++) sorted[i] = lengths[i]
    for (i = 1; i <= n; i++)
        for (j = i+1; j <= n; j++)
            if (sorted[j] > sorted[i]) {
                tmp = sorted[i]; sorted[i] = sorted[j]; sorted[j] = tmp
            }

    # Walk through sorted contigs, accumulating length
    cumsum = 0
    for (i = 1; i <= n; i++) {
        cumsum += sorted[i]
        # N50: the contig where cumulative sum crosses 50% of total
        if (n50 == 0 && cumsum >= total / 2) {
            n50 = sorted[i]    # This contig's length is the N50
            l50 = i            # We needed i contigs to get here
        }
        # N90: same but at 90%
        if (n90 == 0 && cumsum >= total * 0.9) {
            n90 = sorted[i]
            l90 = i
        }
    }

    # Print results
    printf "N50: %d bp\n", n50
    printf "L50: %d\n", l50
}' assembly.fa
```

There is no external tool needed — it's pure text processing on a FASTA file.

---

## Part 3: BUSCO — Biological Completeness

### 3.1 The Problem Contiguity Doesn't Solve

You can have an assembly with great N50 but that's missing 10% of the genes.
N50 tells you nothing about whether the *right* sequences are there — only
how long the pieces are. BUSCO fills this gap.

### 3.2 What BUSCO Stands For

**B**enchmarking **U**niversal **S**ingle-**C**opy **O**rthologs.

Every word matters:

- **Benchmarking**: It's a quality metric, not a correction tool
- **Universal**: The genes it looks for should exist in virtually ALL species of a given lineage
- **Single-Copy**: Each gene should appear exactly ONCE in the genome (not in multi-gene families)
- **Orthologs**: Genes inherited from a common ancestor (not paralogs from duplication)

### 3.3 The Core Idea

Certain genes are so essential for life that they're conserved across
entire taxonomic groups. For example, ribosomal proteins, DNA polymerase
subunits, and basic metabolic enzymes are found in virtually every insect.

BUSCO maintains curated databases of these genes for different lineages:

- **Diptera odb10**: 3,285 genes expected in all flies (what the CiFi paper used)
- **Mammalia odb10**: ~9,226 genes expected in all mammals
- **Eukaryota odb10**: 255 genes expected in all eukaryotes
- And many more

### 3.4 How BUSCO Works Step by Step

```
Your assembly (FASTA)
        │
        ▼
  ┌─────────────────────┐
  │  1. Gene Prediction  │  BUSCO runs a gene finder (AUGUSTUS or
  │     (ab initio)      │  MetaEuk) on your assembly to predict
  │                      │  where genes might be
  └──────────┬──────────┘
             │
             ▼
  ┌─────────────────────┐
  │  2. Protein Search   │  It searches for the expected conserved
  │     (HMMER/miniprot) │  proteins using profile hidden Markov
  │                      │  models from OrthoDB
  └──────────┬──────────┘
             │
             ▼
  ┌─────────────────────┐
  │  3. Classification   │  Each expected gene is classified as:
  │                      │
  │   C = Complete       │  Found, full-length, single copy
  │   D = Duplicated     │  Found, but in 2+ copies
  │   F = Fragmented     │  Found, but partial/broken
  │   M = Missing        │  Not found at all
  └─────────────────────┘
```

### 3.5 Interpreting BUSCO Results

A typical BUSCO output looks like:

```
C:96.7%[S:95.1%,D:1.6%],F:1.2%,M:2.1%,n:3285
```

Breaking this down:

| Symbol | Value | Meaning |
|--------|-------|---------|
| C | 96.7% | Complete: 3,176 of 3,285 Diptera genes found in full |
| S | 95.1% | Single-copy: found exactly once (good) |
| D | 1.6% | Duplicated: found 2+ times (may indicate false duplications) |
| F | 1.2% | Fragmented: partially present (assembly broke the gene) |
| M | 2.1% | Missing: not found at all |
| n | 3,285 | Total genes searched for |

### 3.6 What Do the Numbers Mean?

- **>95% Complete** for a well-studied lineage = excellent assembly
- **90–95%** = good, typical for complex or heterozygous genomes
- **<90%** = something may be wrong (insufficient data, high repeat content, or the lineage database is a poor fit)
- **High Duplicated (>5%)** = may indicate the assembler failed to merge haplotypes and kept both copies as separate contigs
- **High Fragmented** = contigs are breaking in the middle of genes

The CiFi paper reported **96.7–99.3%** complete BUSCO for their *C. capitata* assembly, which is excellent.

### 3.7 Important Caveats

BUSCO has real limitations you should be aware of:

1. **Lineage choice matters**: Using "eukaryota" (255 genes) vs "diptera" (3,285 genes)
   gives very different granularity. Always use the most specific lineage available.

2. **BUSCO can underestimate completeness on large genomes**: The T2T human genome
   (essentially perfect) only scores ~95.7% in BUSCO because the gene finder
   has trouble with some human gene structures. When run on the annotated proteins
   directly, it scores 99.2%.

3. **BUSCO only checks conserved genes**: If your organism has unique genes not in
   the database, BUSCO cannot assess them.

4. **Version matters**: Always report which BUSCO version and which lineage database
   you used. Results from BUSCO v3 and v5 are not directly comparable.

### 3.8 How to Run BUSCO

```bash
# Basic command
busco -i assembly.fa \
      -l diptera_odb10 \
      -o busco_output \
      -m genome \
      --cpu 16

# Flags:
#   -i  : input FASTA file (your assembly)
#   -l  : lineage dataset (diptera_odb10 for flies)
#   -o  : output directory name
#   -m  : mode (genome, transcriptome, or proteins)
#   --cpu : number of threads
```

---

## Part 4: Quality Value (QV) — Base-Level Accuracy

### 4.1 The Problem Neither N50 Nor BUSCO Solve

N50 tells you about contiguity. BUSCO tells you about gene completeness.
Neither tells you whether the individual A/T/G/C bases in your assembly
are **correct**. A contig could be 100 Mbp long, contain all expected genes,
and still have thousands of single-base errors scattered throughout.

### 4.2 What Is QV?

**QV** (Quality Value) measures the probability that any given base in your
assembly is wrong, expressed on the **Phred scale**:

```
QV = -10 × log₁₀(error_rate)
```

This is the same scale used for sequencing quality scores:

| QV  | Error Rate | Accuracy | Meaning |
|-----|-----------|----------|---------|
| Q20 | 1 in 100 | 99% | Low quality |
| Q30 | 1 in 1,000 | 99.9% | Good |
| Q40 | 1 in 10,000 | 99.99% | Very good |
| Q50 | 1 in 100,000 | 99.999% | Excellent |
| Q60 | 1 in 1,000,000 | 99.9999% | Near-perfect |

The CiFi paper reported **QV 41.4–58.1** for the *C. capitata* assembly,
meaning accuracy between 99.993% and 99.9998% — very high.

### 4.3 How K-mer-Based QV Works (the Key Insight)

Both **YAK** and **Merqury** estimate QV without a reference genome. The
logic is beautifully simple:

**Step 1: Build a k-mer database from your reads**

A "k-mer" is a substring of length k. For example, with k=4, the
sequence `ACGTACG` contains these k-mers:

```
ACGT, CGTA, GTAC, TACG
```

If you count all k-mers in your HiFi reads, you get the "truth" — what
k-mers should exist in the genome.

**Step 2: Count k-mers in your assembly**

Do the same for your assembled contigs.

**Step 3: Compare**

- If a k-mer exists in both reads and assembly → probably correct
- If a k-mer is in the assembly but NOT in the reads → probably a base error
  (because the error created a new k-mer that doesn't exist in nature)
- If a k-mer is in the reads but NOT in the assembly → missing sequence

**Step 4: Calculate error rate**

```
error_rate = (# k-mers only in assembly) / (total k-mers in assembly)
QV = -10 × log₁₀(error_rate)
```

### 4.4 Why Does a Single Base Error Affect K K-mers?

This is the key insight that makes k-mer QV work. Consider k=5:

```
Correct:    ...ACGTACGTA...
K-mers:     ACGTA CGTAC GTACG TACGT
                                     (4 k-mers passing through this region)

With error at position 4 (T→G):
Erroneous:  ...ACGGACGTA...
K-mers:     ACGGA CGGAC GGACG GACGT
                                     (4 new k-mers that don't exist in reads!)
```

A single base error corrupts up to k different k-mers. This amplification
effect is what makes the method sensitive enough to detect rare errors.

### 4.5 YAK vs Merqury

Both tools do essentially the same thing, with a key difference:

**Merqury** (by Arang Rhie, NHGRI):
- The most widely used tool
- Needs a separate k-mer counting step (using Meryl)
- Can evaluate haplotype phasing accuracy with trio data
- Produces nice visualizations (spectrum plots)

**YAK** (by Heng Li, the author of hifiasm):
- Designed specifically for HiFi data
- Has a correction for overestimated QV at very high read depth
- When read depth is very high, some rare k-mers in the genome
  might still be missing from reads by chance; YAK models this
- Used by the CiFi paper for QV assessment

Both are reference-free — you don't need a "truth" genome to compare against.

### 4.6 How to Run YAK

```bash
# Step 1: Count k-mers from HiFi reads (k=31 is standard)
yak count -b37 -t16 -o reads.yak reads.fastq.gz

# Step 2: Evaluate assembly QV
yak qv reads.yak assembly.fa

# Output looks like:
#   CC  QV  raw_QV  adjusted_QV
#   QV  hap1  52.3  48.7  52.3
```

### 4.7 How to Run Merqury

```bash
# Step 1: Build k-mer database with Meryl
meryl count k=21 output reads.meryl reads.fastq.gz

# Step 2: Run Merqury
merqury reads.meryl assembly.fa output_prefix

# Outputs include:
#   output_prefix.qv               — QV per assembly
#   output_prefix.completeness.stats — k-mer completeness
#   output_prefix.*.spectra-cn.*.png — spectrum plots
```

---

## Part 5: How These Metrics Work Together

Here's a framework for thinking about all three dimensions:

```
                        ┌──────────────┐
                        │  Your        │
                        │  Assembly    │
                        └──────┬───────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
     ┌────────────────┐ ┌──────────────┐ ┌──────────────┐
     │  CONTIGUITY    │ │ COMPLETENESS │ │  ACCURACY    │
     │  N50, L50, etc │ │ BUSCO        │ │  QV (YAK)    │
     │                │ │              │ │              │
     │  "How long     │ │ "Are all the │ │ "Are the     │
     │   are the      │ │  genes       │ │  bases       │
     │   pieces?"     │ │  there?"     │ │  correct?"   │
     └────────────────┘ └──────────────┘ └──────────────┘
     
     Can be great        Can be great      Can be great
     even with           even with         even with
     base errors         fragmentation     missing genes
```

A real-world example of why you need all three:

| Scenario | N50 | BUSCO | QV | Problem |
|----------|-----|-------|-----|---------|
| Great assembly | High | >96% | >Q40 | None! |
| Fragmented | Low | >96% | >Q40 | Contigs broke at repeats, but genes are intact |
| Missing regions | High | <85% | >Q40 | Large chunks of genome absent |
| Base errors | High | >96% | <Q30 | Contigs are long but full of mismatches |
| Haplotype merging | High | High D% | >Q40 | Both alleles kept as separate contigs |

---

## Part 6: A Visual Walkthrough with Real Numbers

Let's trace through what you'll see when your *C. capitata* assembly completes.

### Expected Result (from the CiFi paper)

The fruit fly genome has 2n = 12 (6 chromosome pairs: 5 autosomes + X/Y).
Each haplotype should be ~600 Mbp.

**Ideal haplotype assembly:**

```
Chromosome 1:  ~140 Mbp   ████████████████
Chromosome 2:  ~120 Mbp   █████████████
Chromosome 3:  ~110 Mbp   ████████████
Chromosome 4:  ~100 Mbp   ███████████
Chromosome 5:   ~80 Mbp   █████████
Chromosome X:   ~50 Mbp   ██████
                ─────────
Total:         ~600 Mbp
N50:           ~110 Mbp (Chromosome 3 pushes past 50%)
L50:           3
```

**What you'll actually get:**

Instead of 6 perfect chromosomes, you'll get dozens to hundreds of contigs.
Some chromosomes will be in one piece, others fragmented at complex repeats.
The CiFi-phased assembly (Mode 2) should be significantly more contiguous
than the HiFi-only assembly (Mode 1), because the CiFi data provides
long-range information that connects contigs across gaps.

---

## Part 7: Other Tools Worth Knowing

### QUAST

Computes contiguity statistics (like our awk script does) but also performs
reference-based evaluation if you provide a reference genome. Generates
nice HTML reports. Available at: https://github.com/ablab/quast

### Bandage

A GUI tool for visualizing assembly graphs (GFA files). You can see how
contigs connect, identify bubbles (haplotype variants), and spot potential
misassemblies. Available at: https://rrwick.github.io/Bandage/

### Compleasm

A recent reimplementation of BUSCO that uses miniprot for protein-to-genome
alignment instead of AUGUSTUS/MetaEuk. It's faster and can give higher
completeness scores on large genomes where BUSCO's gene finder struggles.

---

## Part 8: Quick Reference — Putting It All Together

When you report an assembly in a paper or to your professor, include at minimum:

```
Assembly: C. capitata (haplotype 1, CiFi-phased)
Total size:   XXX Mbp
# Contigs:    XXX
N50:          XXX Mbp
L50:          XXX
BUSCO:        C:XX.X%[S:XX.X%,D:X.X%],F:X.X%,M:X.X%,n:3285 (diptera_odb10)
QV:           XX.X (YAK)
Assembler:    hifiasm v0.24.0 with --dual-scaf
Input data:   PacBio HiFi WGS + CiFi HindIII (same individual)
```

This gives a complete picture: the pieces are long (N50), the genes are
there (BUSCO), and the bases are correct (QV).
