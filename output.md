# Hifiasm Output Files Explained

## Source of this information

The file naming conventions come from the **hifiasm GitHub README**
(https://github.com/chhylp123/hifiasm), not the CiFi paper. The README
documents what each output suffix means. The CiFi paper only reports
final assembly statistics — it doesn't explain hifiasm's output files.

---

## The Assembly Pipeline (what hifiasm does internally)

Hifiasm builds the assembly in stages, and each stage produces its own output:

```
Raw HiFi reads
      │
      ▼
  ┌──────────┐
  │ r_utg    │  "Raw unitigs" — initial overlap graph, unprocessed.
  │ (1.9G)   │  Contains ALL sequences including both haplotypes
  └────┬─────┘  tangled together + duplications.
       │
       ▼
  ┌──────────┐
  │ p_utg    │  "Processed unitigs" — graph cleaned up, bubbles
  │ (1.4G)   │  (haplotype variants) partially resolved.
  └────┬─────┘  Still contains some redundancy.
       │
       ▼
  ┌──────────┐
  │ p_ctg    │  "Primary contigs" — the consensus assembly.
  │ (580M)   │  One best path through the graph.
  └────┬─────┘  Haplotypes merged into one sequence.
       │
       ├──────────────┐
       ▼              ▼
  ┌──────────┐   ┌──────────┐
  │ hap1     │   │ hap2     │  Haplotype-separated contigs.
  │ p_ctg    │   │ p_ctg    │  Mom vs Dad (or arbitrary hap1/hap2).
  │ (551M)   │   │ (456M)   │
  └──────────┘   └──────────┘
```

---

## File Types (the three suffixes)

Every assembly stage produces three versions of the same file:

| Suffix | Contains | Use case |
|--------|----------|----------|
| `.gfa` | Full GFA with sequences | The actual assembly — convert to FASTA for downstream tools |
| `.noseq.gfa` | GFA graph structure only (no sequences) | For visualization in Bandage (small enough to download) |
| `.lowQ.bed` | Genomic coordinates of low-confidence regions | QC — regions where hifiasm is uncertain about the bases |

---

## HiFi-only Assembly Files

These files have the `.bp.` infix (for "bipartition" — hifiasm's
internal haplotype separation without external phasing data).

### Assembly stages (GFA files)

| File | Size | What it is |
|------|------|------------|
| `ccap_hifi_only.bp.r_utg.gfa` | 1.9G | **Raw unitigs** — the earliest assembly stage. All sequences from both haplotypes, not yet separated. Largest because nothing has been collapsed yet. |
| `ccap_hifi_only.bp.p_utg.gfa` | 1.4G | **Processed unitigs** — cleaned up graph with bubbles (haplotype differences) partially resolved. Smaller because redundancy is being removed. |
| `ccap_hifi_only.bp.p_ctg.gfa` | 581M | **Primary contigs** — the main assembly. Both haplotypes merged into a single consensus. This is what you'd use if you don't care about phasing. ~600 Mbp expected. |
| `ccap_hifi_only.bp.hap1.p_ctg.gfa` | 551M | **Haplotype 1 contigs** — one parent's version of the genome. Without Hi-C/CiFi data, this phasing is partial and approximate. |
| `ccap_hifi_only.bp.hap2.p_ctg.gfa` | 456M | **Haplotype 2 contigs** — the other parent. Note it's smaller (456M vs 551M) — without CiFi data, hifiasm can't fully separate the haplotypes, so one ends up with more sequence than the other. |

### Graph-only files (for visualization)

| File | Size | What it is |
|------|------|------------|
| `ccap_hifi_only.bp.r_utg.noseq.gfa` | 26M | Raw unitig graph structure (no sequences) |
| `ccap_hifi_only.bp.p_utg.noseq.gfa` | 20M | Processed unitig graph structure |
| `ccap_hifi_only.bp.p_ctg.noseq.gfa` | 9.6M | Primary contig graph structure |
| `ccap_hifi_only.bp.hap1.p_ctg.noseq.gfa` | 9.0M | Haplotype 1 graph structure |
| `ccap_hifi_only.bp.hap2.p_ctg.noseq.gfa` | 7.3M | Haplotype 2 graph structure |

### Low-quality region files

| File | Size | What it is |
|------|------|------------|
| `ccap_hifi_only.bp.r_utg.lowQ.bed` | 5.8M | Low-confidence regions in raw unitigs |
| `ccap_hifi_only.bp.p_utg.lowQ.bed` | 4.4M | Low-confidence regions in processed unitigs |
| `ccap_hifi_only.bp.p_ctg.lowQ.bed` | 1.2M | Low-confidence regions in primary contigs |
| `ccap_hifi_only.bp.hap1.p_ctg.lowQ.bed` | 1.3M | Low-confidence regions in haplotype 1 |
| `ccap_hifi_only.bp.hap2.p_ctg.lowQ.bed` | 1.3M | Low-confidence regions in haplotype 2 |

### Binary index files (reusable caches)

| File | Size | What it is |
|------|------|------------|
| `ccap_hifi_only.ec.bin` | 8.9G | **Error-corrected reads** — binary index of all reads after hifiasm's internal error correction. Reusable: if you re-run hifiasm with different parameters, it skips the error correction step. |
| `ccap_hifi_only.ovlp.source.bin` | 11G | **Overlap index (source)** — all read-to-read overlaps computed during assembly. The most expensive computation. Also reusable on re-runs. |
| `ccap_hifi_only.ovlp.reverse.bin` | 2.3G | **Overlap index (reverse)** — reverse complement overlap index. Companion to source.bin. |

### Log file

| File | Size | What it is |
|------|------|------------|
| `hifiasm_hifi_only.log` | 65K | **Runtime log** — messages, statistics, warnings printed during the run. Check this if something looks wrong. |

---

## CiFi-phased Assembly Files

These files have the `.hic.` infix (hifiasm's naming convention when
Hi-C or CiFi data is provided via `--h1`/`--h2`).

### Assembly stages (GFA files)

| File | Size | What it is |
|------|------|------------|
| `ccap_cifi_phased.hic.r_utg.gfa` | 1.9G | **Raw unitigs** — same as HiFi-only (CiFi data only affects the phasing step, not the initial graph construction). |
| `ccap_cifi_phased.hic.p_utg.gfa` | 1.4G | **Processed unitigs** — same as HiFi-only for the same reason. |
| `ccap_cifi_phased.hic.p_ctg.gfa` | 580M | **Primary contigs** — consensus assembly. |
| `ccap_cifi_phased.hic.hap1.p_ctg.gfa` | 587M | **Haplotype 1 (CiFi-phased)** — THIS IS THE MAIN RESULT. Note it's larger than the HiFi-only hap1 (587M vs 551M) because CiFi data helped assign more contigs to this haplotype correctly. |
| `ccap_cifi_phased.hic.hap2.p_ctg.gfa` | 437M | **Haplotype 2 (CiFi-phased)** — the other haplotype. Together with hap1, these should total ~1,000–1,200 Mbp (the full diploid genome). |

### CiFi-specific binary files

| File | Size | What it is |
|------|------|------------|
| `ccap_cifi_phased.hic.lk.bin` | 167M | **CiFi link binary** — stores which CiFi read pairs link which contigs together. This IS the phasing information — which contigs belong to the same chromosome and haplotype. |
| `ccap_cifi_phased.hic.tlb.bin` | 17G | **CiFi table binary** — full CiFi read alignment table. Large because it stores every CiFi read's mapping position against every contig. |

### Other files

The remaining files (`.ec.bin`, `.ovlp.*.bin`, `.noseq.gfa`, `.lowQ.bed`,
`.log`) serve the same purposes as in the HiFi-only assembly.

---

## Which Files Actually Matter?

For your results and comparison with the paper, you only need **six GFA files**:

```
HiFi-only:
  ★ ccap_hifi_only.bp.p_ctg.gfa      → primary assembly (both haplotypes merged)
  ★ ccap_hifi_only.bp.hap1.p_ctg.gfa → haplotype 1 (partial phasing)
  ★ ccap_hifi_only.bp.hap2.p_ctg.gfa → haplotype 2 (partial phasing)

CiFi-phased:
  ★ ccap_cifi_phased.hic.p_ctg.gfa      → primary assembly
  ★ ccap_cifi_phased.hic.hap1.p_ctg.gfa → haplotype 1 (CiFi-phased) ← MAIN RESULT
  ★ ccap_cifi_phased.hic.hap2.p_ctg.gfa → haplotype 2 (CiFi-phased) ← MAIN RESULT
```

The `.bin` files are reusable caches — you can delete them to save ~50G
of disk space once you're happy with the results.

The `.noseq.gfa` files are only needed if you want to visualize the graph
in Bandage.

The `.lowQ.bed` files are for advanced QC — identifying uncertain regions.

---

# CiFi2PE Output Files Explained

## Source of this information

The file descriptions come from **reading the actual CiFi2PE source code**
(`cifi2pe_full_length_args.py` in https://github.com/sheinasim-USDA/CiFi2PE).
I traced every output file to the specific lines in the Python code that
produce them.

---


### The paired-end FASTQ files (the main output)

| File | What it is |
|------|------------|
| `ccap_cifi_HiC_R1.fastq.gz` | **Forward reads** — one fragment from each pair. Every line corresponds to one line in R2. This is what you feed to hifiasm via `--h1`. |
| `ccap_cifi_HiC_R2.fastq.gz` | **Reverse reads** — the other fragment from each pair. Paired 1:1 with R1 (line N in R1 is paired with line N in R2). This is what you feed to hifiasm via `--h2`. |

These two files together mimic **standard Hi-C paired-end data**. Each R1/R2
pair represents two genomic fragments that were physically close in the
nucleus. Hifiasm uses this proximity information to figure out:

- Which contigs belong to the **same chromosome** (scaffolding)
- Which contigs belong to the **same haplotype** (phasing)

**Why are these files larger than the input?** The input is one FASTQ file
(5.8G compressed). But a single CiFi read with 10 fragments produces
10×9/2 = 45 pairs. So the output is a combinatorial expansion — many more
records than the input, though each record is shorter (individual fragments
instead of full concatemer reads).

### The fragment statistics files

| File | What it is |
|------|------------|
| `ccap_cifi_fraglens.txt` | **Fragment lengths** — one number per line, each being the length (in bp) of a single fragment after in silico digestion. Tells you the distribution of fragment sizes produced by HindIII cutting. The CiFi paper reported a median of 1,893 bp for HindIII. |
| `ccap_cifi_fragcounts.txt` | **Fragment counts per read** — one number per line, each being the number of fragments a single CiFi read was cut into. Tells you how many HindIII sites were present per read. The CiFi paper reported a median of 2 fragments per read for HindIII. |

### The histogram plots

| File | What it is |
|------|------------|
| `ccap_cifi_fraglenhist.png` | **Fragment length histogram** — a bar chart showing the distribution of fragment lengths. Lets you visually check if the distribution matches the paper's Figure 1D (median ~1,893 bp for HindIII, range ~350 bp to 10,000 bp). |
| `ccap_cifi_fragcounthist.png` | **Fragment count histogram** — a bar chart showing how many fragments each read was split into. Should match the paper's Figure 1C (median 2 fragments for HindIII, with a tail extending to higher numbers). |

### The log file

| File | What it is |
|------|------------|
| `cifi2pe.log` | **Runtime log** — the standard output from the CiFi2PE run. Contains progress messages and any warnings or errors. Check this if the output files look wrong or are missing. |

---

## Which Files Actually Matter?

```
Essential (feed to hifiasm):
  ★ ccap_cifi_HiC_R1.fastq.gz  → hifiasm --h1
  ★ ccap_cifi_HiC_R2.fastq.gz  → hifiasm --h2

QC (compare with the paper's Figure 1C and 1D):
  ○ ccap_cifi_fraglenhist.png   → should show median ~1,893 bp
  ○ ccap_cifi_fragcounthist.png → should show median ~2 fragments/read

Reference data (useful for custom analysis):
  ○ ccap_cifi_fraglens.txt
  ○ ccap_cifi_fragcounts.txt
```

Only the R1 and R2 FASTQ files are needed for the assembly. Everything
else is for quality control and understanding the data.

---
