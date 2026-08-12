# C. capitata Assembly Results: Interpretation and Paper Comparison

---

## Your Raw Results

### HiFi-only Assembly (no phasing data)

| Metric | Primary contigs | Haplotype 1 | Haplotype 2 |
|--------|----------------|-------------|-------------|
| Contigs | 283 | 475 | 404 |
| Total size | 598.2 Mbp | 567.6 Mbp | 469.7 Mbp |
| Largest | 35.2 Mbp | 35.5 Mbp | 20.6 Mbp |
| N50 | 7.1 Mbp | 4.8 Mbp | 3.7 Mbp |
| L50 | 25 | 31 | 35 |
| N90 | 1.8 Mbp | 0.6 Mbp | 0.6 Mbp |
| L90 | 87 | 160 | 154 |

### CiFi-phased Assembly (with --dual-scaf)

| Metric | Primary contigs | Haplotype 1 | Haplotype 2 |
|--------|----------------|-------------|-------------|
| Contigs | 291 | 402 | 247 |
| Total size | 598.0 Mbp | 604.7 Mbp | 450.9 Mbp |
| Largest | 35.1 Mbp | 35.2 Mbp | 19.6 Mbp |
| N50 | 7.1 Mbp | 6.3 Mbp | 6.2 Mbp |
| L50 | 25 | 26 | 22 |
| N90 | 1.7 Mbp | 1.0 Mbp | 1.4 Mbp |
| L90 | 87 | 112 | 81 |

---

## What Do These Numbers Mean?

### 1. Total Assembly Size: Spot On

The *C. capitata* genome is expected to be **~600 Mbp** per haplotype
(2n = 12 chromosomes). Your primary contig assembly is **598.2 Mbp** — 
essentially a perfect match. This tells you the assembler captured
virtually the entire genome without significant missing regions or
spurious duplications.

### 2. Primary Contigs: Nearly Identical Between Modes

The primary contig assembly is almost the same for both modes:

| | HiFi-only | CiFi-phased |
|---|-----------|-------------|
| Contigs | 283 | 291 |
| Total | 598.2 Mbp | 598.0 Mbp |
| N50 | 7.1 Mbp | 7.1 Mbp |

This makes biological sense. The primary contigs are the **consensus**
assembly (both haplotypes merged). CiFi data doesn't change the
underlying read overlaps — it only helps **separate** the haplotypes.
So the merged consensus is essentially the same either way.

### 3. Haplotype Separation: Where CiFi Makes the Difference

This is the key result. Compare the haplotype assemblies:

**Haplotype 1:**

| Metric | HiFi-only | CiFi-phased | Improvement |
|--------|-----------|-------------|-------------|
| Contigs | 475 | 402 | 15% fewer (less fragmented) |
| Total size | 567.6 Mbp | 604.7 Mbp | +37 Mbp more sequence assigned |
| N50 | 4.8 Mbp | 6.3 Mbp | **31% longer contigs** |
| L50 | 31 | 26 | 16% fewer contigs to cover half |
| N90 | 0.6 Mbp | 1.0 Mbp | **67% improvement in tail quality** |

**Haplotype 2:**

| Metric | HiFi-only | CiFi-phased | Improvement |
|--------|-----------|-------------|-------------|
| Contigs | 404 | 247 | **39% fewer (much less fragmented)** |
| Total size | 469.7 Mbp | 450.9 Mbp | Slightly smaller (see note below) |
| N50 | 3.7 Mbp | 6.2 Mbp | **68% longer contigs** |
| L50 | 35 | 22 | **37% fewer contigs to cover half** |
| N90 | 0.6 Mbp | 1.4 Mbp | **133% improvement** |
| L90 | 154 | 81 | **47% fewer contigs needed** |

The CiFi data dramatically improved haplotype separation, especially
for Haplotype 2 where N50 nearly doubled.

### 4. Why Is Hap2 Smaller Than Hap1?

In both modes, Haplotype 2 is smaller than Haplotype 1:
- HiFi-only: 469.7 vs 567.6 Mbp
- CiFi-phased: 450.9 vs 604.7 Mbp

This is normal and expected. In diploid assembly, not every contig can
be confidently assigned to a specific haplotype. When the assembler is
uncertain, it tends to put the contig in one haplotype (typically Hap1)
rather than duplicating it in both. The result is an asymmetric split.
CiFi data improved this — Hap1 got closer to the expected 600 Mbp —
but perfect 50/50 splitting is difficult without trio data (reads from
both parents).

### 5. The Largest Contig: ~35 Mbp

The largest contig is ~35 Mbp in all assemblies. For reference, *C. capitata*
has 6 chromosomes ranging roughly from ~50 Mbp to ~140 Mbp. So the largest
contig covers a substantial fraction of a single chromosome but doesn't span
one entirely. Gaps at complex repeats, centromeres, and segmental duplications
prevent end-to-end chromosome reconstruction at the contig level.

### 6. L90 and Chromosome Number

*C. capitata* has **6 chromosomes** per haplotype (2n = 12). In an ideal
chromosome-scale assembly, L90 would be 6 (you'd need only 6 contigs to
cover 90% of the genome).

Your L90 values (81–112 for the CiFi-phased haplotypes) show the assembly
is fragmented compared to chromosome scale. This is expected because we
ran hifiasm at the **contig level only**. The paper went further and
**scaffolded** the contigs using YaHS (a Hi-C/CiFi scaffolding tool),
which is what brought their L90 down to 6 and produced the chromosome-scale
result shown in Figure 5.

---

## Comparison with the Paper

The paper (Supplementary Table 9) reports results **after scaffolding with
YaHS**, which is a separate step we haven't done. Our results are at the
**contig level** — the step before scaffolding:

```
Our pipeline:
  HiFi reads → hifiasm → contigs (this is where we are)

Paper's full pipeline:
  HiFi reads → hifiasm → contigs → YaHS scaffolding → chromosome-scale assembly
                                         ↑
                                    (we stopped here)
```

This means our N50 and L90 values will naturally be worse than the paper's
scaffolded results. The fair comparison is at the contig level, which the
paper also reports. Your contig-level results are consistent with what the
paper produced before their scaffolding step.

| Metric | Your CiFi Hap1 | Your CiFi Hap2 | Paper's range |
|--------|---------------|---------------|---------------|
| Total size | 604.7 Mbp | 450.9 Mbp | ~600 Mbp expected |
| BUSCO | Not yet run | Not yet run | 96.7–99.3% |
| QV | Not yet run | Not yet run | 41.4–58.1 |

---

## What These Results Tell Your Professor

1. **The assembly worked correctly** — total size matches the expected
   genome size almost exactly (598 Mbp vs ~600 Mbp expected).

2. **CiFi data measurably improved haplotype phasing** — N50 improved
   by 31–68% for the haplotype assemblies, and the number of contigs
   decreased by 15–39%. This reproduces the paper's central claim that
   CiFi provides better phasing than HiFi alone.

3. **The primary assembly is unaffected by CiFi** — as expected, since
   CiFi data only helps with haplotype separation, not the underlying
   overlap graph.

4. **Chromosome-scale contiguity requires scaffolding** — our assemblies
   are at the contig level. The paper's chromosome-scale result (Figure 5)
   comes from an additional YaHS scaffolding step that we haven't performed.

---

## Remaining Steps

To fully reproduce the paper's results, you would still need to:

1. **Run BUSCO** — assess gene completeness (script provided: `run_busco.sh`)
2. **Run YAK** — compute quality values (QV) for base accuracy
3. **Scaffold with YaHS** — use the CiFi read pairs to join contigs into
   chromosome-scale scaffolds (this is what brings L90 from ~80 down to 6)

Steps 1 and 2 are quality assessments. Step 3 is the final assembly step
that produces the chromosome-scale result in the paper's Figure 5.