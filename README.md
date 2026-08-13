# CiFi + Hifiasm: Reproducing Phased Diploid Assembly of *Ceratitis capitata*

[![Organism](https://img.shields.io/badge/Organism-Ceratitis%20capitata-blue)]()
[![Genome Size](https://img.shields.io/badge/Genome%20Size-~600%20Mbp-green)]()
[![Chromosomes](https://img.shields.io/badge/Chromosomes-2n%20%3D%2012-orange)]()
[![Status](https://img.shields.io/badge/Status-In%20Progress-yellow)]()

This repository contains scripts and documentation for reproducing the **de novo phased diploid genome assembly** of a single Mediterranean fruit fly (*Ceratitis capitata*) from the CiFi paper (McGinty et al., 2025). The project uses **CiFi-phased PacBio HiFi reads** assembled with **hifiasm** and scaffolded with **YaHS** to produce chromosome-scale haplotype-resolved assemblies.

---

## Table of Contents
 
- [Background](#background)
- [Datasets](#datasets)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Pipeline Steps](#pipeline-steps)
- [Results](#results)
- [Remaining Work](#remaining-work)
- [References](#references)
- [Acknowledgments](#acknowledgments)
---

## Background

**CiFi** (Chromosome conformation capture with HiFi sequencing) is a method that couples chromatin conformation capture (3C) with PacBio HiFi long-read sequencing. It enables analysis of genomic interactions across repetitive regions with low-input requirements (~60,000 cells), producing multi-kilobasepair reads containing multiple interacting concatenated segments.

This project reproduces the *C. capitata* assembly from the CiFi paper because it is the only dataset where **hifiasm** was used, and it exercises hifiasm's most advanced feature — `--dual-scaf` for CiFi-phased assembly. The medfly genome (~600 Mbp, 2n = 12 chromosomes) is also small enough to be feasible on an institutional HPC cluster.

All analyses were performed on the **NSCC HPC cluster** (SLURM scheduler).

---

## Datasets

All sequencing data is from **BioProject [PRJEB83708](https://www.ebi.ac.uk/ena/browser/view/PRJEB83708)** (European Nucleotide Archive). Only the Medfly HiFi WGS reads and CiFi HindIII reads were used. 


## Repository Structure

```
├── README.md                   
├── scripts/
│   ├── run_cifi2pe.py          # CiFi2PE paired-end conversion
│   ├── 02_run_hifiasm.sh       # Hifiasm: HiFi-only + CiFi-phased assembly
│   ├── process_cifi.sh         # GFA → FASTA conversion
│   ├── assembly_stats.sh       # N50/L50/N90/L90 computation
│   ├── run_busco.sh            # BUSCO quality assessment
│   ├── run_yahs.sh             # YaHS scaffolding (BWA method)
│   ├── run_yahs_hap1.sh        # YaHS for hap1 only (BWA)
│   ├── run_yahs_hap2.sh        # YaHS for hap2 only (BWA)
│   ├── hap1_stats.sh           # Scaffold stats for hap1
│   ├── run_yahs_porec_full.sh  # Full paper method (Pore-C + filter_bam + YaHS)
│   └── run_porec_mapping.sh    # Pore-C mapping only (for validation)
└── docs/
    └── project_status.md       # Detailed project status and notes
```

> **Note:** Raw sequencing data and large assembly outputs are not included in this repository due to file size.

---

## Dependencies

### Software

| Tool | Version Used | Purpose | Installation |
|------|-------------|---------|--------------|
| [hifiasm](https://github.com/chhylp123/hifiasm) | Current GitHub (≥ v0.24.0) | Genome assembly | Compiled from source |
| [CiFi2PE](https://github.com/mydennislab/CiFi) | Latest | CiFi → paired-end conversion | Cloned + pip install |
| [YaHS](https://github.com/c-zhou/yahs) | 1.2.2 | Hi-C scaffolding | Module or compiled |
| [BWA](https://github.com/lh3/bwa) | 0.7.19 | Short-read alignment (initial method) | Compiled from source |
| [minimap2](https://github.com/lh3/minimap2) | ≥ 2.28 | Long-read alignment (paper method) | Compiled from source |
| [samtools](https://github.com/samtools/samtools) | 1.23 | BAM/FASTA manipulation | Module |
| [BUSCO](https://busco.ezlab.org/) | 6.0.0 | Assembly completeness assessment | Module |
| [pore-c-py](https://github.com/epi2me-labs/pore-c-py) | Latest | Pore-C read processing | pip install |
| [filter_bam](https://github.com/anika-nur/2024-sep-mapqfilter) | Latest | MAPQ-based BAM filtering | Compiled from source (C++) |

### HPC Environment

- **Cluster:** NSCC (SLURM scheduler)
- **Available modules:** `samtools/1.23`, `sratoolkit/3.4.1`, `busco`, `yahs`, `anaconda3`

---

## Pipeline Steps

### 1. CiFi2PE: Paired-End Conversion

**Script:** `CiFi2PE.py`

A Python script that clones the CiFi2PE tool and processes the raw CiFi HindIII concatemer reads (`ERR14654081.fastq.gz`) into paired-end FASTQs. CiFi2PE digests each long read at HindIII restriction sites, then generates all pairwise combinations of the resulting fragments as synthetic R1/R2 read pairs. These pairs encode chromatin proximity information and are used downstream by hifiasm for contact-based phasing.

**Output:** `ccap_cifi_HiC_R1.fastq.gz` and `ccap_cifi_HiC_R2.fastq.gz`, plus QC files (fragment length distribution, fragment count histogram).

### 2. Hifiasm Assembly

**Script:** `hifiasm.sh`

Runs hifiasm in two modes:

1. **HiFi-only mode:** Standard assembly using only the WGS reads, producing primary and haplotype-resolved contigs.
2. **CiFi-phased mode (`--dual-scaf`):** Uses CiFi paired-end reads (`--h1`/`--h2`) for contact-based phasing, producing improved haplotype-resolved contigs.

All outputs are in GFA format.

### 3. GFA to FASTA Conversion

**Script:** `process_cifi.sh`
 
Converts all hifiasm GFA output files to FASTA format for use in downstream tools (BUSCO, YaHS, samtools). Uses an awk one-liner that extracts sequence lines (S-lines) from the GFA graph format.

### 4. Assembly Statistics
 
**Script:** `assembly_stats.sh`
 
Computes standard genome assembly metrics — total size, contig/scaffold count, N50, L50, N90, L90, and largest sequence — for all FASTA assemblies. The script reads FASTA files and calculates cumulative length-sorted statistics.


### 3. BUSCO Quality Assessment

**Script:** `quality_scaffolding.sh`

Evaluates assembly completeness against the `diptera_odb10` lineage dataset (3,285 conserved genes). Run in genome mode on all four haplotype assemblies. Also reports the 10 longest scaffolds in each assembly for quick inspection.


### 4. YaHS Scaffolding

**Scripts:** `run_yahs.sh`

Original script that runs the full BWA + YaHS pipeline for both haplotypes. 

**Method:** CiFi2PE paired-end reads are mapped to haplotype contigs with `bwa mem -5SPM` (Hi-C mode flags), duplicates are marked with `samtools markdup`, and the resulting BAM is passed to YaHS for scaffolding with `-e AAGCTT` (HindIII cut site).

[Still in progress. Looking at the YaHS-related scripts right now might be confusing, since I am still working on this part. The README will be updated once I'm done.]

---

## Results (So far)

### Assembly Statistics

Contig-level metrics from hifiasm:

| Metric | HiFi Primary | HiFi Hap1 | HiFi Hap2 | CiFi Primary | CiFi Hap1 | CiFi Hap2 |
|--------|:-----------:|:---------:|:---------:|:-----------:|:---------:|:---------:|
| Contigs | 283 | 475 | 404 | 291 | 402 | 247 |
| Total (Mbp) | 598.2 | 567.6 | 469.7 | 598.0 | 604.7 | 450.9 |
| N50 (Mbp) | 7.1 | 4.8 | 3.7 | 7.1 | 6.3 | 6.2 |
| L50 | 25 | 31 | 35 | 25 | 26 | 22 |
| N90 (Mbp) | 1.8 | 0.6 | 0.6 | 1.7 | 1.0 | 1.4 |
| L90 | 87 | 160 | 154 | 87 | 112 | 81 |

### BUSCO Completeness

BUSCO v6.0.0, `diptera_odb10` (3,285 genes), genome mode:

| Assembly | Complete | Single Copy | Duplicated | Fragmented | Missing |
|----------|:-------:|:-----------:|:----------:|:----------:|:-------:|
| HiFi-only Hap1 | 98.6% | 97.2% | 1.4% | 0.1% | 1.3% |
| HiFi-only Hap2 | 98.0% | 96.9% | 1.0% | 0.1% | 2.0% |
| CiFi-phased Hap1 | **99.7%** | 97.1% | 2.6% | 0.0% | 0.3% |
| CiFi-phased Hap2 | 97.3% | 96.8% | 0.5% | 0.1% | 2.6% |
| **Paper range** | **96.7–99.3%** | — | — | — | — |

Our BUSCO scores are within or above the range reported in the paper.

### Scaffold Statistics

YaHS scaffold results using the BWA mapping method, compared to the paper's reported values:

| Metric | Hap1 (BWA) | Hap2 (BWA) | Paper Hap1 | Paper Hap2 |
|--------|:---------:|:---------:|:---------:|:---------:|
| Scaffolds | 409 | 124 | 164 | 86 |
| N50 (Mbp) | 82.2 | 75.2 | 98.7 | 81.4 |
| L50 | 4 | 3 | 3 | 3 |
| N90 (Mbp) | 1.6 | 39.2 | 72.6 | 77.9 |
| L90 | 19 | 6 | 6 | 5 |
| Largest (Mbp) | 112.5 | 87.6 | 112.9 | 112.4 |

**Key finding:** Hap2 closely matches the paper. Hap1 shows L90 = 19 vs. the paper's L90 = 6. This discrepancy is attributed to using BWA (designed for short Illumina reads) instead of the paper's minimap2-based Pore-C pipeline (optimized for PacBio long reads). See [Remaining Work](#remaining-work) for the planned resolution.

---

## Remaining Work

- [ ] **Run the Pore-C + filter_bam + YaHS pipeline** — Replace BWA with the paper's exact scaffolding method (pore-c-py digest → minimap2 `-x map-hifi` → filter_bam → YaHS) to determine whether this resolves the hap1 L90 discrepancy. 
- [ ] **Compare BWA vs. Pore-C scaffolding results** — Side-by-side evaluation of scaffold metrics from both mapping strategies.
- [ ] **YAK quality value (QV) assessment** — The paper reports adjusted QV of 41.4–58.1; this has not yet been computed for this project.

---

## References

### Primary Paper

McGinty, S. P., Kaya, G., Sim, S. B., Makunin, A., Corpuz, R. L., Quail, M. A., Abuelanin, M., Lawniczak, M. K. N., Geib, S. M., Korlach, J., & Dennis, M. Y. (2025). CiFi: Accurate long-read chromosome conformation capture with low-input requirements. *Nature Communications*, 17, 215. https://doi.org/10.1038/s41467-025-66918-y

### Software and Repositories

| Tool | Citation / Repository |
|------|-----------------------|
| hifiasm | Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*, 18, 170–175. [GitHub](https://github.com/chhylp123/hifiasm) |
| YaHS | Zhou, C., McCarthy, S. A., & Durbin, R. (2023). YaHS: yet another Hi-C scaffolding tool. *Bioinformatics*, 39(1), btac808. [GitHub](https://github.com/c-zhou/yahs) |
| BUSCO | Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Molecular Biology and Evolution*, 38(10), 4647–4654. [Website](https://busco.ezlab.org/) |
| BWA | Li, H. & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler aligner. *Bioinformatics*, 25(14), 1754–1760. [GitHub](https://github.com/lh3/bwa) |
| minimap2 | Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. [GitHub](https://github.com/lh3/minimap2) |
| CiFi2PE / CiFi pipeline | Dennis Lab, UC Davis. [GitHub](https://github.com/mydennislab/CiFi) |
| wf-pore-c | Oxford Nanopore Technologies / EPI2ME Labs. [GitHub](https://github.com/epi2me-labs/wf-pore-c) |
| pore-c-py | Oxford Nanopore Technologies / EPI2ME Labs. [GitHub](https://github.com/epi2me-labs/pore-c-py) |
| filter_bam | Dennis Lab MAPQ filter tool. [GitHub](https://github.com/anika-nur/2024-sep-mapqfilter) |
| samtools | Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. [GitHub](https://github.com/samtools/samtools) |

### Data

- **BioProject:** [PRJEB83708](https://www.ebi.ac.uk/ena/browser/view/PRJEB83708) (European Nucleotide Archive)

---

## Acknowledgments

This work was conducted at Colby College under the guidance of Dr. David R. Angelini.

**Dr. David R. Angelini, Ph.D.**
Professor and Associate Chair, Department of Biology
Colby College, Waterville, Maine, USA
[Faculty Profile](https://www.colby.edu/people/people-directory/dave-angelini/) · [Lab Website](https://web.colby.edu/aphanotus/) · [ResearchGate](https://www.researchgate.net/profile/David-Angelini)

Computations were performed on the NSCC HPC cluster. We thank the authors of the CiFi paper and all open-source tool developers whose software made this analysis possible.
