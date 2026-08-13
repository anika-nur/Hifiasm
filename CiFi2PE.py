#!/usr/bin/env python3
"""
CiFi2PE.py — Clone and run CiFi2PE on the C. capitata CiFi HindIII reads.

What this script does:
    1. Clones the CiFi2PE fork from my github
    2. Installs Python dependencies (biopython, matplotlib, numpy)
    3. Decompresses the CiFi FASTQ (required by BioPython)
    4. Runs cifi2pe_full_length_args.py
    5. Compresses the output R1/R2 files
"""

import subprocess
import os
import sys
import gzip
import shutil



DATA_DIR   = "/export/groups/drangeli/CiFi.test"
CIFI_GZ    = os.path.join(DATA_DIR, "ERR14654081.fastq.gz")   # CiFi HindIII reads 
CIFI_FQ    = os.path.join(DATA_DIR, "ERR14654081.fastq")      
TOOLS_DIR  = os.path.join(DATA_DIR, "tools")
OUTPUT_DIR = os.path.join(DATA_DIR, "cifi_paired")
REPO_URL   = "https://github.com/anika-nur/CiFi2PE.git"
ENZYME     = "HindIII"   # The restriction enzyme used for C. capitata CiFi library
OUTPUT_PREFIX = os.path.join(OUTPUT_DIR, "ccap_cifi")


def run_cmd(cmd, description):
    """Run a shell command, print what we're doing, and check for errors."""
    print(f"\n>>> {description}")
    print(f"    $ {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"    ERROR: Command failed with exit code {result.returncode}")
        sys.exit(1)
    print(f"    [OK]")


def main():
    print("=" * 65)
    print(" CiFi2PE Runner — Converting CiFi reads to paired-end format")
    print(" Input:  ERR14654081.fastq.gz (CiFi HindIII, C. capitata)")
    print(" Enzyme: HindIII (recognition site: AAGCTT)")
    print("=" * 65)

    # output directories
    os.makedirs(TOOLS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\n[OK] Directories ready")

    # Cloning the CiFi2PE repository
    # CiFi2PE contains two scripts:
    #   - cifi2pe_full_length_args.py: keeps full-length fragments
    #   - cifi2pe_150bp_args.py: trims fragments to 150bp
    #
    # I am using the full-length version because hifiasm can handle long reads

    repo_dir = os.path.join(TOOLS_DIR, "CiFi2PE")
    if not os.path.isdir(repo_dir):
        run_cmd(f"git clone {REPO_URL} {repo_dir}",
                "Cloning CiFi2PE repository")
    else:
        print(f"\n[OK] CiFi2PE already cloned at {repo_dir}")

    script_path = os.path.join(repo_dir, "cifi2pe_full_length_args.py")
    if not os.path.isfile(script_path):
        print(f"ERROR: {script_path} not found!")
        print("Contents of repo:")
        for f in os.listdir(repo_dir):
            print(f"  {f}")
        sys.exit(1)

    # Installing Python dependencies
    #
    # CiFi2PE depends on:
    #   - BioPython: for reading FASTQ files (SeqIO) and for the
    #     Restriction module that knows enzyme cut sites
    #   - matplotlib + numpy: for generating QC histograms
    # ─────────────────────────────────────────────────────────────
    run_cmd("pip install --user --break-system-packages biopython matplotlib numpy 2>/dev/null || "
            "pip install --user biopython matplotlib numpy",
            "Installing Python dependencies (biopython, matplotlib, numpy)")

   
    # Decompressing the CiFi FASTQ
    if not os.path.isfile(CIFI_GZ):
        print(f"ERROR: Input file not found: {CIFI_GZ}")
        sys.exit(1)

    if not os.path.isfile(CIFI_FQ):
        run_cmd(f"gunzip -k {CIFI_GZ}",
                f"Decompressing {os.path.basename(CIFI_GZ)} (this takes a few minutes)")
    else:
        print(f"\n[OK] Uncompressed FASTQ already exists: {CIFI_FQ}")

    # Show file size
    size_gb = os.path.getsize(CIFI_FQ) / (1024**3)
    print(f"    Uncompressed size: {size_gb:.1f} GB")

    # ─────────────────────────────────────────────────────────────
    # STEP 5: Run CiFi2PE
    #
    #   For each HiFi read in the FASTQ file:
    #     1. Use BioPython Restriction.HindIII.catalyse(read.seq)
    #        to cut the read at every AAGCTT site.
    #
    #     2. If the read has ≤3 fragments, then skip
    #        (Too few fragments to generate meaningful pairs.)
    #
    #     3. If >3 fragments, generate all pairwise combinations.
    #        For example, if a read has fragments [A, B, C, D], produce:
    #          Pair 1: R1=A, R2=B
    #          Pair 2: R1=A, R2=C
    #          Pair 3: R1=A, R2=D
    #          Pair 4: R1=B, R2=C
    #          Pair 5: R1=B, R2=D
    #          Pair 6: R1=C, R2=D
    #
    #
    #     4. Write R1s to {prefix}_HiC_R1.fastq
    #        Write R2s to {prefix}_HiC_R2.fastq
    #        (with quality scores sliced to match each fragment)
    #
    #     5. Also record fragment lengths and counts for QC histograms.
    #
    # OUTPUT FILES:
    #   {prefix}_HiC_R1.fastq      → forward reads for hifiasm --h1
    #   {prefix}_HiC_R2.fastq      → reverse reads for hifiasm --h2
    #   {prefix}_fraglens.txt       → one fragment length per line
    #   {prefix}_fragcounts.txt     → fragments per read
    #   {prefix}_fraglenhist.png    → histogram of fragment lengths
    #   {prefix}_fragcounthist.png  → histogram of fragment counts
    # ─────────────────────────────────────────────────────────────
    cmd = (f"python3 {script_path} "
           f"{CIFI_FQ} "
           f"{ENZYME} "
           f"--out {OUTPUT_PREFIX}")

    run_cmd(cmd, "Running CiFi2PE")

    # Compressing output and verify
    # ─────────────────────────────────────────────────────────────
    r1 = OUTPUT_PREFIX + "_HiC_R1.fastq"
    r2 = OUTPUT_PREFIX + "_HiC_R2.fastq"

    for fq in [r1, r2]:
        if os.path.isfile(fq):
            run_cmd(f"gzip {fq}", f"Compressing {os.path.basename(fq)}")
        elif os.path.isfile(fq + ".gz"):
            print(f"[OK] {os.path.basename(fq)}.gz already exists")
        else:
            print(f"WARNING: Expected output not found: {fq}")


    # Cleaning up the uncompressed FASTQ
    if os.path.isfile(CIFI_FQ) and os.path.isfile(CIFI_GZ):
        os.remove(CIFI_FQ)
        print(f"\n[OK] Removed uncompressed FASTQ")


    print("\n" + "=" * 65)
    print(" DONE! Output files:")
    print("=" * 65)
    for f in sorted(os.listdir(OUTPUT_DIR)):
        full = os.path.join(OUTPUT_DIR, f)
        if os.path.isfile(full):
            size = os.path.getsize(full)
            if size > 1e9:
                print(f"  {f}  ({size/1e9:.1f} GB)")
            elif size > 1e6:
                print(f"  {f}  ({size/1e6:.1f} MB)")
            else:
                print(f"  {f}  ({size/1e3:.1f} KB)")


if __name__ == "__main__":
    main()