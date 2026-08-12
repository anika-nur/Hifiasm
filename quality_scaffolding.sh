#!/bin/bash
#SBATCH --job-name=quality_scaffolding
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal

#======================================================================
# quality_scaffolding.sh
#
# PURPOSE:
#   1. Extract the 10 longest contigs/scaffolds from each assembly
#   2. Run BUSCO on the CiFi-phased haplotype assemblies
#      using the diptera_odb10 database (3,285 genes)
#
# WHAT THE PAPER USED:
#   "Assessment of genome completeness was performed using BUSCO v5.7.1
#    to identify the presence of the 3285 genes in the Diptera v10
#    database (diptera_odb10)"
#======================================================================

set -euo pipefail

ASSEMBLY="/export/groups/drangeli/CiFi.test/assembly"
QC="/export/groups/drangeli/CiFi.test/qc"
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "${QC}"

echo "=========================================="
echo " Date: $(date)"
echo "=========================================="

# Scaffolds

echo ""
echo "=========================================="
echo " 10 Longest Contigs/Scaffolds"
echo "=========================================="

longest_contigs() {
    local fa="$1"
    local label="$2"

    if [ ! -f "$fa" ]; then
        echo "  [SKIP] ${label} — file not found: $(basename $fa)"
        return
    fi

    echo ""
    echo "=== ${label} ==="
    echo "    Rank    Contig Name                    Length (bp)      Length (Mbp)"
    echo "    ----    ------------                   -----------      ------------"

    awk '
    /^>/ {
        if (n > 0) { names[n] = name; lengths[n] = len }
        n++
        name = substr($1, 2)  # remove the ">"
        len = 0
        next
    }
    { len += length($0) }
    END {
        if (n > 0) { names[n] = name; lengths[n] = len }

        # Sort descending by length
        for (i = 1; i <= n; i++) {
            for (j = i+1; j <= n; j++) {
                if (lengths[j] > lengths[i]) {
                    tmp = lengths[i]; lengths[i] = lengths[j]; lengths[j] = tmp
                    tmps = names[i]; names[i] = names[j]; names[j] = tmps
                }
            }
        }

        # Print top 10
        top = (n < 10) ? n : 10
        for (i = 1; i <= top; i++) {
            printf "    %-7d %-30s %12d     %8.2f\n", i, names[i], lengths[i], lengths[i]/1e6
        }
        printf "\n    Total contigs: %d\n", n
    }' "$fa"
}

# HiFi-only
echo ""
echo "── HiFi-only Assembly ──"
longest_contigs "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.p_ctg.fa"      "Primary contigs"
longest_contigs "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap1.p_ctg.fa" "Haplotype 1"
longest_contigs "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap2.p_ctg.fa" "Haplotype 2"

# CiFi-phased
echo ""
echo "── CiFi-phased Assembly ──"
longest_contigs "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.p_ctg.fa"      "Primary contigs"
longest_contigs "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap1.p_ctg.fa" "Haplotype 1"
longest_contigs "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap2.p_ctg.fa" "Haplotype 2"



# BUSCO


echo ""
echo "=========================================="
echo " BUSCO Analysis"
echo "=========================================="



echo "BUSCO version: $(busco --version 2>&1)"
echo "Database: diptera_odb10 (3,285 genes)"
echo "Mode: genome"
echo ""

# ── Run BUSCO on each haplotype assembly ──
# The paper ran BUSCO on the phased, chromosome-scale scaffolds.
# I run on all 4 haplotype assemblies for comparison.

run_busco() {
    local fa="$1"
    local outname="$2"
    local label="$3"

    if [ ! -f "$fa" ]; then
        echo "  [SKIP] ${label} — file not found"
        return
    fi

    echo ">>> Running BUSCO: ${label}"
    echo "    Input:  $(basename $fa)"
    echo "    Output: ${QC}/${outname}"

    busco \
        -i "$fa" \
        -l diptera_odb10 \
        -o "$outname" \
        --out_path "$QC" \
        -m genome \
        --cpu ${THREADS} \
        -f \
        2>&1 | tail -20

    # Print the short summary 
    SUMMARY=$(find "${QC}/${outname}" -name "short_summary*.txt" 2>/dev/null | head -1)
    if [ -n "$SUMMARY" ] && [ -f "$SUMMARY" ]; then
        echo ""
        echo "    --- BUSCO Summary ---"
        grep -E "^[[:space:]]*(C:|S:|D:|F:|M:|Complete)" "$SUMMARY" || cat "$SUMMARY"
        echo ""
    fi
}

cd "${QC}"

# HiFi-only haplotypes
run_busco "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap1.p_ctg.fa" \
    "busco_hifi_hap1" "HiFi-only Haplotype 1"

run_busco "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap2.p_ctg.fa" \
    "busco_hifi_hap2" "HiFi-only Haplotype 2"

# CiFi-phased haplotypes
run_busco "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap1.p_ctg.fa" \
    "busco_cifi_hap1" "CiFi-phased Haplotype 1"

run_busco "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap2.p_ctg.fa" \
    "busco_cifi_hap2" "CiFi-phased Haplotype 2"

echo ""
echo "=========================================="
echo " Compare with Paper"
echo "=========================================="
echo ""
echo "All done: $(date)"