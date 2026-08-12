#!/bin/bash
#SBATCH --job-name=process_cifi
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --partition=normal

#======================================================================
# process_cifi.sh
#
# PURPOSE:
#   Convert hifiasm GFA output to FASTA 
#   for both the HiFi-only and CiFi-phased assemblies.
#======================================================================


set -euo pipefail

ASSEMBLY="/export/groups/drangeli/CiFi.test/assembly"

echo "=========================================="
echo " GFA to FASTA Conversion"
echo " Date: $(date)"
echo "=========================================="


process() {
    local gfa="$1"
    local label="$2"

    if [ ! -f "$gfa" ]; then
        echo "  [SKIP] ${label} — not found"
        return
    fi

    local fa="${gfa%.gfa}.fa"
    awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"

    echo ""
    echo "=== ${label} ==="
    echo "    Created: $(basename $fa) ($(du -h $fa | cut -f1))"
}

echo ""
echo "══════════════════════════════════════"
echo " HiFi-only Assembly"
echo "══════════════════════════════════════"

process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.p_ctg.gfa"     "Primary contigs"
process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap1.p_ctg.gfa" "Haplotype 1"
process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap2.p_ctg.gfa" "Haplotype 2"

echo ""
echo "══════════════════════════════════════"
echo " CiFi-phased Assembly"
echo "══════════════════════════════════════"

process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.p_ctg.gfa"     "Primary contigs"
process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap1.p_ctg.gfa"   "Haplotype 1"
process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap2.p_ctg.gfa"   "Haplotype 2"

echo ""
echo "Conversion complete."
echo "Done: $(date)"