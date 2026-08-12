#!/bin/bash
#SBATCH --job-name=assembly_stats
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH --partition=normal

#======================================================================
# 03_assembly_stats.sh
#
# PURPOSE:
#   Compute assembly statistics (contig count, total size, N50, L50, N90, L90)
#   for both the HiFi-only and CiFi-phased assemblies.
#======================================================================

set -euo pipefail

ASSEMBLY="/export/groups/drangeli/CiFi.test/assembly"

echo "=========================================="
echo " Assembly Statistics"
echo " Date: $(date)"
echo "=========================================="

process() {
    local fa="$1"
    local label="$2"

    if [ ! -f "$fa" ]; then
        echo "  [SKIP] ${label} — not found"
        return
    fi

    echo ""
    echo "=== ${label} ==="
    echo "    File: $(basename $fa) ($(du -h $fa | cut -f1))"

    awk '
    /^>/ {
        if (n>0) { L[n]=len; tot+=len; if(len>mx) mx=len }
        n++; len=0; next
    }
    { len+=length($0) }
    END {
        if (n>0) { L[n]=len; tot+=len; if(len>mx) mx=len }
        for(i=1;i<=n;i++) s[i]=L[i]
        for(i=1;i<=n;i++) for(j=i+1;j<=n;j++) if(s[j]>s[i]){t=s[i];s[i]=s[j];s[j]=t}
        c=0
        for(i=1;i<=n;i++){
            c+=s[i]
            if(!n50&&c>=tot/2){n50=s[i];l50=i}
            if(!n90&&c>=tot*0.9){n90=s[i];l90=i}
        }
        printf "    Contigs:    %d\n",n
        printf "    Total size: %d bp (%.1f Mbp)\n",tot,tot/1e6
        printf "    Largest:    %d bp (%.1f Mbp)\n",mx,mx/1e6
        printf "    N50:        %d bp (%.1f Mbp)\n",n50,n50/1e6
        printf "    L50:        %d\n",l50
        printf "    N90:        %d bp (%.1f Mbp)\n",n90,n90/1e6
        printf "    L90:        %d\n",l90
    }' "$fa"
}

echo ""
echo "══════════════════════════════════════"
echo " HiFi-only Assembly"
echo "══════════════════════════════════════"

process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.p_ctg.fa" "Primary contigs"
process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap1.p_ctg.fa" "Haplotype 1"
process "${ASSEMBLY}/hifi_only/ccap_hifi_only.bp.hap2.p_ctg.fa" "Haplotype 2"

echo ""
echo "══════════════════════════════════════"
echo " CiFi-phased Assembly"
echo "══════════════════════════════════════"

process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.p_ctg.fa" "Primary contigs"
process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap1.p_ctg.fa" "Haplotype 1"
process "${ASSEMBLY}/cifi_phased/ccap_cifi_phased.hic.hap2.p_ctg.fa" "Haplotype 2"

echo ""
echo "══════════════════════════════════════"
echo " Compare with Paper Results"
echo "══════════════════════════════════════"
echo ""
echo "  Paper (Supplementary Table 9):"
echo "    Genome:  ~600 Mbp (2n=12, 5 autosomes + X/Y)"
echo "    BUSCO:   96.7–99.3% complete (Diptera odb10)"
echo "    Adj. QV: 41.4–58.1"
echo ""
echo "Done: $(date)"