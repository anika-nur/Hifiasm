#!/bin/bash
#SBATCH --job-name=yahs_final
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=normal

set -euo pipefail

module load yahs
module load samtools

SCAFFOLD="/export/groups/drangeli/CiFi.test/scaffolding"
ENZYME="AAGCTT"

echo "=========================================="
echo " Running YaHS (mapping already done)"
echo " Start: $(date)"
echo "=========================================="

# ── Haplotype 1 ──
echo ""
echo ">>> Scaffolding Haplotype 1..."
cd "${SCAFFOLD}/hap1"
yahs contigs.fa mapped_cifi.bam -e ${ENZYME} -o ccap_hap1
echo "[OK] Hap1 done: $(date)"

# ── Haplotype 2 ──
echo ""
echo ">>> Scaffolding Haplotype 2..."
cd "${SCAFFOLD}/hap2"
yahs contigs.fa mapped_cifi.bam -e ${ENZYME} -o ccap_hap2
echo "[OK] Hap2 done: $(date)"

echo ""
echo "=========================================="
echo " Scaffold Statistics"
echo "=========================================="

for HAP in hap1 hap2; do
    FA="${SCAFFOLD}/${HAP}/ccap_${HAP}_scaffolds_final.fa"
    if [ ! -f "$FA" ]; then
        echo "[SKIP] ${HAP} — output not found"
        continue
    fi

    samtools faidx "$FA"

    echo ""
    echo "=== ${HAP} ==="
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
        printf "    Scaffolds:  %d\n",n
        printf "    Total size: %d bp (%.1f Mbp)\n",tot,tot/1e6
        printf "    Largest:    %d bp (%.1f Mbp)\n",mx,mx/1e6
        printf "    N50:        %d bp (%.1f Mbp)\n",n50,n50/1e6
        printf "    L50:        %d\n",l50
        printf "    N90:        %d bp (%.1f Mbp)\n",n90,n90/1e6
        printf "    L90:        %d\n",l90
    }' "$FA"

    echo ""
    echo "    Top 10 scaffolds:"
    awk '/^>/{if(n)print name,len; name=substr($1,2); len=0; n++; next}
         {len+=length($0)}
         END{print name,len}' "$FA" \
    | sort -k2 -nr | head -10 \
    | awk '{printf "      %-30s %12d bp (%6.1f Mbp)\n", $1, $2, $2/1e6}'
done

echo ""
echo "Expected: L90 ≈ 6 (one scaffold per chromosome)"
echo "Done: $(date)"