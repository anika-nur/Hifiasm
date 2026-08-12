#!/bin/bash
#SBATCH --job-name=yahs_porec
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal

#======================================================================
# yahs_porec.sh
#
# From the paper,
#   Step 1: wf-pore-c Nextflow with --paired_end
#           (digest → minimap2 → annotate → pair, all handled internally)
#   Step 2: filter_bam (Dennis Lab MAPQ filter)
#           (retains only pairs where BOTH mates have MAPQ ≥ 1)
#   Step 3: YaHS scaffolding
#
# DEPENDENCIES: nextflow, samtools, yahs
#======================================================================

set -eo pipefail


source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pore-c-py

DATADIR="/export/groups/drangeli/CiFi.test"
CIFI_RAW="${DATADIR}/ERR14654081.fastq.gz"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
ENZYME="HindIII"
ENZYME_SEQ="AAGCTT"
MAPQ_THRESHOLD=1

FILTER_BAM="$(command -v filter_bam)"

echo "=========================================="
echo " YaHS with Paper's Full Pipeline"
echo " Step 1: wf-pore-c Nextflow (--paired_end)"
echo " Step 2: filter_bam (MAPQ >= ${MAPQ_THRESHOLD})"
echo " Step 3: YaHS scaffolding"
echo " Start:  $(date)"
echo "=========================================="

echo ""
echo "Active Conda environment: ${CONDA_DEFAULT_ENV:-unknown}"
echo "nextflow:  $(command -v nextflow)"
echo "samtools:  $(command -v samtools)"
echo "yahs:      $(command -v yahs)"
echo "filter_bam: ${FILTER_BAM}"

# ------------------------------------------------------------
# Verify required programs
# ------------------------------------------------------------

for PROGRAM in nextflow samtools yahs filter_bam; do
    if ! command -v "${PROGRAM}" > /dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} was not found in PATH."
        echo ""
        echo "Loaded modules:"
        module list 2>&1
        exit 1
    fi
done

if [ ! -x "${FILTER_BAM}" ]; then
    echo "ERROR: filter_bam was found but is not executable:"
    echo "       ${FILTER_BAM}"
    exit 1
fi

# ── Convert CiFi FASTQ to BAM (wf-pore-c expects BAM input) ──
CIFI_BAM="${DATADIR}/scaffolding_porec/cifi_input.bam"
mkdir -p "${DATADIR}/scaffolding_porec"
if [ ! -f "${CIFI_BAM}" ]; then
    echo ">>> Converting CiFi FASTQ to BAM..."
    samtools import -0 "${CIFI_RAW}" -o "${CIFI_BAM}"
fi

# ═══════════════════════════════════════════════════════
#  Function: full pipeline for one haplotype
# ═══════════════════════════════════════════════════════
scaffold_hap() {
    local HAP_FA="$1"
    local OUTDIR="$2"
    local PREFIX="$3"
    local LABEL="$4"

    echo ""
    echo "=========================================="
    echo " ${LABEL}"
    echo " Start: $(date)"
    echo "=========================================="

    mkdir -p "${OUTDIR}"
    cd "${OUTDIR}"

    # Copy contigs as reference
    cp "${HAP_FA}" contigs.fa
    samtools faidx contigs.fa

    # ── Step 1: wf-pore-c Nextflow pipeline ──
    # This handles: digest → minimap2 → annotate → pair
    # The --paired_end flag is what generates the paired-end BAM
    # --minimap2_settings '-x map-hifi' matches the paper's modification
    echo ""
    echo ">>> Step 1: Running wf-pore-c Nextflow pipeline..."
    echo "    (digest → minimap2 -x map-hifi → annotate → pair)"

    nextflow run epi2me-labs/wf-pore-c \
        -work-dir "${OUTDIR}/workflow" \
        --bam "${CIFI_BAM}" \
        --ref "${OUTDIR}/contigs.fa" \
        --cutter "${ENZYME}" \
        --out_dir "${OUTDIR}/results" \
        --threads ${THREADS} \
        --minimap2_settings '-x map-hifi' \
        --paired_end \
        --sample "${PREFIX}"

    echo "[OK] wf-pore-c done: $(date)"

    # The paired-end BAM is in: results/paired_end/
    PAIRED_BAM=$(find "${OUTDIR}/results" -name "*.paired_end.bam" | head -1)

    if [ -z "${PAIRED_BAM}" ] || [ ! -f "${PAIRED_BAM}" ]; then
        echo "ERROR: Paired-end BAM not found in ${OUTDIR}/results/"
        echo "Contents:"
        find "${OUTDIR}/results" -name "*.bam" -ls
        return
    fi
    echo "    Paired BAM: ${PAIRED_BAM} ($(du -h ${PAIRED_BAM} | cut -f1))"

    # ── Step 2: filter_bam (Dennis Lab MAPQ filter) ──
    # From the paper: "A filtering tool was developed to retain only
    # paired-end reads where both mates are mapped, exceed a MAPQ
    # threshold, and display proper pairing flags"
    # Tutorial: run with MAPQ threshold of 1
    echo ""
    echo ">>> Step 2: Filtering paired BAM (MAPQ >= ${MAPQ_THRESHOLD})..."
    echo "    ${FILTER_BAM} ${PAIRED_BAM} filtered.bam ${MAPQ_THRESHOLD} ${THREADS}"

    ${FILTER_BAM} "${PAIRED_BAM}" "${OUTDIR}/filtered.bam" ${MAPQ_THRESHOLD} ${THREADS}

    echo "[OK] Filtering done: $(date)"
    echo "    Before: $(samtools view -c ${PAIRED_BAM}) reads"
    echo "    After:  $(samtools view -c ${OUTDIR}/filtered.bam) reads"

    # Sort by name for YaHS
    samtools sort -n -@ 4 -o "${OUTDIR}/filtered.nsort.bam" "${OUTDIR}/filtered.bam"

    # ── Step 3: YaHS scaffolding ──
    echo ""
    echo ">>> Step 3: Running YaHS..."
    yahs contigs.fa filtered.nsort.bam -e ${ENZYME_SEQ} -o ${PREFIX}
    echo "[OK] YaHS done: $(date)"

    # ── Stats ──
    FA="${OUTDIR}/${PREFIX}_scaffolds_final.fa"
    if [ -f "$FA" ]; then
        samtools faidx "$FA"
        echo ""
        echo "=== Scaffold Stats: ${LABEL} ==="
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
    fi
}

# ═══════════════════════════════════════════════════════
#  Run both haplotypes
# ═══════════════════════════════════════════════════════
scaffold_hap \
    "${DATADIR}/assembly/cifi_phased/ccap_cifi_phased.hic.hap1.p_ctg.fa" \
    "${DATADIR}/scaffolding_porec/hap1" \
    "ccap_hap1" \
    "Haplotype 1 (Pore-C + filter_bam)"

scaffold_hap \
    "${DATADIR}/assembly/cifi_phased/ccap_cifi_phased.hic.hap2.p_ctg.fa" \
    "${DATADIR}/scaffolding_porec/hap2" \
    "ccap_hap2" \
    "Haplotype 2 (Pore-C + filter_bam)"

# ═══════════════════════════════════════════════════════
echo ""
echo "=========================================="
echo " Comparison"
echo "=========================================="
echo ""
echo " Paper targets (HindIII):"
echo "   Hap1: L90=6, N90=72.6 Mbp, N50=98.7 Mbp"
echo "   Hap2: L90=5, N90=77.9 Mbp, N50=81.4 Mbp"
echo ""
echo " BWA results (first attempt):"
echo "   Hap1: L90=19, N90=1.6 Mbp"
echo "   Hap2: L90=6,  N90=39.2 Mbp"
echo ""
echo "=========================================="
echo "Done: $(date)"