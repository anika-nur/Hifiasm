#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=normal

#======================================================================
# hifiasm.sh
#
# PURPOSE:
#   Run two hifiasm assmblies:
#     1: Uses only PacBio HiFi whole-genome sequencing reads (baseline) 
#     2: Uses the HiFi reads plus CiFi-derived paired-end data.
#   Comparing both shows exactly what CiFi data adds.
#
# INPUT FILES:
#   HiFi WGS:   /export/groups/drangeli/CiFi.test/ERR14654084.fastq.gz 
#   CiFi R1/R2: /export/groups/drangeli/CiFi.test/cifi_paired/ccap_cifi_HiC_R1.fastq.gz
#               /export/groups/drangeli/CiFi.test/cifi_paired/ccap_cifi_HiC_R2.fastq.gz
#======================================================================

set -euo pipefail

DATADIR="/export/groups/drangeli/CiFi.test"
ASSEMBLY="${DATADIR}/assembly"
TOOLS="${DATADIR}/tools"
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Input files
HIFI="${DATADIR}/ERR14654084.fastq.gz"
CIFI_R1="${DATADIR}/cifi_paired/ccap_cifi_HiC_R1.fastq.gz"
CIFI_R2="${DATADIR}/cifi_paired/ccap_cifi_HiC_R2.fastq.gz"

mkdir -p "${ASSEMBLY}/hifi_only" "${ASSEMBLY}/cifi_phased"

echo "=========================================="
echo " Hifiasm Assembly: Ceratitis capitata"
echo " Threads: ${THREADS}"
echo " Start:   $(date)"
echo "=========================================="

# Verifying hifiasm
module load hifiasm
echo "Using hifiasm version:"
hifiasm --version


#  HiFi-only assembly 
echo "=========================================="
echo " HiFi-only assembly"
echo "=========================================="

# Verify HiFi input 
if [ ! -f "${HIFI}" ]; then
    echo "ERROR: HiFi WGS file not found: ${HIFI}"
    exit 1
fi
echo "HiFi WGS: ${HIFI} ($(du -h ${HIFI} | cut -f1))"


cd "${ASSEMBLY}/hifi_only"

hifiasm \
    -o ccap_hifi_only \
    -t ${THREADS} \
    "${HIFI}" \
    2>&1 | tee hifiasm_hifi_only.log

echo "HiFi-only done: $(date)"


#  CiFi-phased with --dual-scaf
echo "=========================================="
echo " CiFi-phased assembly (--dual-scaf)"
echo "=========================================="

# Check CiFi paired-end files
for file in "${CIFI_R1}" "${CIFI_R2}"; do
    if [ ! -f "$file" ]; then
        echo "WARNING: Missing CiFi file: $file"
        echo "Skipping Mode 2. HiFi-only assembly already completed."
        echo "Done: $(date)"
        exit 0
    fi
done

echo "CiFi R1: ${CIFI_R1} ($(du -h ${CIFI_R1} | cut -f1))"
echo "CiFi R2: ${CIFI_R2} ($(du -h ${CIFI_R2} | cut -f1))"

cd "${ASSEMBLY}/cifi_phased"

hifiasm \
    -o ccap_cifi_phased \
    -t ${THREADS} \
    --dual-scaf \
    --h1 "${CIFI_R1}" \
    --h2 "${CIFI_R2}" \
    "${HIFI}" \
    2>&1 | tee hifiasm_cifi_phased.log

echo "CiFi-phased done: $(date)"

# outputs 
echo ""
echo "=========================================="
echo " Output files"
echo "=========================================="
echo ""
echo "--- HiFi-only ---"
ls -lhS "${ASSEMBLY}/hifi_only/ccap_hifi_only"*.gfa 2>/dev/null || echo "  (none)"
echo ""
echo "--- CiFi-phased ---"
ls -lhS "${ASSEMBLY}/cifi_phased/ccap_cifi_phased"*.gfa 2>/dev/null || echo "  (none)"
echo ""
echo "Finished at: $(date)"