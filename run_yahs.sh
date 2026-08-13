#!/bin/bash
#SBATCH --job-name=yahs_scaffold
#SBATCH --output=/export/groups/drangeli/CiFi.test/logs/%x_%A.out
#SBATCH --error=/export/groups/drangeli/CiFi.test/logs/%x_%A.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal

#======================================================================
# run_yahs.sh
#
# PURPOSE:
#   Scaffold the CiFi-phased hifiasm contigs into chromosome-scale
#   sequences using YaHS, reproducing Figure 5 of the CiFi paper.
#
# WHAT THE PAPER DID:
#   "The contig assemblies were scaffolded with the matching individual
#    CiFi reads using the Pore-C bioinformatic pipeline for
#    desegmentation, mapping, and read pairing and YAHS v.1.2a.2 for
#    scaffolding using the resulting read pairs."
#
# THE PIPELINE HAS 3 STEPS:
#
#   Step 1: Index the contig FASTA
#   Step 2: Map the CiFi paired-end reads to the contigs
#     and then sort by name, fix mates, sort by coord, mark duplicates,
#     and final sort by name (what YaHS prefers).
#   Step 3: RUN YaHS
# EXPECTED RESULT:
#   L90 should drop from ~80-112 (contigs) to ~6 (chromosomes)
#   because C. capitata has 2n=12 (6 chromosomes per haplotype).
#
# DEPENDENCIES:
#   bwa, samtools, yahs (all need to be available via module or PATH)
#======================================================================

set -euo pipefail

DATADIR="/export/groups/drangeli/CiFi.test"
ASM="${DATADIR}/assembly/cifi_phased"
SCAFFOLD="${DATADIR}/scaffolding"
TOOLS="${DATADIR}/tools"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# CiFi paired-end reads (output from CiFi2PE)
CIFI_R1="${DATADIR}/cifi_paired/ccap_cifi_HiC_R1.fastq.gz"
CIFI_R2="${DATADIR}/cifi_paired/ccap_cifi_HiC_R2.fastq.gz"

# HindIII restriction enzyme cut site
ENZYME="AAGCTT"

mkdir -p "${SCAFFOLD}/hap1" "${SCAFFOLD}/hap2"

echo "=========================================="
echo " YaHS Scaffolding: C. capitata"
echo " Threads:  ${THREADS}"
echo " Enzyme:   HindIII (${ENZYME})"
echo " Start:    $(date)"
echo "=========================================="

# Load modules
module load samtools/1.23 2>/dev/null || true

# ── Check/install bwa ──
if [ -f "${TOOLS}/bwa/bwa" ]; then
    export PATH="${TOOLS}/bwa:$PATH"
    echo "bwa: ${TOOLS}/bwa/bwa (already compiled)"
elif command -v bwa &> /dev/null; then
    echo "bwa: $(which bwa)"
else
    echo ">>> bwa not found. Compiling from source..."
    mkdir -p "${TOOLS}"
    cd "${TOOLS}"
    git clone https://github.com/lh3/bwa.git
    cd bwa && make -j${THREADS}
    export PATH="${TOOLS}/bwa:$PATH"
    
fi
YAHS_BIN="${TOOLS}/yahs/yahs"

# Verify CiFi reads exist
for f in "${CIFI_R1}" "${CIFI_R2}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: CiFi reads not found: $f"
        exit 1
    fi
done

#  Function: scaffold one haplotype
scaffold_hap() {
    local CONTIGS="$1"     # path to contig FASTA
    local OUTDIR="$2"      # output directory
    local LABEL="$3"       # label for logging
    local PREFIX="$4"      # output prefix for YaHS

    if [ ! -f "${CONTIGS}" ]; then
        echo "[SKIP] ${LABEL} — contig file not found: ${CONTIGS}"
        return
    fi

    echo ""
    echo "=========================================="
    echo " Scaffolding: ${LABEL}"
    echo " Contigs: $(basename ${CONTIGS})"
    echo "=========================================="

    cd "${OUTDIR}"

    cp "${CONTIGS}" "${OUTDIR}/contigs.fa"

    # Step 1: Index the contigs 
    echo ""
    echo ">>> Step 1: Indexing contigs..."

    # bwa index: creates the FM-index that bwa mem needs to map reads
    bwa index contigs.fa

    # samtools faidx: creates .fai index (contig names + lengths)
    samtools faidx contigs.fa

    echo "[OK] Indexing complete"

    # Step 2: Map CiFi reads to contigs
    echo ""
    echo ">>> Step 2: Mapping CiFi reads to contigs..."
    echo "    bwa mem -5SPM -t ${THREADS} contigs.fa R1 R2"
    echo "    → sort by name → fixmate → sort by coord → markdup → sort by name"
    echo ""

    # This is the standard Hi-C mapping pipeline:
    #
    # bwa mem -5SPM        Map with Hi-C flags
    #   | samtools view    Convert SAM to BAM
    #   | samtools sort -n Sort by read NAME (for fixmate)
    #   | samtools fixmate Fill in mate information
    #   | samtools sort    Sort by COORDINATE (for markdup)
    #   | samtools markdup Remove PCR duplicates
    #   | samtools sort -n Sort by NAME again (YaHS preference)
    #   > final.bam

    bwa mem -5SPM -t ${THREADS} contigs.fa "${CIFI_R1}" "${CIFI_R2}" \
        | samtools view -@ 4 -buS - \
        | samtools sort -@ 4 -n -T tmp_nsort -O bam - \
        | samtools fixmate -mr - - \
        | samtools sort -@ 4 -T tmp_csort -O bam - \
        | samtools markdup -r -s - - 2> markdup_stats.txt \
        | samtools sort -@ 4 -n -T tmp_nsort2 -O bam - \
        > mapped_cifi.bam

    echo "[OK] Mapping complete"
    echo "    Duplicate stats:"
    cat markdup_stats.txt

    # Step 3: Run YaHS
    echo ""
    echo ">>> Step 3: Running YaHS scaffolding..."
    echo "    ${YAHS_BIN} contigs.fa mapped_cifi.bam -e ${ENZYME} -o ${PREFIX}"
    echo ""

    ${YAHS_BIN} contigs.fa mapped_cifi.bam \
        -e "${ENZYME}" \
        -o "${PREFIX}" \
        2>&1 | tee yahs.log

    echo "[OK] YaHS complete"

    # Check output 
    # YaHS produces:
    #   {prefix}_scaffolds_final.fa  — the scaffolded FASTA 
    #   {prefix}_scaffolds_final.agp — AGP file showing how contigs are ordered
    #   {prefix}.bin                 — binary contact file
    #   {prefix}_scaffolds_final.fa.fai — index of scaffolded FASTA

    SCAFF_FA="${OUTDIR}/${PREFIX}_scaffolds_final.fa"
    if [ -f "${SCAFF_FA}" ]; then
        echo ""
        echo "=== Scaffold Stats: ${LABEL} ==="

        # Index the scaffold FASTA
        samtools faidx "${SCAFF_FA}"

        # Count scaffolds and compute basic stats
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
        }' "${SCAFF_FA}"

        echo ""
        echo "    Top 10 scaffolds:"
        awk '/^>/{if(n)print name,len; name=substr($1,2); len=0; n++; next}
             {len+=length($0)}
             END{print name,len}' "${SCAFF_FA}" \
        | sort -k2 -nr | head -10 \
        | awk '{printf "      %-30s %12d bp (%6.1f Mbp)\n", $1, $2, $2/1e6}'
    else
        echo "WARNING: Scaffold output not found: ${SCAFF_FA}"
        echo "Check yahs.log for errors."
    fi

    # Clean up temp files 
    rm -f tmp_nsort* tmp_csort* tmp_nsort2*
}

#  Scaffold both haplotypes

scaffold_hap \
    "${ASM}/ccap_cifi_phased.hic.hap1.p_ctg.fa" \
    "${SCAFFOLD}/hap1" \
    "CiFi-phased Haplotype 1" \
    "ccap_hap1"

scaffold_hap \
    "${ASM}/ccap_cifi_phased.hic.hap2.p_ctg.fa" \
    "${SCAFFOLD}/hap2" \
    "CiFi-phased Haplotype 2" \
    "ccap_hap2"

#  Final summary
echo ""
echo "=========================================="
echo " Scaffolding Complete"
echo "=========================================="
echo ""
echo " Output files:"
echo "   ${SCAFFOLD}/hap1/ccap_hap1_scaffolds_final.fa"
echo "   ${SCAFFOLD}/hap2/ccap_hap2_scaffolds_final.fa"
echo ""
echo " Expected result (paper Figure 5):"
echo "   6 chromosome-scale scaffolds per haplotype"
echo "   L90 ≈ 6 (down from ~80-112 at contig level)"
echo "   Chr5-Y translocation visible between chromosomes"
echo ""
echo "Done: $(date)"