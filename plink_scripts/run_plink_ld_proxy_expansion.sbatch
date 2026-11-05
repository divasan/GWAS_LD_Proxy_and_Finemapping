#!/bin/bash
#SBATCH --job-name=LDproxy_1020
#SBATCH --time=36:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=plink_LDproxy_1020_%j.out
#SBATCH --error=plink_LDproxy_1020_%j.err

set -euo pipefail

# Load PLINK 1.9
# Update this line if your cluster uses a different module name.
module load plink/1.90b7.2

command -v plink >/dev/null 2>&1 || {
    echo "ERROR: PLINK not found in PATH"
    exit 1
}

plink --version


# --------------------------------------------------
# Input files
# --------------------------------------------------

# Locate the repository directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Sentinel SNPs from Shrine et al. (2023), Supplementary Table 6
# The file should contain one rsID per line.
SNP_LIST="${REPO_DIR}/data/1020_rsID_plink_input.txt"

# 1000 Genomes Phase 3 EUR reference in PLINK binary format
# One BED/BIM/FAM file set per chromosome (1-22)
#
# Set REF_DIR when submitting the job.
REF_DIR="${REF_DIR:?Please set REF_DIR to your 1000 Genomes EUR reference directory}"
REF_BASENAME="${REF_BASENAME:-1000G.EUR.QC}"

# LD parameters used in the original analysis
R2_THRESH="0.8"
WINDOW_KB="1000"

OUTTAG="LDproxy_1000G_EUR_r2ge${R2_THRESH}_win${WINDOW_KB}kb"


# --------------------------------------------------
# Output directories
# --------------------------------------------------

# Save the outputs under the repository's results directory.
OUTPUT_BASE="${REPO_DIR}/results"

RUNSTAMP="$(date +%Y%m%d_%H%M%S)"
JOB_ID="${SLURM_JOB_ID:-manual}"

OUT_DIR="${OUTPUT_BASE}/${OUTTAG}__${RUNSTAMP}__job${JOB_ID}"

mkdir -p "${OUT_DIR}"

# Use temporary scratch space for the PLINK intermediate files.
# WORK_BASE can be set to a preferred scratch directory.
WORK_BASE="${WORK_BASE:-${SLURM_TMPDIR:-${TMPDIR:-${OUT_DIR}}}}"
WORKDIR="${WORK_BASE}/${OUTTAG}__job${JOB_ID}"

mkdir -p "${WORKDIR}"

echo "SNP list: ${SNP_LIST}"
echo "Reference directory: ${REF_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "Working directory: ${WORKDIR}"
echo

# Check that the sentinel input exists.
if [[ ! -s "${SNP_LIST}" ]]; then
    echo "ERROR: Sentinel SNP list not found or empty: ${SNP_LIST}"
    exit 1
fi

cd "${WORKDIR}"


# --------------------------------------------------
# Prepare the sentinel SNP list
# --------------------------------------------------

# Copy the original SNP list and remove empty lines and
# duplicate IDs, if any.
cp "${SNP_LIST}" sentinels_raw.txt

grep -v '^[[:space:]]*$' sentinels_raw.txt \
    | sort -u > sentinels.txt

echo "Unique input IDs: $(wc -l < sentinels.txt | awk '{print $1}')" \
    | tee "${OUT_DIR}/counts.txt"


# --------------------------------------------------
# Find sentinel SNPs in the EUR reference
# --------------------------------------------------

: > all_pairs.ld
: > matched_all.txt
: > warnings.txt

# Process chromosomes separately because the reference
# files are stored as chromosome-specific PLINK datasets.
for chr in $(seq 1 22); do

    REF_PREFIX="${REF_DIR}/${REF_BASENAME}.${chr}"

    # Skip chromosomes with missing reference files.
    if [[ ! -f "${REF_PREFIX}.bed" ||
          ! -f "${REF_PREFIX}.bim" ||
          ! -f "${REF_PREFIX}.fam" ]]; then

        echo "WARNING: Missing reference for chr${chr}: ${REF_PREFIX}.*" \
            | tee -a warnings.txt
        continue
    fi

    echo "=== chr${chr} ==="

    # Identify sentinel SNPs present in the reference.
    plink \
        --bfile "${REF_PREFIX}" \
        --extract sentinels.txt \
        --make-just-bim \
        --out "chr${chr}_match" \
        >/dev/null

    if [[ ! -f "chr${chr}_match.bim" ]]; then

        echo "WARNING: chr${chr}_match.bim not created; skipping chr${chr}" \
            | tee -a warnings.txt
        continue
    fi

    cut -f2 "chr${chr}_match.bim" \
        | sort -u > "chr${chr}_matched.txt"

    nmatch=$(wc -l < "chr${chr}_matched.txt")

    if [[ "${nmatch}" -eq 0 ]]; then
        echo "chr${chr}: 0 matched IDs"
        continue
    fi

    cat "chr${chr}_matched.txt" >> matched_all.txt


    # --------------------------------------------------
    # Calculate LD proxies
    # --------------------------------------------------

    # Calculate pairwise LD between each matched sentinel
    # and nearby variants within a 1 Mb window.
    #
    # Retain SNP pairs with r2 >= 0.8.
    # --ld-window 999999 increases the pair-count limit
    # so dense regions are not restricted by the default.

    plink \
        --bfile "${REF_PREFIX}" \
        --r2 \
        --ld-snp-list "chr${chr}_matched.txt" \
        --ld-window-kb "${WINDOW_KB}" \
        --ld-window 999999 \
        --ld-window-r2 "${R2_THRESH}" \
        --out "chr${chr}_ld" \
        >/dev/null


    # Append chromosome-level LD pairs into one file.
    # Keep the header only from the first chromosome.
    if [[ -f "chr${chr}_ld.ld" ]]; then

        if [[ ! -s all_pairs.ld ]]; then
            cat "chr${chr}_ld.ld" > all_pairs.ld
        else
            tail -n +2 "chr${chr}_ld.ld" >> all_pairs.ld
        fi

    else

        echo "WARNING: chr${chr}_ld.ld missing for chr${chr}" \
            | tee -a warnings.txt
    fi

done


# --------------------------------------------------
# Matched and missing sentinel SNPs
# --------------------------------------------------

sort -u matched_all.txt > matched_sentinels_in_reference.txt

comm -23 \
    sentinels.txt \
    matched_sentinels_in_reference.txt \
    > missing_sentinels_not_in_reference.txt

echo "Matched sentinel IDs: $(wc -l < matched_sentinels_in_reference.txt | awk '{print $1}')" \
    | tee -a "${OUT_DIR}/counts.txt"

echo "Missing sentinel IDs: $(wc -l < missing_sentinels_not_in_reference.txt | awk '{print $1}')" \
    | tee -a "${OUT_DIR}/counts.txt"


# --------------------------------------------------
# Extract LD proxies
# --------------------------------------------------

# The PLINK 1.9 LD output has these columns:
#
# CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
#
# SNP_B (column 6) contains the linked variants.

if [[ ! -s all_pairs.ld ]]; then
    echo "ERROR: all_pairs.ld is empty. Check warnings.txt" >&2
    cp warnings.txt "${OUT_DIR}/" || true
    exit 1
fi

awk 'NR>1 {print $6}' all_pairs.ld \
    | sort -u > "proxies_r2ge${R2_THRESH}.txt"

# Combine matched sentinels with the LD proxy SNPs.
cat matched_sentinels_in_reference.txt \
    "proxies_r2ge${R2_THRESH}.txt" \
    | sort -u > "sentinels_plus_proxies_r2ge${R2_THRESH}.txt"


echo "Unique proxies r2>=${R2_THRESH}: $(wc -l < "proxies_r2ge${R2_THRESH}.txt" | awk '{print $1}')" \
    | tee -a "${OUT_DIR}/counts.txt"

echo "Union (sentinels+proxies): $(wc -l < "sentinels_plus_proxies_r2ge${R2_THRESH}.txt" | awk '{print $1}')" \
    | tee -a "${OUT_DIR}/counts.txt"


# --------------------------------------------------
# Copy final results
# --------------------------------------------------

cp -v all_pairs.ld \
    "${OUT_DIR}/${OUTTAG}.allchr.ld"

cp -v "proxies_r2ge${R2_THRESH}.txt" "${OUT_DIR}/"

cp -v "sentinels_plus_proxies_r2ge${R2_THRESH}.txt" "${OUT_DIR}/"

cp -v matched_sentinels_in_reference.txt "${OUT_DIR}/"

cp -v missing_sentinels_not_in_reference.txt "${OUT_DIR}/"

cp -v warnings.txt "${OUT_DIR}/"

echo
echo "LD proxy expansion completed."
echo "Results saved in: ${OUT_DIR}"
