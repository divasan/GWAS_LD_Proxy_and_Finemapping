#!/bin/bash
#SBATCH --job-name=ld_susie_all_loci
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=ld_susie_%j.out
#SBATCH --error=ld_susie_%j.err

set -euo pipefail

module load plink/1.90b7.2

# Input GWAS summary statistics with hg19 coordinates
SUMSTATS="data/icgcukb_b38_070825.with_hg19_coords.tsv"

# 1000 Genomes Phase 3 EUR reference
REF_DIR="${REF_DIR:?Please set REF_DIR}"
REF_BASENAME="${REF_BASENAME:-1000G.EUR.QC}"

# Input loci from Step 1
OUTDIR="results/susie_inputs"
LOCI="${OUTDIR}/loci.bed"

mkdir -p "${OUTDIR}/loci_sumstats"
mkdir -p "${OUTDIR}/loci_snplists"
mkdir -p "${OUTDIR}/loci_ld"

# Check the required input files
[[ -s "${SUMSTATS}" ]] || { echo "Missing GWAS summary statistics"; exit 1; }
[[ -s "${LOCI}" ]] || { echo "Missing loci.bed"; exit 1; }

# Process each merged 1 Mb locus
while read -r chr start end locus_id; do

    base="${chr}_${start}_${end}"

    echo "Processing ${base}"

    SS="${OUTDIR}/loci_sumstats/${base}.tsv"
    REF_PREFIX="${REF_DIR}/${REF_BASENAME}.${chr}"
    SNP_PREFIX="${OUTDIR}/loci_snplists/${base}"
    LD_PREFIX="${OUTDIR}/loci_ld/${base}"

    # Extract GWAS variants within the locus
    # Columns 17 and 18 contain the hg19 coordinates
    awk -v c="$chr" -v s="$start" -v e="$end" \
        'BEGIN {FS=OFS="\t"}
        NR==1 {
            print "chr","pos","rsid","a1","a2","beta","se","n"
            next
        }
        $17==c && $18>=s && $18<=e {
            print $17,$18,$4,$5,$6,$8,$9,$13
        }' "${SUMSTATS}" > "${SS}"

    # Extract the variant IDs for this locus
    tail -n +2 "${SS}" | cut -f3 | sort -u \
        > "${SNP_PREFIX}.input_snps"

    if [[ ! -s "${SNP_PREFIX}.input_snps" ]]; then
        echo "No variants found for ${base}"
        continue
    fi

    # Match the GWAS variants to the EUR reference
    # Save a locus-specific PLINK dataset
    plink \
        --bfile "${REF_PREFIX}" \
        --extract "${SNP_PREFIX}.input_snps" \
        --make-bed \
        --out "${SNP_PREFIX}" \
        >/dev/null

    nsnps=$(wc -l < "${SNP_PREFIX}.bim")

    if [[ "${nsnps}" -lt 2 ]]; then
        echo "Skipping ${base}: fewer than 2 matched variants"
        continue
    fi

    # Save the matched SNP IDs
    cut -f2 "${SNP_PREFIX}.bim" \
        > "${SNP_PREFIX}.snplist"

    # Calculate the signed LD correlation matrix
    # SuSiE uses correlation (r), not r-squared
    plink \
        --bfile "${SNP_PREFIX}" \
        --r square \
        --out "${LD_PREFIX}" \
        >/dev/null

done < "${LOCI}"

echo "Finished preparing SuSiE inputs."
