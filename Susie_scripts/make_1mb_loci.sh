#!/bin/bash
set -euo pipefail

# Input files
SENTINELS="results/plink_output/matched_sentinels_in_reference.txt"

# 1000 Genomes EUR reference
REF_DIR="${REF_DIR:?Please set REF_DIR}"
REF_BASENAME="${REF_BASENAME:-1000G.EUR.QC}"

# Output directory
OUTDIR="results/susie_inputs"
mkdir -p "${OUTDIR}"

WINDOW_BP=1000000

# Get chromosome positions of the matched sentinel SNPs
: > "${OUTDIR}/sentinels_chrpos.tsv"

for chr in $(seq 1 22); do

    BIM="${REF_DIR}/${REF_BASENAME}.${chr}.bim"

    awk 'NR==FNR {ids[$1]=1; next}
         ($2 in ids) {print $1, $4, $2}' \
        "${SENTINELS}" "${BIM}" \
        >> "${OUTDIR}/sentinels_chrpos.tsv"

done

# Create +/- 1 Mb windows around each sentinel SNP
awk -v W="${WINDOW_BP}" 'BEGIN{OFS="\t"}
{
    start=$2-W
    if(start<1) start=1
    print $1,start,$2+W
}' "${OUTDIR}/sentinels_chrpos.tsv" \
| sort -k1,1n -k2,2n \
| awk 'BEGIN{OFS="\t"}
    NR==1 {c=$1;s=$2;e=$3;next}
    {
        if($1==c && $2<=e) {
            if($3>e) e=$3
        } else {
            print c,s,e
            c=$1;s=$2;e=$3
        }
    }
    END {if(NR>0) print c,s,e}
' > "${OUTDIR}/merged_loci.tsv"

# Add a locus identifier
awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,"chr"$1"_"$2"_"$3
}' "${OUTDIR}/merged_loci.tsv" > "${OUTDIR}/loci.bed"

echo "Number of loci: $(wc -l < "${OUTDIR}/loci.bed")"
echo "Saved: ${OUTDIR}/loci.bed"
