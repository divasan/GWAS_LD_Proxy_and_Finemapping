#!/bin/bash
#SBATCH --job-name=susie_all_loci
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --output=susie_%j.out
#SBATCH --error=susie_%j.err

set -euo pipefail

# Load R
module load R

# Run SuSiE fine-mapping for all prepared loci
Rscript Susie_scripts/run_susie_all_loci.R

echo "SuSiE fine-mapping completed."
