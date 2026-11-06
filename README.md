# GWAS LD Proxy Analysis and Statistical Fine-Mapping

## Overview

Genome-wide association studies (GWAS) identify genetic variants associated with complex traits, but sentinel variants may not necessarily be the functional variants underlying these associations. Due to linkage disequilibrium (LD), multiple variants within a genomic region can be correlated with the same GWAS signal. Identifying these correlated variants provides a broader set of candidates for downstream functional genomic analyses.

We performed LD proxy expansion of lung-function-associated sentinel variants reported by Shrine et al. (2023), using PLINK 1.9 and the European ancestry (EUR) subset of the 1000 Genomes Project Phase 3 reference panel. Variants in strong LD (r² ≥ 0.8) within a 1 Mb window of each sentinel variant were retained.

## Input Data

`data/1020_rsID_plink_input.txt` contains 1,020 lung-function-associated sentinel variant identifiers from Supplementary Table 6 of Shrine et al. (2023), reconstructed from the matched and unmatched sentinel lists of the original PLINK analysis.

## Methods

LD proxy identification was performed using chromosome-specific PLINK binary genotype files (`.bed`, `.bim`, `.fam`) from the 1000 Genomes Phase 3 EUR reference panel.

| Parameter | Setting |
|---|---|
| Software | PLINK 1.9 |
| LD reference | 1000 Genomes Phase 3 EUR |
| Genome assembly | GRCh37/hg19 |
| LD threshold | r² ≥ 0.8 |
| Genomic window | 1,000 kb |
| Chromosomes | 1–22 |

The analysis script is available at `plink_scripts/run_plink_ld_proxy_expansion.sbatch`.

Reference resources: [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/resources) | [1000 Genomes Project](https://doi.org/10.1038/nature15393)

## Results

Of the 1,020 GWAS sentinel variants, 879 were represented in the EUR reference panel and used for LD expansion.

| Metric | Result |
|---|---:|
| Input GWAS sentinel variants | 1,020 |
| Sentinels matched to the EUR reference | 879 |
| Unmatched sentinel variants | 141 |
| Unique sentinel and LD proxy variants | 27,332 |
| Median LD-linked variants per sentinel | 13 |
| 95th percentile | 107 |
| 99th percentile | 227 |
| Maximum LD-linked variants for a single sentinel | 2,716 |

The figure below shows the distribution of LD-linked variant counts across the 879 reference-matched GWAS sentinel SNPs, ranked from the largest to the smallest LD neighborhood.

![LD expansion of lung-function-associated sentinel SNPs](figures/ld_proxy_expansion.png)

LD neighborhood sizes varied substantially across GWAS signals. The median sentinel SNP was associated with 13 LD-linked variants, indicating that most GWAS signals expanded to relatively small sets of correlated variants. At the 95th percentile, a sentinel SNP had 107 LD-linked variants, meaning approximately 95% of sentinel SNPs had
