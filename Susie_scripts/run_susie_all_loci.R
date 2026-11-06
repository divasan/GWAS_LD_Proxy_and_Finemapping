library(data.table)
library(susieR)

# Input and output directories
OUTDIR <- "results/susie_inputs"

sumstats_dir <- file.path(OUTDIR, "loci_sumstats")
snplist_dir  <- file.path(OUTDIR, "loci_snplists")
ld_dir       <- file.path(OUTDIR, "loci_ld")

rds_dir     <- file.path(OUTDIR, "loci_susie_rds")
pip_dir     <- file.path(OUTDIR, "loci_susie_pip_tables")
summary_dir <- file.path(OUTDIR, "loci_susie_summary_tables")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pip_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(sumstats_dir, pattern = "\\.tsv$", full.names = TRUE)

for (f in files) {

    base <- sub("\\.tsv$", "", basename(f))

    bim_file <- file.path(snplist_dir, paste0(base, ".bim"))
    ld_file  <- file.path(ld_dir, paste0(base, ".ld"))

    if (!file.exists(bim_file) || !file.exists(ld_file)) {
        cat("Skipping", base, "- missing reference files\n")
        next
    }

    cat("Processing", base, "\n")

    # Read GWAS summary statistics and reference variant information
    ss <- unique(fread(f), by = "rsid")
    bim <- fread(bim_file, header = FALSE)
    Rmat <- as.matrix(fread(ld_file, header = FALSE))

    setnames(bim, c("chr_ref", "rsid", "cm", "pos_ref",
                    "ref_a1", "ref_a2"))

    # Check that the LD matrix matches the reference variants
    if (nrow(Rmat) != nrow(bim) || ncol(Rmat) != nrow(bim)) {
        cat("Skipping", base, "- LD matrix dimension mismatch\n")
        next
    }

    # Match GWAS variants to the reference SNP order
    ss2 <- ss[match(bim$rsid, ss$rsid)]

    # Check allele alignment
    direct <- ss2$a1 == bim$ref_a1 & ss2$a2 == bim$ref_a2
    reverse <- ss2$a1 == bim$ref_a2 & ss2$a2 == bim$ref_a1

    valid <- !is.na(direct) & !is.na(reverse) &
             (direct | reverse) &
             is.finite(ss2$beta) &
             is.finite(ss2$se) & ss2$se > 0 &
             is.finite(ss2$n) & ss2$n > 0

    # Keep usable variants and subset the LD matrix in the same order
    idx <- which(valid)

    if (length(idx) < 2) {
        cat("Skipping", base, "- insufficient matched SNPs\n")
        next
    }

    ss2 <- ss2[idx]
    Rmat <- Rmat[idx, idx, drop = FALSE]

    # Calculate GWAS Z-scores
    z <- ss2$beta / ss2$se

    # Reverse the sign where effect alleles are swapped
    z[reverse[idx]] <- -z[reverse[idx]]

    # Use the reported GWAS sample size
    if (length(unique(ss2$n)) != 1) {
        cat("Skipping", base, "- inconsistent sample sizes\n")
        next
    }

    n <- ss2$n[1]

    # Run SuSiE-RSS
    fit <- tryCatch(
        susie_rss(
            z = z,
            R = Rmat,
            n = n,
            L = 10,
            estimate_residual_variance = FALSE
        ),
        error = function(e) {
            cat("SuSiE failed:", base, conditionMessage(e), "\n")
            NULL
        }
    )

    if (is.null(fit)) next

    # Save the fitted SuSiE model
    saveRDS(fit, file.path(rds_dir, paste0(base, ".rds")))

    # Extract PIPs and credible sets
    pip_dt <- copy(ss2)
    pip_dt[, locus := base]
    pip_dt[, pip := fit$pip]
    pip_dt[, credible_set := NA_character_]

    cs <- fit$sets$cs

    if (!is.null(cs)) {
        for (i in seq_along(cs)) {
            cs_idx <- cs[[i]]
            pip_dt[cs_idx, credible_set :=
                ifelse(
                    is.na(credible_set),
                    paste0("CS", i),
                    paste0(credible_set, ";CS", i)
                )
            ]
        }
    }

    fwrite(
        pip_dt[, .(locus, rsid, chr, pos, a1, a2,
                   beta, se, n, pip, credible_set)],
        file.path(pip_dir, paste0(base, ".pip.tsv")),
        sep = "\t"
    )

    # Summarize the fine-mapping results for each locus
    top_idx <- which.max(pip_dt$pip)

    summary_dt <- data.table(
        locus = base,
        n_snps = nrow(pip_dt),
        n_credible_sets = if (is.null(cs)) 0L else length(cs),
        top_rsid = pip_dt$rsid[top_idx],
        top_pip = pip_dt$pip[top_idx]
    )

    fwrite(
        summary_dt,
        file.path(summary_dir, paste0(base, ".summary.tsv")),
        sep = "\t"
    )
}

cat("SuSiE fine-mapping completed.\n")
