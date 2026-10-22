#!/usr/bin/env Rscript
# Step 0 — SoupX ambient RNA correction
# Uses load10X() — reads outs/ folder including CellRanger clustering info,
# which makes autoEstCont() more accurate than manual SoupChannel().
#
# Input change from previous version:
#   --matrix_path  path to filtered_feature_bc_matrix directory
#                  (outs_dir is derived as dirname(matrix_path))
#   This replaces the old --outs_dir argument, which was constructed from
#   base_dir + sample_id.  Explicit paths in the samplesheet give more
#   flexibility (STARsolo, different CellRanger output structures, etc.)
#
# Output: standard 10x MTX dir → read with sc.read_10x_mtx() in Python.

suppressPackageStartupMessages({ library(SoupX); library(Matrix) })

# ── Parse arguments ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
p    <- list()
for (i in seq(1, length(args) - 1, 2)) p[[sub("^--", "", args[i])]] <- args[i + 1]

sample_id    <- p$sample_id
matrix_path  <- p$matrix_path    # path to filtered_feature_bc_matrix
out_mtx      <- p$outdir
report       <- p$report

# Derive outs_dir from matrix_path:
#   .../outs/filtered_feature_bc_matrix  →  .../outs
# For STARsolo or other structures the same logic applies:
#   .../Solo.out/GeneFull/filtered  →  .../Solo.out/GeneFull
outs_dir <- dirname(matrix_path)

cat(sprintf("\n[SoupX] %s\n", sample_id))
cat(sprintf("  matrix_path : %s\n", matrix_path))
cat(sprintf("  outs_dir    : %s\n", outs_dir))

dir.create(out_mtx, recursive = TRUE, showWarnings = FALSE)

# ── Helper: save corrected matrix as gzipped 10x MTX ─────────────────────────
save_mtx <- function(mat, outdir) {
    writeMM(mat, file.path(outdir, "matrix.mtx"))
    system(paste("gzip -f", file.path(outdir, "matrix.mtx")))
    writeLines(colnames(mat), file.path(outdir, "barcodes.tsv"))
    system(paste("gzip -f", file.path(outdir, "barcodes.tsv")))
    write.table(
        data.frame(V1 = rownames(mat), V2 = rownames(mat), V3 = "Gene Expression"),
        file.path(outdir, "features.tsv"),
        sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
    )
    system(paste("gzip -f", file.path(outdir, "features.tsv")))
    cat(sprintf("  Saved MTX → %s/\n", outdir))
}

# ── Helper: fallback — copy filtered matrix uncorrected ──────────────────────
fallback <- function(reason, matrix_path, out_mtx, report, sample_id) {
    cat(sprintf("  FALLBACK (%s) — using uncorrected filtered matrix\n", reason))
    filt <- Read10X(matrix_path)
    save_mtx(filt, out_mtx)
    write.csv(
        data.frame(sample_id = sample_id, rho = NA,
                   status = paste0("fallback: ", reason)),
        report, row.names = FALSE
    )
}

# ── Validate paths ────────────────────────────────────────────────────────────
if (!dir.exists(matrix_path)) {
    fallback("matrix_path not found", matrix_path, out_mtx, report, sample_id)
    quit(save = "no")
}

raw_dir <- file.path(outs_dir, "raw_feature_bc_matrix")

if (!dir.exists(raw_dir)) {
    cat("  WARNING: raw_feature_bc_matrix not found — skipping SoupX, using filtered\n")
    filt <- Read10X(matrix_path)
    save_mtx(filt, out_mtx)
    write.csv(
        data.frame(sample_id = sample_id, rho = NA, status = "no_raw"),
        report, row.names = FALSE
    )
    quit(save = "no")
}

# ── SoupX via load10X() ───────────────────────────────────────────────────────
# load10X() reads both filtered + raw matrices and CellRanger clustering,
# which improves contamination fraction estimation.
tryCatch({
    cat("  Running load10X()...\n")
    sc <- load10X(outs_dir)

    cat("  Running autoEstCont()...\n")
    sc <- autoEstCont(sc, doPlot = FALSE, verbose = FALSE)

    rho <- sc$metaData$rho[1]
    cat(sprintf("  Contamination fraction (rho): %.4f (%.1f%%)\n", rho, rho * 100))

    cat("  Running adjustCounts()...\n")
    corrected <- adjustCounts(sc, roundToInt = TRUE)

    save_mtx(corrected, out_mtx)

    # Y-gene counts before/after correction (mLOY QC)
    y_genes   <- c("AMELY", "RPS4Y1", "RPS4Y2", "TBL1Y", "DDX3Y", "UTY",
                   "EIF1AY", "KDM5D", "ZFY", "USP9Y", "TMSB4Y", "SRY",
                   "PCDH11Y", "NLGN4Y")
    y_present <- intersect(y_genes, rownames(corrected))
    y_raw     <- if (length(y_present) > 0) sum(sc$toc[y_present, ]) else NA
    y_corr    <- if (length(y_present) > 0) sum(corrected[y_present, ]) else NA

    write.csv(data.frame(
        sample_id          = sample_id,
        rho                = round(rho, 5),
        rho_pct            = round(rho * 100, 2),
        n_cells            = ncol(corrected),
        n_genes            = nrow(corrected),
        y_counts_raw       = y_raw,
        y_counts_corrected = y_corr,
        status             = "ok"
    ), report, row.names = FALSE)

    cat(sprintf("  Y-gene counts: raw=%s → corrected=%s\n",
                ifelse(is.na(y_raw), "NA", y_raw),
                ifelse(is.na(y_corr), "NA", y_corr)))
    cat("  SoupX done!\n")

}, error = function(e) {
    cat(sprintf("  SoupX error: %s\n", e$message))
    fallback(e$message, matrix_path, out_mtx, report, sample_id)
})
