#!/usr/bin/env python3
"""
qc_per_sample.py
──────────────────────────────────────────────────────────────────────────────
Per-sample QC: ambient correction ingestion, scrublet doublet removal,
MT% filter, adaptive gene-count filter, Y-gene annotation.

--corrected accepts three types, auto-detected by the script:
  • .h5 file    → CellBender output  (sc.read_10x_h5)
  • directory   → SoupX MTX output   (sc.read_10x_mtx)
  • directory   → raw filtered mtx   (fallback / skip mode)

Metadata is received as a JSON file (all samplesheet columns).
Raw counts are stored in adata.layers["counts"] for scVI downstream.
"""
import argparse, json, warnings
from pathlib import Path

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# ── QC thresholds ─────────────────────────────────────────────────────────────
QC = {
    "min_genes"  : 500,   # hard floor
    "max_mt_pct" : 2.0,   # low for snRNA-seq (nuclei have little cytoplasmic MT)
    "min_nuclei" : 50,    # samples below this after QC are excluded
}
NUMERIC_COLS = {"braak", "age", "pmi"}
Y_GENES = [
    "AMELY","NLGN4Y","PCDH11Y","SRY","RPS4Y1","TBL1Y",
    "DDX3Y","UTY","EIF1AY","KDM5D","ZFY","USP9Y","TMSB4Y",
]

# ── Args ──────────────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--sample_id", required=True)
p.add_argument("--corrected", required=True,
               help="Ambient-corrected output: .h5 (CellBender), directory (SoupX), "
                    "or filtered_feature_bc_matrix dir (skip mode)")
p.add_argument("--filtered",  required=True,
               help="Raw filtered_feature_bc_matrix — fallback if corrected not found")
p.add_argument("--meta_json", required=True,
               help="JSON file with all samplesheet metadata columns")
args = p.parse_args()

sample_id = args.sample_id
print(f"\n[QC] {sample_id}")

# ── Load metadata ─────────────────────────────────────────────────────────────
with open(args.meta_json) as fh:
    meta = json.load(fh)

# ── Load count matrix — auto-detect source type ───────────────────────────────
corrected = Path(args.corrected)
filtered  = Path(args.filtered)

if corrected.suffix == ".h5" and corrected.exists():
    # ── CellBender h5 ──────────────────────────────────────────────────────────
    adata = sc.read_10x_h5(str(corrected))
    print("  CellBender h5 loaded")
elif corrected.is_dir() and corrected.exists():
    try:
        files = list(corrected.iterdir())
    except PermissionError:
        files = []
    if files:
        adata = sc.read_10x_mtx(corrected, var_names="gene_symbols", cache=False)
        label = "SoupX MTX" if "soupx" in corrected.name.lower() else "ambient-corrected MTX"
        print(f"  {label} loaded")
    else:
        # Empty dir — fall through to filtered
        adata = sc.read_10x_mtx(filtered, var_names="gene_symbols", cache=False)
        print("  WARNING: corrected dir empty — using raw filtered counts")
else:
    # skip mode OR corrected path not found
    if filtered.exists():
        adata = sc.read_10x_mtx(filtered, var_names="gene_symbols", cache=False)
        print("  Using raw filtered counts (no ambient correction)")
    else:
        print(f"  ERROR: no count matrix found at {corrected} or {filtered}")
        exit(1)

adata.var_names_make_unique()
adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

# Store raw integer counts for scVI
adata.layers["counts"] = adata.X.copy()

# ── Attach samplesheet metadata to obs ───────────────────────────────────────
SKIP_COLS = {"matrix_path", "raw_matrix_path"}
for col, val in meta.items():
    if col in SKIP_COLS:
        continue
    if col in NUMERIC_COLS:
        try:    adata.obs[col] = float(val)
        except: adata.obs[col] = np.nan
    else:
        v = str(val).strip() if val not in (None, "", "nan", "NA", "NaN") else "NA"
        adata.obs[col] = v

# ── QC metrics ───────────────────────────────────────────────────────────────
adata.var["mt"]   = adata.var_names.str.upper().str.startswith("MT-")
adata.var["rb"]   = adata.var_names.str.upper().str.match(r"^RP[SL]\d+$")
adata.var["chrY"] = adata.var_names.isin(Y_GENES)

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "rb", "chrY"],
    percent_top=None, log1p=False, inplace=True,
)

sc.pp.scrublet(adata, random_state=42)
n_doublets = int(adata.obs["predicted_doublet"].sum())

# ── Adaptive gene-count filter ────────────────────────────────────────────────
mean_g = adata.obs["n_genes_by_counts"].mean()
std_g  = adata.obs["n_genes_by_counts"].std()

# Individual masks — to report per-filter cell loss
f_doublet  = ~adata.obs["predicted_doublet"]
f_min_g    = adata.obs["n_genes_by_counts"] >= QC["min_genes"]
f_adap_lo  = adata.obs["n_genes_by_counts"] >= mean_g - 2 * std_g
f_adap_hi  = adata.obs["n_genes_by_counts"] <= mean_g + 2 * std_g
f_mt       = adata.obs["pct_counts_mt"]     <= QC["max_mt_pct"]

keep = f_doublet & f_min_g & f_adap_lo & f_adap_hi & f_mt

# Count cells lost to each filter independently (not cumulative)
n0 = adata.n_obs
lost_doublet  = int((~f_doublet).sum())
lost_min_g    = int((~f_min_g).sum())
lost_adap_lo  = int((~f_adap_lo).sum())
lost_adap_hi  = int((~f_adap_hi).sum())
lost_mt       = int((~f_mt).sum())

adata = adata[keep].copy()
print(f"  {n0:,} → {adata.n_obs:,} nuclei "
      f"({100*adata.n_obs/n0:.1f}% kept) | doublets: {n_doublets} | "
      f"MT>2%: {lost_mt} | low_genes: {lost_min_g} | "
      f"adaptive_low: {lost_adap_lo} | adaptive_high: {lost_adap_hi}")

if adata.n_obs < QC["min_nuclei"]:
    print(f"  EXCLUDED — {adata.n_obs} nuclei (threshold: {QC['min_nuclei']})")
    with open(f"{sample_id}.EXCLUDED", "w") as fh:
        fh.write(
            f"sample_id: {sample_id}\n"
            f"reason: too_few_nuclei\n"
            f"nuclei_after_qc: {adata.n_obs}\n"
            f"threshold: {QC['min_nuclei']}\n"
            f"nuclei_before_qc: {n0}\n"
            f"pct_kept: {100*adata.n_obs/n0:.1f}\n"
            f"doublets_removed: {lost_doublet}\n"
            f"lost_to_min_genes: {lost_min_g}\n"
            f"lost_to_adaptive_low: {lost_adap_lo}\n"
            f"lost_to_adaptive_high: {lost_adap_hi}\n"
            f"lost_to_mt_pct: {lost_mt}\n"
        )
    exit(0)

adata.write_h5ad(f"{sample_id}.h5ad")
print(f"  Saved: {sample_id}.h5ad")
