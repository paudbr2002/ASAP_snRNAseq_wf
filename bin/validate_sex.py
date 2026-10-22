#!/usr/bin/env python3
"""
validate_sex.py
────────────────────────────────────────────────────────────────────────────────
Run AFTER QC_PER_SAMPLE, BEFORE INTEGRATE.

Flags patients with wrong sex annotation:
  - Males   with XIST_high  → probably female
  - Females with chrY_high  → probably male

Outputs:
  1. validated_samples.csv  — clean sample list for INTEGRATE
  2. sex_corrections.csv    — patients whose sex was corrected (for xlsx update)
  3. excluded_samples.csv   — samples excluded from integration
  4. heatmap_chrY_XIST.png  — per-patient expression heatmap

Usage:
    python validate_sex.py
        --h5ad_dir  <dir with per-sample .h5ad files>
        --samplesheet <original ASAP_samplesheet.xlsx>
        --outdir    <output directory>
"""
import argparse, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# ── Args ──────────────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--h5ad_dir",    required=True)
p.add_argument("--samplesheet", required=True)
p.add_argument("--outdir",      default=".")
args = p.parse_args()

OUT = Path(args.outdir); OUT.mkdir(parents=True, exist_ok=True)

Y_GENES   = ["RPS4Y1","DDX3Y","UTY","EIF1AY","KDM5D","ZFY","USP9Y","TMSB4Y"]
X_GENE    = ["XIST"]

# ── Load all per-sample h5ads (light — only need obs + a few genes) ───────────
h5ads = sorted(Path(args.h5ad_dir).glob("*.h5ad"))
print(f"\n-- Loading {len(h5ads)} sample h5ads --")

frames = []
for f in h5ads:
    ad_full = sc.read_h5ad(f, backed='r')
    genes   = [g for g in Y_GENES + X_GENE if g in ad_full.var_names]
    ad      = ad_full[:, genes].to_memory()
    ad_full.file.close(); del ad_full

    arr_y = ad[:, [g for g in Y_GENES if g in ad.var_names]].X
    arr_x = ad[:, [g for g in X_GENE  if g in ad.var_names]].X
    if hasattr(arr_y, "toarray"): arr_y = arr_y.toarray()
    if hasattr(arr_x, "toarray"): arr_x = arr_x.toarray()

    tmp = ad.obs[["sample_id","patient","sex"]].copy()
    tmp["chrY_score"] = np.asarray(arr_y).sum(axis=1)
    tmp["XIST_score"] = np.asarray(arr_x).sum(axis=1) if arr_x.size else 0.0
    frames.append(tmp)
    print(f"  [{f.stem}] {len(tmp):,} nuclei")

obs = pd.concat(frames, ignore_index=True)

# ── Aggregate per patient ─────────────────────────────────────────────────────
per_patient = (obs.groupby(["patient","sex"])
               .agg(chrY_score=("chrY_score","mean"),
                    XIST_score=("XIST_score","mean"),
                    n_nuclei  =("sample_id","count"))
               .reset_index())

# ── Thresholds ────────────────────────────────────────────────────────────────
m_mask = per_patient["sex"].str.strip().isin(["M","male"])
f_mask = per_patient["sex"].str.strip().isin(["F","female"])

xist_high_thr = per_patient.loc[f_mask, "XIST_score"].quantile(0.50) if f_mask.sum() else np.inf
chry_high_thr = per_patient.loc[m_mask, "chrY_score"].quantile(0.75) if m_mask.sum() else np.inf
xist_low_thr  = per_patient.loc[f_mask, "XIST_score"].quantile(0.25) if f_mask.sum() else 0.0

print(f"\n  Thresholds:")
print(f"    XIST high (males)   > {xist_high_thr:.3f}  (Q90 females)")
print(f"    chrY high (females) > {chry_high_thr:.3f}  (Q75 males)")
print(f"    XIST low  (females) < {xist_low_thr:.3f}   (Q25 females)")

def flag_patient(row):
    sex = str(row["sex"]).strip()
    if sex in ("M","male"):
        if row["XIST_score"] > xist_high_thr:
            return "XIST_high→misannotated", "F"
    elif sex in ("F","female"):
        if row["chrY_score"] > chry_high_thr and row["XIST_score"] < xist_low_thr:
            return "chrY_high+XIST_low→misannotated", "M"
    return "ok", row["sex"]

per_patient[["flag","corrected_sex"]] = per_patient.apply(
    flag_patient, axis=1, result_type="expand")

flagged   = per_patient[per_patient["flag"] != "ok"]
print(f"\n  Flagged patients: {len(flagged)}")
if len(flagged):
    print(flagged[["patient","sex","chrY_score","XIST_score","flag"]].to_string(index=False))

# ── Build outputs ─────────────────────────────────────────────────────────────
bad_patients = set(flagged["patient"].tolist())

# samples to exclude
excl_samples = obs.loc[obs["patient"].isin(bad_patients), "sample_id"].unique().tolist()

# all sample → patient map (from obs)
all_samples = obs[["sample_id","patient"]].drop_duplicates()
valid_samples = all_samples.loc[~all_samples["patient"].isin(bad_patients), "sample_id"].tolist()

print(f"\n  Excluded samples : {len(excl_samples)}  ({excl_samples})")
print(f"  Valid samples    : {len(valid_samples)}")

pd.DataFrame({"sample_id": valid_samples}).to_csv(
    OUT / "validated_samples.csv", index=False)
pd.DataFrame({"sample_id": excl_samples}).to_csv(
    OUT / "excluded_samples.csv", index=False)

# sex corrections for xlsx update
corrections = flagged[["patient","sex","corrected_sex","flag"]].copy()
corrections.rename(columns={"sex":"original_sex"}, inplace=True)
corrections.to_csv(OUT / "sex_corrections.csv", index=False)
print(f"\n  sex_corrections.csv  ({len(corrections)} patients to fix in xlsx)")

# ── Update xlsx ────────────────────────────────────────────────────────────────
try:
    import openpyxl
    df_sheet = pd.read_excel(args.samplesheet)
    df_sheet.columns = [c.strip().lower() for c in df_sheet.columns]

    corr_map = dict(zip(corrections["patient"], corrections["corrected_sex"]))
    patient_col = "patient" if "patient" in df_sheet.columns else None

    if patient_col:
        before = df_sheet["sex"].copy()
        df_sheet["sex"] = df_sheet.apply(
            lambda r: corr_map.get(str(r[patient_col]).strip(), r["sex"]), axis=1)
        n_changed = (before != df_sheet["sex"]).sum()
        out_xlsx = OUT / "ASAP_samplesheet_sexcorrected.xlsx"
        df_sheet.to_excel(out_xlsx, index=False)
        print(f"  Updated xlsx: {n_changed} rows changed → {out_xlsx.name}")
    else:
        print("  WARNING: 'patient' column not found in samplesheet — xlsx not updated")
except ImportError:
    print("  WARNING: openpyxl not available — xlsx not updated (csv corrections saved)")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT — Heatmap per patient (chrY genes + XIST)
# ══════════════════════════════════════════════════════════════════════════════
pal = {"M":"#4C72B0","male":"#4C72B0","F":"#DD8452","female":"#DD8452"}

# per-patient mean per gene
gene_frames = []
for f in h5ads:
    ad_full = sc.read_h5ad(f, backed='r')
    genes   = [g for g in Y_GENES + X_GENE if g in ad_full.var_names]
    ad      = ad_full[:, genes].to_memory()
    ad_full.file.close(); del ad_full

    gdf = pd.DataFrame(index=ad.obs_names)
    for g in genes:
        arr = ad[:, g].X
        if hasattr(arr, "toarray"): arr = arr.toarray()
        gdf[g] = np.asarray(arr).flatten()
    gdf["patient"] = ad.obs["patient"].values
    gene_frames.append(gdf)

gdf_all   = pd.concat(gene_frames, ignore_index=True)
gene_cols = [g for g in Y_GENES + X_GENE if g in gdf_all.columns]
heat_df   = (gdf_all.groupby("patient")[gene_cols]
             .mean()
             .fillna(0).replace([np.inf,-np.inf], 0))

# row annotations
patient_meta = per_patient.set_index("patient")[["sex","flag"]]
row_colors   = pd.DataFrame({
    "sex"    : patient_meta["sex"].map(lambda s: pal.get(str(s).strip(),"grey")),
    "flagged": patient_meta["flag"].map(lambda f: "red" if f != "ok" else "whitesmoke"),
}, index=heat_df.index)

fig_h = max(8, len(heat_df) * 0.22)
fig_w = max(6, len(gene_cols) * 0.9)

g = sns.clustermap(
    heat_df,
    row_colors      = row_colors,
    col_cluster     = False,
    row_cluster     = True,
    cmap            = "YlOrRd",
    figsize         = (fig_w, fig_h),
    xticklabels     = True,
    yticklabels     = True,
    dendrogram_ratio= 0.08,
    cbar_pos        = (0.02, 0.82, 0.025, 0.12),
    linewidths      = 0.3,
    linecolor       = "white",
)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_title("chrY genes + XIST — mean log-norm expression per patient\n"
                        "(blue = M annotated, orange = F annotated, red stripe = flagged)",
                        fontsize=9, pad=10)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=8, rotation=45, ha="right")
g.savefig(OUT / "heatmap_chrY_XIST_per_patient.png", dpi=150, bbox_inches="tight")
plt.close()
print("\n  heatmap_chrY_XIST_per_patient.png")

print(f"\nDONE → {OUT}/")
