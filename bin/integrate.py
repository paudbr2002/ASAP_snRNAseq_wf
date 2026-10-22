#!/usr/bin/env python3
"""
integrate.py — Harmony-based integration (CPU)
──────────────────────────────────────────────────────────────────────────────
Load per-sample h5ads (already filtered to verified samples by the Nextflow
channel), apply optional sex filter, then:
  Merge → normalize → HVG → PCA → Harmony → PAGA → UMAP → Leiden

sex_filter  : "all" | "M" | "F"
  Applied after merge.  The raw h5ads staged here may contain both sexes;
  this parameter restricts the integration to one sex.
  Use "all" to integrate all available donors (default).
"""
import argparse, os, warnings

# ── Thread config before any numba/scanpy import ─────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument("--h5ad_dir",   default=".")
_p.add_argument("--outdir",     default=".")
_p.add_argument("--sex_filter", default="all",
                choices=["all", "M", "F"],
                help="Integrate 'all' donors, males ('M') or females ('F') only")
args = _p.parse_args()

n_cores = int(os.environ.get("SLURM_CPUS_PER_TASK",
              os.environ.get("SLURM_NTASKS", 8)))
for v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
          "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"]:
    if v not in os.environ:
        os.environ[v] = str(n_cores)

# ── Imports ───────────────────────────────────────────────────────────────────
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

try:
    import harmonypy as hm
except ImportError:
    raise ImportError("harmonypy not found — pip install harmonypy in the container")

warnings.filterwarnings("ignore")
sc.settings.n_jobs    = n_cores
sc.settings.verbosity = 1

PLOTS  = Path(args.outdir) / "plots";  PLOTS.mkdir(parents=True, exist_ok=True)
TABLES = Path(args.outdir) / "tables"; TABLES.mkdir(parents=True, exist_ok=True)

Y_GENES = ["AMELY","NLGN4Y","PCDH11Y","SRY","RPS4Y1","TBL1Y","DDX3Y","UTY",
           "EIF1AY","KDM5D","ZFY","USP9Y","TMSB4Y"]

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD + MERGE
# ══════════════════════════════════════════════════════════════════════════════
h5ads = sorted(Path(args.h5ad_dir).glob("*.h5ad"))
print(f"\n-- Loading {len(h5ads)} h5ad files --")
adatas = {}
for f in h5ads:
    ad = sc.read_h5ad(f)
    adatas[f.stem] = ad
    print(f"  [{f.stem}]  {ad.n_obs:,} nuclei")

print("\n-- Merging (join=outer) --")
adata = sc.concat(list(adatas.values()), join="outer", merge="same")
adata.obs_names_make_unique()
sc.pp.filter_genes(adata, min_cells=3)
print(f"  Total: {adata.n_obs:,} nuclei × {adata.n_vars:,} genes")

# ── Sex filter ────────────────────────────────────────────────────────────────
if args.sex_filter != "all":
    if "sex" not in adata.obs.columns:
        print(f"  WARNING: 'sex' column not found in obs — sex_filter ignored")
    else:
        sex_keep = adata.obs["sex"].str.strip().isin(
            ["M", "male"] if args.sex_filter == "M" else ["F", "female"]
        )
        n_before = adata.n_obs
        adata = adata[sex_keep].copy()
        print(f"  Sex filter ({args.sex_filter}): {n_before:,} → {adata.n_obs:,} nuclei")
        if adata.n_obs == 0:
            raise ValueError(
                f"No nuclei remain after sex_filter='{args.sex_filter}'. "
                "Check that the 'sex' column in your samplesheet uses 'M'/'F'."
            )

# ══════════════════════════════════════════════════════════════════════════════
#  NORMALISE + HVG
# ══════════════════════════════════════════════════════════════════════════════
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, batch_key="seqnum",
    flavor="seurat", subset=False,
)

# ══════════════════════════════════════════════════════════════════════════════
#  PCA  (on HVGs only)
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- PCA + Harmony --")
hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(hvg, max_value=10)
sc.tl.pca(hvg, n_comps=30, svd_solver="arpack", random_state=42)
adata.obsm["X_pca"] = hvg.obsm["X_pca"]
adata.uns["pca"]    = hvg.uns["pca"]
del hvg

# ── Pre-Harmony UMAP (stored for comparison) ──────────────────────────────────
sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=20, random_state=42)
sc.tl.leiden(adata, resolution=0.1, key_added="leiden_pre",
             directed=False, flavor="igraph", random_state=42)
sc.tl.paga(adata, groups="leiden_pre")
sc.pl.paga(adata, show=False)
sc.tl.umap(adata, init_pos="paga", random_state=42)
adata.obsm["X_umap_pre_harmony"] = adata.obsm["X_umap"].copy()

# ── Harmony (direct harmonypy call — avoids scanpy wrapper shape bug) ─────────
keys = [k for k in ["seqnum", "sex", "patient"] if k in adata.obs.columns]
if args.sex_filter != "all" and "sex" in keys:
    keys.remove("sex")   # constant after filtering — would cause Harmony to error
print(f"  Harmony keys: {keys}")

pca_mat = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float64)
ho = hm.run_harmony(pca_mat, adata.obs, keys, random_state=42)

Z = np.asarray(ho.Z_corr)
print(f"  Harmony Z_corr shape: {Z.shape}")
if   Z.ndim == 2 and Z.shape[0] == adata.n_obs: adata.obsm["X_pca_harmony"] = Z
elif Z.ndim == 2 and Z.shape[1] == adata.n_obs: adata.obsm["X_pca_harmony"] = Z.T
else: raise ValueError(f"Unexpected harmony shape {Z.shape} for n_obs={adata.n_obs}")

# ── Post-Harmony UMAP + Leiden ─────────────────────────────────────────────────
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=20,
                random_state=42)
sc.tl.leiden(adata, resolution=0.1, key_added="leiden",
             directed=False, flavor="igraph", random_state=42)
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, show=False)
sc.tl.umap(adata, init_pos="paga", random_state=42)
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- Generating plots --")

def gene_score(ad, genes):
    present = [g for g in genes if g in ad.var_names]
    if not present:
        return pd.Series(0.0, index=ad.obs_names)
    return pd.Series(
        np.asarray(ad[:, present].X.sum(axis=1)).flatten(),
        index=ad.obs_names,
    )

adata.obs["chrY_score"] = gene_score(adata, Y_GENES)
adata.obs["XIST_score"] = gene_score(adata, ["XIST"])

# 1. XIST vs chrY scatter
fig, ax = plt.subplots(figsize=(6, 5))
palette = {"M": "#4C72B0", "male": "#4C72B0", "F": "#DD8452", "female": "#DD8452"}
if "sex" in adata.obs.columns:
    for sex, grp in adata.obs.groupby("sex"):
        c = palette.get(str(sex).strip(), "grey")
        ax.scatter(grp["chrY_score"], grp["XIST_score"],
                   s=1, alpha=0.2, color=c, label=sex, rasterized=True)
    ax.legend(markerscale=6, title="metadata sex", framealpha=0.7)
ax.set_xlabel("chrY score"); ax.set_ylabel("XIST score")
ax.set_title("Sex validation: XIST vs chrY")
fig.tight_layout()
fig.savefig(PLOTS / "sex_validation_XIST_vs_chrY.png", dpi=150, bbox_inches="tight")
plt.close()

# 2. chrY score per male sample
male_mask = adata.obs.get("sex", pd.Series(dtype=str)).str.strip().isin(["M","male"])
if male_mask.sum() > 0:
    male_obs = adata.obs.loc[male_mask].copy()
    order = (male_obs.groupby("sample_id")["chrY_score"]
             .median().sort_values().index.tolist())
    fig, ax = plt.subplots(figsize=(max(6, len(order) * 0.35), 5))
    sns.boxplot(data=male_obs, x="sample_id", y="chrY_score",
                order=order, ax=ax, flierprops={"marker": ".", "markersize": 2},
                color="#4C72B0")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
    ax.set_title("chrY score per male sample (low = candidate LOY)")
    fig.tight_layout()
    fig.savefig(PLOTS / "loy_chrY_score_per_sample.png", dpi=150, bbox_inches="tight")
    plt.close()

# 3. UMAPs pre/post Harmony
import shutil
umap_keys = [k for k in ["seqnum", "condition", "sex", "leiden_pre", "leiden"]
             if k in adata.obs.columns]

adata.obsm["_umap_bak"]  = adata.obsm["X_umap"].copy()
adata.obsm["X_umap"]     = adata.obsm["X_umap_pre_harmony"]
for key in umap_keys:
    sc.pl.umap(adata, color=key, show=False,
               title=f"pre-Harmony — {key}", save=f"_pre_{key}.png")

adata.obsm["X_umap"] = adata.obsm.pop("_umap_bak")
for key in umap_keys:
    sc.pl.umap(adata, color=key, show=False,
               title=f"post-Harmony — {key}", save=f"_post_{key}.png")

for png in Path(sc.settings.figdir).glob("umap_*.png"):
    shutil.move(str(png), str(PLOTS / png.name))
print(f"  UMAPs: pre/post-Harmony × {len(umap_keys)} keys")

# 4. QC summary table
qc_cols = [c for c in ["sample_id","condition","sex","seqnum","patient",
                        "braak","age","pmi"] if c in adata.obs.columns]
summary = (
    adata.obs.groupby("sample_id")
    .agg(
        n_nuclei     = ("sample_id", "count"),
        median_genes = ("n_genes_by_counts", "median"),
        pct_mt_mean  = ("pct_counts_mt",     "mean"),
        chrY_median  = ("chrY_score",         "median"),
        XIST_median  = ("XIST_score",         "median"),
    )
    .reset_index()
)
meta_tbl = adata.obs[qc_cols].drop_duplicates("sample_id").set_index("sample_id")
summary  = summary.set_index("sample_id").join(meta_tbl).reset_index()
summary.to_csv(TABLES / "integration_summary.csv", index=False)

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE
# ══════════════════════════════════════════════════════════════════════════════
out_h5 = Path(args.outdir) / "ASAP_integrated.h5ad"
adata.write_h5ad(out_h5)
print(f"\nDONE (Harmony): {adata.n_obs:,} nuclei × {adata.n_vars:,} genes → {out_h5}")
