#!/usr/bin/env python3
"""
integrate_scvi.py — scVI-tools GPU integration
──────────────────────────────────────────────────────────────────────────────
Drop-in replacement for integrate.py that uses a variational autoencoder
(scVI, Lopez et al. 2018) for batch correction instead of Harmony.

scVI advantages for this project
  • Probabilistic model (negative binomial) is better suited to raw UMI counts
    than Harmony's PCA-residual approach.
  • Learns a non-linear latent space, which can resolve subtle cell-type
    distinctions better than linear PCA.
  • GPU acceleration makes the runtime comparable to Harmony on CPU for
    datasets of this size (~100-300k nuclei).
  • The trained model can be reused for reference mapping of future datasets.

Batch correction
  batch_key="seqnum"  — corrects sequencing batch effects
  continuous_covariate_keys=["age","pmi"]  — regresses out donor covariates
    that are not biologically interesting (added if columns present)
  sex is intentionally NOT a batch key — it is the biological variable
    of interest; use sex_filter to subset instead.

Outputs
  ASAP_integrated.h5ad  — annotated object with X_scVI, X_umap, leiden in obs
  scvi_model/           — trained model (can be reloaded with scvi.model.SCVI.load)
  plots/                — same plots as the Harmony script
  tables/               — integration_summary.csv

Usage
  integrate_scvi.py --h5ad_dir . --outdir . --sex_filter M
"""
import argparse, os, shutil, warnings

# ── Thread config before any numba/scanpy import ─────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument("--h5ad_dir",   default=".")
_p.add_argument("--outdir",     default=".")
_p.add_argument("--sex_filter", default="all",
                choices=["all", "M", "F"])
_p.add_argument("--n_latent",   type=int, default=30)
_p.add_argument("--n_top_genes",type=int, default=3000)
_p.add_argument("--max_epochs", type=int, default=400)
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
    import scvi
except ImportError:
    raise ImportError(
        "scvi-tools not found.  Make sure you are running inside "
        "scvitools_1.5.1-gpu.sif (or install: pip install scvi-tools)."
    )

warnings.filterwarnings("ignore")
sc.settings.n_jobs    = n_cores
sc.settings.verbosity = 1
scvi.settings.seed    = 42

PLOTS  = Path(args.outdir) / "plots";  PLOTS.mkdir(parents=True, exist_ok=True)
TABLES = Path(args.outdir) / "tables"; TABLES.mkdir(parents=True, exist_ok=True)
MODEL  = Path(args.outdir) / "scvi_model"

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
        print("  WARNING: 'sex' column not found — sex_filter ignored")
    else:
        sex_keep = adata.obs["sex"].str.strip().isin(
            ["M","male"] if args.sex_filter == "M" else ["F","female"]
        )
        n_before = adata.n_obs
        adata = adata[sex_keep].copy()
        print(f"  Sex filter ({args.sex_filter}): {n_before:,} → {adata.n_obs:,} nuclei")
        if adata.n_obs == 0:
            raise ValueError(
                f"No nuclei remain after sex_filter='{args.sex_filter}'. "
                "Check the 'sex' column in your samplesheet."
            )

# ── Sanitize counts for scVI ──────────────────────────────────────────────────
# scVI needs non-negative integer-like counts (NB distribution).
# SoupX and CellBender return float counts that can be slightly negative
# after ambient correction — log(neg) or log(0) in the NB encoder = NaN.
# Fix: round to integers and clip negatives to 0.
import scipy.sparse as _sp

if "counts" not in adata.layers:
    print("  WARNING: 'counts' layer missing — using adata.X as raw counts")
    adata.layers["counts"] = adata.X.copy()

X = adata.layers["counts"]
if _sp.issparse(X):
    X.data = np.round(X.data).clip(0).astype(np.float32)
    X.eliminate_zeros()
else:
    adata.layers["counts"] = np.round(X).clip(0).astype(np.float32)

n_neg = int((X.data < 0).sum()) if _sp.issparse(X) else int((adata.layers["counts"] < 0).sum())
print(f"  Counts sanitised — negative values clipped: {n_neg}")

# ── Disable TF32 on Ampere GPUs ───────────────────────────────────────────────
# Ampere (A40, A100) uses TF32 (19-bit mantissa) by default for matmul.
# This trades precision for speed and can cause NaN in the scVI encoder
# on some datasets. Disable it to use full fp32.
import torch as _torch
_torch.backends.cuda.matmul.allow_tf32 = False
_torch.backends.cudnn.allow_tf32      = False
print("  TF32 disabled (full fp32 matmul)")

# ══════════════════════════════════════════════════════════════════════════════
#  HVG SELECTION  (seurat_v3 on raw counts)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n-- Selecting top {args.n_top_genes} HVGs (seurat_v3, batch-aware) --")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes   = args.n_top_genes,
    flavor        = "seurat_v3",   # uses raw counts (variance-stabilizing)
    batch_key     = "seqnum" if "seqnum" in adata.obs.columns else None,
    layer         = "counts",
    subset        = True,          # subset adata to HVGs in-place
)
print(f"  Retained {adata.n_vars:,} HVGs")

# Log-normalise for UMAP/plotting (scVI will use the raw counts layer)
adata.layers["lognorm"] = sc.pp.normalize_total(adata, target_sum=1e4,
                                                  layer="counts", copy=True).X
sc.pp.log1p(adata, layer="lognorm")

# ══════════════════════════════════════════════════════════════════════════════
#  scVI SETUP + TRAINING
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- Setting up scVI model --")

# Continuous covariates: regress out donor-level confounders that are not
# biologically interesting.  Only add if the column is present and non-constant.
cont_keys = []
for col in ["age", "pmi"]:
    if col in adata.obs.columns:
        vals = pd.to_numeric(adata.obs[col], errors="coerce")
        n_nan = vals.isna().sum()
        if n_nan > 0:
            # NaN covariates propagate directly through linear layers → NaN latent
            # Fill with median so the column is still usable as a covariate
            fill = vals.median()
            adata.obs[col] = vals.fillna(fill)
            print(f"  Covariate '{col}': filled {n_nan} NaN with median {fill:.1f}")
        if vals.dropna().nunique() > 1:
            cont_keys.append(col)
if cont_keys:
    print(f"  Continuous covariates: {cont_keys}")
else:
    print("  No continuous covariates (age/pmi missing or constant)")

batch_key = "seqnum" if "seqnum" in adata.obs.columns else None
if batch_key and adata.obs[batch_key].nunique() < 2:
    print(f"  WARNING: batch_key '{batch_key}' has only one unique value — "
          "scVI will run without batch correction")
    batch_key = None

scvi.model.SCVI.setup_anndata(
    adata,
    layer                      = "counts",
    batch_key                  = batch_key,
    continuous_covariate_keys  = cont_keys if cont_keys else None,
)

model = scvi.model.SCVI(
    adata,
    n_latent         = args.n_latent,
    n_layers         = 1,     # 2→1: simpler encoder is more stable
    n_hidden         = 128,   # 256→128: reduces risk of exploding activations
    gene_likelihood  = "nb",
    dropout_rate     = 0.1,
)
print(f"  Model: n_latent={args.n_latent}, n_layers=1, n_hidden=128, gene_likelihood=nb")
print(model)

print(f"\n-- Training (max_epochs={args.max_epochs}, early_stopping=True) --")
# NaN explosion fixes:
#   precision="32"        — disables fp16/bf16 mixed precision which causes NaN
#                           explosions on some datasets; fp32 is slower but stable
#   gradient_clip_val=1.0 — clips exploding gradients before they produce NaN
#   lr=1e-4               — conservative learning rate (default 1e-3 can diverge)
#   batch_size=512        — smaller batches improve gradient estimate quality
model.train(
    max_epochs              = args.max_epochs,
    early_stopping          = True,
    early_stopping_patience = 20,
    batch_size              = 512,
    plan_kwargs             = {"lr": 1e-4},
    precision               = "32",        # passed directly to Lightning Trainer
    gradient_clip_val       = 1.0,         # idem — not wrapped in trainer_kwargs
)
print("  Training complete.")

# ── Latent representation ─────────────────────────────────────────────────────
adata.obsm["X_scVI"] = model.get_latent_representation()

# Save the trained model for future reference-mapping
model.save(str(MODEL), overwrite=True)
print(f"  Model saved → {MODEL}/")

# ══════════════════════════════════════════════════════════════════════════════
#  NEIGHBOURS + UMAP + LEIDEN  (on scVI latent space)
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- Neighbours / UMAP / Leiden on X_scVI --")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, random_state=42)
sc.tl.leiden(adata, resolution=0.1, key_added="leiden",
             directed=False, flavor="igraph", random_state=42)
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, show=False)
sc.tl.umap(adata, init_pos="paga", random_state=42)
print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOTS  (same set as the Harmony script)
# ══════════════════════════════════════════════════════════════════════════════
print("\n-- Generating plots --")

def gene_score(ad, genes, layer="lognorm"):
    present = [g for g in genes if g in ad.var_names]
    if not present:
        return pd.Series(0.0, index=ad.obs_names)
    if layer and layer in ad.layers:
        X = ad[:, present].layers[layer]
    else:
        X = ad[:, present].X
    return pd.Series(np.asarray(X.sum(axis=1)).flatten(), index=ad.obs_names)

adata.obs["chrY_score"] = gene_score(adata, Y_GENES)
adata.obs["XIST_score"] = gene_score(adata, ["XIST"])

# 1. XIST vs chrY scatter
palette = {"M":"#4C72B0","male":"#4C72B0","F":"#DD8452","female":"#DD8452"}
fig, ax = plt.subplots(figsize=(6, 5))
if "sex" in adata.obs.columns:
    for sex, grp in adata.obs.groupby("sex"):
        c = palette.get(str(sex).strip(), "grey")
        ax.scatter(grp["chrY_score"], grp["XIST_score"],
                   s=1, alpha=0.2, color=c, label=sex, rasterized=True)
    ax.legend(markerscale=6, title="metadata sex", framealpha=0.7)
ax.set_xlabel("chrY score"); ax.set_ylabel("XIST score")
ax.set_title("Sex validation: XIST vs chrY (post-scVI)")
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
                order=order, ax=ax, flierprops={"marker":".", "markersize":2},
                color="#4C72B0")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
    ax.set_title("chrY score per male sample (low = candidate LOY)")
    fig.tight_layout()
    fig.savefig(PLOTS / "loy_chrY_score_per_sample.png", dpi=150, bbox_inches="tight")
    plt.close()

# 3. UMAPs
umap_keys = [k for k in ["seqnum","condition","sex","leiden"] if k in adata.obs.columns]
for key in umap_keys:
    sc.pl.umap(adata, color=key, show=False,
               title=f"scVI — {key}", save=f"_scvi_{key}.png")
for png in Path(sc.settings.figdir).glob("umap_*.png"):
    shutil.move(str(png), str(PLOTS / png.name))
print(f"  UMAPs: {len(umap_keys)} keys")

# 4. ELBO training history
if hasattr(model, "history") and "elbo_train" in model.history:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(model.history["elbo_train"],   label="train", lw=1.5)
    if "elbo_validation" in model.history:
        ax.plot(model.history["elbo_validation"], label="validation", lw=1.5, ls="--")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO")
    ax.set_title("scVI training history"); ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "scvi_training_history.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  scvi_training_history.png")

# 5. QC summary table
qc_cols = [c for c in ["sample_id","condition","sex","seqnum","patient",
                        "braak","age","pmi"] if c in adata.obs.columns]
summary = (
    adata.obs.groupby("sample_id")
    .agg(
        n_nuclei     = ("sample_id",         "count"),
        median_genes = ("n_genes_by_counts",  "median"),
        pct_mt_mean  = ("pct_counts_mt",      "mean"),
        chrY_median  = ("chrY_score",          "median"),
        XIST_median  = ("XIST_score",          "median"),
    )
    .reset_index()
)
meta_tbl = adata.obs[qc_cols].drop_duplicates("sample_id").set_index("sample_id")
summary  = summary.set_index("sample_id").join(meta_tbl).reset_index()
summary.to_csv(TABLES / "integration_summary.csv", index=False)
print("  integration_summary.csv")

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE
# ══════════════════════════════════════════════════════════════════════════════
out_h5 = Path(args.outdir) / "ASAP_integrated.h5ad"
adata.write_h5ad(out_h5)
print(f"\nDONE (scVI): {adata.n_obs:,} nuclei × {adata.n_vars:,} genes → {out_h5}")
