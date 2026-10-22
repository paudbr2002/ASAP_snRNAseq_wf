# ASAP snRNA-seq Pipeline

Single-nucleus RNA-seq processing pipeline for the integration of ASAP's snRNAseq cohorts.  
---

## Overview

```
<img src="pipeline_diagram.pdf" alt="Descripción de la imagen">
```

---

## Quick start

```bash
# Default: SoupX + Harmony, all samples
bash run.sh

# Resume after failure
bash run.sh -resume

# Skip ambient correction
bash run.sh --skip-ambient

# scVI integration (requires GPU node)
bash run.sh --scvi

# Males or females only
bash run.sh --males
bash run.sh --females

# Combine flags
bash run.sh --scvi --males -resume
```

---

## Parameters

| Parameter | Default | Options |
|---|---|---|
| `ambient_tool` | `soupx` | `soupx` · `skip` |
| `integration_tool` | `harmony` | `harmony` · `scvi` |
| `sex_filter` | `all` | `all` · `M` · `F` |
| `extra_filter` | *(none)* | `"col:val1,val2;col2:val3"` |
| `max_parallel` | `8` | any int |

Edit `EXTRA_FILTER` directly in `run.sh` to subset by region, condition, etc.

---

## Pipeline steps

### 1 · PREPARE_SAMPLES
Reads the Excel samplesheet, drops `hypoxia == 1`, applies `extra_filter`, and validates that every `matrix_path` exists on disk.  
**Output:** `samples_valid.csv`, `samples_missing.txt`

### 2 · SOUPX *(CPU · rocker_sif · 4 cpu · 24 GB · 2h)*
Ambient RNA correction via [SoupX](https://github.com/constantAmateur/SoupX).  
`load10X()` → `autoEstCont()` → `adjustCounts(roundToInt=TRUE)`.  
Tracks Y-gene counts before/after correction (`y_counts_raw`, `y_counts_corrected`) for mLOY QC.  
Falls back to uncorrected matrix if `raw_feature_bc_matrix` is missing.  
**Output:** `{sid}_soupx/` (MTX), `{sid}_soupx_report.csv`

### 3 · QC_PER_SAMPLE
Per-sample quality control with scanpy + scrublet.

| Filter | Threshold |
|---|---|
| Doublets | scrublet (predicted_doublet) |
| min genes | 500 |
| max MT% | 2.0% |
| Adaptive | mean ± 2 SD on n_genes |
| min nuclei | 50 (else → `.EXCLUDED`) |

Stores raw counts in `adata.layers["counts"]` before normalisation (required by scVI).  
**Output:** `{sample_id}.h5ad` or `{sample_id}.EXCLUDED`

### 4 · VALIDATE_SEX *(CPU · scrnaseq_sif · 2 cpu · 64 GB · 30 min)*
Computes XIST expression and a chrY gene score per donor using 8 Y-genes:  
`RPS4Y1 · DDX3Y · UTY · EIF1AY · KDM5D · ZFY · USP9Y · TMSB4Y`

| Flag | Condition |
|---|---|
| ♂ mislabeled | XIST > Q50(females) |
| ♀ mislabeled | chrY > Q75(males) **and** XIST < Q25(females) |

Flagged donors are excluded before integration.  
**Output:** `validated_samples.csv`, `excluded_samples.csv`, `heatmap_chrY_XIST_per_patient.png`

### 5a · INTEGRATE_HARMONY *(CPU · scrnaseq_sif · 16 cpu · 100 GB · 12h · sens-xl)*
`normalize_total(1e6)` → `log1p` → `HVG(3000)` → `PCA(30)` → Harmony → UMAP → Leiden.  
Batch key: `seqnum`. Sex added as correction key when `sex_filter=all`.

### 5b · INTEGRATE_SCVI *(GPU ⚡ · scvitools_sif · 8 cpu · 64 GB · 6h · sens-gpu)*
scVI variational autoencoder on raw counts (`layers["counts"]`).  
`n_latent=30 · n_layers=1 · n_hidden=128 · NegBinom · dropout=0.1`  
`max_epochs=400 · early_stopping(patience=20) · batch_size=512`  
TF32 disabled (prevents NaN on Ampere GPUs).  
**Output:** `ASAP_integrated.h5ad` + `scvi_model/`

---

## Output structure

```
results_ASAP/
├── tables/
│   └── samples_valid.csv
├── ambient_reports/
│   └── {sid}_soupx_report.csv
├── per_sample_qc/
│   ├── {sid}.h5ad
│   └── {sid}.EXCLUDED
├── sex_validation/
│   ├── validated_samples.csv
│   ├── excluded_samples.csv
│   └── heatmap_chrY_XIST_per_patient.png
├── integration/
│   ├── ASAP_integrated.h5ad
│   ├── scvi_model/          ← only with --scvi
│   ├── plots/
│   └── tables/
├── pipeline_timeline.html
├── pipeline_report.html
└── pipeline_dag.html
```

---

## Containers

| Container | Used by |
|---|---|
| `scrnaseq_sif` | PREPARE_SAMPLES, QC_PER_SAMPLE, VALIDATE_SEX, HARMONY |
| `rocker_sif` | SOUPX |
| `scvitools_sif` | INTEGRATE_SCVI |

Paths defined in `nextflow.config`. Singularity is run with `--nv` so GPU containers can access the host NVIDIA driver.

---

## Notes

- **mLOY study:** SoupX is strongly recommended. Sample ASAP106 (♂, rho=41.7%) had 1 340 raw chrY counts corrected to 321 — without SoupX it would have passed sex QC incorrectly.  
- The `extra_filter` syntax is `"column:val1,val2;column2:val3"` and is applied on top of the built-in hypoxia filter.  
- Only `INTEGRATE_SCVI` requires a GPU node. All other steps run on CPU.
