#!/usr/bin/env python3
"""
prepare_samples.py
──────────────────────────────────────────────────────────────────────────────
Reads the ASAP samplesheet, applies filters, validates paths, and derives
raw_matrix_path for CellBender (raw_feature_bc_matrix).

Filters applied (in order)
  1. Hypoxia filter  — rows where hypoxia == 1 are dropped (always active).
  2. Extra filters   — optional CLI argument.
     Syntax: "col1:val1,val2;col2:val3"  (case-insensitive)

Derived columns added to the output CSV
  raw_matrix_path — path to raw_feature_bc_matrix (.h5 preferred, then dir).
                    Required by CellBender; empty string if not found.
"""
import argparse
from pathlib import Path
import pandas as pd


# ── Args ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--samplesheet")
parser.add_argument("--out_csv")
parser.add_argument("--out_missing")
parser.add_argument("--extra_filter", default="",
    help="Semicolon-separated col:val1,val2 pairs. Example: 'region:SN,putamen;condition:PD'")
args = parser.parse_args()


# ── Helper: parse extra_filter string ────────────────────────────────────────
def parse_filter(s: str) -> dict:
    if not s or not s.strip():
        return {}
    out = {}
    for part in s.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        col, vals = part.split(":", 1)
        out[col.strip().lower()] = [v.strip().lower() for v in vals.split(",")]
    return out


# ── Helper: find raw matrix path for CellBender ───────────────────────────────
def get_bam_path(matrix_path_str: str) -> str:
    """
    Derive BAM path for Souporcell from filtered_feature_bc_matrix path.
    Prefers possorted_genome_bam.bam (always indexed by CellRanger).
    Falls back to cellsorted version.
    Returns empty string if not found (Souporcell will fail with clear error).
    """
    outs_dir = Path(matrix_path_str).parent
    for name in [
        "possorted_genome_bam.bam",            # standard CellRanger, always indexed
        "cellsorted_possorted_genome_bam.bam",  # re-sorted by barcode
    ]:
        p = outs_dir / name
        if p.exists():
            return str(p)
    return ""


def get_raw_matrix_path(matrix_path_str: str) -> str:
    """
    Derive raw_feature_bc_matrix path from filtered_feature_bc_matrix.
    CellBender prefers the HDF5 file; falls back to the MTX directory.
    Returns empty string if neither is found (will trigger a clear error
    when CellBender tries to run).
    """
    outs_dir = Path(matrix_path_str).parent   # .../outs/
    h5  = outs_dir / "raw_feature_bc_matrix.h5"
    mtx = outs_dir / "raw_feature_bc_matrix"
    if h5.exists():  return str(h5)
    if mtx.exists(): return str(mtx)
    return ""


# ── Load samplesheet ──────────────────────────────────────────────────────────
df = pd.read_excel(args.samplesheet)
df.columns = [c.strip().lower() for c in df.columns]
df.rename(columns={"sample": "sample_id", "dx": "condition"}, inplace=True)

if "sample_id" not in df.columns:
    raise ValueError("Samplesheet must have a 'sample_id' (or 'sample') column.")
if "matrix_path" not in df.columns:
    raise ValueError(
        "Samplesheet must have an explicit 'matrix_path' column pointing to "
        "each sample's filtered_feature_bc_matrix directory."
    )

df["sample_id"]   = df["sample_id"].astype(str).str.strip()
df["matrix_path"] = df["matrix_path"].astype(str).str.strip()

# ── 1. Hypoxia filter ─────────────────────────────────────────────────────────
if "hypoxia" in df.columns:
    n = len(df)
    df = df[df["hypoxia"].isna() | (df["hypoxia"] == 0)].reset_index(drop=True)
    print(f"Hypoxia filter : {n} → {len(df)} samples")
else:
    print("Hypoxia filter : column not found, skipping")

# ── 2. Extra filters ──────────────────────────────────────────────────────────
for col, allowed in parse_filter(args.extra_filter).items():
    if col not in df.columns:
        print(f"  WARNING: extra_filter column '{col}' not in samplesheet — skipped")
        continue
    n = len(df)
    df = df[df[col].astype(str).str.strip().str.lower().isin(allowed)].reset_index(drop=True)
    print(f"  Extra filter [{col} ∈ {allowed}] : {n} → {len(df)} samples")

# ── Standardise condition labels ──────────────────────────────────────────────
if "condition" in df.columns:
    _map = {"ctl": "Control", "ctrl": "Control", "control": "Control", "pd": "PD"}
    df["condition"] = df["condition"].astype(str).str.lower().map(_map).fillna(df["condition"])

# ── 3. Derive raw_matrix_path (CellBender) and bam_path (Souporcell) ─────────
df["raw_matrix_path"] = df["matrix_path"].apply(get_raw_matrix_path)
n_no_raw = (df["raw_matrix_path"] == "").sum()
if n_no_raw:
    print(f"  WARNING: {n_no_raw} samples have no raw_feature_bc_matrix "
          f"(CellBender will fail for those; SoupX / skip unaffected)")

df["bam_path"] = df["matrix_path"].apply(get_bam_path)
n_no_bam = (df["bam_path"] == "").sum()
if n_no_bam:
    print(f"  WARNING: {n_no_bam} samples have no BAM file "
          f"(Souporcell will fail for those)")

# ── 4. Validate filtered matrix paths ────────────────────────────────────────
found, missing = [], []
for _, row in df.iterrows():
    sid = row["sample_id"]
    mp  = Path(row["matrix_path"])
    if mp.exists():
        found.append(sid)
    else:
        missing.append(sid)
        print(f"  MISSING matrix_path: {mp}  [{sid}]")

valid = df[df["sample_id"].isin(found)].copy()
valid.to_csv(args.out_csv, index=False)
with open(args.out_missing, "w") as fh:
    fh.write("\n".join(missing))

print(f"\nValid samples  : {len(found)}")
print(f"Missing paths  : {len(missing)}")
