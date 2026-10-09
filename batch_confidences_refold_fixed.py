#!/usr/bin/env python3
"""
batch_confidences_refold.py

Adapted from filter_binders.py for the MEGF8 refold batch. See module
docstring history for the flat-directory / per-output rfd_run-derivation
background.

RMSD METHOD (updated):
  Binder-only self-RMSD (aligning chain A onto chain A directly) doesn't
  capture what actually matters here: whether the binder is positioned
  correctly RELATIVE TO the target. It's also not directly usable as-is,
  because the RF3 output's target (chain B) is C-terminally truncated
  (abbreviated target) relative to the RFD backbone's target (chain B,
  full-length) - same start, RF3's is just shorter.

  So instead:
    1. Take the RFD backbone's target CA atoms (chain B, full-length),
       truncate to the RF3 model's target length (first N - the
       truncation is C-terminal, so this is a straightforward N-terminal
       matching slice, not a sequence-alignment problem).
    2. Superimpose RF3's target CA onto that truncated backbone target
       CA - this solves for the rigid-body transform that puts RF3's
       target into the RFD reference frame.
    3. Apply that SAME transform (not a new fit) to RF3's binder chain
       (chain A).
    4. RMSD between the (untouched) RFD backbone binder CA and the
       transformed RF3 binder CA is the reported "ca_rmsd_A" - this
       reflects binder placement error given a fixed target frame,
       rather than binder shape error alone.

  ASSUMPTION: binder is chain A and target is chain B in BOTH the RFD
  backbone and the RF3 output. This hasn't been independently confirmed
  for the RF3 side's target chain letter - if wrong, this will surface
  as either a "no CA atoms found for chain B" error or a wildly-off CA
  count mismatch (not a subtle few-residue discrepancy), so it should
  fail loudly rather than silently give a wrong number.
"""

import argparse
import csv
import gzip
import io
import json
import os
import time
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

if hasattr(pdbx, "CIFFile"):
    PDBX_FILE_CLS = pdbx.CIFFile
else:
    PDBX_FILE_CLS = pdbx.PDBxFile


# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--rf3_outputs_dir", required=True,
                     help="Flat directory containing one subdirectory per "
                          "RF3 refold output (replaces FOLD_RUNS_DIR)")
args = parser.parse_args()

RF3_OUTPUTS_DIR = Path(args.rf3_outputs_dir)
RFD_RUNS_DIR = Path("/mmfs1/gscratch/stf/igem/all_out/all_bookkeep/rfd_runs")

REFOLD_DIR = Path("/mmfs1/gscratch/stf/igem/all_out/megf8_refold")
OUT_CSV = REFOLD_DIR / "batch_results_confidences_refold_posewise.csv"
CACHE_FILE = REFOLD_DIR / ".batch_results_confidences_refold_posewise_cache.json"

PROGRESS_EVERY = 50
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"  # assumed same letter in both RFD and RF3 output - unverified on RF3 side

FIELDNAMES = [
    "rf3_output", "run_id", "backbone_file", "rfd_run",
    "min_pae_cross", "ca_rmsd_A", "n_ca_aligned", "n_target_ca_aligned",
    "status", "rf3_model_file",
]


# --------------------------------------------------------------------------
# run_id extraction (same repeat-match logic as megf8_refold_jsoner.py)
# --------------------------------------------------------------------------

def get_run_id(rf3_folder_name):
    after_rfd = rf3_folder_name.split("rfd_", 1)[1]
    marker = "megf8_"
    start = 0
    while True:
        idx = after_rfd.find(marker, start)
        if idx == -1:
            raise ValueError(f"Could not determine run_id for: {rf3_folder_name}")
        candidate = after_rfd[:idx].rstrip("_")
        remainder = after_rfd[idx + len(marker):]
        if remainder == candidate or remainder.startswith(candidate + "_"):
            return candidate
        start = idx + len(marker)


# --------------------------------------------------------------------------
# mmCIF loading
# --------------------------------------------------------------------------

def _read_cif_file(path):
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            text = f.read()
        return PDBX_FILE_CLS.read(io.StringIO(text))
    return PDBX_FILE_CLS.read(path)


def load_structure(path):
    """Load the full parsed structure (all chains, all atoms) from a .cif
    or .cif.gz path."""
    cif_file = _read_cif_file(path)
    return pdbx.get_structure(cif_file, model=1)


def get_ca(structure, chain_id):
    """Return the CA atoms of a single chain from an already-loaded
    structure."""
    mask = (
        (structure.atom_name == "CA")
        & struc.filter_amino_acids(structure)
        & (structure.chain_id == chain_id)
    )
    return structure[mask]


# --------------------------------------------------------------------------
# minPAE
# --------------------------------------------------------------------------

def cross_chain_min_pae_from_summary(summary_json_path):
    with open(summary_json_path) as f:
        summary = json.load(f)

    matrix = summary.get("chain_pair_pae_min")
    if matrix is None:
        raise ValueError("chain_pair_pae_min missing from summary_confidences.json")

    values = [v for row in matrix for v in row if v is not None]
    if len(values) == 0:
        raise ValueError("chain_pair_pae_min has no non-null entries")
    if len(values) > 1:
        raise ValueError(
            f"chain_pair_pae_min has {len(values)} non-null entries "
            f"(>2 chains?) - cannot unambiguously determine binder-vs-target minPAE"
        )
    return float(values[0])


# --------------------------------------------------------------------------
# Target-superimposed binder-pose RMSD
# --------------------------------------------------------------------------

def target_superimposed_binder_rmsd(backbone_binder_ca, backbone_target_ca,
                                     model_binder_ca, model_target_ca):
    """Superimpose model target onto (truncated) backbone target, apply
    that transform to the model binder, and return RMSD against the
    backbone binder (untouched, reference frame). Returns
    (rmsd, n_binder_aligned, n_target_aligned)."""

    if len(model_target_ca) == 0:
        raise ValueError("no target CA atoms found in RF3 model")
    if len(backbone_target_ca) == 0:
        raise ValueError("no target CA atoms found in RFD backbone")

    n_target = len(model_target_ca)
    if n_target > len(backbone_target_ca):
        raise ValueError(
            f"RF3 target ({n_target} CA) longer than RFD backbone target "
            f"({len(backbone_target_ca)} CA) - truncation assumption violated"
        )
    backbone_target_ca_matched = backbone_target_ca[:n_target]

    _, transform = struc.superimpose(backbone_target_ca_matched, model_target_ca)
    fitted_binder_ca = transform.apply(model_binder_ca)

    if len(backbone_binder_ca) != len(fitted_binder_ca):
        raise ValueError(
            f"binder CA count mismatch: RFD {len(backbone_binder_ca)} vs "
            f"RF3 {len(fitted_binder_ca)}"
        )
    if len(backbone_binder_ca) == 0:
        raise ValueError("no binder CA atoms found")

    rmsd = struc.rmsd(backbone_binder_ca, fitted_binder_ca)
    return float(rmsd), len(backbone_binder_ca), n_target


# --------------------------------------------------------------------------
# Backbone discovery and matching (scoped per-output via derived rfd_run)
# --------------------------------------------------------------------------

def build_backbone_index_for_run(rfd_run_dir):
    index = {}
    for f in sorted(rfd_run_dir.glob("*.cif.gz")):
        stem = f.name[:-3]  # strip ".gz", keep ".cif"
        index[stem] = f
    return index


def match_backbone(rf3_folder_name, backbone_index):
    matches = [stem for stem in backbone_index if rf3_folder_name.startswith(stem)]
    if not matches:
        return None
    best = max(matches, key=len)
    return backbone_index[best]


# --------------------------------------------------------------------------
# RF3 output discovery (flat directory - no fold_run nesting)
# --------------------------------------------------------------------------

def find_rf3_outputs(rf3_outputs_dir):
    return sorted(p for p in rf3_outputs_dir.iterdir() if p.is_dir())


def find_rf3_files(folder):
    model_matches = list(folder.glob("*_model.cif"))
    summary_matches = list(folder.glob("*_summary_confidences.json"))
    model_cif = model_matches[0] if len(model_matches) == 1 else None
    summary_json = summary_matches[0] if len(summary_matches) == 1 else None
    return model_cif, summary_json


# --------------------------------------------------------------------------
# Cache
# --------------------------------------------------------------------------

def load_cache():
    if not CACHE_FILE.exists():
        return {}
    try:
        with open(CACHE_FILE) as f:
            return json.load(f)
    except Exception as e:
        print(f"WARNING: could not read cache file ({e}); starting with empty cache")
        return {}


def save_cache(cache):
    tmp_path = str(CACHE_FILE) + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump(cache, f)
    os.replace(tmp_path, CACHE_FILE)
    print(f"Cache written to {CACHE_FILE} ({len(cache)} cached entries)")


def mtime_signature(model_cif, summary_json):
    try:
        return max(os.path.getmtime(model_cif), os.path.getmtime(summary_json))
    except OSError:
        return None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def process():
    if not RF3_OUTPUTS_DIR.is_dir():
        raise SystemExit(f"RF3 outputs dir not found: {RF3_OUTPUTS_DIR}")
    if not RFD_RUNS_DIR.is_dir():
        raise SystemExit(f"RFD runs dir not found: {RFD_RUNS_DIR}")

    rf3_outputs = find_rf3_outputs(RF3_OUTPUTS_DIR)
    if not rf3_outputs:
        raise SystemExit(f"No RF3 output subdirectories found under {RF3_OUTPUTS_DIR}")

    os.makedirs(REFOLD_DIR, exist_ok=True)

    cache = load_cache()
    new_cache = {}
    backbone_index_cache = {}  # rfd_run_name -> {stem: path}
    backbone_ca_cache = {}     # backbone_path (str) -> {"binder": ca, "target": ca}

    rows = []
    n_cached = 0
    n_computed = 0
    total = len(rf3_outputs)

    try:
        for i, folder in enumerate(rf3_outputs, start=1):
            cache_key = folder.name
            model_cif, summary_json = find_rf3_files(folder)
            sig = mtime_signature(model_cif, summary_json) if (model_cif and summary_json) else None

            cached_entry = cache.get(cache_key)
            if cached_entry is not None and sig is not None and cached_entry.get("mtime") == sig:
                rows.append(cached_entry["row"])
                new_cache[cache_key] = cached_entry
                n_cached += 1
            else:
                row = compute_row(folder, model_cif, summary_json,
                                   backbone_index_cache, backbone_ca_cache)
                rows.append(row)
                n_computed += 1
                if row["status"] == "ok" and sig is not None:
                    new_cache[cache_key] = {"mtime": sig, "row": row}

            if i % PROGRESS_EVERY == 0 or i == total:
                print(f"[{i}/{total}] processed ({n_cached} cached, {n_computed} computed)")
    finally:
        save_cache(new_cache)

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUT_CSV}")
    print(f"({n_cached} reused from cache, {n_computed} freshly computed)")


def compute_row(folder, model_cif, summary_json, backbone_index_cache, backbone_ca_cache):
    row = {
        "rf3_output": folder.name,
        "run_id": "",
        "backbone_file": "",
        "rfd_run": "",
        "min_pae_cross": "",
        "ca_rmsd_A": "",
        "n_ca_aligned": "",
        "n_target_ca_aligned": "",
        "status": "",
        "rf3_model_file": "",
    }
    issues = []

    if model_cif is None:
        issues.append("model.cif missing or ambiguous")
    else:
        row["rf3_model_file"] = model_cif.name

    if summary_json is None:
        issues.append("summary_confidences.json missing or ambiguous")

    # --- minPAE ---
    if summary_json is not None:
        try:
            row["min_pae_cross"] = round(cross_chain_min_pae_from_summary(summary_json), 4)
        except Exception as e:
            issues.append(f"minPAE: {e}")
            print(f"ERROR (minPAE) {folder.name}: {e}")

    # --- derive originating run_id, then scoped backbone lookup ---
    try:
        run_id = get_run_id(folder.name)
        row["run_id"] = run_id
        rfd_run_name = f"megf8_rfd_{run_id}"
        row["rfd_run"] = rfd_run_name
    except ValueError as e:
        issues.append(str(e))
        print(f"ERROR (run_id) {folder.name}: {e}")
        row["status"] = "; ".join(issues)
        return row

    if rfd_run_name not in backbone_index_cache:
        rfd_run_dir = RFD_RUNS_DIR / rfd_run_name
        if rfd_run_dir.is_dir():
            backbone_index_cache[rfd_run_name] = build_backbone_index_for_run(rfd_run_dir)
        else:
            backbone_index_cache[rfd_run_name] = {}
            print(f"WARNING: expected rfd_run dir not found: {rfd_run_dir}")

    backbone_index = backbone_index_cache[rfd_run_name]
    backbone_path = match_backbone(folder.name, backbone_index)

    if backbone_path is None:
        issues.append(f"no matching backbone found in {rfd_run_name}")
        print(f"WARNING: no matching backbone for {folder.name} (looked in {rfd_run_name})")
    elif model_cif is not None:
        row["backbone_file"] = backbone_path.name
        try:
            bb_key = str(backbone_path)
            if bb_key not in backbone_ca_cache:
                bb_structure = load_structure(backbone_path)
                backbone_ca_cache[bb_key] = {
                    "binder": get_ca(bb_structure, BINDER_CHAIN),
                    "target": get_ca(bb_structure, TARGET_CHAIN),
                }
            backbone_binder_ca = backbone_ca_cache[bb_key]["binder"]
            backbone_target_ca = backbone_ca_cache[bb_key]["target"]

            model_structure = load_structure(model_cif)
            model_binder_ca = get_ca(model_structure, BINDER_CHAIN)
            model_target_ca = get_ca(model_structure, TARGET_CHAIN)

            rmsd, n_binder_aligned, n_target_aligned = target_superimposed_binder_rmsd(
                backbone_binder_ca, backbone_target_ca,
                model_binder_ca, model_target_ca,
            )
            row["ca_rmsd_A"] = round(rmsd, 3)
            row["n_ca_aligned"] = n_binder_aligned
            row["n_target_ca_aligned"] = n_target_aligned
        except Exception as e:
            issues.append(f"RMSD: {e}")
            print(f"ERROR (RMSD) {folder.name}: {e}")

    row["status"] = "ok" if not issues else "; ".join(issues)
    return row


if __name__ == "__main__":
    start = time.time()
    process()
    print(f"Done in {time.time() - start:.1f}s")
