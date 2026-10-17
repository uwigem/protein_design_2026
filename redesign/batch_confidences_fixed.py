#!/usr/bin/env python3
"""
filter_binders.py

================================================================================
PURPOSE
================================================================================
For every RoseTTAFold3 (RF3) folded output in a large, ever-growing batch,
compute:
  1. cross-chain minimum PAE (binder vs. target),
  2. mean pLDDT of the binder chain (chain A) from the RF3 model, and
  3. CA-RMSD (chain A only) of the RF3 prediction superimposed onto its
     RFdiffusion3 (RFD3) backbone
and write ONE row per RF3 output to a single CSV. Meant to be rerun
repeatedly as new outputs are added, so already-computed rows are cached
and never recomputed unless their underlying files change.

================================================================================
ASSUMPTIONS / PREFERENCES (as specified in conversation, do not change
without re-confirming with the user)
================================================================================
- No CLI arguments. All paths are hardcoded absolute paths (see CONFIG below).

- DIRECTORY STRUCTURE (two levels of nesting under all_out/all_bookkeep/):

    all_out/all_bookkeep/rfd_runs/
        <target>_rfd_<name>/              e.g. mosmo_rfd_testrun2/
            <backbone>.cif.gz             one gzipped file per backbone
            ...
        <target>_rfd_<other_name>/
            ...

    all_out/all_bookkeep/fold_runs/
        <target>_fold_<name>/             e.g. mosmo_fold_testrun2/
            <rf3_output_folder>/          one subdir per RF3 output (5 per
                <prefix>_model.cif        backbone: different b/d indices)
                <prefix>_confidences.json
                <prefix>_summary_confidences.json
                <prefix>_ranking_scores.csv
                seed-0_sample-0/          (unused - duplicate per-sample copy)
            ...
        <target>_fold_<other_name>/
            ...

  "target" is either "mosmo" or "megf8". "name" is a unique identifier for
  one run through the whole pipeline (one diffusion run producing many
  backbones + sequences, forward-folded into many RF3 outputs).

  IMPORTANT: the backbone for any RF3 output inside
  "fold_runs/<target>_fold_<name>/" ALWAYS lives under the correspondingly
  named "rfd_runs/<target>_rfd_<name>/" directory (i.e. same target, same
  name, only "fold"/"rfd" differs). So backbone lookup is SCOPED per fold
  run rather than global: for each fold_run directory, only its one
  matching rfd_run directory is indexed/searched (lazily, and cached in
  memory per rfd_run so it's only built once even though 5+ RF3 outputs
  per backbone will need it repeatedly). Within that one rfd_run directory,
  matching is by filename prefix (the backbone's filename, including its
  trailing ".cif", is a literal string-prefix of the RF3 output folder's
  name - e.g. backbone "..._model_0.cif.gz" matches RF3 output folder
  "..._model_0.cif_b0_d0").

- Backbones are ALWAYS gzip-compressed (.cif.gz) - never plain .cif. They
  are decompressed IN MEMORY only (via gzip + io.StringIO) when actually
  read for RMSD; nothing is ever written back to rfd_runs/ on disk. This
  avoids write-permission issues in a shared/read-only directory and avoids
  unnecessary disk I/O for files that are only read, never re-read outside
  this script's in-run cache.

- Each RF3 output folder's "*_model.cif" and "*_summary_confidences.json"
  are the only files read - "*_confidences.json" (full PAE matrix) and
  "*_ranking_scores.csv" are not used, and "seed-0_sample-0/" (duplicate
  per-sample copy) is ignored.

- minPAE is read directly from "*_summary_confidences.json"'s
  "chain_pair_pae_min" field: an NxN upper-triangular matrix (null on/below
  the diagonal). All designs here have exactly 2 chains (binder=chain A,
  target=MEGF8 or MOSMO), so there is exactly one non-null entry, which IS
  the cross-chain minPAE - taken directly (flatten, drop nulls, take the
  min). If a design ever has >2 chains this is flagged in "status" rather
  than guessed at, since a flat min would not unambiguously give the
  binder-vs-target value in that case.

- pLDDT is NOT read from "*_confidences.json" or "*_summary_confidences.json"
  (per the "only *_model.cif and *_summary_confidences.json are read" rule
  above, and summary_confidences.json does not carry a usable per-residue
  or binder-only pLDDT field anyway). Instead, RF3 writes per-atom pLDDT
  into the mmCIF B_iso_or_equiv (B-factor) column of "*_model.cif" itself
  - confirmed directly against a real RF3 model.cif. In these files that
  column is a 0-1 fraction (e.g. 0.726, 0.738, ...) rather than the
  conventional 0-100 pLDDT scale, so mean_plddt_A is rescaled x100 to
  match the standard 0-100 convention (e.g. AlphaFold DB). mean_plddt_A
  is the mean of that rescaled per-atom value over chain A (binder) CA
  atoms in the RF3 model - reusing the same chain-A CA parse of model.cif
  that's already done for RMSD, so no extra file is read and model.cif is
  parsed only once per RF3 output (b_factor is requested as an extra
  field on that same parse). This is computed independently of whether a
  matching backbone is found, since pLDDT only depends on model.cif, not
  on the RFD3 backbone.

- RMSD: CA-only Kabsch superposition of the RF3 model onto its matched RFD3
  backbone, via biotite, restricted to chain A (the binder) on BOTH sides.
  This matters because the RFD3 backbone typically contains only the binder
  chain while the RF3 model contains the full binder+target complex - using
  all CA atoms on both sides (rather than chain-A-only on both sides)
  previously produced a spurious "CA count mismatch" on every single row.
  Chain-A length is constant across the 5 RF3 outputs sharing a backbone,
  so backbone chain-A CA coordinates are parsed once per backbone per run
  and cached in memory (not on disk) for reuse across those 5 outputs.

- No thresholds / no pass-fail filtering - raw values only, user filters
  and sorts the CSV themselves after download.
- No ranking / no grouping by backbone - flat CSV, one row per RF3 output.
- ALL RF3 outputs found on disk appear in the CSV, even if data is missing
  or broken (missing model.cif, no backbone match, malformed JSON, etc.),
  with a free-text "status" column explaining what went wrong (status is
  "ok" only when minPAE, pLDDT, and RMSD were all computed cleanly with no
  issues) - this is meant to make troubleshooting the batch straightforward.
- Progress is printed to stdout every 50 RF3 outputs processed; per-file
  error detail is also printed to stdout (not the CSV, beyond the short
  status string).

- OUTPUT CSV: /mmfs1/gscratch/stf/igem/all_out/batch_results_confidences.csv
  Rewritten each run, but any row whose underlying files haven't changed is
  pulled straight from the cache rather than recomputed - see CACHING.

- CACHE FILE: /mmfs1/gscratch/stf/igem/all_out/.batch_results_confidences_cache.json
  A small JSON file recording, for every RF3 output previously computed
  with status == "ok", its row data plus an mtime signature (max of its
  model.cif and summary_confidences.json mtimes) at computation time.
    * cache key = RF3 output's path relative to fold_runs/, i.e.
      "<target>_fold_<name>/<output_folder>" (leaf folder names alone are
      not guaranteed unique across different runs)
    * if a key is cached AND its current mtime signature matches -> reuse
      the cached row, skip recomputation entirely
    * otherwise (new output, changed files, or previous status != "ok")
      -> recompute from scratch
  Rows with any error/missing-data status are NEVER cached as final - they
  are retried every run, since the cause may get fixed later. Cache
  entries for RF3 outputs no longer present on disk are pruned. The cache
  is saved after every run (even if the run is interrupted partway by an
  unexpected error, via a try/finally) so progress on a large batch is
  never silently lost.
================================================================================
"""

import csv
import gzip
import io
import json
import os
import time
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

# biotite >=0.39 renamed PDBxFile -> CIFFile. Support both.
if hasattr(pdbx, "CIFFile"):
    PDBX_FILE_CLS = pdbx.CIFFile
else:
    PDBX_FILE_CLS = pdbx.PDBxFile


# --------------------------------------------------------------------------
# CONFIG - hardcoded absolute paths, no CLI arguments
# --------------------------------------------------------------------------

ALL_OUT_DIR = Path("/mmfs1/gscratch/stf/igem/all_out")
BOOKKEEP_DIR = ALL_OUT_DIR / "all_bookkeep"
RFD_RUNS_DIR = BOOKKEEP_DIR / "rfd_runs"     # contains <target>_rfd_<name>/ subdirs
FOLD_RUNS_DIR = BOOKKEEP_DIR / "fold_runs"   # contains <target>_fold_<name>/ subdirs

OUT_CSV = ALL_OUT_DIR / "batch_results_confidences_plDDT_fixed.csv"
CACHE_FILE = ALL_OUT_DIR / ".batch_results_confidences_cache_plDDT.json"

PROGRESS_EVERY = 50

BINDER_CHAIN = "A"  # RMSD and pLDDT are both computed on chain A only, in
                     # both the RFD3 backbone (which typically contains only
                     # the binder chain) and the RF3 model (which contains
                     # the full binder+target complex) - restricting both
                     # sides to chain A avoids comparing mismatched atom
                     # counts, and keeps pLDDT specific to the binder.

FIELDNAMES = [
    "rf3_output", "fold_run", "backbone_file", "rfd_run",
    "min_pae_cross", "mean_plddt_A", "ca_rmsd_A", "n_ca_aligned",
    "status", "rf3_model_file",
]


# --------------------------------------------------------------------------
# mmCIF loading (backbones are always .cif.gz, decompressed in memory only)
# --------------------------------------------------------------------------

def _read_cif_file(path):
    """Return a parsed PDBx/CIFFile object from a .cif or .cif.gz path.
    .cif.gz is decompressed in memory (io.StringIO) - nothing is written
    to disk.
    """
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            text = f.read()
        return PDBX_FILE_CLS.read(io.StringIO(text))
    return PDBX_FILE_CLS.read(path)


def load_ca(path, chain_id=None, extra_fields=None):
    """Parse a .cif/.cif.gz file and return its CA atoms (optionally
    restricted to one chain). extra_fields is passed straight through to
    biotite's get_structure - pass extra_fields=["b_factor"] to also pull
    back per-atom B-factor (used for pLDDT, since RF3 writes per-atom
    pLDDT into the B-factor column of *_model.cif).
    """
    cif_file = _read_cif_file(path)
    structure = pdbx.get_structure(cif_file, model=1, extra_fields=extra_fields)
    mask = (structure.atom_name == "CA") & struc.filter_amino_acids(structure)
    if chain_id is not None:
        mask = mask & (structure.chain_id == chain_id)
    return structure[mask]


def ca_rmsd(ref_ca, mob_ca):
    n = min(len(ref_ca), len(mob_ca))
    if n == 0:
        raise ValueError("no chain-A CA atoms found in one or both structures")
    fitted_ca, _ = struc.superimpose(ref_ca[:n], mob_ca[:n])
    rmsd = struc.rmsd(ref_ca[:n], fitted_ca)
    return float(rmsd), n, (len(ref_ca) != len(mob_ca))


def mean_plddt(ca_atoms):
    """Mean per-atom pLDDT (RF3 stores this in the B-factor column, as a
    0-1 fraction rather than the conventional 0-100 scale) over a set of
    chain-A CA atoms from an RF3 model.cif. Rescaled x100 here so the
    output matches the standard 0-100 pLDDT convention used elsewhere
    (e.g. AlphaFold DB).
    """
    if len(ca_atoms) == 0:
        raise ValueError("no chain-A CA atoms found in RF3 model for pLDDT")
    if "b_factor" not in ca_atoms.get_annotation_categories():
        raise ValueError("model.cif has no B-factor column to read pLDDT from")
    return float(ca_atoms.b_factor.mean()) * 100.0


# --------------------------------------------------------------------------
# minPAE straight from summary_confidences.json
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
        # more than 2 chains - flat min is not guaranteed to be the
        # binder-vs-target value. Flag rather than silently guess.
        raise ValueError(
            f"chain_pair_pae_min has {len(values)} non-null entries "
            f"(>2 chains?) - cannot unambiguously determine binder-vs-target minPAE"
        )

    return float(values[0])


# --------------------------------------------------------------------------
# Backbone discovery and matching (scoped per rfd_run directory)
# --------------------------------------------------------------------------

def rfd_run_name_for_fold_run(fold_run_name):
    """'<target>_fold_<name>' -> '<target>_rfd_<name>'."""
    return fold_run_name.replace("_fold_", "_rfd_", 1).partition("_redesign")[0]


def build_backbone_index_for_run(rfd_run_dir):
    """Map backbone stem (filename ending in .cif, .gz stripped) -> full
    .cif.gz path, for every backbone in a single rfd_run directory.
    """
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
# RF3 output discovery
# --------------------------------------------------------------------------

def find_rf3_outputs(fold_runs_dir):
    """Yield (fold_run_name, output_folder_path) for every RF3 output
    subdirectory under every "<target>_fold_<name>/" run directory.
    """
    fold_run_dirs = sorted(p for p in fold_runs_dir.iterdir() if p.is_dir())
    for run_dir in fold_run_dirs:
        for output_dir in sorted(p for p in run_dir.iterdir() if p.is_dir()):
            yield run_dir.name, output_dir


def find_rf3_files(folder):
    """Return (model_cif_path, summary_json_path) or (None, None) if absent."""
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
    if not RFD_RUNS_DIR.is_dir():
        raise SystemExit(f"RFD runs dir not found: {RFD_RUNS_DIR}")
    if not FOLD_RUNS_DIR.is_dir():
        raise SystemExit(f"Fold runs dir not found: {FOLD_RUNS_DIR}")

    rf3_outputs = list(find_rf3_outputs(FOLD_RUNS_DIR))
    if not rf3_outputs:
        raise SystemExit(f"No RF3 output subdirectories found under {FOLD_RUNS_DIR}")

    cache = load_cache()
    new_cache = {}
    backbone_index_cache = {}  # rfd_run_name -> {stem: path} (built lazily, once per run)
    backbone_ca_cache = {}     # backbone_path (str) -> chain-A CA AtomArray

    rows = []
    n_cached = 0
    n_computed = 0
    total = len(rf3_outputs)

    try:
        for i, (fold_run_name, folder) in enumerate(rf3_outputs, start=1):
            cache_key = f"{fold_run_name}/{folder.name}"
            model_cif, summary_json = find_rf3_files(folder)
            sig = mtime_signature(model_cif, summary_json) if (model_cif and summary_json) else None

            cached_entry = cache.get(cache_key)
            if cached_entry is not None and sig is not None and cached_entry.get("mtime") == sig:
                rows.append(cached_entry["row"])
                new_cache[cache_key] = cached_entry
                n_cached += 1
            else:
                row = compute_row(fold_run_name, folder, model_cif, summary_json,
                                   backbone_index_cache, backbone_ca_cache)
                rows.append(row)
                n_computed += 1
                if row["status"] == "ok" and sig is not None:
                    new_cache[cache_key] = {"mtime": sig, "row": row}

            if i % PROGRESS_EVERY == 0 or i == total:
                print(f"[{i}/{total}] processed ({n_cached} cached, {n_computed} computed)")
    finally:
        # always persist whatever progress was made, even if the loop above
        # raised partway through a large batch
        save_cache(new_cache)

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {OUT_CSV}")
    print(f"({n_cached} reused from cache, {n_computed} freshly computed)")


def compute_row(fold_run_name, folder, model_cif, summary_json, backbone_index_cache, backbone_ca_cache):
    row = {
        "rf3_output": folder.name,
        "fold_run": fold_run_name,
        "backbone_file": "",
        "rfd_run": "",
        "min_pae_cross": "",
        "mean_plddt_A": "",
        "ca_rmsd_A": "",
        "n_ca_aligned": "",
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
            print(f"ERROR (minPAE) {fold_run_name}/{folder.name}: {e}")

    # --- chain-A CA parse of the RF3 model (with B-factor for pLDDT) ---
    # Parsed once here (if model_cif exists) and reused below for RMSD, so
    # model.cif is never read/parsed twice for the same RF3 output.
    mob_ca = None
    if model_cif is not None:
        try:
            mob_ca = load_ca(model_cif, chain_id=BINDER_CHAIN, extra_fields=["b_factor"])
            row["mean_plddt_A"] = round(mean_plddt(mob_ca), 3)
        except Exception as e:
            issues.append(f"pLDDT: {e}")
            print(f"ERROR (pLDDT) {fold_run_name}/{folder.name}: {e}")

    # --- scoped backbone lookup: only the one rfd_run dir that must contain it ---
    rfd_run_name = rfd_run_name_for_fold_run(fold_run_name)
    row["rfd_run"] = rfd_run_name

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
        print(f"WARNING: no matching backbone for {fold_run_name}/{folder.name} "
              f"(looked in {rfd_run_name})")
    else:
        row["backbone_file"] = backbone_path.name
        if model_cif is not None:
            try:
                bb_key = str(backbone_path)
                if bb_key not in backbone_ca_cache:
                    backbone_ca_cache[bb_key] = load_ca(backbone_path, chain_id=None) #BINDER_CHAIN
                ref_ca = backbone_ca_cache[bb_key]
                # reuse the chain-A CA parse from the pLDDT step above rather
                # than re-parsing model.cif; only fall back to a fresh parse
                # if that earlier parse failed for some reason.
                # if mob_ca is None:
                mob_ca = load_ca(model_cif, chain_id=None) #BINDER_CHAIN
                rmsd, n_aligned, mismatch = ca_rmsd(ref_ca, mob_ca)
                row["ca_rmsd_A"] = round(rmsd, 3)
                row["n_ca_aligned"] = n_aligned
                if mismatch:
                    issues.append(f"CA count mismatch, truncated to {n_aligned}")
            except Exception as e:
                issues.append(f"RMSD: {e}")
                print(f"ERROR (RMSD) {fold_run_name}/{folder.name}: {e}")

    row["status"] = "ok" if not issues else "; ".join(issues)
    return row


if __name__ == "__main__":
    start = time.time()
    process()
    print(f"Done in {time.time() - start:.1f}s")


