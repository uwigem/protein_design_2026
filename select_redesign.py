#!/usr/bin/env python3
"""
select_redesign.py

Pick the RF3 outputs worth a second pass, then hand them to redesign.sh.

Walks a directory of RF3 outputs, scores each one on the same two metrics the
batch_confidences_*.py scripts report -- cross-chain minPAE and chain-A
CA-RMSD against the originating RFD3 backbone -- keeps everything at or below
the thresholds you give, and invokes scripts/redesign.sh on the survivors.

    "equal or better" means minPAE <= --min_pae AND ca_rmsd <= --rms.

Optional exclusive floors turn each threshold into a band, which is how you
target the near-misses -- designs already at 2/2 have little room to improve,
while the ones just outside are where a redesign pass actually pays:

    --min_pae_floor 2 --min_pae 6   ->  2 < minPAE < 6
    --rms_floor 2     --rms 3       ->  2 < RMSD   < 3

Metrics come from one of two places:

  --confidences_csv PATH   join against an existing batch_results_confidences
                           csv on the RF3 output folder name. Fast, and uses
                           numbers you have already looked at.

  (default)                recompute from disk: minPAE out of each output's
                           *_summary_confidences.json, RMSD by finding the
                           output's backbone under rfd_runs/<target>_rfd_<run_id>/
                           and superimposing chain A. Self-contained, slower.

Outputs are grouped by target, so a mixed directory (passing_refold/, which
holds both mosmo and megf8 outputs) launches one redesign run per target
rather than failing. Both use the same --name, and since redesign.sh's
directories are prefixed by target they never collide.

Usage:
    python3 scripts/select_redesign.py \
        --input_dir all_out/all_bookkeep/fold_runs/mosmo_fold_eva_0823 \
        --min_pae 2.0 --rms 2.0 \
        --name eva_0823

    # preview the selection without launching anything
    python3 scripts/select_redesign.py ... --dry_run
"""

import argparse
import csv
import gzip
import io
import json
import os
import subprocess
import sys
from pathlib import Path

# biotite is imported lazily: it is only needed to RECOMPUTE RMSD from
# structures. With --confidences_csv the metrics come from a csv and this
# script then runs on a stock python with no third-party deps at all -- which
# matters on klone, where the default python has neither biotite nor numpy.
struc = None
pdbx = None
PDBX_FILE_CLS = None


def _require_biotite():
    global struc, pdbx, PDBX_FILE_CLS
    if struc is not None:
        return
    try:
        import biotite.structure as _struc
        import biotite.structure.io.pdbx as _pdbx
    except ImportError as e:
        raise SystemExit(
            f"Error: recomputing RMSD needs biotite, which is not importable ({e}).\n"
            "Either pass --confidences_csv to read metrics from a csv instead, or "
            "run this script inside foundry.sif, which has biotite:\n"
            "  apptainer exec --bind $PWD:$PWD --pwd $PWD containers/foundry.sif \\\n"
            "    /app/foundry/.venv/bin/python scripts/select_redesign.py ... --dry_run\n"
            "Use --dry_run in the container: launching redesign.sh from inside it "
            "would nest apptainer. The --report it writes is itself a valid "
            "--confidences_csv for a second, real run outside the container."
        )
    struc, pdbx = _struc, _pdbx
    PDBX_FILE_CLS = pdbx.CIFFile if hasattr(pdbx, "CIFFile") else pdbx.PDBxFile


BINDER_CHAIN = "A"
TARGETS = ("mosmo", "megf8")

FIELDNAMES = [
    "rf3_output", "target", "run_id", "rfd_run", "backbone_file",
    "min_pae_cross", "ca_rmsd_A", "n_ca_aligned", "selected", "status",
    "rf3_model_file",
]


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------

def find_rf3_models(input_dir):
    """Every RF3 model.cif under input_dir, flat or nested.

    fold_runs/<run>/<output>/ and the flat passing_refold/<output>/ layouts
    both work. The per-sample "seed-N_sample-M/" subdirectories are skipped:
    they are duplicate copies of the parent's model, and counting them would
    double every design.
    """
    models = []
    for path in sorted(Path(input_dir).rglob("*_model.cif")):
        if path.parent.name.startswith("seed-"):
            continue
        models.append(path)
    return models


def summary_for(model_cif):
    matches = list(model_cif.parent.glob("*_summary_confidences.json"))
    return matches[0] if len(matches) == 1 else None


def output_name(model_cif):
    """The RF3 output's identity: its containing folder name."""
    return model_cif.parent.name


def detect_target(name, override=None):
    if override:
        return override
    lowered = name.lower()
    for target in TARGETS:
        if lowered.startswith(target):
            return target
    return None


# --------------------------------------------------------------------------
# minPAE  (same contract as batch_confidences_refold_mosmo.py)
# --------------------------------------------------------------------------

def cross_chain_min_pae(summary_json_path):
    with open(summary_json_path) as f:
        summary = json.load(f)

    matrix = summary.get("chain_pair_pae_min")
    if matrix is None:
        raise ValueError("chain_pair_pae_min missing from summary_confidences.json")

    values = [v for row in matrix for v in row if v is not None]
    if not values:
        # A 1x1 all-null matrix means this was folded as a monomer (no target
        # chain), e.g. the passing_refold outputs from passing_jsoner.py.
        # There is no interface to score and nothing to extract a binder from.
        raise ValueError(
            "chain_pair_pae_min has no non-null entries — this looks like a "
            "single-chain (monomer) refold, not a cofolded complex"
        )
    if len(values) > 1:
        raise ValueError(
            f"chain_pair_pae_min has {len(values)} non-null entries (>2 chains?) "
            "- cannot unambiguously determine binder-vs-target minPAE"
        )
    return float(values[0])


# --------------------------------------------------------------------------
# RMSD against the originating RFD3 backbone
# --------------------------------------------------------------------------

def get_run_id(rf3_folder_name, target):
    """Recover the pipeline run name from an RF3 output folder name.

    Folder names look like "<target>_rfd_{run_id}_{target}_{run_id}_..." and
    the run_id can itself contain the target token, so the split point is the
    one where the prefix text repeats after the next "<target>_". Same logic
    as mosmo_refold_jsoner.py.
    """
    if "rfd_" not in rf3_folder_name:
        raise ValueError(f"no 'rfd_' in name: {rf3_folder_name}")
    after_rfd = rf3_folder_name.split("rfd_", 1)[1]
    marker = f"{target}_"
    start = 0
    while True:
        idx = after_rfd.find(marker, start)
        if idx == -1:
            raise ValueError(f"could not determine run_id for: {rf3_folder_name}")
        candidate = after_rfd[:idx].rstrip("_")
        remainder = after_rfd[idx + len(marker):]
        if remainder == candidate or remainder.startswith(candidate + "_"):
            return candidate
        start = idx + len(marker)


def read_cif(path):
    _require_biotite()
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            text = f.read()
        return PDBX_FILE_CLS.read(io.StringIO(text))
    return PDBX_FILE_CLS.read(path)


def load_ca(path, chain_id=None):
    _require_biotite()
    structure = pdbx.get_structure(read_cif(path), model=1)
    mask = (structure.atom_name == "CA") & struc.filter_amino_acids(structure)
    if chain_id is not None:
        mask = mask & (structure.chain_id == chain_id)
    return structure[mask]


def ca_rmsd(ref_ca, mob_ca):
    _require_biotite()
    n = min(len(ref_ca), len(mob_ca))
    if n == 0:
        raise ValueError("no chain-A CA atoms in one or both structures")
    fitted, _ = struc.superimpose(ref_ca[:n], mob_ca[:n])
    return float(struc.rmsd(ref_ca[:n], fitted)), n, (len(ref_ca) != len(mob_ca))


def build_backbone_index(rfd_run_dir):
    index = {}
    for f in sorted(rfd_run_dir.glob("*.cif.gz")):
        index[f.name[:-3]] = f          # strip ".gz", keep ".cif"
    for f in sorted(rfd_run_dir.glob("*.cif")):
        index.setdefault(f.name, f)
    return index


def match_backbone(rf3_folder_name, backbone_index):
    matches = [stem for stem in backbone_index if rf3_folder_name.startswith(stem)]
    if not matches:
        return None
    return backbone_index[max(matches, key=len)]


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def in_band(min_pae, rmsd, args):
    """Whether one output's metrics fall inside the requested window.

    Ceilings are inclusive (<=, "equal or better"); floors are exclusive (>),
    matching how a band like "2 < minPAE < 6" is normally written. A missing
    RMSD only reaches here under --allow_missing_rmsd, in which case the RMSD
    half of the window is skipped rather than assumed to pass.
    """
    if min_pae > args.min_pae:
        return False
    if args.min_pae_floor is not None and min_pae <= args.min_pae_floor:
        return False
    if rmsd is not None:
        if rmsd > args.rms:
            return False
        if args.rms_floor is not None and rmsd <= args.rms_floor:
            return False
    return True


def describe_band(args):
    pae = (f"{args.min_pae_floor} < minPAE <= {args.min_pae}"
           if args.min_pae_floor is not None else f"minPAE <= {args.min_pae}")
    rms = (f"{args.rms_floor} < RMSD <= {args.rms}"
           if args.rms_floor is not None else f"RMSD <= {args.rms}")
    return f"{pae} and {rms}"


def score_from_csv(rows_by_name, name):
    row = rows_by_name.get(name)
    if row is None:
        return None, None, "", f"not present in --confidences_csv"

    pae_str = (row.get("min_pae_cross") or "").strip()
    rmsd_str = (row.get("ca_rmsd_A") or "").strip()
    if not pae_str or not rmsd_str:
        return None, None, row.get("n_ca_aligned", ""), (
            f"blank metric in csv (status: {row.get('status', '')})"
        )
    try:
        return float(pae_str), float(rmsd_str), row.get("n_ca_aligned", ""), "ok"
    except ValueError:
        return None, None, "", "unparseable metric in csv"


def score_from_disk(model_cif, name, target, rfd_runs_dir, caches):
    index_cache, ca_cache = caches
    issues = []
    min_pae = None
    rmsd = None
    n_aligned = ""
    rfd_run = ""
    backbone_name = ""

    summary_json = summary_for(model_cif)
    if summary_json is None:
        issues.append("summary_confidences.json missing or ambiguous")
    else:
        try:
            min_pae = round(cross_chain_min_pae(summary_json), 4)
        except Exception as e:
            issues.append(f"minPAE: {e}")

    try:
        run_id = get_run_id(name, target)
    except ValueError as e:
        issues.append(str(e))
        return min_pae, rmsd, n_aligned, "", "", "; ".join(issues)

    rfd_run = f"{target}_rfd_{run_id}"
    if rfd_run not in index_cache:
        rfd_run_dir = Path(rfd_runs_dir) / rfd_run
        index_cache[rfd_run] = (
            build_backbone_index(rfd_run_dir) if rfd_run_dir.is_dir() else {}
        )
        if not rfd_run_dir.is_dir():
            print(f"WARNING: expected rfd_run dir not found: {rfd_run_dir}")

    backbone = match_backbone(name, index_cache[rfd_run])
    if backbone is None:
        issues.append(f"no matching backbone in {rfd_run}")
    else:
        backbone_name = backbone.name
        try:
            key = str(backbone)
            if key not in ca_cache:
                ca_cache[key] = load_ca(backbone, chain_id=BINDER_CHAIN)
            value, n_aligned, mismatch = ca_rmsd(
                ca_cache[key], load_ca(model_cif, chain_id=BINDER_CHAIN)
            )
            rmsd = round(value, 3)
            if mismatch:
                issues.append(f"CA count mismatch, truncated to {n_aligned}")
        except Exception as e:
            issues.append(f"RMSD: {e}")

    return min_pae, rmsd, n_aligned, rfd_run, backbone_name, ("ok" if not issues else "; ".join(issues))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Select RF3 outputs by minPAE/RMSD and redesign them"
    )
    parser.add_argument("--input_dir", required=True,
                        help="directory of RF3 outputs (searched recursively)")
    parser.add_argument("--min_pae", type=float, required=True,
                        help="keep outputs with cross-chain minPAE <= this")
    parser.add_argument("--rms", type=float, required=True,
                        help="keep outputs with chain-A CA-RMSD <= this")
    parser.add_argument("--min_pae_floor", type=float, default=None,
                        help="also require minPAE > this (exclusive). Use with "
                             "--min_pae to select a band of near-misses rather "
                             "than everything below a ceiling.")
    parser.add_argument("--rms_floor", type=float, default=None,
                        help="also require CA-RMSD > this (exclusive)")
    parser.add_argument("--name", required=True,
                        help="run name handed to redesign.sh")

    parser.add_argument("--confidences_csv", default=None,
                        help="read metrics from this csv instead of recomputing "
                             "(joined on the rf3_output column)")
    parser.add_argument("--rfd_runs_dir", default="./all_out/all_bookkeep/rfd_runs",
                        help="where <target>_rfd_<run_id>/ backbone dirs live, for "
                             "recomputing RMSD (default: %(default)s)")
    parser.add_argument("--target", choices=TARGETS, default=None,
                        help="force the target instead of reading it off each "
                             "output's name prefix")
    parser.add_argument("--allow_missing_rmsd", action="store_true",
                        help="select on minPAE alone when the backbone for an "
                             "output cannot be found (default: drop it)")

    parser.add_argument("--report", default=None,
                        help="write the full scored table here as CSV "
                             "(every output, selected or not)")
    parser.add_argument("--dry_run", action="store_true",
                        help="score and report, but do not launch redesign.sh")

    parser.add_argument("--redesign_sh", default="./scripts/redesign.sh",
                        help="path to redesign.sh (default: %(default)s)")
    parser.add_argument("--list_dir", default="./all_out/all_bookkeep/redesign_runs",
                        help="where the per-target cif lists handed to redesign.sh "
                             "are written (default: %(default)s)")
    parser.add_argument("--rfdiffusion", action="store_true",
                        help="pass --rfdiffusion through to redesign.sh")
    parser.add_argument("--mpnn_context", choices=["binder", "complex"], default=None,
                        help="pass --mpnn-context through to redesign.sh: what "
                             "ProteinMPNN sees (binder alone, or binder+target)")
    parser.add_argument("--partial_noise", type=float, default=None,
                        help="pass --partial-noise through to redesign.sh "
                             "(Angstroms of noise for partial diffusion)")
    parser.add_argument("--batch_size", type=int, default=None,
                        help="pass --batch-size through to redesign.sh")
    parser.add_argument("--hotspots", default=None,
                        help="pass --hotspots through to redesign.sh")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        raise SystemExit(f"Error: --input_dir '{input_dir}' is not a directory")

    models = find_rf3_models(input_dir)
    if not models:
        raise SystemExit(f"Error: no *_model.cif files found under {input_dir}")

    print(f"Found {len(models)} RF3 output(s) under {input_dir}")

    rows_by_name = {}
    if args.confidences_csv:
        with open(args.confidences_csv, newline="") as f:
            for row in csv.DictReader(f):
                rows_by_name[row["rf3_output"]] = row
        print(f"Loaded {len(rows_by_name)} rows from {args.confidences_csv}")

    caches = ({}, {})
    rows = []
    selected_by_target = {}

    for i, model_cif in enumerate(models, start=1):
        name = output_name(model_cif)
        target = detect_target(name, args.target)

        row = {
            "rf3_output": name,
            "target": target or "",
            "run_id": "",
            "rfd_run": "",
            "backbone_file": "",
            "min_pae_cross": "",
            "ca_rmsd_A": "",
            "n_ca_aligned": "",
            "selected": "no",
            "status": "",
            "rf3_model_file": model_cif.name,
        }

        if target is None:
            row["status"] = ("cannot tell mosmo from megf8 by name; "
                             "pass --target to force one")
            rows.append(row)
            continue

        if args.confidences_csv:
            min_pae, rmsd, n_aligned, status = score_from_csv(rows_by_name, name)
            rfd_run = rows_by_name.get(name, {}).get("rfd_run", "")
            backbone = rows_by_name.get(name, {}).get("backbone_file", "")
        else:
            min_pae, rmsd, n_aligned, rfd_run, backbone, status = score_from_disk(
                model_cif, name, target, args.rfd_runs_dir, caches
            )

        row.update({
            "rfd_run": rfd_run,
            "backbone_file": backbone,
            "min_pae_cross": "" if min_pae is None else min_pae,
            "ca_rmsd_A": "" if rmsd is None else rmsd,
            "n_ca_aligned": n_aligned,
            "status": status,
        })
        try:
            row["run_id"] = get_run_id(name, target)
        except ValueError:
            pass

        if min_pae is None:
            row["status"] = (row["status"] or "") + "; no minPAE — cannot select"
        elif rmsd is None and not args.allow_missing_rmsd:
            row["status"] = (row["status"] or "") + "; no RMSD — cannot select"
        elif in_band(min_pae, rmsd, args):
            row["selected"] = "yes"
            selected_by_target.setdefault(target, []).append(str(model_cif))

        rows.append(row)

        if i % 50 == 0 or i == len(models):
            n_sel = sum(len(v) for v in selected_by_target.values())
            print(f"[{i}/{len(models)}] scored ({n_sel} selected so far)")

    if args.report:
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        with open(args.report, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDNAMES)
            w.writeheader()
            w.writerows(rows)
        print(f"Report -> {args.report}")

    total_selected = sum(len(v) for v in selected_by_target.values())
    print(f"\n{total_selected}/{len(models)} output(s) pass {describe_band(args)}")
    for target, paths in sorted(selected_by_target.items()):
        print(f"  {target}: {len(paths)}")

    unscored = [r for r in rows if r["selected"] == "no" and not r["min_pae_cross"]]
    if unscored:
        print(f"({len(unscored)} output(s) could not be scored — see the status "
              f"column{' in ' + args.report if args.report else '; pass --report to inspect'})")

    if total_selected == 0:
        print("Nothing to redesign.")
        return

    if args.dry_run:
        print("\n[dry-run] would launch:")

    if not args.dry_run:
        os.makedirs(args.list_dir, exist_ok=True)
    failures = []

    for target, paths in sorted(selected_by_target.items()):
        list_path = Path(args.list_dir) / f"{target}_redesign_{args.name}_inputs.txt"
        if list_path.exists() and not args.dry_run:
            raise SystemExit(
                f"Error: '{list_path}' already exists. Choose a different --name "
                "or remove it — redesign.sh would refuse the collision anyway."
            )
        if not args.dry_run:
            with open(list_path, "w") as f:
                f.write("\n".join(paths) + "\n")

        cmd = ["bash", args.redesign_sh,
               "--name", args.name,
               f"--{target}",
               "--cif-list", str(list_path)]
        if args.rfdiffusion:
            cmd.append("--rfdiffusion")
        if args.mpnn_context:
            cmd += ["--mpnn-context", args.mpnn_context]
        if args.partial_noise is not None:
            cmd += ["--partial-noise", str(args.partial_noise)]
        if args.batch_size is not None:
            cmd += ["--batch-size", str(args.batch_size)]
        if args.hotspots:
            cmd += ["--hotspots", args.hotspots]

        printable = " ".join(cmd)
        if args.dry_run:
            print(f"  {printable}   # {len(paths)} cif(s)")
            continue

        print(f"\n== Redesigning {len(paths)} {target} binder(s) ==")
        print(f"   {printable}")
        # redesign.sh writes straight to the terminal; flush first so our own
        # buffered output does not end up printed after its.
        sys.stdout.flush()
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"ERROR: redesign.sh failed for {target} "
                  f"(exit {result.returncode})", file=sys.stderr)
            failures.append(target)

    if failures:
        raise SystemExit(f"redesign.sh failed for: {', '.join(failures)}")


if __name__ == "__main__":
    main()
