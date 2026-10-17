#!/usr/bin/env python3
"""
python_1.py — build the RFdiffusion input json for a MOSMO binder design run.

Usage:
    python python_1.py --name NAME --hotspots "A39,A47,A56" \
        --output_json all_out/all_bookkeep/all_json/mosmo_rfd_NAME.json

Writes the json to --output_json (parent directory must already exist —
this script does not create it), e.g.:

{
    "mosmo_<NAME>": {
        "dialect": 2,
        "infer_ori_strategy": "hotspots",
        "input": "/igem/targets/MOSMO.pdb",
        "contig": "40-120,/0,A2-150",
        "select_hotspots": {
            "A39": "CA",
            "A47": "CA",
            "A56": "CA"
        },
        "is_non_loopy": true
    }
}
"""
import argparse
import json
import os

TARGET = "mosmo"

# ---------------------------------------------------------------------------
# TODO: confirm these fixed config values for MOSMO runs
# ---------------------------------------------------------------------------
INPUT_PDB = "/igem/targets/MOSMO.pdb"
CONTIG = "40-120,/0,A2-150"
DIALECT = 2
INFER_ORI_STRATEGY = "hotspots"
IS_NON_LOOPY = True
HOTSPOT_ATOM = "CA"  # TODO: confirm all hotspots use CA

# Fallback used only if --output_json is not given, kept for backward
# compatibility with any existing manual invocations of this script.
DEFAULT_OUTPUT_JSON_TEMPLATE = "all_out/all_bookkeep/all_json/{target}_rfd_{name}.json"


def parse_hotspots(hotspots_str):
    """Turn "A39,A47,A56" into {"A39": "CA", "A47": "CA", "A56": "CA"}."""
    hotspots = {}
    for token in hotspots_str.split(","):
        token = token.strip()
        if not token:
            continue
        hotspots[token] = HOTSPOT_ATOM
    return hotspots


def main():
    parser = argparse.ArgumentParser(description="Build MOSMO RFdiffusion input json")
    parser.add_argument("--name", required=True, help="base run name")
    parser.add_argument(
        "--hotspots",
        required=True,
        help='comma-separated hotspot residues, e.g. "A39,A47,A56"',
    )
    parser.add_argument(
        "--output_json",
        required=False,
        default=None,
        help=(
            "path to write the RFdiffusion input json to. Parent directory "
            "must already exist — it is not created by this script. "
            f"Defaults to '{DEFAULT_OUTPUT_JSON_TEMPLATE}' if omitted."
        ),
    )
    args = parser.parse_args()

    key = f"{TARGET}_{args.name}"
    out_path = args.output_json or DEFAULT_OUTPUT_JSON_TEMPLATE.format(
        target=TARGET, name=args.name
    )

    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        raise SystemExit(
            f"Error: output directory '{out_dir}' does not exist. "
            "This script does not create it — create it first or point "
            "--output_json at an existing directory."
        )

    payload = {
        key: {
            "dialect": DIALECT,
            "infer_ori_strategy": INFER_ORI_STRATEGY,
            "input": INPUT_PDB,
            "contig": CONTIG,
            "select_hotspots": parse_hotspots(args.hotspots),
            "is_non_loopy": IS_NON_LOOPY,
        }
    }

    with open(out_path, "w") as f:
        json.dump(payload, f, indent=4)

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
