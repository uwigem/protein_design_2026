#!/usr/bin/env python3
"""
extract_binder.py

Pull the binder out of a cofolded RF3 (or RFD3) mmCIF so it can be fed back
into the design pipeline.

Every structure this repo produces downstream of rfd3.sh is a two-chain
complex: chain A is the binder, chain B is the target (MOSMO or MEGF8). This
script reads such a file and writes out either

  * the binder chain alone            (default -- what ProteinMPNN gets), or
  * binder + target, chains preserved (--keep_target -- what RFD3 partial
    diffusion gets, since it needs the target to hold the binding pose)

Output filenames are the input's name with a suffix inserted before the
extension, so provenance survives the whole rest of the pipeline:

    mosmo_rfd_eva_0823_..._model_8.cif_b1_d0_model.cif
        -> mosmo_rfd_eva_0823_..._model_8.cif_b1_d0_redesign.cif

The trailing "_model" that RF3 appends is stripped first (see --name_from),
so the redesign is named after the RF3 *output*, not its model file.

Usage:
    python3 scripts/extract_binder.py --out_dir DIR [options] CIF [CIF ...]
    python3 scripts/extract_binder.py --out_dir DIR --cif_list list.txt
"""

import argparse
import csv
import gzip
import io
import os
import sys
from pathlib import Path

import numpy as np

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

# biotite >=0.39 renamed PDBxFile -> CIFFile. Support both, same as the
# batch_confidences_*.py scripts do.
if hasattr(pdbx, "CIFFile"):
    PDBX_FILE_CLS = pdbx.CIFFile
else:
    PDBX_FILE_CLS = pdbx.PDBxFile


def read_cif(path):
    """Read a .cif or .cif.gz. Gzipped files are decompressed in memory only."""
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            text = f.read()
        return PDBX_FILE_CLS.read(io.StringIO(text))
    return PDBX_FILE_CLS.read(path)


def load_structure(path):
    return pdbx.get_structure(read_cif(path), model=1)


def base_name(path, name_from):
    """Derive the output stem from an input path.

    "auto"  strip .gz/.cif, then drop a trailing "_model" if present. For an
            RF3 output "<folder>/<folder>_model.cif" this reproduces the
            RF3 output folder name exactly.
    "file"  strip .gz/.cif only.
    "dir"   use the containing directory's name.
    """
    p = Path(path)
    if name_from == "dir":
        return p.parent.name

    stem = p.name
    if stem.endswith(".gz"):
        stem = stem[: -len(".gz")]
    if stem.endswith(".cif"):
        stem = stem[: -len(".cif")]

    if name_from == "auto" and stem.endswith("_model"):
        stem = stem[: -len("_model")]
    return stem


def write_cif(atoms, out_path, block_name, compress=False):
    """Write an mmCIF, optionally gzipped.

    ProteinMPNN's wrapper (ProteinMPNN/scripts/mpnn.sh) iterates over
    "$INPUT_DIR"/*.cif.gz, so anything destined for it MUST be gzipped or it
    is silently invisible. rfd3 backbones are gzipped for the same reason.
    """
    out_file = PDBX_FILE_CLS()
    pdbx.set_structure(out_file, atoms, data_block=block_name)
    if compress:
        buf = io.StringIO()
        out_file.write(buf)
        with gzip.open(str(out_path), "wt") as f:
            f.write(buf.getvalue())
    else:
        out_file.write(str(out_path))


def extract(path, args):
    """Returns (out_name, atoms) or raises."""
    atoms = load_structure(path)
    atoms = atoms[struc.filter_amino_acids(atoms)]

    present = sorted(set(atoms.chain_id))
    if args.binder_chain not in present:
        raise ValueError(
            f"binder chain '{args.binder_chain}' not found (chains present: "
            f"{', '.join(present) or 'none'})"
        )

    if args.keep_target:
        if args.target_chain not in present:
            raise ValueError(
                f"--keep_target given but target chain '{args.target_chain}' "
                f"not found (chains present: {', '.join(present)})"
            )
        wanted = (args.binder_chain, args.target_chain)
    else:
        wanted = (args.binder_chain,)

    selected = atoms[np.isin(atoms.chain_id, wanted)]

    if selected.array_length() == 0:
        raise ValueError("no atoms left after chain selection")

    return selected


def main():
    parser = argparse.ArgumentParser(
        description="Extract the binder (optionally binder+target) from cofolded mmCIFs"
    )
    parser.add_argument("cifs", nargs="*", help="input .cif / .cif.gz files")
    parser.add_argument("--cif_list", default=None,
                        help="plain text file, one input cif path per line "
                             "(combined with any positional CIFs)")
    parser.add_argument("--out_dir", required=True,
                        help="directory to write extracted structures into (created if missing)")
    parser.add_argument("--suffix", default="_redesign",
                        help="inserted before .cif in every output name (default: _redesign)")
    parser.add_argument("--binder_chain", default="A",
                        help="chain id of the binder in the input (default: A)")
    parser.add_argument("--target_chain", default="B",
                        help="chain id of the target in the input (default: B)")
    parser.add_argument("--keep_target", action="store_true",
                        help="write binder AND target instead of the binder alone")
    parser.add_argument("--name_from", choices=["auto", "file", "dir"], default="auto",
                        help="how to derive the output stem (default: auto)")
    parser.add_argument("--manifest", default=None,
                        help="optional CSV recording source -> output for each file")
    parser.add_argument("--gzip", action="store_true",
                        help="write .cif.gz instead of .cif. REQUIRED for anything "
                             "feeding mpnn.sh, which only globs *.cif.gz.")
    parser.add_argument("--overwrite", action="store_true",
                        help="allow overwriting existing outputs (default: refuse)")
    args = parser.parse_args()

    inputs = list(args.cifs)
    if args.cif_list:
        with open(args.cif_list) as f:
            inputs.extend(l.strip() for l in f if l.strip() and not l.startswith("#"))

    if not inputs:
        parser.error("no input CIFs given (pass paths positionally or via --cif_list)")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    n_ok = 0
    n_err = 0

    for path in inputs:
        stem = base_name(path, args.name_from)
        ext = ".cif.gz" if args.gzip else ".cif"
        out_path = out_dir / f"{stem}{args.suffix}{ext}"

        if out_path.exists() and not args.overwrite:
            print(f"[skip] {out_path.name} already exists (use --overwrite to replace)")
            rows.append({"source": path, "output": str(out_path), "status": "exists"})
            continue

        try:
            selected = extract(path, args)
            block = out_path.name[: -len(".cif.gz")] if args.gzip else out_path.stem
            write_cif(selected, out_path, block, compress=args.gzip)
        except Exception as e:
            print(f"[error] {path}: {e}", file=sys.stderr)
            rows.append({"source": path, "output": "", "status": f"error: {e}"})
            n_err += 1
            continue

        chains = sorted(set(selected.chain_id))
        n_res = struc.get_residue_count(selected)
        print(f"[ok]   {Path(path).name} -> {out_path.name} "
              f"(chains {'+'.join(chains)}, {n_res} residues)")
        rows.append({"source": path, "output": str(out_path), "status": "ok"})
        n_ok += 1

    if args.manifest:
        os.makedirs(os.path.dirname(args.manifest) or ".", exist_ok=True)
        with open(args.manifest, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["source", "output", "status"])
            w.writeheader()
            w.writerows(rows)
        print(f"Manifest -> {args.manifest}")

    print(f"\nExtracted {n_ok}/{len(inputs)} into {out_dir}"
          + (f" ({n_err} failed)" if n_err else ""))

    if n_ok == 0:
        raise SystemExit("nothing was extracted")


if __name__ == "__main__":
    main()
