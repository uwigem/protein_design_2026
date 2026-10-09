#!/usr/bin/env bash
#
# redesign.sh
#
# Takes cofolded RF3 output .cif files and pushes their binders back through
# the design pipeline, on the theory that a second pass sharpens them.
#
# This is batch_design.sh's sibling. batch_design.sh starts from a target and
# invents binders; redesign.sh starts from binders that already exist and
# improves them. Everything downstream of the first step is identical, so the
# same rfd3.sh / mpnn.sh / <target>_jsoner.py / fold.sh wrappers are reused
# verbatim and the outputs slot into the same bookkeeping tree.
#
# Two modes:
#
#   default (ProteinMPNN, complex context)
#     1. extract_binder.py --keep_target -> binder+target .cif per input
#     2. mpnn.sh (--designed_chain A)    -> re-designed binder sequences (.fa)
#     3. <target>_jsoner.py              -> <target>_fold_<name>_redesign.json
#                                          (trims the target off at the junction)
#     4. fold.sh                         -> refolded complexes
#
#   ProteinMPNN sees the target by default because it has to: designing the
#   binder in isolation (--mpnn-context binder) has no interface to satisfy,
#   so MPNN optimises the monomer and collapses toward poly-alanine. Measured
#   on two aug18_1 binders, binder-only took minPAE 4.15/4.60 -> 17-20 and
#   CA-RMSD 2.1/2.8 -> 5-12 A. Use "binder" only if you want that behaviour.
#
#   --rfdiffusion (partial diffusion first, to lightly correct the backbone)
#     1. extract_binder.py --keep_target  -> binder+target .cif per input
#     2. redesign_json_rfd.py             -> partial-diffusion rfd json
#     3. rfd3.sh                          -> corrected backbones
#     4. extract_binder.py                -> binder-only again, for MPNN
#     5. mpnn.sh / jsoner / fold.sh       -> as above
#
# Every artifact is named "<original>_redesign" (or lives in a directory
# suffixed "_redesign"), so a redesign is always traceable to its parent and
# never collides with the run it came from.
#
# Run from the igem root (/mmfs1/gscratch/stf/igem), same as batch_design.sh.

set -euo pipefail

on_error() {
    echo "Error: redesign aborted — command on line ${1} exited non-zero." >&2
    exit 1
}
trap 'on_error $LINENO' ERR

# ---------------------------------------------------------------------------
# EDIT ME per run — script locations
# ---------------------------------------------------------------------------
PY_EXTRACT="./scripts/extract_binder.py"             # cif -> binder (or binder+target)
PY_RFD_JSON="./scripts/redesign_json_rfd.py"         # partial-diffusion rfd json
PYTHON_3="./RosettaFold/scripts/mosmo_jsoner.py"        # mosmo RF3-json-prep step
PYTHON_4="./RosettaFold/scripts/megf8_jsoner.py"        # megf8 RF3-json-prep step
BASH_1="./RFDiffusion/scripts/rfd3.sh"                  # partial diffusion
BASH_2="./ProteinMPNN/scripts/mpnn.sh"                  # run ProteinMPNN over a design dir
BASH_3="./RosettaFold/scripts/fold.sh"                  # run RF3 fold over fold json

# ---------------------------------------------------------------------------
# EDIT ME per run — relative base directories for each stage's I/O.
# These match batch_design.sh so redesigns land in the same tree, with
# "_redesign" appended to each run directory name.
# ---------------------------------------------------------------------------
REDESIGN_IN_DIR="./all_out/all_bookkeep/redesign_runs"  # extracted binders (created if missing)
RFD_JSON_DIR="./all_out/all_bookkeep/all_json"          # where the rfd json is written/read
DESIGN_OUT_DIR="./all_out/all_bookkeep/rfd_runs"        # rfd3.sh --out_dir root
MPNN_OUT_DIR="./all_out/all_bookkeep/mpnn_runs"         # mpnn.sh --output_dir root
FOLD_JSON_DIR="./all_out/all_bookkeep/all_json"         # where the fold json is written/read
FOLD_OUT_DIR="./all_out/all_bookkeep/fold_runs"         # fold.sh --out_dir root

SUFFIX="_redesign"

# ---------------------------------------------------------------------------
# Python for this script's own helpers.
#
# klone's default python3 has neither biotite nor numpy, but foundry.sif ships
# biotite 1.4.0 / numpy 2.4.2, so the helpers run inside the container by
# default and nothing needs installing. The igem root is bound to itself (not
# to /igem) so host paths — including absolute paths to input CIFs anywhere
# under it — resolve identically inside and out.
#
# Override with REDESIGN_PYTHON=python3 if you have a host env with biotite.
# ---------------------------------------------------------------------------
IGEM_ROOT="$(pwd -P)"
if [[ -n "${REDESIGN_PYTHON:-}" ]]; then
    read -r -a PY_RUN <<< "$REDESIGN_PYTHON"
else
    PY_RUN=(apptainer exec
            --bind "${IGEM_ROOT}:${IGEM_ROOT}"
            --pwd "${IGEM_ROOT}"
            ./containers/foundry.sif
            /app/foundry/.venv/bin/python)
fi

usage() {
    cat <<EOF
Usage: $(basename "$0") --name NAME (--mosmo | --megf8) [--rfdiffusion] CIF [CIF ...]
       $(basename "$0") --name NAME (--mosmo | --megf8) --cif-list FILE

  --name NAME            base name used to construct output filenames/dirs
  --mosmo                target MOSMO for this run
  --megf8                target MEGF8 for this run
  --cif-list FILE        plain text file, one input .cif path per line
                         (combined with any CIFs given positionally)

  --rfdiffusion          run RFD3 partial diffusion before ProteinMPNN
                         (default is ProteinMPNN directly)
  --partial-noise A      Angstroms of noise for partial diffusion (default: 5).
                         Lower = lighter correction; rfd3 recommends <= 15.
                         --rfdiffusion only.
  --batch-size N         diffusion_batch_size for rfd3.sh (default: 10).
                         --rfdiffusion only.
  --hotspots "B41,B47"   fixed hotspots for every entry. --rfdiffusion only;
                         if omitted, each complex's own interface is used.

  --mpnn-context C       what ProteinMPNN sees: "complex" (default, binder +
                         target, designing chain A in the target's presence)
                         or "binder" (the binder chain alone). Binder-only
                         cannot see the interface, so it optimises the monomer
                         and destroys binding -- see the note at the top of
                         this script. Ignored with --rfdiffusion for the
                         pre-diffusion step.
  --binder-chain C       binder chain id in the input cifs (default: A)
  --target-chain C       target chain id in the input cifs (default: B)
  --dry-run              build inputs and print the commands, run nothing expensive

Exactly one of --mosmo / --megf8 must be given.
Run from the igem root.
EOF
    exit 1
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
NAME=""
TARGET=""          # "mosmo" or "megf8"
CIF_LIST=""
USE_RFD=0
PARTIAL_NOISE=5
BATCH_SIZE=10
HOTSPOTS=""
BINDER_CHAIN="A"
TARGET_CHAIN="B"
MPNN_CONTEXT="complex"
DRY_RUN=0
CIFS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)          NAME="$2"; shift 2 ;;
        --cif-list)      CIF_LIST="$2"; shift 2 ;;
        --mosmo)         TARGET="mosmo"; shift ;;
        --megf8)         TARGET="megf8"; shift ;;
        --rfdiffusion)   USE_RFD=1; shift ;;
        --partial-noise) PARTIAL_NOISE="$2"; shift 2 ;;
        --batch-size)    BATCH_SIZE="$2"; shift 2 ;;
        --hotspots)      HOTSPOTS="$2"; shift 2 ;;
        --mpnn-context)  MPNN_CONTEXT="$2"; shift 2 ;;
        --binder-chain)  BINDER_CHAIN="$2"; shift 2 ;;
        --target-chain)  TARGET_CHAIN="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       usage ;;
        -*)              echo "Unknown argument: $1" >&2; usage ;;
        *)               CIFS+=("$1"); shift ;;
    esac
done

if [[ -z "$NAME" || -z "$TARGET" ]]; then
    echo "Error: missing required argument(s)." >&2
    usage
fi

if [[ ${#CIFS[@]} -eq 0 && -z "$CIF_LIST" ]]; then
    echo "Error: no input CIFs given (pass paths positionally or via --cif-list)." >&2
    usage
fi

if [[ -n "$CIF_LIST" && ! -f "$CIF_LIST" ]]; then
    echo "Error: --cif-list '$CIF_LIST' not found." >&2
    exit 1
fi

if ! [[ "$BATCH_SIZE" =~ ^[0-9]+$ ]]; then
    echo "Error: --batch-size must be a positive integer, got '$BATCH_SIZE'" >&2
    exit 1
fi

# rfd3 declares partial_t as Optional[float] with ge=0.0, so accept decimals.
if ! [[ "$PARTIAL_NOISE" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: --partial-noise must be a non-negative number, got '$PARTIAL_NOISE'" >&2
    exit 1
fi

if [[ $USE_RFD -eq 0 && -n "$HOTSPOTS" ]]; then
    echo "Warning: --hotspots is ignored without --rfdiffusion." >&2
fi

if [[ "$MPNN_CONTEXT" != "binder" && "$MPNN_CONTEXT" != "complex" ]]; then
    echo "Error: --mpnn-context must be 'binder' or 'complex', got '$MPNN_CONTEXT'" >&2
    exit 1
fi

MODE="proteinmpnn"
[[ $USE_RFD -eq 1 ]] && MODE="rfdiffusion + proteinmpnn"

n_positional=${#CIFS[@]}
n_listed=0
[[ -n "$CIF_LIST" ]] && n_listed=$(grep -cve '^[[:space:]]*$' -e '^#' "$CIF_LIST" || true)

echo "== Redesign config =="
echo "  name:       $NAME"
echo "  target:     $TARGET"
echo "  mode:       $MODE"
echo "  mpnn sees:  $MPNN_CONTEXT"
echo "  inputs:     $((n_positional + n_listed)) cif(s)"
if [[ $USE_RFD -eq 1 ]]; then
    echo "  partial_t:  $PARTIAL_NOISE A of noise"
    echo "  batch_size: $BATCH_SIZE"
    echo "  hotspots:   ${HOTSPOTS:-<from each interface>}"
fi
echo "====================="

# ---------------------------------------------------------------------------
# Preflight: compute every output path up front and refuse to run if any of
# them already exist. Same contract as batch_design.sh — a naming collision
# must fail before any GPU time is burned, not after.
# ---------------------------------------------------------------------------
extract_dir="${REDESIGN_IN_DIR}/${TARGET}_redesign_${NAME}"
rfd_json="${RFD_JSON_DIR}/${TARGET}_rfd_${NAME}${SUFFIX}.json"
design_dir="${DESIGN_OUT_DIR}/${TARGET}_rfd_${NAME}${SUFFIX}"
mpnn_in_dir="${extract_dir}_mpnn_in"        # only used in --rfdiffusion mode
manifest="${extract_dir}_manifest.csv"
mpnn_dir="${MPNN_OUT_DIR}/${TARGET}_mpnn_${NAME}${SUFFIX}"
fold_json="${FOLD_JSON_DIR}/${TARGET}_fold_${NAME}${SUFFIX}.json"
fold_dir="${FOLD_OUT_DIR}/${TARGET}_fold_${NAME}${SUFFIX}"

conflicts=()
[[ -e "$extract_dir" ]] && conflicts+=("$extract_dir")
[[ -e "$manifest" ]]    && conflicts+=("$manifest")
[[ -e "$mpnn_dir" ]]    && conflicts+=("$mpnn_dir")
[[ -e "$fold_json" ]]   && conflicts+=("$fold_json")
[[ -e "$fold_dir" ]]    && conflicts+=("$fold_dir")
if [[ $USE_RFD -eq 1 ]]; then
    [[ -e "$rfd_json" ]]    && conflicts+=("$rfd_json")
    [[ -e "$design_dir" ]]  && conflicts+=("$design_dir")
    [[ -e "$mpnn_in_dir" ]] && conflicts+=("$mpnn_in_dir")
fi

if [[ ${#conflicts[@]} -gt 0 ]]; then
    echo "Error: refusing to run — the following already exist:" >&2
    printf '  %s\n' "${conflicts[@]}" >&2
    echo "Choose a different --name, or remove/rename the conflicting path(s)." >&2
    exit 1
fi

# The stage roots are shared with batch_design.sh and must already exist — if
# one is missing, the EDIT ME block above is wrong and silently creating it
# would scatter outputs. REDESIGN_IN_DIR is the exception: it is new to this
# script, so it gets created on first use.
for parent in "$RFD_JSON_DIR" "$DESIGN_OUT_DIR" "$MPNN_OUT_DIR" "$FOLD_JSON_DIR" "$FOLD_OUT_DIR"; do
    if [[ ! -d "$parent" ]]; then
        echo "Error: expected parent directory '$parent' does not exist. Not creating it — check the EDIT ME paths at the top of this script." >&2
        exit 1
    fi
done
mkdir -p "$REDESIGN_IN_DIR"

if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] would write:"
    echo "  extract dir: $extract_dir"
    [[ $USE_RFD -eq 1 ]] && echo "  rfd json:    $rfd_json"
    [[ $USE_RFD -eq 1 ]] && echo "  design dir:  $design_dir"
    echo "  mpnn dir:    $mpnn_dir"
    echo "  fold json:   $fold_json"
    echo "  fold dir:    $fold_dir"
    exit 0
fi

mkdir "$extract_dir"
mkdir "$mpnn_dir"
mkdir "$fold_dir"
[[ $USE_RFD -eq 1 ]] && mkdir "$design_dir"
echo "-> reserved run directories: $extract_dir, $mpnn_dir, $fold_dir"

if [[ "$TARGET" == "mosmo" ]]; then
    jsoner="$PYTHON_3"
else
    jsoner="$PYTHON_4"
fi

# ---------------------------------------------------------------------------
# Step 1: pull the binder out of each cofolded input.
#
# In ProteinMPNN mode the binder is extracted alone. In RFdiffusion mode the
# target comes along too, because partial diffusion has to see the target to
# keep the binding pose it is meant to be correcting.
# ---------------------------------------------------------------------------
extract_args=(--out_dir "$extract_dir"
              --suffix "$SUFFIX"
              --binder_chain "$BINDER_CHAIN"
              --target_chain "$TARGET_CHAIN"
              --manifest "$manifest")
if [[ $USE_RFD -eq 1 ]]; then
    # RFD mode: plain .cif complexes; rfd3 reads either, and the json builder
    # is simpler against uncompressed input.
    extract_args+=(--keep_target)
else
    # MPNN mode: these go straight to mpnn.sh, which globs *.cif.gz ONLY.
    extract_args+=(--gzip)
    # With complex context the target rides along and mpnn.sh's default
    # --designed_chain A redesigns only the binder. The fasta then comes back
    # as binder+target concatenated, which the jsoner trims at the junction.
    [[ "$MPNN_CONTEXT" == "complex" ]] && extract_args+=(--keep_target)
fi
[[ -n "$CIF_LIST" ]] && extract_args+=(--cif_list "$CIF_LIST")

"${PY_RUN[@]}" "$PY_EXTRACT" "${extract_args[@]}" ${CIFS[@]+"${CIFS[@]}"}

shopt -s nullglob
extracted=("$extract_dir"/*.cif "$extract_dir"/*.cif.gz)
shopt -u nullglob
if [[ ${#extracted[@]} -eq 0 ]]; then
    echo "Error: extract_binder.py reported success but '$extract_dir' has no .cif/.cif.gz files." >&2
    exit 1
fi
echo "-> extracted binders: $extract_dir"

# ---------------------------------------------------------------------------
# Steps 2-4 (--rfdiffusion only): partial diffusion, then re-extract the
# corrected binder so ProteinMPNN sees the same binder-only input it would
# have seen in the default mode.
# ---------------------------------------------------------------------------
if [[ $USE_RFD -eq 1 ]]; then
    rfd_json_args=(--input_dir "$extract_dir"
                   --output_json "$rfd_json"
                   --partial_noise "$PARTIAL_NOISE"
                   --binder_chain "$BINDER_CHAIN"
                   --target_chain "$TARGET_CHAIN")
    [[ -n "$HOTSPOTS" ]] && rfd_json_args+=(--hotspots "$HOTSPOTS")

    "${PY_RUN[@]}" "$PY_RFD_JSON" "${rfd_json_args[@]}"

    if [[ ! -f "$rfd_json" ]]; then
        echo "Error: expected RFdiffusion json '$rfd_json' not found." >&2
        exit 1
    fi
    echo "-> rfd json: $rfd_json"

    bash "$BASH_1" --out_dir "$design_dir" --inputs "$rfd_json" --batch_size "$BATCH_SIZE"

    if [[ -z "$(ls -A "$design_dir" 2>/dev/null)" ]]; then
        echo "Error: rfd3.sh reported success but '$design_dir' is empty." >&2
        exit 1
    fi
    echo "-> design dir: $design_dir"

    # rfd3 writes binder+target complexes; strip back to the binder so MPNN
    # gets a monomer, matching the default mode. Backbones come out as either
    # .cif or .cif.gz depending on the rfd3 build, so collect both.
    shopt -s nullglob
    backbones=("$design_dir"/*.cif "$design_dir"/*.cif.gz)
    shopt -u nullglob

    if [[ ${#backbones[@]} -eq 0 ]]; then
        echo "Error: no .cif/.cif.gz backbones found in '$design_dir'." >&2
        exit 1
    fi

    mkdir "$mpnn_in_dir"
    "${PY_RUN[@]}" "$PY_EXTRACT" \
        --out_dir "$mpnn_in_dir" \
        --suffix "" \
        --name_from file \
        --binder_chain "$BINDER_CHAIN" \
        --target_chain "$TARGET_CHAIN" \
        --gzip \
        --manifest "${mpnn_in_dir}_manifest.csv" \
        "${backbones[@]}"

    echo "-> mpnn input dir: $mpnn_in_dir (${#backbones[@]} backbones)"
else
    mpnn_in_dir="$extract_dir"
fi

# ---------------------------------------------------------------------------
# Step: mpnn.sh --input_dir <binder cifs> --output_dir <mpnn_dir>
# (designed_chain defaults to A, and the inputs are binder-only here, so it is
# intentionally never passed — same as batch_design.sh)
# ---------------------------------------------------------------------------
bash "$BASH_2" --input_dir "$mpnn_in_dir" --output_dir "$mpnn_dir"

if [[ -z "$(ls -A "$mpnn_dir" 2>/dev/null)" ]]; then
    echo "Error: mpnn.sh reported success but '$mpnn_dir' is empty." >&2
    exit 1
fi
echo "-> mpnn dir: $mpnn_dir"

# ---------------------------------------------------------------------------
# Step: mosmo_jsoner.py / megf8_jsoner.py --fasta_dir <mpnn_dir> --output_json <fold_json>
#
# The jsoner re-attaches the target as a {"path": ...} component, so folding
# the redesigned binder back against the real target works even though MPNN
# only ever saw the binder. Its junction-trimming branch is a no-op on
# binder-only fastas, which is exactly what we want here.
# ---------------------------------------------------------------------------
fasta_dir="$mpnn_dir"

"${PY_RUN[@]}" "$jsoner" --fasta_dir "$fasta_dir" --output_json "$fold_json"

if [[ ! -f "$fold_json" ]]; then
    echo "Error: expected fold-prep output '$fold_json' not found." >&2
    exit 1
fi
echo "-> fold json: $fold_json"

# ---------------------------------------------------------------------------
# Step: fold.sh --out_dir <fold_dir> --inputs <fold_json>
# ---------------------------------------------------------------------------
bash "$BASH_3" --out_dir "$fold_dir" --inputs "$fold_json"

if [[ -z "$(ls -A "$fold_dir" 2>/dev/null)" ]]; then
    echo "Error: fold.sh reported success but '$fold_dir' is empty." >&2
    exit 1
fi
echo "-> fold dir: $fold_dir"

echo "== Redesign complete =="
echo "  extract dir: $extract_dir"
if [[ $USE_RFD -eq 1 ]]; then
echo "  rfd json:    $rfd_json"
echo "  design dir:  $design_dir"
fi
echo "  mpnn dir:    $mpnn_dir"
echo "  fold json:   $fold_json"
echo "  fold dir:    $fold_dir"
