#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Orchestrates the binder design pipeline:
#   1. RFdiffusion step (python_1.py / python_2.py) -> writes <target>_rfd_<name>.json
#   2. rfd3.sh: post-process rfd json into a directory of designs
#   3. mpnn.sh: run ProteinMPNN over that directory -> <target>_mpnn_<name> dir
#   4. RF3 json-prep step (mosmo_jsoner.py / megf8_jsoner.py) -> <target>_fold_<name>.json
#   5. fold.sh: run RF3 fold over that json -> <target>_fold_<name> dir
#
# where <target> is "mosmo" or "megf8" depending on the --mosmo/--megf8 flag.
#
# NOTE: --batch-size is only consumed by rfd3.sh (Step 2's post-processing),
# not by the RFdiffusion generation step itself (Step 1).

set -euo pipefail

on_error() {
    echo "Error: pipeline aborted — command on line ${1} exited non-zero." >&2
    exit 1
}
trap 'on_error $LINENO' ERR

# ---------------------------------------------------------------------------
# EDIT ME per run — script locations
# ---------------------------------------------------------------------------
PYTHON_1="./RFDiffusion/scripts/mosmo_json_rfd.py"      # mosmo RFdiffusion step
PYTHON_2="./RFDiffusion/scripts/megf8_json_rfd.py"      # megf8 RFdiffusion step
PYTHON_3="./RosettaFold/scripts/mosmo_jsoner.py"        # mosmo RF3-json-prep step
PYTHON_4="./RosettaFold/scripts/megf8_jsoner.py"        # megf8 RF3-json-prep step
BASH_1="./RFDiffusion/scripts/rfd3.sh"                  # post-process rfd json -> design dir
BASH_2="./ProteinMPNN/scripts/mpnn.sh"                  # run ProteinMPNN over design dir
BASH_3="./RosettaFold/scripts/fold.sh"                  # run RF3 fold over fold json

# ---------------------------------------------------------------------------
# EDIT ME per run — relative base directories for each stage's I/O.
# Actual dir/file names are still built as <base>/<target>_<stage>_<name>
# further down, so you only need to change the roots here.
# ---------------------------------------------------------------------------
RFD_JSON_DIR="./all_out/all_bookkeep/all_json"    # where <target>_rfd_<name>.json is written/read
DESIGN_OUT_DIR="./all_out/all_bookkeep/rfd_runs"  # rfd3.sh --out_dir root
MPNN_OUT_DIR="./all_out/all_bookkeep/mpnn_runs"   # mpnn.sh --output_dir root
FOLD_JSON_DIR="./all_out/all_bookkeep/all_json"   # where <target>_fold_<name>.json is written/read
FOLD_OUT_DIR="./all_out/all_bookkeep/fold_runs"   # fold.sh --out_dir root

usage() {
    cat <<EOF
Usage: $(basename "$0") --name NAME --batch-size BATCH_SIZE --hotspots HOTSPOTS (--mosmo | --megf8)

  --name NAME            base name used to construct output filenames/dirs
  --batch-size BATCH     integer batch size passed to rfd3.sh (post-processing step)
  --hotspots HOTSPOTS    hotspot residue spec, e.g. "A35,A50,A123"
  --mosmo                target MOSMO for this run
  --megf8                target MEGF8 for this run

Exactly one of --mosmo / --megf8 must be given.
EOF
    exit 1
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
NAME=""
BATCH_SIZE=""
HOTSPOTS=""
TARGET=""   # "mosmo" or "megf8"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)
            NAME="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --hotspots)
            HOTSPOTS="$2"
            shift 2
            ;;
        --mosmo)
            TARGET="mosmo"
            shift
            ;;
        --megf8)
            TARGET="megf8"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            ;;
    esac
done

if [[ -z "$NAME" || -z "$BATCH_SIZE" || -z "$HOTSPOTS" || -z "$TARGET" ]]; then
    echo "Error: missing required argument(s)." >&2
    usage
fi

if ! [[ "$BATCH_SIZE" =~ ^[0-9]+$ ]]; then
    echo "Error: --batch-size must be a positive integer, got '$BATCH_SIZE'" >&2
    exit 1
fi

echo "== Run config =="
echo "  name:       $NAME"
echo "  batch_size: $BATCH_SIZE"
echo "  hotspots:   $HOTSPOTS"
echo "  target:     $TARGET"
echo "================="

# ---------------------------------------------------------------------------
# Preflight: compute every output path up front and refuse to run if any of
# them already exist. This must happen BEFORE any expensive step (RFdiffusion,
# ProteinMPNN, RF3) so a naming collision fails immediately instead of after
# burning compute. It also reserves the three run directories right away.
# ---------------------------------------------------------------------------
rfd_json="${RFD_JSON_DIR}/${TARGET}_rfd_${NAME}.json"
design_dir="${DESIGN_OUT_DIR}/${TARGET}_rfd_${NAME}"
mpnn_dir="${MPNN_OUT_DIR}/${TARGET}_mpnn_${NAME}"
fold_json="${FOLD_JSON_DIR}/${TARGET}_fold_${NAME}.json"
fold_dir="${FOLD_OUT_DIR}/${TARGET}_fold_${NAME}"

conflicts=()
[[ -e "$rfd_json" ]]   && conflicts+=("$rfd_json")
[[ -e "$design_dir" ]] && conflicts+=("$design_dir")
[[ -e "$mpnn_dir" ]]   && conflicts+=("$mpnn_dir")
[[ -e "$fold_json" ]]  && conflicts+=("$fold_json")
[[ -e "$fold_dir" ]]   && conflicts+=("$fold_dir")

if [[ ${#conflicts[@]} -gt 0 ]]; then
    echo "Error: refusing to run — the following already exist:" >&2
    printf '  %s\n' "${conflicts[@]}" >&2
    echo "Choose a different --name, or remove/rename the conflicting path(s)." >&2
    exit 1
fi

for parent in "$RFD_JSON_DIR" "$DESIGN_OUT_DIR" "$MPNN_OUT_DIR" "$FOLD_JSON_DIR" "$FOLD_OUT_DIR"; do
    if [[ ! -d "$parent" ]]; then
        echo "Error: expected parent directory '$parent' does not exist. Not creating it — check the EDIT ME paths at the top of this script." >&2
        exit 1
    fi
done

mkdir "$design_dir"
mkdir "$mpnn_dir"
mkdir "$fold_dir"
echo "-> reserved run directories: $design_dir, $mpnn_dir, $fold_dir"

# ---------------------------------------------------------------------------
# Step 1: RFdiffusion — produces <target>_rfd_<name>.json
# ---------------------------------------------------------------------------
if [[ "$TARGET" == "mosmo" ]]; then
    python3 "$PYTHON_1" --name "$NAME" --hotspots "$HOTSPOTS" --output_json "$rfd_json"
    jsoner="$PYTHON_3"
else
    python3 "$PYTHON_2" --name "$NAME" --hotspots "$HOTSPOTS" --output_json "$rfd_json"
    jsoner="$PYTHON_4"
fi

if [[ ! -f "$rfd_json" ]]; then
    echo "Error: expected RFdiffusion output '$rfd_json' not found." >&2
    exit 1
fi
echo "-> RFdiffusion json: $rfd_json"

# ---------------------------------------------------------------------------
# Step 2: rfd3.sh --out_dir <design_dir> --inputs <rfd_json> --batch_size <n>
# ---------------------------------------------------------------------------
bash "$BASH_1" --out_dir "$design_dir" --inputs "$rfd_json" --batch_size "$BATCH_SIZE"

if [[ -z "$(ls -A "$design_dir" 2>/dev/null)" ]]; then
    echo "Error: rfd3.sh reported success but '$design_dir' is empty." >&2
    exit 1
fi
echo "-> design dir: $design_dir"

# ---------------------------------------------------------------------------
# Step 3: mpnn.sh --input_dir <design_dir> --output_dir <mpnn_dir>
# (designed_chain always defaults to A when this script runs end to end,
# so it's intentionally never passed here)
# ---------------------------------------------------------------------------
bash "$BASH_2" --input_dir "$design_dir" --output_dir "$mpnn_dir"

if [[ -z "$(ls -A "$mpnn_dir" 2>/dev/null)" ]]; then
    echo "Error: mpnn.sh reported success but '$mpnn_dir' is empty." >&2
    exit 1
fi
echo "-> mpnn dir: $mpnn_dir"

# ---------------------------------------------------------------------------
# Step 4: mosmo_jsoner.py/megf8_jsoner.py --fasta_dir <fasta_dir> --output_json <fold_json>
# ---------------------------------------------------------------------------
# TODO: confirm whether mpnn.sh writes fasta files directly into
# --output_dir, or into a subdirectory (e.g. "$mpnn_dir/seqs"). Adjust
# fasta_dir below if it's the latter.
fasta_dir="$mpnn_dir"

python3 "$jsoner" --fasta_dir "$fasta_dir" --output_json "$fold_json"

if [[ ! -f "$fold_json" ]]; then
    echo "Error: expected fold-prep output '$fold_json' not found." >&2
    exit 1
fi
echo "-> fold json: $fold_json"

# ---------------------------------------------------------------------------
# Step 5: fold.sh --out_dir <fold_dir> --inputs <fold_json>
# ---------------------------------------------------------------------------
bash "$BASH_3" --out_dir "$fold_dir" --inputs "$fold_json"

if [[ -z "$(ls -A "$fold_dir" 2>/dev/null)" ]]; then
    echo "Error: fold.sh reported success but '$fold_dir' is empty." >&2
    exit 1
fi
echo "-> fold dir: $fold_dir"

echo "== Pipeline complete =="
echo "  rfd json:   $rfd_json"
echo "  design dir: $design_dir"
echo "  mpnn dir:   $mpnn_dir"
echo "  fold json:  $fold_json"
echo "  fold dir:   $fold_dir"
