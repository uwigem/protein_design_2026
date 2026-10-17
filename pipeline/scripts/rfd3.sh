#!/bin/bash

# defaults
OUT_DIR=""
INPUTS=""
BATCH_SIZE=10

# parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --inputs)
            INPUTS="$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# validate
if [[ -z "$OUT_DIR" || -z "$INPUTS" ]]; then
    echo "Usage: bash eliza_rfd3.sh --out_dir <path> --inputs <path> [--batch_size <int>]"
    exit 1
fi

apptainer exec --nv \
  --bind .:/igem \
  ./containers/foundry.sif \
  /app/foundry/.venv/bin/python /app/foundry/models/rfd3/src/rfd3/cli.py design \
  out_dir=/igem/${OUT_DIR} \
  inputs=/igem/${INPUTS} \
  ckpt_path=rfd3 \
  diffusion_batch_size=${BATCH_SIZE}
