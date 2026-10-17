#!/bin/bash

# defaults
INPUTS=""
OUT_DIR=""

# parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --inputs)
            INPUTS="$2"
            shift 2
            ;;
        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# validate
if [[ -z "$INPUTS" || -z "$OUT_DIR" ]]; then
    echo "Usage: bash eliza_rf3.sh --inputs <path> --out_dir <path>"
    exit 1
fi

export DISABLE_CUEQUIVARIANCE=1

apptainer exec --nv \
  --bind .:/igem \
  ./containers/foundry_fixed.sif \
  rf3 fold \
  inputs=/igem/${INPUTS} \
  out_dir=/igem/${OUT_DIR} \
  ckpt_path=/weights/rf3_foundry_01_24_latest_remapped.ckpt \
  diffusion_batch_size=1 \
  num_steps=50 \
  early_stopping_plddt_threshold=0.5
