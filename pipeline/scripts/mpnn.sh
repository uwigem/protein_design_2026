#!/bin/bash

# defaults
INPUT_DIR=""
OUTPUT_DIR=""
DESIGNED_CHAIN="A"

# parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --designed_chain)
            DESIGNED_CHAIN="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# validate
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: bash eliza_mpnn.sh --input_dir <path> --output_dir <path> [--designed_chain <chain>]"
    exit 1
fi

for cif in ${INPUT_DIR}/*.cif.gz; do
  apptainer exec \
    --bind .:/igem \
    ./containers/foundry.sif \
    mpnn \
    --structure_path /igem/${cif} \
    --out_directory /igem/${OUTPUT_DIR} \
    --checkpoint_path /weights/proteinmpnn_v_48_020.pt \
    --model_type protein_mpnn \
    --is_legacy_weights True \
    --designed_chains ${DESIGNED_CHAIN} \
    --number_of_batches 3
done
