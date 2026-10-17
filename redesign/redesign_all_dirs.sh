#!/usr/bin/env bash
# Loops through all generated binders in the directory given below.
# Auto runs the redesign script to quickly redesign all binders.

FOLD_RUNS="all_out/all_bookkeep/fold_runs"
MEGF8_CSV="all_out/batch_results_confidences_plDDT_fixed.csv"
MOSMO_CSV="all_out/tmtest_batch_results_confidences_fixed.csv"

counter=2

for dir in "$FOLD_RUNS"/*/; do
    sub=$(basename "$dir")

    if [[ "$sub" == *mosmo* ]]; then
        csv="$MOSMO_CSV"
    else
        csv="$MEGF8_CSV"
    fi

    echo "=== $sub (csv: $csv) ==="
    python3 scripts/select_redesign.py \
        --input_dir "$FOLD_RUNS/$sub" \
        --min_pae_floor 0 --min_pae 3 \
        --rms_floor 0 --rms 6 \
        --name "defne_$counter" \
        --confidences_csv "$csv" \
        --report "scripts/_logs/${sub}_preview.csv"
    echo

    counter=$((counter + 1))
done