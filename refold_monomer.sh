set -euo pipefail

on_error() {
    echo "Error: pipeline aborted  ~@~T command on line ${1} exited non-zero." >&2
    exit 1
}
trap 'on_error $LINENO' ERR

# ---------------------------------------------------------------------------

jsoner="./scripts/monomer_jsoner.py"                            # path to json generation script

fold_script="./RosettaFold/scripts/fold.sh"                     # path to RosettaFold script

confidences_csv="./all_out/batch_results_confidences_plDDT_fixed.csv"       # where binder confidences are written/read

mpnn_out_dir="./all_out/all_bookkeep/mpnn_runs"                 # where mpnn outputs are written/read

fold_json_dir="./all_out/all_bookkeep/all_json"                 # where refold_<name>.json is written/read

fold_out_dir="./all_out/monomer_refold"  # where the refolded binders are written


fold_json="${fold_json_dir}/monomer_refold.json"                # json being written to for folding

python3 "$jsoner" --confidences_csv "$confidences_csv" --mpnn_dir "$mpnn_out_dir" --output_json "$fold_json"

if [[ ! -f "$fold_json" ]]; then
    echo "Error: expected fold-prep output '$fold_json' not found." >&2
    exit 1
fi
echo "-> fold json: $fold_json"

bash "$fold_script" --out_dir "$fold_out_dir" --inputs "$fold_json"

if [[ -z "$(ls -A "$fold_out_dir" 2>/dev/null)" ]]; then
    echo "Error: fold.sh reported success but '$fold_out_dir' is empty." >&2
    exit 1
fi
echo "-> fold dir: $fold_out_dir"
