# Protein Design 2026

## Overview
This repository contains the design and validation workflow used by the Washington iGEM Protein Design Team to generate novel minibinder proteins against target structures such as MOSMO and MEGF8. The project is organized around a repeatable batch-design pipeline: generate backbones with RFdiffusion, refine sequences with ProteinMPNN, prepare RF3 inputs, fold the complexes, and then evaluate the resulting binders.

The active workflow is centered in the `pipeline/` directory and is meant to be run from the project root. The pipeline is automated, it takes inital user input and runs through the 3 steps. The `archive/` folder contains earlier exploratory scripts and one-off utilities that are kept for historical context, while the `redesign/` folder contains a second-pass optimization loop for improving promising binder candidates after an initial fold.

### Typical design pipeline

```mermaid
flowchart LR
    A[Target structure + hotspot residues] --> B[RFdiffusion JSON builder]
    B --> C[RFdiffusion backbone generation]
    C --> D[rfd3.sh: convert designs to structure set]
    D --> E[ProteinMPNN sequence redesign]
    E --> F[RF3 input JSON prep]
    F --> G[fold.sh: RosettaFold structure prediction]
    G --> H[Evaluate binders and select promising hits]
    H --> I[redesign/ second-pass optimization]
    I --> J[Refined binder candidates]
```

### Repository layout

- `archive/`: historical scripts used during earlier project iterations. Useful for reference, but not the main active workflow.
- `pipeline/`: primary automated binder-generation pipeline and the scripting layer that orchestrates each stage.
- `redesign/`: post-processing and redesign utilities for refining already-folded binder candidates.

### Quick start
The main end-to-end batch generation script is `pipeline/batch_design.sh`. It expects a target name, batch size, hotspot specification, and either `--mosmo` or `--megf8`.

Example:

```bash
bash pipeline/batch_design.sh \
  --name my_run \
  --batch-size 16 \
  --hotspots "A35,A50,A123" \
  --mosmo
```

This workflow generates a set of RFdiffusion candidates, converts them into design directories, runs ProteinMPNN, prepares RosettaFold inputs, and folds the final complexes.

The directories of the input and output files can be edited at the top of `batch_design.sh` to match any file structure.

## Archive
The `archive/` directory contains older scripts that were used during earlier generations of the project and may still be useful for debugging or reproducing legacy results. These files are not part of the current primary workflow and should be treated as historical reference material.

Current contents:

- `archive/megf8_refold_jsoner.py`: legacy JSON preparation utility for MEGF8 refolding workflows.
- `archive/monomer_jsoner.py`: earlier script for monomer preparation and input conversion.
- `archive/passing_sequence_extraction.py`: helper for extracting or filtering sequence candidates that passed an earlier stage.

### How to use this folder
- Use it for reference when checking old assumptions or reproducing a previous analysis.
- Do not rely on it for the active design pipeline unless you are intentionally reproducing legacy analysis.
- If a task is not clearly documented elsewhere, inspect the script headers before using it in production runs.

## Pipeline
The `pipeline/` directory is the heart of the active binder design workflow. It contains the orchestration script and the helper scripts that generate the RFdiffusion input JSON files used by the external structure-generation tools.

### Folder structure

- `pipeline/batch_design.sh`: main batch pipeline controller. It runs the full five-step design process in order.
- `pipeline/scripts/`: helper scripts used by the orchestration layer.
  - `megf8_json_rfd.py`: builds the RFdiffusion input JSON for the MEGF8 target.
  - `mosmo_json_rfd.py`: builds the RFdiffusion input JSON for the MOSMO target.

### Pipeline purpose and operating logic
The pipeline follows this pattern:

1. Build a target-specific RFdiffusion JSON input.
2. Generate backbone candidates using RFdiffusion.
3. Post-process the diffusion outputs into a design directory.
4. Run ProteinMPNN to redesign sequences for the generated backbone.
5. Prepare a RosettaFold-compatible JSON file.
6. Run the final fold step to produce folded complexes for evaluation.

The wrapper script enforces a clear naming convention for all major artifacts:

- `<target>_rfd_<name>.json`: diffusion input or output JSON
- `<target>_rfd_<name>/`: RFdiffusion-generated design directory
- `<target>_mpnn_<name>/`: ProteinMPNN output directory
- `<target>_fold_<name>.json`: fold-prep JSON
- `<target>_fold_<name>/`: final folded outputs

### Running the pipeline
From the repository root, run:

```bash
bash pipeline/batch_design.sh --name <run_name> --batch-size <N> --hotspots "A35,A50,A123" --mosmo
```

or:

```bash
bash pipeline/batch_design.sh --name <run_name> --batch-size <N> --hotspots "B35,B50,B123" --megf8
```

Notes:

- Use a unique `--name` for every run to avoid collisions.
- `--hotspots` should be supplied as a comma-separated list of target residues. Use chain `A` for MOSMO and chain `B` for MEGF8.
- The script expects the relevant external tooling (RFdiffusion, ProteinMPNN, RosettaFold) to be available in the surrounding project environment, as referenced by the wrapper paths.

## Redesign
The `redesign/` directory is for iterative improvement of already-folded candidates. Instead of starting from an unbound target and inventing a binder from scratch, the redesign workflow takes promising complex structures and refines them to improve interface compatibility, binder geometry, or sequence quality.

### Current files and purpose

- `redesign/scripts/redesign.sh`: the main second-pass design workflow. It refines binder complexes after initial folds and can operate in either ProteinMPNN-only mode or a partial RFdiffusion + ProteinMPNN mode.
- `redesign/scripts/extract_binder.py`: extracts the binder chain (and optionally the target chain) from an RF3 or RFD3 structure so it can be re-used downstream.
- `redesign/select_redesign.py`: selects the best redesign candidates from a set of generated outputs. Main script to run the redesign.
- `redesign/scripts/compile_success_confidences.py`: aggregates confidence metrics from redesign or folding runs to identify successful candidates.
- `redesign/batch_confidences_fixed.py`: batch-processing script for confidence scoring and evaluation across many outputs. Can be run separately to generate a CSV with minPAE, plDDT, and RMSD thresholds calculated.
- `redesign/scripts/tmtest_batch_confidences_fixed.py`: evaluates candidate structures with TM-score-based checks. Used explicitly for MOSMO redesign.
- `redesign/scripts/refold_monomer.sh`: helper for monomer refolding tasks or monomer-focused evaluation.
- `redesign/redesign_all_dirs.sh`: helper script that loops through all generated binders and runs the redesign script on each one.

### Redesign workflow
The redesign loop is designed to answer a simple question, once a candidate binder has already been folded into a complex, how can we improve it without restarting the whole design search from scratch?

Typical operation:

1. Take a folded binder-target complex.
2. Extract the binder and, if needed, the target context.
3. Re-run ProteinMPNN or a light RFdiffusion correction step.
4. Rebuild the fold-prep JSON and refold the complex.
5. Compare metrics and select the strongest candidate.

This is especially useful when the first-generation design is close to correct but needs a second pass to tighten binding geometry, sequence quality, or structural confidence.

### Example use
```bash
python3 scripts/select_redesign.py \
  --input_dir all_out/all_bookkeep/fold_runs/[[fold_run_folder]] \
  --min_pae_floor 2 --min_pae 6 --rms_floor 2 --rms 3 \
  --name [[run_name]] \
  --confidences_csv all_out/batch_results_confidences.csv \
  --report scripts/_logs/[[run_name]]_selection.csv"
```

The limits listed for min_pae and rms (RMSD) are suggestions that ensure the redesigned binders have the maximum potential for success. However, all metrics can be adjusted as needed.

---

This repository is intended to be used as a modular protein-design workflow. Archive for legacy context, pipeline for initial binder generation, and redesign for refinement of promising designs. This process is applicable to a wide variety of cilliary targets and is intended to be adapted for generalized therapuetic design.
