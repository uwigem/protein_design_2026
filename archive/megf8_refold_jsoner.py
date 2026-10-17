import json
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rf3_names_list", required=True,
                     help="Plain text file, one RF3 output name to refold per line (rewritten by user each run)")
parser.add_argument("--mpnn_out_dir", default="./all_out/all_bookkeep/mpnn_runs",
                     help="Base MPNN output dir (default matches MPNN_OUT_DIR; assumes script is run from igem/)")
parser.add_argument("--target_pdb", required=True,
                     help="Path/name of the abbreviated MEGF8 target PDB")
parser.add_argument("--cache_file", default="all_out/megf8_refold/refold_cache.txt",
                     help="Plain text cache of already-processed RF3 names")
args = parser.parse_args()

# output json is named after the input list: NAME.txt -> megf8_refold_NAME.json
list_basename = os.path.splitext(os.path.basename(args.rf3_names_list))[0]
output_json = os.path.join("all_out/megf8_refold", f"megf8_refold_{list_basename}.json")

MEGF8_JUNCTION = "ESFHGSPLGGQQCYRLISVEQE"


def get_run_id(rf3_name):
    # run ID = token(s) between "rfd_" and the next "megf8_"
    after_rfd = rf3_name.split("rfd_", 1)[1]
    run_id = after_rfd.split("megf8_", 1)[0].rstrip("_")
    return run_id


def get_fasta_path(rf3_name, mpnn_out_dir):
    run_id = get_run_id(rf3_name)
    mpnn_dir = os.path.join(mpnn_out_dir, f"megf8_mpnn_{run_id}")
    fasta_name = rf3_name.split(".cif", 1)[0] + ".cif.fa"
    return os.path.join(mpnn_dir, fasta_name)


def parse_fasta(fasta_path):
    """Returns dict: header_name -> sequence"""
    records = {}
    with open(fasta_path) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
    for i in range(0, len(lines), 2):
        header = lines[i].lstrip(">").split(",")[0].strip()
        seq = lines[i + 1]
        records[header] = seq
    return records


# --- load input list ---
with open(args.rf3_names_list) as f:
    requested_names = [l.strip() for l in f.readlines() if l.strip()]

# --- load cache ---
if os.path.exists(args.cache_file):
    with open(args.cache_file) as f:
        cache = set(l.strip() for l in f.readlines() if l.strip())
else:
    cache = set()

to_process = [n for n in requested_names if n not in cache]
skipped_cached = len(requested_names) - len(to_process)

# --- group by fasta path ---
by_fasta = {}
for name in to_process:
    fasta_path = get_fasta_path(name, args.mpnn_out_dir)
    by_fasta.setdefault(fasta_path, []).append(name)

entries = []
newly_cached = []
errors = []

for fasta_path, names in by_fasta.items():
    if not os.path.exists(fasta_path):
        errors.append((None, f"FASTA not found: {fasta_path}", names))
        continue

    records = parse_fasta(fasta_path)

    for name in names:
        if name not in records:
            errors.append((name, f"header not found in {fasta_path}", None))
            continue

        seq = records[name]

        if MEGF8_JUNCTION not in seq:
            errors.append((name, "junction motif not found in sequence", None))
            continue

        seq = seq[:seq.index(MEGF8_JUNCTION)]

        entries.append({
            "name": name + "_abbrev",
            "components": [
                {
                    "seq": seq,
                    "chain_id": "A"
                },
                {
                    "path": args.target_pdb
                }
            ],
            "template_selection": ["B"]
        })
        newly_cached.append(name)

# --- write output json (this run's new entries only) ---
os.makedirs(os.path.dirname(output_json), exist_ok=True)
with open(output_json, "w") as f:
    json.dump(entries, f, indent=4)

# --- update cache (successful entries only) ---
cache.update(newly_cached)
os.makedirs(os.path.dirname(args.cache_file), exist_ok=True)
with open(args.cache_file, "w") as f:
    for name in sorted(cache):
        f.write(name + "\n")

# --- report ---
print(f"Output JSON: {output_json}")
print(f"Requested: {len(requested_names)}")
print(f"Skipped (already cached): {skipped_cached}")
print(f"Newly processed: {len(newly_cached)}")
print(f"Errors: {len(errors)}")
for name, reason, group in errors:
    if group is not None:
        print(f"  [FASTA missing] {reason} (affects: {', '.join(group)})")
    else:
        print(f"  [{name}] {reason}")
