import json
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_dir", required=True, help="Directory containing .fa files from ProteinMPNN")
parser.add_argument("--output_json", required=True, help="Path to output JSON file")
args = parser.parse_args()

target_pdb = "/igem/targets/MOSMO_chainb.pdb"
mosmo_start = "DKLTIISGCLFLAA"

entries = []

for fa_file in sorted(glob.glob(f"{args.fasta_dir}/*.fa")):
    with open(fa_file) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    for i in range(0, len(lines), 2):
        header = lines[i].lstrip(">")
        header = header.split(",")[0].strip()
        seq = lines[i+1]

        if mosmo_start in seq:
            seq = seq[:seq.index(mosmo_start)]

        entries.append({
            "name": header,
            "components": [
                {
                    "seq": seq,
                    "chain_id": "A"
                },
                {
                    "path": target_pdb
                }
            ],
            "template_selection": ["B"]
        })

with open(args.output_json, "w") as f:
    json.dump(entries, f, indent=4)

print(f"Wrote {len(entries)} entries to {args.output_json}")
