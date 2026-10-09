import json
import argparse
import csv
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--confidences_csv", required=True, help="Path to input confidences csv")
parser.add_argument("--mpnn_dir", required=True, help="Directory containing ProteinMPNN old results and their .fa files")
parser.add_argument("--output_json", required=True, help="Path to output JSON file")
args = parser.parse_args()

target_pdb = "/igem/targets/MOSMO_chainb.pdb"
mosmo_start = "DKLTIISGCLFLAA"

output_dir = Path("/mmfs1/gscratch/stf/igem/all_out/monomer_refold")

entries = []


passing_rows = []

with open(args.confidences_csv, mode='r', newline='', encoding='utf-8') as file:
    reader = csv.DictReader(file)
    
    for row in reader:
        # Filter condition (e.g., keep only rows where Age > 25)
        try:
            minPAE = float(row['min_pae_cross'])
            rmsd = float(row['ca_rmsd_A'])
            plDDT = float(row['mean_plddt_A'])
            file_path = output_dir / row['rf3_output']
            if minPAE <= 2.0 and rmsd <= 2.5 and plDDT >= 80.0 and not file_path.is_dir():
                passing_rows.append(row)
        except ValueError:
            print(f"Error reading csv for row {row['rf3_output']}")

# 2. Use the remaining data
for row in passing_rows:
    fa_file = f"{args.mpnn_dir}/{row['fold_run'][:6]}mpnn{row['fold_run'][10:]}/{row['rf3_output'][:-6]}.fa"
    with open(fa_file) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    for i in range(0, len(lines), 2):
        header = lines[i].lstrip(">")
        header = header.split(",")[0].strip()
        seq = lines[i+1]

        if mosmo_start in seq:
            seq = seq[:seq.index(mosmo_start)]

        if row['rf3_output'] in header:
            entries.append({
                "name": header,
                "components": [
                    {
                        "seq": seq,
                        "chain_id": "A"
                    }
                ]
            })

with open(args.output_json, "w") as f:
    json.dump(entries, f, indent=4)

print(f"Wrote {len(entries)} entries to {args.output_json}")
