import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--list_csv", required=True, help="Path to all fold confidences csv")
parser.add_argument("--output_csv", required=True, help="Path to output csv")
args = parser.parse_args()

mpnn_outputs = Path("/mmfs1/gscratch/stf/igem/all_out/all_bookkeep/mpnn_runs")

passing_binders = pd.read_csv(args.list_csv)
output = pd.DataFrame()

blacklist = [
    "megf8_rfd_Trevor_Set2_megf8_Trevor_Set2_0_model_35.cif_b0_d0",
    "megf8_rfd_eva0810_2_megf8_eva0810_2_0_model_10.cif_b1_d0_redesign.cif_b0_d0",
    "megf8_rfd_eliza81_megf8_eliza81_0_model_2.cif_b4_d0",
    "megf8_rfd_eliza92_megf8_eliza92_0_model_12.cif_b0_d0",
    "megf8_rfd_eva0810_2_megf8_eva0810_2_0_model_8.cif_b2_d0_redesign.cif_b2_d0"
]

target_pdb = "/igem/targets/MOSMO_chainb.pdb"
mosmo_start = "DKLTIISGCLFLAA"
megf8_junction = "ESFHGSPLGGQQCYRLISVEQE"

def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

sequence_column = []
name_column = []
target_column = []

copy_counter = 0
error_counter = 0
blacklist_counter = 0

for row in passing_binders.itertuples():
    if row.rf3_output not in blacklist:
        current_fa = row.rf3_output.rpartition("cif_b")[0] + "cif.fa"
        mpnn_folder = row.fold_run.replace("_fold_", "_mpnn_", 1)
        mpnn_path = mpnn_outputs / mpnn_folder / current_fa

        with open(mpnn_path) as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]

        copied = False
        for i in range(0, len(lines), 2):
            header = lines[i].lstrip(">")
            header = header.split(",")[0].strip()
            seq = lines[i+1]
            target = lines[i+1]

            if mosmo_start in seq:
                target = seq[seq.index(mosmo_start):]
                seq = seq[:seq.index(mosmo_start)]

            if megf8_junction in seq:
                target = seq[seq.index(megf8_junction):]
                seq = seq[:seq.index(megf8_junction)]

            if row.rf3_output in header:
                sequence_column.append(seq)
                target_column.append(target)
                if("redesign" in row.rf3_output):
                    count = 0
                    while f"{row.rf3_output}_{count}" in name_column:
                        count = count+1
                    name_column.append(f"{row.rf3_output}_{count}")
                else:
                    name_column.append(f"{row.rf3_output}")
                copied = True
        
        if copied:
            print(f"Copied {row.rf3_output} sequence to csv")
            copy_counter = copy_counter + 1
        else:
            print(f"ERROR: Problem copying {row.rf3_output} to destination")
            error_counter = error_counter + 1
    else:
        print(f"Blacklist: {row.rf3_output}")
        blacklist_counter = blacklist_counter + 1

output['Name'] = name_column
output['Seq'] = sequence_column
output['Target'] = target_column

output.to_csv(args.output_csv, index=False)

print(f"Binder Summary - Successfully Copied: {copy_counter}, Errors: {error_counter}, Blacklisted: {blacklist_counter}")
