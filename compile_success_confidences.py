import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--confidences_csv", required=True, help="Path to all fold confidences csv")
parser.add_argument("--monomer_csv", required=True, help="Path to monomer confidences csv")
parser.add_argument("--megf8_csv", required=True, help="Path to monomer confidences csv")
parser.add_argument("--tm_csv", required=True, help="Path to tm filter csv")
parser.add_argument("--output_csv", required=True, help="Path to output csv")
args = parser.parse_args()

all_confidences = pd.read_csv(args.confidences_csv)
monomer_confidences = pd.read_csv(args.monomer_csv)
tm_confidences = pd.read_csv(args.tm_csv)
megf8_confidences = pd.read_csv(args.megf8_csv)

min_pae_threshold = 2.0
rmsd_threshold = 2.0
plddt_threshold = 80.0

filtered_confidences = all_confidences[(all_confidences['min_pae_cross'] <= min_pae_threshold) & 
                                       (all_confidences['ca_rmsd_A'] <= rmsd_threshold) & 
                                       (all_confidences['mean_plddt_A'] >= plddt_threshold)]

monomer_column = []

for row in filtered_confidences.itertuples():
    try:
        monomer_column.append(monomer_confidences.loc[monomer_confidences['rf3_output'] == row.rf3_output, 'ca_rmsd_A'].values[0])
    except:
        print(f"No monomer rmsd for {row.rf3_output}")

filtered_confidences['monomer_ca_rmsd'] = monomer_column


tm_column = []

for row in filtered_confidences.itertuples():
    try:
        tm_column.append(tm_confidences.loc[tm_confidences['rf3_output'] == row.rf3_output, 'tm_passing'].values[0])
        # passed = tm_confidences.loc[tm_confidences['cif_path'].str.partition("/")[0] == row.rf3_output, 'cif_path'].values[0]
        # tm_column.append("pass")
    except:
        print(f"Failed reading tm filter for {row.rf3_output}")
        # tm_column.append("fail")

filtered_confidences['tm_filter'] = tm_column


megf8_refold_column = []

filtered_megf8_confidences = all_confidences[(all_confidences['min_pae_cross'] <= min_pae_threshold) & (all_confidences['ca_rmsd_A'] <= rmsd_threshold)]

for row in filtered_confidences.itertuples():
    try:
        # megf8_refold_column.append(megf8_confidences.loc[megf8_confidences['rf3_output'] == row.rf3_output, 'tm_passing'].values[0])
        if "megf8" in row.rf3_output:
            passed = filtered_megf8_confidences.loc[filtered_megf8_confidences['rf3_output'].str.partition("_abbrev")[0] == row.rf3_output, 'rf3_output'].values[0]
            megf8_refold_column.append(True)
        else:
            megf8_refold_column.append(False)
    except:
        # print(f"Failed reading tm filter for {row.rf3_output}")
        megf8_refold_column.append(False)

filtered_confidences['megf8_filter'] = megf8_refold_column

master_list = filtered_confidences[(filtered_confidences['monomer_ca_rmsd'] <= rmsd_threshold) & 
                                   (((filtered_confidences['rf3_output'].str.contains("mosmo", na=False)) & (filtered_confidences['tm_filter'])) |
                                    ((filtered_confidences['rf3_output'].str.contains("megf8", na=False)) & (filtered_confidences['megf8_filter'])))]

master_list.to_csv(args.output_csv, index=False)

