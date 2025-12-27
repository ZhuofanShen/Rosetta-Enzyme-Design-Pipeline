import argparse
import json
import os
import csv
from pyrosetta import *

def load_scores_from_fasc(file_name):
    decoy_scores_str = "["
    with open(file_name) as fasc:
        for line in fasc:
            decoy_scores_str += line[:-1]
            decoy_scores_str += ","
    decoy_scores_str = decoy_scores_str[:-1] + "]"
    return json.loads(decoy_scores_str)

def extract_n_decoys(scores, n=1, is_reversed=False):
    scores = sorted(scores, key=lambda entry: entry["total_score"], reverse=is_reversed)
    stable_scores = {}
    for score in scores[:n]:
        stable_scores[str(score["decoy"])] = (
            score["total_score"]
            - score["res_type_constraint"]
            - score["coordinate_constraint"]
        )
    return stable_scores

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str)
parser.add_argument('-fold', type=str)
parser.add_argument("-sub", "--substrate", type=str,
                    default="substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-stereo", "--stereoisomer", type=str, required=True,
                    choices=["1R2R", "1S2S", "1R2S", "1S2R"])
parser.add_argument("-clean", "--clean_decoys", action="store_true")
args = parser.parse_args()

substrate = args.substrate.strip("/").split("/")[-1]
stereoisomer_abbr = args.stereoisomer[1] + args.stereoisomer[3]

# ---- build header ----
header = ["PDB_ID", "Classification", "WT"]
for rot in range(1, 5):
    for ester in ['+', '-']:
        for score_type in ["TS", "ddG"]:
            header.append(f"{stereoisomer_abbr}{rot}{ester}_{score_type}")

rows = [header]

# ---- main loop ----
for pdb in filter(
    lambda x: not args.fold or
    os.path.isfile(os.path.join(args.directory, x, args.fold + '.fold')),
    sorted(os.listdir(args.directory))
):
    fold_categories = []
    pdb_dir = os.path.join(args.directory, pdb)

    for fold in filter(lambda x: x.endswith(".fold"), os.listdir(pdb_dir)):
        fold_categories.append(fold.strip(".fold").lower())

    with open(os.path.join(pdb_dir, pdb + '_relaxed.pdb'), 'r') as pf:
        for line in pf:
            if line.startswith('pose'):
                scores = line[5:-1].split(' ')
    apo_score = float(scores[-1]) - float(scores[-13])

    row = [pdb, ";".join(fold_categories), apo_score]

    for open_site in ['distal', 'proximal']:
        pdb_sub_path = os.path.join(pdb_dir, f"{substrate}_{open_site}")
        if not os.path.isdir(pdb_sub_path):
            row.extend([None] * 16)
            continue

        for rot in range(1, 5):
            for ester in ['+', '-']:
                stereo_state = f"{pdb}_{args.stereoisomer}-rot{rot}{ester}"
                design_path = os.path.join(pdb_sub_path, "FastDesign", stereo_state)

                fasc_file = next(
                    (f for f in os.listdir(design_path) if f.endswith(".fasc")),
                    None
                )

                if not fasc_file:
                    row.extend([None, None])
                    continue

                ssd_scores = load_scores_from_fasc(
                    os.path.join(design_path, fasc_file)
                )

                decoy_name, score = next(
                    iter(extract_n_decoys(ssd_scores).items())
                )

                decoy_path = os.path.join(design_path, decoy_name)

                if not os.path.isfile(decoy_path):
                    decoy_name = list(filter(lambda x: x.endswith(".pdb"), os.listdir(design_path)))[0]
                    decoy_path = os.path.join(design_path, decoy_name)

                if args.clean_decoys:
                    for f in os.listdir(design_path):
                        if (
                            f.endswith(".pdb")
                            and f != decoy_name
                            and not f.endswith(".fasc")
                            and not f.endswith(".sh")
                        ):
                            os.remove(os.path.join(design_path, f))

                with open(decoy_path, 'r') as pf:
                    ts_info = pf.readlines()[-5].strip().split(' ')
                ts_score = float(ts_info[-1]) - float(ts_info[-2])

                score_diff = score - apo_score

                row.extend([round(ts_score, 2), score_diff])

    rows.append(row)

# ---- write CSV ----
if not args.fold:
    args.fold = "default"

out_csv = f"{args.directory}_{args.fold}_{args.stereoisomer}.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(rows)

print(f"Saved CSV to: {out_csv}")
