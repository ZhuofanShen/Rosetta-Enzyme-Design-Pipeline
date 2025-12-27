import argparse
import json
import os
import csv
from pyrosetta import *

def load_scores_from_fasc(file_name):
    decoy_scores_str = "["
    with open(file_name) as fasc:
        for line in fasc:
            decoy_scores_str += line.rstrip() + ","
    decoy_scores_str = decoy_scores_str[:-1] + "]"
    return json.loads(decoy_scores_str)

def extract_n_decoys(scores, n=1, is_reversed=False):
    scores = sorted(scores, key=lambda e: e["total_score"], reverse=is_reversed)
    best = scores[:n]
    return {
        str(s["decoy"]): (
            s["total_score"]
            - s["res_type_constraint"]
            - s["coordinate_constraint"]
        )
        for s in best
    }

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=str)
parser.add_argument("-f", "--fold", type=str)
parser.add_argument(
    "-sub", "--substrate", type=str,
    default="substrates/cyclopropanation_styrene_EDA"
)
parser.add_argument(
    "-stereo", "--stereoisomer",
    required=True,
    choices=["1R2R", "1S2S", "1R2S", "1S2R"]
)
parser.add_argument("-clean", "--clean_decoys", action="store_true")
args = parser.parse_args()

substrate = args.substrate.strip("/").split("/")[-1]
stereo_abbr = args.stereoisomer[1] + args.stereoisomer[3]

# ---------------- header ----------------
header = ["PDB_ID", "Site", "Classification", "WT"]
for rot in range(1, 5):
    for ester in ["+", "-"]:
        for score_type in ["TS", "ddG"]:
            header.append(f"{stereo_abbr}{rot}{ester}_{score_type}")

rows = [header]

# ---------------- main loop ----------------
for pdb in sorted(os.listdir(args.directory)):
    pdb_dir = os.path.join(args.directory, pdb)
    if not os.path.isdir(pdb_dir):
        continue

    if args.fold:
        if not os.path.isfile(os.path.join(pdb_dir, args.fold + ".fold")):
            continue

    # fold categories
    fold_categories = [
        f.replace(".fold", "").lower()
        for f in os.listdir(pdb_dir)
        if f.endswith(".fold")
    ]

    # WT score
    with open(os.path.join(pdb_dir, pdb + "_relaxed.pdb")) as pf:
        for line in pf:
            if line.startswith("pose"):
                scores = line[5:].split()
                break
    apo_score = float(scores[-1]) - float(scores[-13])

    # ---- one row per site ----
    for site in ["distal", "proximal"]:
        pdb_sub_path = os.path.join(
            pdb_dir, f"{substrate}_{site}"
        )
        if not os.path.isdir(pdb_sub_path):
            continue

        row = [
            pdb,
            site,
            ";".join(fold_categories),
            apo_score,
        ]

        for rot in range(1, 5):
            for ester in ["p", "m"]:
                stereo_state = f"{pdb}_{args.stereoisomer}-rot{rot}{ester}"
                design_path = os.path.join(
                    pdb_sub_path, "FastDesign", stereo_state
                )

                if not os.path.isdir(design_path):
                    row.extend(["", ""])
                    continue

                fasc_file = next(
                    (f for f in os.listdir(design_path) if f.endswith(".fasc")),
                    None
                )
                if not fasc_file:
                    row.extend(["", ""])
                    continue

                scores = load_scores_from_fasc(
                    os.path.join(design_path, fasc_file)
                )
                decoy, score = next(iter(extract_n_decoys(scores).items()))

                decoy_path = os.path.join(design_path, decoy)

                if not os.path.isfile(decoy_path):
                    decoy_name = list(filter(lambda x: x.endswith(".pdb"), os.listdir(design_path)))[0]
                    decoy_path = os.path.join(design_path, decoy_name)

                if args.clean_decoys:
                    for f in os.listdir(design_path):
                        if f.endswith(".pdb") and f != decoy:
                            os.remove(os.path.join(design_path, f))

                with open(decoy_path) as pf:
                    ts_info = pf.readlines()[-5].split()
                ts_score = float(ts_info[-1]) - float(ts_info[-2])
                ddg = score - apo_score

                row.extend([round(ts_score, 2), ddg])

        rows.append(row)

# ---------------- write CSV ----------------
if not args.fold:
    args.fold = "default"

out_csv = f"{args.directory}_{args.fold}_{args.stereoisomer}.csv"
with open(out_csv, "w", newline="") as f:
    csv.writer(f).writerows(rows)

print(f"Saved CSV: {out_csv}")
