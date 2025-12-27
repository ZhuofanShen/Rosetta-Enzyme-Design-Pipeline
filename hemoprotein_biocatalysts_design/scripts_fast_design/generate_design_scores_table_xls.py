import argparse
import json
import os
from pyrosetta import *
import xlwt

def load_scores_from_fasc(file_name):
    decoy_scores_str = "["
    with open(file_name) as fasc:
        for line in fasc:
            decoy_scores_str += line[:-1]
            decoy_scores_str += ","
    decoy_scores_str = decoy_scores_str[:-1] + "]"
    return json.loads(decoy_scores_str)

def extract_n_decoys(scores, n=1, is_reversed=False):
    # sort decoys considering residue type constraint
    scores = sorted(scores, key=lambda entry: entry["total_score"], reverse=is_reversed)
    stable_scores = dict()
    for score in scores[0:n]:
        # record without considering residue type constraint
        stable_scores.update({str(score["decoy"]): score["total_score"] - \
            score["res_type_constraint"] - score["coordinate_constraint"]})
    return stable_scores

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str)
parser.add_argument('-fold', type=str)
parser.add_argument("-sub", "--substrate", type=str, default="substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-stereo", "--stereoisomer", type=str, required=True, choices=["1R2R", "1S2S", "1R2S", "1S2R"])
parser.add_argument("-clean", "--clean_decoys", action="store_true")
args = parser.parse_args()

substrate = args.substrate.strip("/").split("/")[-1]

workbook = xlwt.Workbook(encoding="ascii")
design_sheet = workbook.add_sheet("design")
red_style = xlwt.XFStyle()
red_font = xlwt.Font()
red_font.colour_index = 2
red_style.font = red_font
orange_style = xlwt.XFStyle()
orange_font = xlwt.Font()
orange_font.colour_index = 52
orange_style.font = orange_font

design_sheet.write(0, 0, "PDB_ID")
design_sheet.write(0, 1, "Classification")
design_sheet.write(0, 2, "WT")
stereoisomer_abbr = args.stereoisomer[1] + args.stereoisomer[3]
col = 3
for rot in range(1, 5):
    for ester in ['+', '-']:
        for score_type in ["TS", "ddG"]:
            design_sheet.write(0, col, stereoisomer_abbr + str(rot) + ester + "_" + score_type)
            col += 1
row = 1
for pdb in filter(lambda x: not args.fold or \
                  os.path.isfile(os.path.join(args.directory, x, args.fold + '.fold')), \
                  sorted(os.listdir(args.directory))):
    fold_categories = list()
    for fold in filter(lambda x: x.endswith(".fold"), os.listdir(os.path.join(args.directory, pdb))):
        fold_categories.append(fold.strip(".fold").lower())
    with open(os.path.join(args.directory, pdb, pdb + '_relaxed.pdb'), 'r') as pf:
        for line in pf:
            if line.startswith('pose'):
                scores = line[5:-1].split(' ')
    apo_score = float(scores[-1]) - float(scores[-13])
    design_sheet.write(row, 0, pdb)
    design_sheet.write(row, 1, ";".join(fold_categories))
    design_sheet.write(row, 2, apo_score)
    col = 3
    for open_coordination_site in ['distal', 'proximal']:
        pdb_sub_path = os.path.join(args.directory, pdb, substrate + "_" + open_coordination_site)
        if not os.path.isdir(pdb_sub_path):
            continue
        for rot in range(1, 5):
            for ester in ['+', '-']:
                stereo_conf_state = pdb + "_" + args.stereoisomer + "-rot" + str(rot) + ester
                design_path = os.path.join(pdb_sub_path, "FastDesign", stereo_conf_state)
                fasc_file = None
                for f in os.listdir(design_path):
                    if f.endswith(".fasc"):
                        fasc_file = f
                        break
                if not fasc_file:
                    col += 2
                    continue
                ssd_scores = load_scores_from_fasc(os.path.join(design_path, fasc_file))
                for name, score in extract_n_decoys(ssd_scores).items():
                    break
                decoy_path = None
                for decoy in filter(lambda x: x.endswith(".pdb"), os.listdir(design_path)):
                    if decoy == name:
                        break
                    elif args.clean_decoys and not decoy.endswith('.fasc') and not decoy.endswith('.sh'):
                        os.remove(decoy_path)
                decoy_path = os.path.join(design_path, decoy)
                with open(decoy_path, 'r') as pf:
                    ts_score_info = pf.readlines()[-5][:-1].split(' ')
                ts_score = float(ts_score_info[-1]) - float(ts_score_info[-2])
                if ts_score >= 0:
                    design_sheet.write(row, col, str(round(ts_score, 2)))
                elif ts_score >= -2:
                    design_sheet.write(row, col, str(round(ts_score, 2)), orange_style)
                else:
                    design_sheet.write(row, col, str(round(ts_score, 2)), red_style)
                # design_sheet.write(row, col, name[:-4].split('_')[-1])
                score_diff = score - apo_score
                design_sheet.write(row, col + 1, score_diff)
                # if score_diff > -20:
                #     design_sheet.write(row, col + 1, score_diff)
                # elif score_diff > -25:
                #     design_sheet.write(row, col + 1, score_diff, orange_style)
                # else:
                #     design_sheet.write(row, col + 1, score_diff, red_style)
                col += 2
    row += 1
if not args.fold:
    args.fold = "default"
workbook.save(args.directory + '_' + args.fold + '_' + args.stereoisomer + '.xls')
