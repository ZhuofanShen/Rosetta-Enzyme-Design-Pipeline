import argparse
import json
import os
from pyrosetta import *
import xlwt

def load_scores_from_fasc(variant):
    file_names = os.listdir(variant)
    for file_name in file_names:
        if file_name.endswith(".fasc"):
            decoy_scores_str = "["
            with open(variant + "/" + file_name) as fasc:
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
parser.add_argument('dir', type=str)
parser.add_argument('-fold', type=str)
args = parser.parse_args()
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
row = 0
for pdb in filter(lambda x: os.path.isdir(args.dir + '/' + x + '/cyclopropanation_EDA_styrene') \
    and (not args.fold or os.path.isfile(args.dir + '/' + x + '/' + args.fold + '.fold')), sorted(os.listdir(args.dir))):
    with open(args.dir + '/' + pdb + '/' + pdb + '_relax.pdb', 'r') as pf:
        for line in pf:
            if line.startswith('pose'):
                scores = line[5:-1].split(' ')
    apo_score = float(scores[-1]) - float(scores[-13])
    design_sheet.write(row, 0, pdb)
    design_sheet.write(row, 1, apo_score)
    col = 2
    for pos in ['BELOW', 'ABOVE']:
        if os.path.isdir(args.dir + '/' + pdb + '/cyclopropanation_EDA_styrene/' + pos):
            for stereo in ['1R2S']: # for stereo in ['1R2R', '1S2S', '1R2S', '1S2R']:
                for rot in range(1, 5):
                    for ester in ['+', '-']:
                        path = args.dir + '/' + pdb + '/cyclopropanation_EDA_styrene/' + \
                            pos + '/' + pdb + '_' + pos + '/' + pdb + '_' + pos + '_' + stereo + '-rot' + str(rot) + ester
                        if os.path.isfile(path + '/' + pdb + '_' + pos + '_' + stereo + '-rot' + str(rot) + ester + '.fasc'):
                            ssd_scores = load_scores_from_fasc(path)
                            for name, score in extract_n_decoys(ssd_scores).items():
                                break
                            for decoy in os.listdir(path):
                                if decoy == name:
                                    with open(path + '/' + decoy, 'r') as pf:
                                        ts_score_info = pf.readlines()[-5][:-1].split(' ')
                                    ts_score = float(ts_score_info[-1]) - float(ts_score_info[-2])
                                elif not decoy.endswith('.fasc') and not decoy.endswith('.sh') and not decoy.endswith('.in_progress'):
                                    os.remove(path + '/' + decoy)
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
workbook.save(args.dir + '_' + args.fold + '.xls')
