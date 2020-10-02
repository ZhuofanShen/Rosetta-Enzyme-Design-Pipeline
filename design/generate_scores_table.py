#!/usr/bin/python3
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

def write_score_table(directory, params_files):
    init('-extra_res_fa ' + ' '.join(params_files))
    scaffold = directory.split('_')[0]
    workbook = xlwt.Workbook(encoding="ascii")
    score_sheet = workbook.add_sheet("score")
    num_sheet = workbook.add_sheet("num")
    row = 0
    for match in sorted(filter(lambda x: x.startswith('X'), os.listdir(directory))):
        # design
        num_sheet.write(row, 0, match)
        for design in os.listdir(directory + '/' + match):
            if design.startswith(scaffold + '_'):
                break
        design_new_name = scaffold
        pose = pose_from_pdb(directory + '/' + match + '/' + design)
        pdb_info = pose.pdb_info()
        pose2pdb_mapping = dict()
        for point_mutation in design[len(scaffold) + 1:-4].split('_'):
            pdb_index_chain = pdb_info.pose2pdb(int(point_mutation[1:-1])).split(' ')
            point_mutation_pdb_index = pdb_index_chain[1] + point_mutation[0] + pdb_index_chain[0] + point_mutation[-1]
            design_new_name += '_' + point_mutation_pdb_index
            pose2pdb_mapping[point_mutation] = point_mutation_pdb_index
        score_sheet.write(row, 0, design_new_name)
        scores = extract_n_decoys(load_scores_from_fasc(directory + '/' + match + '/design'))
        for name, score in scores.items():
            number = int(name.split("_")[-1][:-4])
            break
        score_sheet.write(row, 1, round(score, 2))
        num_sheet.write(row, 1, number)
        # revert
        col = 2
        for reverted_point_mutation in filter(lambda x: x.startswith('revert_'), os.listdir(directory + '/' + match)):
            num_sheet.write(row, col, reverted_point_mutation)
            reverted_point_mutation_pdb_index = pose2pdb_mapping[reverted_point_mutation[7:]]
            score_sheet.write(row, col, reverted_point_mutation_pdb_index)
            col += 1
            scores = extract_n_decoys(load_scores_from_fasc(directory + '/' + match + '/' + reverted_point_mutation))
            for name, score in scores.items():
                number = int(name.split("_")[-1][:-4])
                break
            score_sheet.write(row, col, round(score, 2))
            num_sheet.write(row, col, number)
            col += 1
        row = row + 1
    workbook.save(directory + '/' + directory + '.xls')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str, help='$scaffold_$substrate')
    parser.add_argument('-params', '--params_files', type=str, nargs='*', required=True, \
        help='including params files for both matching residues and protein ligands/cofactors.')
    args = parser.parse_args()
    write_score_table(args.directory, params_files=args.params_files)
