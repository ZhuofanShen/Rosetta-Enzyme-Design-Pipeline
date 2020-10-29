#!/usr/bin/python3
import argparse
import json
import os
from pyrosetta import *
import shutil

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str, help='$scaffold_$substrate')
    parser.add_argument('-ref', '--reference_pdb', type=str, required=True)
    parser.add_argument('-params', '--params_files', type=str, nargs='*', required=True, \
        help='including params files for both matching residues and protein ligands/cofactors.')
    parser.add_argument('-linker', '--linker_res_name', type=str)
    parser.add_argument('-symm', '--symmetry', action='store_true')
    return parser.parse_args()

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

def read_annotated_sequence(annotated_sequence, symmetric):
    sequence = str()
    annotation_dict = dict() # key of the dict is sequence string index
    is_annotation = False
    index_diff = 0
    for index in range(len(annotated_sequence)):
        if is_annotation:
            index_diff += 1
            if annotated_sequence[index] == ']':
                annotation_dict.update({str(index - index_diff): annotated_sequence[begin_index + 1:index]})
                is_annotation = False
                if symmetric and annotated_sequence[begin_index + 1:index].endswith(':CtermProteinFull'):
                    break
        else:
            if annotated_sequence[index] == '[':
                begin_index = index
                index_diff += 1
                is_annotation = True
            else:
                sequence += annotated_sequence[index]
    return sequence, annotation_dict

def insert_linker_residue_into_reference_sequence(reference_sequence, sequence, annotation_dict, linker_res_name):
    extra_res_index_list = list()
    for index, annotation in annotation_dict.items():
        if annotation[:annotation.find(':')] == linker_res_name:
            extra_res_index_list.append(int(index)) # sequence string index
    for index in sorted(extra_res_index_list):
        reference_sequence = reference_sequence[:index] + sequence[index] + reference_sequence[index:]
    return reference_sequence

def rename_best_decoy(path_to_variant, initial_ref_seq, linker_res_name, symmetric, output_prefix):
    # get the file name of the best decoy
    scores = load_scores_from_fasc(path_to_variant + '/design')
    stable_scores = extract_n_decoys(scores)
    for name, score in stable_scores.items():
        break
    # get the sequence of the best decoy
    pose = pose_from_pdb(path_to_variant + '/design/' + name)
    seq, annotation_dict = read_annotated_sequence(pose.annotated_sequence(), symmetric)
    # insert the linker residue into the original reference pose
    ref_seq = insert_linker_residue_into_reference_sequence(initial_ref_seq, seq, annotation_dict, linker_res_name)
    # compare the sequences and rename the best decoy with the point mutation information
    for index in range(len(ref_seq)):
        if seq[index] != ref_seq[index]:
            output_prefix += '_' + ref_seq[index] + str(index + 1) + seq[index]
    shutil.copyfile(path_to_variant + '/design/' + name, path_to_variant + '/' + output_prefix + '.pdb')

if __name__ == "__main__":
    # initialize pyrosetta
    args = parse_arguments()
    opts = '-extra_res_fa'
    for params_file in args.params_files:
        opts += ' ' + params_file
    init(opts)
    # use as the prefix to renamed the best decoy
    protein_name = args.directory[:args.directory.find('_')]
    # load initial reference sequence
    ref_pose = pose_from_pdb(args.reference_pdb)
    initial_ref_seq, ref_annotation_dict = read_annotated_sequence(ref_pose.annotated_sequence(), args.symmetry)
    # insert the linker residue into the initial reference sequence and make comparison
    for variant in filter(lambda x: x.startswith('X'), os.listdir(args.directory)):
        # rename the best decoy of each variant
        rename_best_decoy(args.directory + '/' + variant, initial_ref_seq, args.linker_res_name, args.symmetry, protein_name)
