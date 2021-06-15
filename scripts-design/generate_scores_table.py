#!/usr/bin/python3
import argparse
import json
import os
from pyrosetta import *
import xlrd, xlwt
from xlutils.copy import copy

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str, help='$scaffold_$substrate')
    parser.add_argument('-s', '--step', type=int, choices=[1, 2], help='step1 or step2')
    parser.add_argument('-dup', '--duplicate', action='store_true', help='Duplicated variant sites in the symmetric chain.')
    parser.add_argument('-mem', '--memory', type=int, default=2000, help='slurm job memory allocation')
    args = parser.parse_args()
    # scaffold substrate
    scaffold_substrate = args.directory.split('_')
    args.scaffold = scaffold_substrate[0] # CPG2-AB
    args.protein = args.scaffold[:args.scaffold.find('-')] # CPG2
    args.substrate = scaffold_substrate[1] # pCaaF-product
    args.ligand = args.substrate[:args.substrate.find('-')] # pCaaF
    # symmetry file
    args.symmetry = None
    for protein_symm in filter(lambda x: x.startswith(args.protein) and x.endswith('.symm'), os.listdir('../../proteins/' + args.protein)):
        args.symmetry = protein_symm
        break
    # pdb file
    for protein_pdb in filter(lambda x: x.startswith(args.protein) and x.endswith('.pdb'), os.listdir('../../proteins/' + args.protein)):
        args.pdb = protein_pdb
        break
    # params files
    params_files = '-extra_res_fa '
    params_files_py = '-params '
    for params in filter(lambda x: x.endswith('.params'), os.listdir('../../ligands/' + args.ligand + '/' + args.substrate)):
        params_files += '../../ligands/' + args.ligand + '/' + args.substrate + '/' +  params + ' '
        if params.endswith('_ligand.params'):
            params_files_py += '../' + params + ' '
        else:
            params_files_py += '../../../../../ligands/' + args.ligand + '/' + args.substrate + '/' + params + ' '
    for params in filter(lambda x: x.endswith('.params'), os.listdir('../../proteins/' + args.protein)):
        params_files += '../../proteins/' + args.protein + '/' + params + ' '
        params_files_py += '../../../../../ligands/' + args.protein + '/' + params + ' '
    return args, params_files, params_files_py

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

def insert_ligand_residue_into_reference_sequence(reference_sequence, sequence, annotation_dict, ligand_res_name):
    extra_res_index_list = list()
    for index, annotation in annotation_dict.items():
        if annotation[:annotation.find(':')] == ligand_res_name:
            extra_res_index_list.append(int(index)) # sequence string index
    for index in sorted(extra_res_index_list):
        reference_sequence = reference_sequence[:index] + sequence[index] + reference_sequence[index:]
    return reference_sequence

def get_best_decoy_and_mutations(path, initial_ref_seq, ligand_res_name, symmetric):
    # get the file name of the best decoy
    scores = load_scores_from_fasc(path)
    best_decoy_scores = extract_n_decoys(scores)
    for best_decoy, best_score in best_decoy_scores.items():
        break
    # get the sequence of the best decoy
    pose = pose_from_pdb(path + '/' + best_decoy)
    seq, annotation_dict = read_annotated_sequence(pose.annotated_sequence(), symmetric)
    # insert the ligand residue into the original reference pose
    ref_seq = insert_ligand_residue_into_reference_sequence(initial_ref_seq, seq, annotation_dict, ligand_res_name)
    # compare the sequences and rename the best decoy with the point mutation information
    design_pose_numbering = list()
    design_pdb_numbering = list()
    for index in range(len(ref_seq)):
        if seq[index] != ref_seq[index]:
            design_pose_numbering.append(ref_seq[index] + str(index + 1) + seq[index])
            [pdb_index, chain] = pose.pdb_info().pose2pdb(index + 1).strip(' ').split(' ')
            design_pdb_numbering.append(chain + ref_seq[index] + pdb_index + seq[index])
    return design_pose_numbering, design_pdb_numbering, best_decoy, best_score

def move_ligand_position(path_to_variant, best_decoy, ligand_res_name):
    with open(path_to_variant + '/design/' + best_decoy, 'r') as pdb:
        ligand_lines = list()
        for line in pdb:
            if ( line.startswith('ATOM') or line.startswith('HETATM') ) and line[17:20] == ligand_res_name:
                ligand_lines.append(line)
    with open(path_to_variant + '/' + path_to_variant.split('/')[-1] + '.pdb', 'r') as pdb:
        lines = pdb.readlines()
    with open(path_to_variant + '/' + path_to_variant.split('/')[-1] + '-ligand.pdb', 'w') as pdb:
        current_ligand_chain = None
        for line in lines:
            if line[17:20] == ligand_res_name:
                if current_ligand_chain != line[21]:
                    current_ligand_chain = line[21]
                    for ligand_line in ligand_lines:
                        if ligand_line[21] == line[21]:
                            pdb.write(ligand_line)
            elif not line.startswith('CONECT'):
                pdb.write(line)

def read_match_res_scores(path_to_best_decoy):
    with open(path_to_best_decoy, 'r') as pdb:
        match_res_score = dict()
        for line in pdb:
            if line[3:7] == ':MP-' and line[:3] not in match_res_score:
                scores = line.split(' ')
                match_res_score[line[:3]] = (str(float(scores[-1]) - float(scores[-2]) - float(scores[-13])), scores[-14], scores[-12], scores[-11])
        return match_res_score

def write_initial_design_scores(arguments, params_files_py):
    # create a new xls workbook
    workbook = xlwt.Workbook(encoding="ascii")
    directory_sheet = workbook.add_sheet("directory")
    initial_design_sheet = workbook.add_sheet("initial_design")
    manual_design_sheet = workbook.add_sheet("manual_design")
    #workbook.add_sheet("greedy_opt_design")
    # get the reference sequence
    for ligand_params in filter(lambda x: x.endswith('_ligand.params'), os.listdir('../../ligands/' + args.ligand + '/' + args.substrate)):
        ligand_res_name = ligand_params[:-14]
    ref_pose = pose_from_pdb('../../proteins/' + arguments.protein + '/' + arguments.pdb)
    initial_ref_seq, ref_annotation_dict = read_annotated_sequence(ref_pose.annotated_sequence(), arguments.symmetry)
    # iterate through all designs
    row = 0
    for variant in sorted(filter(lambda x: x.startswith('X') and not x.endswith('_deprecated') and \
            os.path.isdir(arguments.directory + '/' + x + '/design'), os.listdir(arguments.directory))):
        # get the best decoy and score, and the designed residues
        design_pose_numbering, design_pdb_numbering, best_decoy, best_score = \
            get_best_decoy_and_mutations(arguments.directory + '/' + variant + '/design', initial_ref_seq, ligand_res_name, arguments.symmetry)
        # reformat the best decoy
        move_ligand_position(arguments.directory + '/' + variant, best_decoy, ligand_res_name)
        # write to the score table
        directory_sheet.write(row, 0, variant)
        initial_design_sheet.write(row, 0, arguments.scaffold + '_' + '_'.join(design_pose_numbering))
        initial_design_sheet.write(row, 1, round(best_score, 2))
        initial_design_sheet.write(row, 2, best_decoy[best_decoy.rfind('_') + 1: -4])
        manual_design_sheet.write(row, 0, arguments.scaffold + '_' + '_'.join(design_pdb_numbering))
        match_res_scores = read_match_res_scores(arguments.directory + '/' + variant + '/design/' + best_decoy)
        col = 3
        for match_res in sorted(match_res_scores.keys()):
            for match_res_score_term in match_res_scores[match_res]:
                initial_design_sheet.write(row, col, match_res_score_term)
                col += 1
        # write the relax script
        if len(design_pose_numbering) > 2:
            muts_arg = ' -muts'
            for mutation in filter(lambda x: not x.endswith('X') and not x.endswith('Z'), design_pose_numbering):
                muts_arg += ' ' + mutation[1:-1] + ',' + mutation[-1]
        else:
            muts_arg = ''
        if arguments.duplicate:
            suffix = '_duplicated.cst'
        else:
            suffix = '.cst'
        with open(arguments.directory + '/' + variant + '/manual_design.sh', 'w+') as bash:
            if arguments.symmetry:
                symmetry_arg = ' -symm ../../../../../proteins/' + arguments.protein + '/' + arguments.symmetry
            else:
                symmetry_arg = ''
            bash.write('mkdir manual_design;\ncd manual_design;\nslurmit.py --job ' + variant + ' --mem ' + str(arguments.memory) + 
' --command "python ../../../../../scripts-design/fast_design.py ../' + variant + '-ligand.pdb' + 
symmetry_arg + ' -sf ref2015_cst --score_terms fa_intra_rep_nonprotein:0.545 \
fa_intra_atr_nonprotein:1 ' + params_files_py + '-enzdes_cst ../../../../../ligands/' + 
arguments.ligand + '/' + arguments.substrate + '/' + arguments.substrate + 
suffix + muts_arg + ' -nbh 8.0 -n 50;"\ncd ..;\n')
        with open(arguments.directory + '/' + variant + '/pymol.txt', 'w+') as pymol:
            resi = str()
            for mutation in design_pdb_numbering:
                resi += mutation[2:-1] + '+'
            pymol.write('show sticks, resi ' + resi[:-1] + '; hide (h.);')
        # continue to deal with the next design
        row = row + 1
    # save to a xls file
    workbook.save(arguments.directory + '/' + arguments.directory + '.xls')

def get_relaxed_design_scores(arguments):
    # open the existing xls workbook
    workbook_read = xlrd.open_workbook(arguments.directory + '/' + arguments.directory + '.xls')
    directory_sheet = workbook_read.sheet_by_index(0)
    workbook_write = copy(workbook_read)
    design_sheet = workbook_write.get_sheet(2) # step
    # get the reference sequence
    for ligand_params in filter(lambda x: x.endswith('_ligand.params'), os.listdir('../../ligands/' + args.ligand + '/' + args.substrate)):
        ligand_res_name = ligand_params.strip('_ligand.params')
    ref_pose = pose_from_pdb('../../proteins/' + arguments.protein + '/' + arguments.pdb)
    initial_ref_seq, ref_annotation_dict = read_annotated_sequence(ref_pose.annotated_sequence(), arguments.symmetry)
    # iterate through all designs
    for row in range(directory_sheet.nrows):
        variant = directory_sheet.cell_value(row, 0)
        path = arguments.directory + '/' + variant + '/' + 'manual_design'
        if not os.path.isdir(path):
            continue
        print(variant)
        # get the best decoy and score, and the designed residues
        design_pose_numbering, design_pdb_numbering, best_decoy, best_score = \
            get_best_decoy_and_mutations(path, initial_ref_seq, ligand_res_name, arguments.symmetry)
        # write to the score table
        design_sheet.write(row, 0, arguments.scaffold + '_'  + '_'.join(design_pdb_numbering))
        design_sheet.write(row, 1, str(round(best_score, 2)))
        design_sheet.write(row, 2, best_decoy.split("_")[-1][:-4])
        match_res_scores = read_match_res_scores(arguments.directory + '/' + variant + '/manual_design/' + best_decoy)
        col = 3
        for match_res in sorted(match_res_scores.keys()):
            for match_res_score_term in match_res_scores[match_res]:
                design_sheet.write(row, col, match_res_score_term)
                col += 1
    workbook_write.save(arguments.directory + '/' + arguments.directory + '.xls')

if __name__ == "__main__":
    args, params_files, params_files_py = parse_arguments()
    init(params_files)
    if args.step == 1:
        write_initial_design_scores(args, params_files_py)
    elif args.step == 2:
        get_relaxed_design_scores(args)
