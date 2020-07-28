import argparse
import os
from pymol import *


def collect_output_match_info(directory):
    match_dict = dict()
    matches = os.listdir(directory)
    for match in matches:
        if match.endswith('.pdb'):
            match_info = match.split('_')
            match_info[-1] = match_info[-1][:-4]
            position = match_info[2]
            if not (position in match_dict and int(match_dict[position][-1]) > int(match_info[-1])):
                match_dict[position] = match_info
    return match_dict

def read_enz_des_remark_lines(pdb):
    enzdes_cst_remarks = list()
    with open(pdb, 'r') as p_pdb:
        flag = False
        for line in p_pdb:
            if line.startswith('REMARK 666 MATCH TEMPLATE'):
                if flag == False:
                    flag = True
                enzdes_cst_remarks.append(line)
            elif flag == True:
                break
    return enzdes_cst_remarks

def exclude_d_amino_acids(pdb, enzdes_cst_remarks):
    flag = True
    cmd.load(pdb)
    prefix = '/' + cmd.get_object_list()[0] + '//'
    for remark in enzdes_cst_remarks:
        position_info = list(filter(lambda x : x != '', remark.split(' ')))
        atom_name_prefix = prefix + position_info[-6] + '/' + position_info[-4] + '/'
        dihedral = cmd.get_dihedral(atom_name_prefix + 'N', atom_name_prefix + 'CA', \
            atom_name_prefix + 'CB', atom_name_prefix + 'C')
        if dihedral < 0:
            flag = False
            break
    cmd.delete('*')
    return flag

def get_match_point_mutations(enzdes_cst_remarks):
    point_mutation_dict = dict()
    point_mutation_dict['template'] = enzdes_cst_remarks[0][28:31]
    for enzdes_cst_remark in enzdes_cst_remarks:
        point_mutation_dict[enzdes_cst_remark[49] + str(int(enzdes_cst_remark[56:59]))] = enzdes_cst_remark[51:54]
    return point_mutation_dict

def make_relax_input_files(directory, match_dict, multimer=False):
    for position, match_info_list in match_dict.items():
        match_prefix = '_'.join(match_info_list[:-1])
        match = match_prefix + '_1.pdb'
        remarks = read_enz_des_remark_lines(directory + '/' + match)
        if exclude_d_amino_acids(directory + '/' + match, remarks):
            has_protein_scaffold = False
            for i in range(1, int(match_info_list[-1]) + 1):
                accumulated_objs = len(cmd.get_object_list())
                match = match_prefix + '_' + str(i) + '.pdb'
                cmd.load(directory + '/' + match)
                scaffold_new_name = match[:-4] + '_0001'
                cmd.set_name(cmd.get_object_list()[accumulated_objs], scaffold_new_name)
                if has_protein_scaffold:
                    cmd.delete(scaffold_new_name)
                else:
                    has_protein_scaffold = True
                mutations = get_match_point_mutations(remarks)
                for mutated_res in mutations.values():
                    cmd.remove('chain A & resi 1 & resn ' + mutated_res)
                cmd.alter('resn ' + mutations['template'], 'resi="0"')
            objs = cmd.get_object_list()
            output_name = match_info_list[2]
            if multimer:
                cmd.create('main_chain', 'chain ' + remarks[0][49])
                cmd.create('duplicate_chain', 'chain ' + remarks[1][49])
                cmd.align('duplicate_chain', 'main_chain')
                cmd.save(output_name + '.pdb', 'main_chain')
                cmd.save(output_name + '_duplicate.pdb', 'duplicate_chain')
                cmd.delete('main_chain')
                cmd.delete('duplicate_chain')
                mutation_lines = list()
                with open(output_name + '_duplicate.pdb', 'r') as p_pdb:
                    for line in p_pdb:
                        if line.startswith('ATOM') and line[17:22] == remarks[1][51:54] + \
                            ' ' + remarks[1][49] and int(line[23:26]) == int(remarks[1][56:59]):
                            mutation_lines.append(line[:21] + remarks[0][49] + line[22:])
                os.remove(output_name + '_duplicate.pdb')
                with open(output_name + '.pdb', 'r+') as p_pdb:
                    lines = p_pdb.readlines()
                    p_pdb.seek(0)
                    written_mutation_lines = False
                    for line in lines:
                        if line.startswith('ATOM') and line[21] == remarks[0][49] and \
                            int(line[23:26]) == int(remarks[1][56:59]):
                            if not written_mutation_lines:
                                p_pdb.writelines(mutation_lines)
                                written_mutation_lines = True
                        elif not line.startswith('CONECT'):
                            p_pdb.write(line)
                cmd.load(output_name + '.pdb')
                cmd.save(output_name + '.pdb', output_name)
                cmd.delete('output_name')
            else:
                cmd.save(output_name + '.pdb', objs[0])
            with open(output_name + '.pdb', 'r+') as p_pdb:
                lines = p_pdb.readlines()
                p_pdb.seek(0)
                p_pdb.write('MODEL    1\n')
                p_pdb.write('HEADER                                            27-JUL-20   XXXX              \n')
                p_pdb.writelines(remarks)
                p_pdb.writelines(lines[:-1])
                if_first_model = True
                for i in range(1, len(objs)):
                    cmd.save(objs[i] + '.pdb', objs[i])
                    if if_first_model:
                        if_first_model = False
                    else:
                        p_pdb.write('MODEL    ' + str(i) + '\n')
                    p_pdb.write('HEADER                                            xx-MMM-xx                     \n')
                    with open(objs[i] + '.pdb', 'r') as rotamer:
                        written_ter = False
                        for line in rotamer:
                            if line.startswith('HETATM'):
                                p_pdb.write(line)
                            elif not written_ter:
                                p_pdb.write('TER                                                                             \n')
                                written_ter = True
                            elif line.startswith('CONECT'):
                                p_pdb.write(line)
                    p_pdb.write('ENDMDL\n')
                    os.remove(objs[i] + '.pdb')
                p_pdb.write(lines[-1])
            cmd.delete('*')
        else:
            continue
        break


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_directory', type=str)
    parser.add_argument('-m', '--multimer', action='store_true')
    args = parser.parse_args()
    match_dict = collect_output_match_info(args.input_directory)
    make_relax_input_files(args.input_directory, match_dict, args.multimer)
