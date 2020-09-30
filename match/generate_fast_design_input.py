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

def exclude_d_amino_acids(pdb):
    flag = True
    cmd.load(pdb)
    prefix = '/' + cmd.get_object_list()[0] + '//'
    enzdes_cst_remarks = read_enz_des_remark_lines(pdb)
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
    for enzdes_cst_remark in enzdes_cst_remarks:
        point_mutation_dict[enzdes_cst_remark[49] + str(int(enzdes_cst_remark[56:59]))] = enzdes_cst_remark[51:54]
    return point_mutation_dict

def make_relax_input_files(directory, match_dict, homomeric=False):
    for position, match_info_list in match_dict.items():
        match_prefix = '_'.join(match_info_list[:-1])
        match = match_prefix + '_1.pdb'
        remarks = read_enz_des_remark_lines(directory + '/' + match)
        if exclude_d_amino_acids(directory + '/' + match):
            # Read
            has_protein_scaffold = False
            for i in range(1, int(match_info_list[-1]) + 1):
                total_objs = len(cmd.get_object_list())
                match = match_prefix + '_' + str(i) + '.pdb'
                cmd.load(directory + '/' + match)
                scaffold_new_name = match[:-4] + '_0001'
                cmd.set_name(cmd.get_object_list()[total_objs], scaffold_new_name)
                if has_protein_scaffold:
                    cmd.delete(scaffold_new_name)
                else:
                    has_protein_scaffold = True
            mutations = get_match_point_mutations(remarks)
            for mutated_res in mutations.values():
                cmd.remove('chain A & resi 1 & resn ' + mutated_res)
            objs = cmd.get_object_list()
            cmd.alter(objs[-1], 'chain="A"')
            cmd.alter(objs[-1], 'resi="999"')
            # Output pdb
            output_name = match_info_list[2]
            os.mkdir(output_name)
            os.chdir(output_name)
            if homomeric:
                cmd.create('main_chain_linker', 'chain ' + remarks[0][49])
                cmd.create('duplicate_chain', 'chain ' + remarks[1][49])
                cmd.align('duplicate_chain', 'main_chain_linker')
                cmd.save(output_name + '_main_chain_linker.pdb', 'main_chain_linker')
                cmd.save(output_name + '_duplicate_chain.pdb', 'duplicate_chain')
                cmd.delete('main_chain_linker')
                cmd.delete('duplicate_chain')
                mutated_res_lines = list()
                with open(output_name + '_duplicate_chain.pdb', 'r') as p_pdb:
                    for line in p_pdb:
                        if line.startswith('ATOM') and line[17:22] == remarks[1][51:54] + \
                            ' ' + remarks[1][49] and int(line[23:26]) == int(remarks[1][56:59]):
                            mutated_res_lines.append(line[:21] + remarks[0][49] + line[22:])
                os.remove(output_name + '_duplicate_chain.pdb')
                with open(output_name + '_main_chain_linker.pdb', 'r+') as p_pdb:
                    lines = p_pdb.readlines()
                    p_pdb.seek(0)
                    written_mutated_res_lines = False
                    for line in lines:
                        if line.startswith('ATOM') and line[21] == remarks[0][49] and \
                            int(line[23:26]) == int(remarks[1][56:59]):
                            if not written_mutated_res_lines:
                                p_pdb.writelines(mutated_res_lines)
                                written_mutated_res_lines = True
                        elif not line.startswith('CONECT'):
                            p_pdb.write(line)
                cmd.load(output_name + '_main_chain_linker.pdb')
                try:
                    cmd.save(output_name + '.pdb', output_name + '_main_chain_linker')
                except:
                    print('ERROR:' + match_info_list[2] + ' is not processed correctly!')
                    cmd.delete('*')
                    os.chdir('..')
                    continue
                cmd.delete(output_name + '_main_chain_linker')
                os.remove(output_name + '_main_chain_linker.pdb')
                with open(output_name + '.pdb', 'r+') as p_pdb:
                    lines = p_pdb.readlines()
                    p_pdb.seek(0)
                    for remark in remarks:
                        p_pdb.write(remark[:26] + 'A' + remark[27:33] + '999' + remark[36:])
                    for remark in remarks:
                        p_pdb.write(remark[:26] + 'B' + remark[27:33] + '999' + remark[36:49])
                        if remark[49] == 'A':
                            p_pdb.write('B' + remark[50:61])
                        elif remark[49] == 'B':
                            p_pdb.write('A' + remark[50:61])
                        if remark[61] == '1':
                            p_pdb.write('3' + remark[62:])
                        elif remark[61] == '2':
                            p_pdb.write('4' + remark[62:])
                    p_pdb.writelines(lines)
            else:
                cmd.save(output_name + '.pdb', objs[0] + '|' + objs[-1])
                with open(output_name + '.pdb', 'r+') as p_pdb:
                    lines = p_pdb.readlines()
                    p_pdb.seek(0)
                    for remark in remarks:
                        p_pdb.write(remark[:26] + 'A' + remark[27:33] + '999' + remark[36:])
                    p_pdb.writelines(lines)
            # Output rotlib
            with open(remarks[0][28:31] + '.rotlib.pdb', 'w+') as p_pdb:
                for i in range(1, len(objs) - 1):
                    rotamer_name = objs[i] + '.pdb'
                    cmd.save(rotamer_name, objs[i])
                    with open(rotamer_name, 'r') as rotamer:
                        p_pdb.write('REMARK ' + str(i) + '\n')
                        for line in rotamer:
                            if line.startswith('HETATM'):
                                p_pdb.write(line)
                        p_pdb.write('ENDMDL\n')
                    os.remove(rotamer_name)
            cmd.delete('*')
            os.chdir('..')
            #break


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str)
    parser.add_argument('-homo', '--homomeric', action='store_true')
    args = parser.parse_args()

    match_dict = collect_output_match_info(args.directory)
    make_relax_input_files(args.directory, match_dict, args.homomeric)
