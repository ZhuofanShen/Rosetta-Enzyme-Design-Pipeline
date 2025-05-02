import argparse
import os
from pymol import *
from random import randint
import shutil
import time


def collect_output_match_info(directory):
    match_dict = dict()
    for match in filter(lambda x: x.startswith("UM_") and x.endswith(".pdb"), os.listdir(directory)):
        match_info = match[:-4].split('_')
        position = match_info[2]
        record = match_dict.get(position)
        if not record:
            match_dict[position] = match_info[:-1] + [{match_info[-1]}]
        elif match_info[-1] not in record[-1]:
            record[-1].add(match_info[-1])
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

def convert_d_amino_acids_to_l(match, enzdes_cst_remarks):
    ncaa_2_caa = {"TYZ": "TYR", "CYX": "CYS", "HIX": "HIS", "GLX": "GLU", "ASX": "ASP", "LYX": "LYS", \
                  "SEX": "SER", "THX": "THR", "TYX": "TYR", "MEX": "MET"}
    chain_position_ncaa_dict = dict()
    cmd.load(match)
    prefix = '/' + cmd.get_object_list()[0] + '//'
    for remark in enzdes_cst_remarks:
        position_info = list(filter(lambda x : x != '', remark.split(' ')))
        atom_name_prefix = prefix + position_info[-6] + '/' + position_info[-4] + '/'
        dihedral = cmd.get_dihedral(atom_name_prefix + 'N', atom_name_prefix + 'CA', \
            atom_name_prefix + 'CB', atom_name_prefix + 'C')
        if dihedral < 0:
            chain_position_ncaa_dict[position_info[-6] + "_" + position_info[-4]] = position_info[-5]
            cmd.wizard("mutagenesis")
            cmd.do("refresh_wizard")
            cmd.get_wizard().set_mode("GLY")
            cmd.get_wizard().do_select(atom_name_prefix)
            cmd.get_wizard().apply()
            cmd.get_wizard().set_mode(ncaa_2_caa[position_info[-5]])
            cmd.get_wizard().do_select(atom_name_prefix)
            cmd.get_wizard().apply()
    for chain_position, ncaa in chain_position_ncaa_dict.items():
        chain, position = chain_position.split("_")
        cmd.alter(prefix + chain + '/' + position + '/', "resn='" + ncaa + "'")
        cmd.remove(prefix + "HH")
        cmd.remove(prefix + "HG")
        cmd.remove(prefix + "1HD")
        cmd.remove(prefix + "2HE")
        cmd.remove(prefix + "2HD")
        cmd.remove(prefix + "3HZ")
    if len(chain_position_ncaa_dict) > 0:
        cmd.save(match)
        time.sleep(10)
    cmd.delete('*')

def get_match_point_mutations(enzdes_cst_remarks):
    point_mutation_dict = dict() # {"A58": "CYX", "B61": "TYZ"}
    for enzdes_cst_remark in enzdes_cst_remarks:
        point_mutation_dict[enzdes_cst_remark[49] + str(int(enzdes_cst_remark[56:59]))] = enzdes_cst_remark[51:54]
    return point_mutation_dict

def read_rotlib_from_cloud_pdb(cloud_pdb_prefix, cloud_pdb_list, cutoff=1000):
    rotamer_library = list()
    rotlib_lines = list()
    for i_cloud_pdb in cloud_pdb_list:
        cloud_pdb = cloud_pdb_prefix + '_' + i_cloud_pdb + '.pdb'
        with open(cloud_pdb, "r") as pf:
            ligand_name3 = None
            read_rotamer = False
            rotamer_index = 0
            for line in pf:
                if not ligand_name3 and line.startswith("REMARK 666 MATCH TEMPLATE"):
                    ligand_name3 = line[28:31]
                    # chain_id = line[49]
                elif ligand_name3 and not read_rotamer and line.startswith("HETNAM     " + ligand_name3):
                    read_rotamer = True
                    rotamer_lines = list()
                    rotamer_index += 1
                    # rotlib_lines.append('REMARK ' + str(rotamer_index) + '\n')
                elif ligand_name3 and read_rotamer and line.startswith("HETATM"):
                    # rotamer_lines.append(line[:21] + chain_id + " 999" + line[26:])
                    rotamer_lines.append(line)
                elif ligand_name3 and read_rotamer and line.startswith("TER"):
                    read_rotamer = False
                    # rotamer_lines.append("TER\n")
                    rotamer_library.append(rotamer_lines)
    cutoff = min(len(rotamer_library), cutoff)
    for rotamer_index in range(1, cutoff + 1):
        if rotamer_index < 10:
            rotamer_index = '   ' + str(rotamer_index)
        elif rotamer_index < 100:
            rotamer_index = '  ' + str(rotamer_index)
        elif rotamer_index < 1000:
            rotamer_index = ' ' + str(rotamer_index)
        else:
            rotamer_index = str(rotamer_index)
        rotlib_lines.append('MODEL  ' + rotamer_index + '\n')
        rotamer_lines = rotamer_library.pop(randint(0, len(rotamer_library) - 1))
        for line in rotamer_lines:
            rotlib_lines.append(line[:22] + rotamer_index + line[26:])
        rotlib_lines.append("ENDMDL \n")
    return rotlib_lines

def make_relax_input_files(directory, match_dict, duplicate_match=False, symmetry=False, debug_mode=False):
    os.chdir(directory)
    for position, match_info_list in match_dict.items(): # match_dict: {"X384Z115": ['UM', '2', 'X384Z115', 'CPG2', 'O2beY-CYS', ['1', '2', '3']]}
        match_prefix = '_'.join(match_info_list[:-1]) # rewriting FastRelax input files is allowed if Matcher output files exist
        rotlib_lines = read_rotlib_from_cloud_pdb(match_prefix, match_info_list[-1])
        # with open(rotamer_name, 'w+') as p_pdb:
        #     for i in range(1, len(objs) - 1):
        #         rotamer_name = objs[i] + '.pdb'
        #         cmd.save(rotamer_name, objs[i])
        #         with open(rotamer_name, 'r') as rotamer:
        #             p_pdb.write('REMARK ' + str(i) + '\n')
        #             for line in rotamer:
        #                 if line.startswith('HETATM'):
        #                     p_pdb.write(line)
        #             p_pdb.write('TER\n')
        #         os.remove(rotamer_name)
        if len(rotlib_lines) == 0:
            continue
        match = match_prefix + "_1.pdb"
        remarks = read_enz_des_remark_lines(match)
        convert_d_amino_acids_to_l(match, remarks)
        # Read
        # has_protein_scaffold = False
        # for i in match_info_list[-1]:
        #     total_objs = len(cmd.get_object_list())
        #     match = match_prefix + '_' + i + '.pdb'
        #     cmd.load(os.path.join(directory, match))
        #     scaffold_new_name = match[:-4] + '_0001'
        #     cmd.set_name(cmd.get_object_list()[total_objs], scaffold_new_name)
        #     if has_protein_scaffold:
        #         cmd.delete(scaffold_new_name)
        #     else:
        #         has_protein_scaffold = True
        # Output pdb
        output_name = position
        if duplicate_match or symmetry:
            cmd.load(os.path.join(directory, match))
            scaffold_new_name = match[:-4] + '_0001'
            cmd.set_name(cmd.get_object_list()[0], scaffold_new_name)
            mutations = get_match_point_mutations(remarks)
            for mutated_res in mutations.values():
                cmd.remove('chain A & resi 1 & resn ' + mutated_res)
        if duplicate_match:
            cmd.create('chain_A_linker', 'chain ' + remarks[0][49])
            cmd.create('chain_B', 'chain ' + remarks[-1][49])
            # for chain A
            cmd.create('chain_B_2', 'chain_B')
            # cmd.align('chain_B_2', 'chain_A_linker')
            for mutated_res_index in mutations.keys():
                if mutated_res_index.startswith(remarks[0][49]):
                    chain_A_mutated_res_index = mutated_res_index
                elif mutated_res_index.startswith(remarks[-1][49]):
                    chain_B_mutated_res_index = mutated_res_index
            cmd.pair_fit('/chain_B_2//' + chain_B_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/N', \
                    '/chain_A_linker//' + chain_A_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/N', \
                    '/chain_B_2//' + chain_B_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/CA', \
                    '/chain_A_linker//' + chain_A_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/CA', \
                    '/chain_B_2//' + chain_B_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/C', \
                    '/chain_A_linker//' + chain_A_mutated_res_index[0] + '/' + chain_B_mutated_res_index[1:] + '/C')
            cmd.save(output_name + '_chain_A_linker.pdb', 'chain_A_linker')
            cmd.save(output_name + '_chain_B_2.pdb', 'chain_B_2')
            # for chain B
            if not symmetry:
                cmd.create('chain_A_linker_2', 'chain_A_linker')
                # cmd.align('chain_A_linker_2', 'chain_B')
                cmd.pair_fit('/chain_A_linker_2//' + chain_A_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/N', \
                    '/chain_B//' + chain_B_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/N', \
                    '/chain_A_linker_2//' + chain_A_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/CA', \
                    '/chain_B//' + chain_B_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/CA', \
                    '/chain_A_linker_2//' + chain_A_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/C', \
                    '/chain_B//' + chain_B_mutated_res_index[0] + '/' + chain_A_mutated_res_index[1:] + '/C')
                cmd.save(output_name + '_chain_B.pdb', 'chain_B')
                cmd.save(output_name + '_chain_A_linker_2.pdb', 'chain_A_linker_2')
            cmd.delete('chain_A_linker*')
            cmd.delete('chain_B*')
            #time.sleep(1)
            # for chain A
            chain_B_mutated_res_lines = list()
            with open(output_name + '_chain_B_2.pdb', 'r') as p_pdb:
                for line in p_pdb:
                    if line.startswith('ATOM') and line[17:22] == remarks[-1][51:54] + \
                        ' ' + remarks[-1][49] and int(line[23:26]) == int(remarks[-1][56:59]):
                        chain_B_mutated_res_lines.append(line[:21] + remarks[0][49] + line[22:])
            os.remove(output_name + '_chain_B_2.pdb')
            with open(output_name + '_chain_A_linker.pdb', 'r') as p_pdb:
                lines = p_pdb.readlines()
            os.remove(output_name + '_chain_A_linker.pdb')
            with open(output_name + '_chain_A_dup_match.pdb', 'w') as p_pdb:
                written_mutated_res_lines = False
                for line in lines:
                    if line.startswith('ATOM') and line[21] == remarks[0][49] and \
                        int(line[23:26]) == int(remarks[-1][56:59]):
                        if not written_mutated_res_lines:
                            p_pdb.writelines(chain_B_mutated_res_lines)
                            written_mutated_res_lines = True
                    elif not line.startswith('CONECT'):
                        p_pdb.write(line)
            cmd.load(output_name + '_chain_A_dup_match.pdb')
            # for chain B
            if not symmetry:
                chain_A_mutated_res_lines = list()
                chain_A_linker_lines = list()
                with open(output_name + '_chain_A_linker_2.pdb', 'r') as p_pdb:
                    for line in p_pdb:
                        # mutated residue lines
                        if line.startswith('ATOM') and line[17:20] == remarks[0][51:54] and \
                            line[21] == remarks[0][49] and int(line[23:26]) == int(remarks[0][56:59]):
                            chain_A_mutated_res_lines.append(line[:21] + remarks[-1][49] + line[22:])
                        #linker lines
                        elif line.startswith('HETATM') and line[17:20] == remarks[0][28:31] and \
                            line[21] == remarks[0][49]:
                            chain_A_linker_lines.append(line[:21] + remarks[-1][49] + line[22:])
                os.remove(output_name + '_chain_A_linker_2.pdb')
                with open(output_name + '_chain_B.pdb', 'r') as p_pdb:
                    lines = p_pdb.readlines()
                os.remove(output_name + '_chain_B.pdb')
                with open(output_name + '_chain_B_dup_match.pdb', 'w') as p_pdb:
                    written_mutated_res_lines = False
                    for line in lines:
                        if line.startswith('ATOM') and line[21] == remarks[-1][49] and \
                            int(line[23:26]) == int(remarks[0][56:59]):
                            if not written_mutated_res_lines:
                                p_pdb.writelines(chain_A_mutated_res_lines)
                                written_mutated_res_lines = True
                        elif not line.startswith('CONECT') and not line.startswith('END'):
                            p_pdb.write(line)
                    p_pdb.writelines(chain_A_linker_lines)
                    p_pdb.write('END\n')
                cmd.load(output_name + '_chain_B_dup_match.pdb')
            if symmetry:
                cmd.save(output_name + '.pdb', output_name + '_chain_A_dup_match')
                cmd.delete(output_name + '_chain_A_linker')
                os.remove(output_name + '_chain_A_dup_match.pdb')
            else:
                cmd.save(output_name + '.pdb', output_name + '_chain_A_dup_match | ' + output_name + '_chain_B_dup_match')
                cmd.delete(output_name + '_chain_A_linker')
                cmd.delete(output_name + '_chain_B')
                os.remove(output_name + '_chain_A_dup_match.pdb')
                os.remove(output_name + '_chain_B_dup_match.pdb')
            # for chain A
            header_lines = list()
            for remark in remarks:
                header_lines.append(remark[:26] + remarks[0][49] + remark[27:33] + '999' + remark[36:])
            # for chain B
            for remark in remarks:
                header_lines.append(remark[:26] + remarks[-1][49] + remark[27:33] + '999' + remark[36:49])
                if remark[49] == remarks[0][49]:
                    header_lines.append(remarks[-1][49] + remark[50:61] + str(int(remark[61]) + len(remarks)) + remark[62:])
                elif remark[49] == remarks[-1][49]:
                    header_lines.append(remarks[0][49] + remark[50:61] + str(int(remark[61]) + len(remarks)) + remark[62:])
            #time.sleep(1)
            with open(output_name + '.pdb', 'r+') as p_pdb:
                lines = list(filter(lambda x: x.startswith("ATOM  ") or x.startswith("HETATM"), p_pdb))
                p_pdb.seek(0)
                p_pdb.writelines(header_lines)
                p_pdb.writelines(lines)
                p_pdb.writelines(rotlib_lines)
        else:
            header_lines = list()
            for remark in remarks:
                header_lines.append(remark[:26] + remarks[0][49] + remark[27:33] + '999' + remark[36:])
            if symmetry:
                cmd.create('chain_A_linker', 'chain ' + remarks[0][49])
                cmd.save(output_name + '.pdb', 'chain_A_linker')
                # time.sleep(1)
                with open(output_name + '.pdb', 'r+') as p_pdb:
                    lines = list(filter(lambda x: x.startswith("ATOM  ") or x.startswith("HETATM"), p_pdb))
                    p_pdb.seek(0)
                    p_pdb.writelines(header_lines)
                    p_pdb.writelines(lines)
                    p_pdb.writelines(rotlib_lines)
            else:
                # cmd.save(output_name + '.pdb', objs[0] + '|' + objs[-1])
                lines = list()
                read_ligand = False
                with open(match, "r") as p_match:
                    for line in p_match:
                        if not read_ligand and (line.startswith("ATOM  ") or line.startswith("HETATM")):
                            if line[17:20] == remarks[0][28:31]:
                                # read_ligand = True
                                # lines.append(line[:21] + remarks[0][49] + " 999" + line[26:])
                                break # DO NOT append the first substrate conformer to pdb lines
                            else:
                                lines.append(line)
                        # elif read_ligand:
                        #     if line[17:20] == remarks[0][28:31]:
                        #         lines.append(line[:21] + remarks[0][49] + " 999" + line[26:])
                        #     else:
                        #         break

                        # if not read_ligand and (line.startswith("ATOM  ") or line.startswith("HETATM")):
                        #     lines.append(line)
                        # elif not read_ligand and line.startswith("HETNAM     " + remarks[0][28:31]):
                        #     read_ligand = True
                        # elif read_ligand and line.startswith("HETATM"):
                        #     lines.append(line[:21] + remarks[0][49] + " 999" + line[26:])
                        # elif read_ligand and line.startswith("TER"):
                        #     break
                with open(output_name + ".pdb", 'w') as p_pdb:
                    p_pdb.writelines(header_lines)
                    p_pdb.writelines(lines)
                    p_pdb.writelines(rotlib_lines)
        if duplicate_match or symmetry:
            cmd.delete('*')
        # Copy the ligand params file to the current directory
        # shutil.copyfile(os.path.join('..', ligand_path, remarks[0][28:31] + ".params"), remarks[0][28:31]+ ".params")
        # with open(os.path.join('..', ligand_path, remarks[0][28:31] + ".params"), "r") as p_params:
        #     params_lines = p_params.readlines()[:-1]
        # with open(position + ".params", "w") as p_params:
        #     p_params.writelines(params_lines)
        #     p_params.write("PDB_ROTAMERS " + position + ".rotlib" + "\n")
        if debug_mode:
            break
        else:
            for i_match in match_info_list[-1]:
                os.remove(match_prefix + '_' + i_match + '.pdb')
    if os.path.isfile(directory + ".log") and not debug_mode:
        os.remove(directory + ".log")
    if os.path.isfile(directory + ".err") and not debug_mode:
        os.remove(directory + ".err")
    os.chdir('..')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str)
    parser.add_argument('-dup', '--duplicate_match', action='store_true')
    parser.add_argument('-symm', '--symmetry', action='store_true')
    parser.add_argument('-debug', '--debug_mode', action='store_true')
    args = parser.parse_args()

    directory = args.directory.strip("/")

    match_dict = collect_output_match_info(directory)
    make_relax_input_files(directory, match_dict, args.duplicate_match, args.symmetry, args.debug_mode)
    if len(os.listdir(directory)) == 0:
        shutil.rmtree(directory)
