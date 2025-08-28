import argparse
import json
# import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
# from PIL import Image
from pymol import *
from pyrosetta import *
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, OrResidueSelector, NotResidueSelector, \
    ResidueIndexSelector, ResidueNameSelector, \
    InterGroupInterfaceByVectorSelector
# import seaborn as sns
import shutil


def read_scores_from_fasc(fasc_path, n=1, is_reversed=False):
    decoy_scores_str = "["
    with open(fasc_path) as p_fasc:
        for line in p_fasc:
            decoy_scores_str += line[:-1]
            decoy_scores_str += ","
    decoy_scores_str = decoy_scores_str[:-1] + "]"
    scores = json.loads(decoy_scores_str)
    scores = list(sorted(scores, key=lambda x: x["total_score"], reverse=is_reversed)[:n])
    return scores

def index_boolean_filter(index, boolean):
    if type(boolean) == np.ndarray:
        for bool in boolean:
            if bool:
                return str(index)
    elif boolean:
        return str(index)
    return False

def boolean_vector_to_indices_set(boolean_vector, n_monomers:int=1):
    boolean_vector = np.array(boolean_vector)
    sequence_length = len(boolean_vector) // n_monomers
    boolean_matrix = boolean_vector[:n_monomers*sequence_length]\
            .reshape(n_monomers, sequence_length).transpose()
    indices = set(filter(lambda x: x, map(index_boolean_filter, \
            range(1, len(boolean_matrix) + 1), boolean_matrix)))
    return indices

xlink_name3_dict = {"O2beY": "OBY"}
uaa_name3_dict = {"O2beY": "TYZ"}
nucleophile_name3_dict = {"CYS": "CYX", "HIS": "HIX", "ASP": "ASX", "GLU": "GLX", "LYS": "LYX"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str)
    parser.add_argument('-lig', "--ligand", type=str, default="O2beY-CYS")
    parser.add_argument('--process_variants', type=str, nargs="*", default=list()) # e.g., 1ysb-A_S6_E126
    parser.add_argument('-s', '--step', choices=['design', 'relax'], default='relax')
    parser.add_argument('-o', '--option', choices=['csv', 'mpnn', 'mpnnall', 'barplot'], default='csv')
    parser.add_argument("--ddg_bind", action="store_true")
    parser.add_argument('-od', '--output_directory', type=str, default='mpnn')
    # parser.add_argument('-cat', '--catalytic_residues', type=str, nargs="*", default=["62", "91", "94", "155", "64", "51", "53", "93", "86", "88", "152", "114", "156", "33", "31", "141", "89", "90", "116", "65", "95", "100"])
    parser.add_argument('-cat', '--catalytic_residues', type=str, nargs="*", default=["1", "18", "22", "25", "28", "29", "30", "32", "33", "34", "37", "43", "49", "51", "53", "55", "62", "64", "71", "72", "73", "79", "86", "87", "88", "90", "91", "94", "96", "103", "106", "109", "110", "114", "116", "122", "126", "127", "137", "141", "144", "149", "152", "154", "155", "156", "157", "93", "31", "89", "65", "95", "100"]) # , "93", "31", "89", "65", "95", "100"
    # parser.add_argument('-cat', '--catalytic_residues', type=str, nargs="*", default=["1", "18", "21", "22", "25", "28", "29", "30", "31", "32", "33", "34", "37", "38", "43", "47", "48", "49", "50", "51", "52", "53", "54", "55", "57", "62", "63", "64", "66", "67", "68", "70", "71", "72", "73", "76", "79", "83", "85", "86", "87", "88", "89", "90", "91", "93", "94", "95", "96", "102", "103", "105", "106", "107", "109", "110", "113", "114", "116", "119", "122", "126", "127", "129", "132", "134", "137", "140", "141", "144", "145", "149", "151", "152", "153", "154", "155", "156", "157", "65", "100"]) # , "65", "100"
    parser.add_argument('-rn', '--read_n_decoys', type=int, default=1)
    parser.add_argument('-concat', '--concatenate', choices=['horizontal', 'vertical'], default='vertical')
    args = parser.parse_args()

    uaa_type, nucleophile_type = args.ligand.split("-")

    uaa_name3 = uaa_name3_dict[uaa_type]
    xlink_name3 = xlink_name3_dict[uaa_type]
    nucleophile_name3 = nucleophile_name3_dict[nucleophile_type]

    params_files_list = ["ligands/" + args.ligand + "/" + uaa_name3 + ".params", 
         "ligands/" + args.ligand + "/" + xlink_name3 + ".params", 
         "ligands/" + args.ligand + "/" + nucleophile_name3 + ".params"]
    
    if args.option == "csv" or args.option == "mpnn" or args.option == "mpnnall":
        init("-ex1 -ex2 -ignore_zero_occupancy false -use_input_sc -no_optH true " +
            "-extra_res_fa " + " ".join(params_files_list) + " " +
            " -enzdes:cstfile ligands/" + args.ligand + "/" + args.ligand + "_Relax.cst " +
            "-run:preserve_header")

    for pdb in os.listdir(args.directory):
        pdb_path = os.path.join(args.directory, pdb)
        for pdb_chain_ligand in filter(lambda x: x.startswith(pdb) and x.endswith(args.ligand) and \
                os.path.isdir(os.path.join(pdb_path, x)), os.listdir(pdb_path)):
            pdb_chain, uaa_nucleophile = pdb_chain_ligand.split("_")
            with open(os.path.join("pdb" + pdb_path[3:], pdb_chain + ".pdb"), "r") as pf:
                xtal_pdb_lines = pf.readlines()
            chains = pdb_chain.split("-")[1:]
            uaa_chain = chains[0]
            if len(chains) == 1:
                nucleophile_chain = chains[0]
                no_uaa_chain = "B" if chains[0] == "A" else "A" # dimer, staple within one chain
                no_nucleophile_chain = no_uaa_chain # dimer, staple within one chain
            else:
                nucleophile_chain = chains[1]
                no_uaa_chain = chains[1] # dimer, staple across two chains
                no_nucleophile_chain = chains[0] # dimer, staple across two chains
            pdb_ligand_path = os.path.join(pdb_path, pdb_chain_ligand)
            if args.option == "csv" or args.option == "mpnn" or args.option == "mpnnall":
                df_list = list()
            for variant in filter(lambda x: os.path.isdir(os.path.join(pdb_ligand_path, x)) and x.startswith(pdb_chain), os.listdir(pdb_ligand_path)):
                if_process_current_variant_anyway = False
                if len(args.process_variants) > 0:
                    if variant in args.process_variants:
                        if_process_current_variant_anyway = True
                    else:
                        continue
                variant_path = os.path.join(pdb_ligand_path, variant)
                mutations = variant.split("_")[1:]
                variant_pdb_numbering = None
                for variant_pdb in os.listdir(variant_path):
                    if variant_pdb.startswith("Z") and variant_pdb.endswith(".pdb") \
                            and not variant_pdb.endswith(".MPNN.pdb") and not variant_pdb.endswith("-unbound.pdb"):
                        pdb_positions = variant_pdb.split(".")[0][1:].split("X")
                        uaa_pdb_index = pdb_positions[0]
                        nucleophile_pdb_index = pdb_positions[1]
                        variant_pdb_numbering = pdb + "_" + mutations[0][0] + uaa_pdb_index + "UAA_" + mutations[1][0] + nucleophile_pdb_index
                        break
                if variant_pdb_numbering is None:
                    shutil.rmtree(variant_path)
                    continue
                fasc_list = list()
                if args.step == "design":
                    design_path = os.path.join(variant_path, "fast_design")
                    fasc = os.path.join(design_path, variant + ".sc")
                    if os.path.isfile(fasc):
                        fasc_list.append(fasc)
                elif args.step == "relax":
                    fasc = os.path.join(variant_path, variant + ".sc")
                    if os.path.isfile(fasc):
                        fasc_list.append(fasc)
                    else:
                        for relax_path in filter(lambda x: x.startswith("fast_relax_"), os.listdir(variant_path)):
                            relax_path = os.path.join(variant_path, relax_path)
                            fasc = os.path.join(relax_path, variant + ".sc")
                            if os.path.isfile(fasc):
                                fasc_list.append(fasc)
                
                pass_threshold = False
                lowest_ddg = 1000000
                lowest_fasc = None
                for fasc in fasc_list:
                    path = fasc[:fasc.rfind("/")]
                    wt_fasc = os.path.join(path, pdb_chain + "_apo.sc")
                    if os.path.isfile(wt_fasc):
                        wt_score_dict = read_scores_from_fasc(wt_fasc)[0]
                        wt_total_score = wt_score_dict["total_score"] - wt_score_dict["coordinate_constraint"]
                        wt_unbound_total_score = 0
                        wt_unbound_fasc = os.path.join(path, pdb_chain + "-unbound_apo.sc")
                        if args.ddg_bind and os.path.isfile(wt_unbound_fasc):
                            wt_unbound_score_dict = read_scores_from_fasc(wt_unbound_fasc)[0]
                            wt_unbound_total_score = wt_unbound_score_dict["total_score"] - wt_unbound_score_dict["coordinate_constraint"]
                    dG_bind_wt = wt_total_score - wt_unbound_total_score
                    if os.path.isfile(fasc):
                        score_dict = read_scores_from_fasc(fasc, n=args.read_n_decoys)[0]
                        total_score = score_dict["total_score"] - score_dict["coordinate_constraint"]
                        unbound_total_score = 0
                        unbound_fasc = fasc[:-3] + "-unbound.sc"
                        if args.ddg_bind and os.path.isfile(unbound_fasc):
                            unbound_score_dict = read_scores_from_fasc(unbound_fasc)[0]
                            unbound_total_score = unbound_score_dict["total_score"] - unbound_score_dict["coordinate_constraint"]
                    dG_bind_mut = total_score - unbound_total_score
                    ddg = dG_bind_mut - dG_bind_wt
                    # ddg = (score_dict["total_score"] - score_dict["coordinate_constraint"]) - wt_total_score
                    residue_energy = score_dict["substrates"]
                    residue_cst_score = score_dict["atom_pair_constraint"] + \
                            score_dict["angle_constraint"] + score_dict["dihedral_constraint"]
                    decoy = os.path.join(path, score_dict["decoy"])
                    
                    if if_process_current_variant_anyway or args.step == "relax" or \
                            (args.step == "design" and residue_energy < 3 and residue_cst_score < 3):
                        pass_threshold = True
                        if args.option == "csv" or args.option == "mpnn" or args.option == "mpnnall":
                            var_dict = {"variant": [variant_pdb_numbering], "ddG": [round(ddg, 2)], \
                                    "xlink": [round(residue_energy, 2)], "xlink_cst": [round(residue_cst_score, 2)]}
                            if args.step == "relax":
                                pseudo_wt_pose = pose_from_pdb(os.path.join(variant_path, variant_pdb))
                                pseudo_wt_sequence = pseudo_wt_pose.sequence()
                                var_pose = pose_from_pdb(decoy)
                                var_sequence = var_pose.sequence()
                                uaa_pose_index = var_pose.pdb_info().pdb2pose(uaa_chain, int(uaa_pdb_index))
                                nucleophile_pose_index = var_pose.pdb_info().pdb2pose(nucleophile_chain, int(nucleophile_pdb_index))
                                mutations_chain_pdb_index = list()
                                mutations_pdb_index = list()
                                for i in range(len(var_sequence)):
                                    if var_sequence[i] != pseudo_wt_sequence[i]:
                                        info = pseudo_wt_pose.pdb_info().pose2pdb(i + 1).split(" ")
                                        if i+1 != nucleophile_pose_index and i+1 != uaa_pose_index and \
                                                var_sequence[i] != "X" and var_sequence[i] != "Z":
                                            mutations_chain_pdb_index.append(info[1] + pseudo_wt_sequence[i] + info[0] + var_sequence[i])
                                            mutations_pdb_index.append(pseudo_wt_sequence[i] + info[0] + var_sequence[i])
                                mutations_chain_pdb_index_str = ",".join(mutations_chain_pdb_index)
                                var_dict["mutations"] = [mutations_chain_pdb_index_str]
                                mutations_pdb_index_str = str()
                                if len(mutations_pdb_index) > 0:
                                    mutations_pdb_index_str = "_" + "_".join(mutations_pdb_index)
                                shutil.copy(decoy, os.path.join(pdb_ligand_path, variant_pdb_numbering + "X" + mutations_pdb_index_str + ".pdb"))
                            df_list.append(pd.DataFrame(var_dict))

                        if ddg < lowest_ddg:
                            lowest_fasc = fasc
                            lowest_ddg = ddg
                            lowest_residue_energy = residue_energy
                            lowest_residue_cst_score = residue_cst_score
                            mpnn_input = os.path.join(variant_path, variant_pdb.split(".")[0] + ".MPNN.pdb")
                if (args.option == "mpnn" or args.option == "mpnnall") and lowest_fasc is not None and not os.path.isfile(mpnn_input):
                    uaa_lines = list()
                    nucleophile_lines = list()
                    xlink_lines = list()
                    with open(decoy, "r") as pf:
                        for pdb_line in filter(lambda x: x.startswith("ATOM") or x.startswith("HETATM"), pf):
                            if pdb_line[17:20] == uaa_name3:
                                uaa_lines.append(pdb_line)
                            elif pdb_line[17:20] == nucleophile_name3:
                                nucleophile_lines.append(pdb_line.replace("CYX", "HIX"))
                            elif pdb_line[17:20] == xlink_name3:
                                xlink_lines.append(pdb_line)
                    # <start: two chains, duplicate the staple
                    cmd.load(decoy)
                    cmd.create('UAA_chain', 'chain ' + uaa_chain)
                    cmd.create('no_UAA_chain', 'chain ' + no_uaa_chain)
                    cmd.align('UAA_chain', 'no_UAA_chain')
                    cmd.alter('UAA_chain', "chain='" + no_uaa_chain + "'")
                    cmd.save(mpnn_input[:-3] + "new_UAA.pdb", 'UAA_chain')
                    cmd.create('nucleophile_chain', 'chain ' + nucleophile_chain)
                    cmd.create('no_nucleophile_chain', 'chain ' + no_nucleophile_chain)
                    cmd.align('nucleophile_chain', 'no_nucleophile_chain')
                    cmd.alter('nucleophile_chain', "chain='" + no_nucleophile_chain + "'")
                    cmd.save(mpnn_input[:-3] + "new_nucleophile.pdb", 'nucleophile_chain')
                    cmd.delete('*')
                    new_uaa_lines = list()
                    new_nucleophile_lines = list()
                    new_xlink_lines = list()
                    with open(mpnn_input[:-3] + "new_UAA.pdb", "r") as pf:
                        for pdb_line in filter(lambda x: x.startswith("ATOM") or x.startswith("HETATM"), pf):
                            if pdb_line[17:20] == uaa_name3:
                                new_uaa_lines.append(pdb_line)
                            elif pdb_line[17:20] == xlink_name3:
                                new_xlink_lines.append(pdb_line)
                    with open(mpnn_input[:-3] + "new_nucleophile.pdb", "r") as pf:
                        for pdb_line in filter(lambda x: x.startswith("ATOM") or x.startswith("HETATM"), pf):
                            if pdb_line[17:20] == nucleophile_name3:
                                new_nucleophile_lines.append(pdb_line.replace("CYX", "HIX"))
                    os.remove(mpnn_input[:-3] + "new_UAA.pdb")
                    os.remove(mpnn_input[:-3] + "new_nucleophile.pdb")
                    # /end>
                    with open(mpnn_input, "w") as pf:
                        write_uaa_lines = False
                        write_nucleophile_lines = False
                        write_new_uaa_lines = False
                        write_new_nucleophile_lines = False
                        for pdb_line in filter(lambda x: x.startswith("ATOM  ") or x.startswith("HETATM"), xtal_pdb_lines):
                            if int(pdb_line[22:26]) == int(uaa_pdb_index) and pdb_line[21] == uaa_chain:
                                if not write_uaa_lines:
                                    pf.writelines(uaa_lines)
                                    write_uaa_lines = True
                            elif int(pdb_line[22:26]) == int(nucleophile_pdb_index) and pdb_line[21] == nucleophile_chain:
                                if not write_nucleophile_lines:
                                    pf.writelines(nucleophile_lines)
                                    write_nucleophile_lines = True
                            # <start: two chains, duplicate the staple
                            elif int(pdb_line[22:26]) == int(uaa_pdb_index) and pdb_line[21] == no_uaa_chain:
                                if not write_new_uaa_lines:
                                    pf.writelines(new_uaa_lines)
                                    write_new_uaa_lines = True
                            elif int(pdb_line[22:26]) == int(nucleophile_pdb_index) and pdb_line[21] == no_nucleophile_chain:
                                if not write_new_nucleophile_lines:
                                    pf.writelines(new_nucleophile_lines)
                                    write_new_nucleophile_lines = True
                            # /end>
                            else:
                                pf.write(pdb_line)
                        pf.writelines(xlink_lines)
                        pf.writelines(new_xlink_lines)
                    # with open(decoy, "r") as pf:
                    #     pdb_lines = pf.readlines()
                    # with open(mpnn_input, "w") as pf:
                    #     for pdb_line in pdb_lines:
                    #         pf.write(pdb_line.replace("CYX", "HIX"))

                if (args.option == "mpnn" or args.option == "mpnnall") and (args.step == "relax" or (args.step == "design" and pass_threshold)):
                    redesigned_res_pdb_indices = set()
                    # intf_pdb_indices = {54,57,58,59,60,61,65,68,69,72,73,74,75,76,79,80,92,93,96,97,99,100,101,117,121,125} - {int(uaa_pdb_index), int(nucleophile_pdb_index)} # 62,91,53, 154,156,157,158
                    if args.option == "mpnn":
                        mpnn_path = os.path.join(variant_path, args.output_directory)
                        mpnn_design_sites_json = os.path.join(variant_path, "mpnn_design_sites.json")
                        if os.path.isfile(mpnn_design_sites_json):
                            with open(mpnn_design_sites_json, "r") as pf:
                                redesigned_res_pdb_indices = json.loads(pf.read())
                        else:
                            pose = pose_from_pdb(mpnn_input)
                            xlink_selector = ResidueNameSelector(xlink_name3)
                            uaa_pose_index = pose.pdb_info().pdb2pose(uaa_chain, int(uaa_pdb_index))
                            nucleophile_pose_index = pose.pdb_info().pdb2pose(nucleophile_chain, int(nucleophile_pdb_index))
                            xlink_selector = OrResidueSelector(xlink_selector, ResidueIndexSelector(\
                                str(uaa_pose_index) + "," + str(nucleophile_pose_index)))
                            neighborhood_selection = InterGroupInterfaceByVectorSelector()
                            neighborhood_selection.group1_selector(xlink_selector)
                            neighborhood_selection.group2_selector(NotResidueSelector(xlink_selector))
                            neighborhood_selection = AndResidueSelector(neighborhood_selection, \
                                NotResidueSelector(xlink_selector))
                            neighborhood_vector = neighborhood_selection.apply(pose)
                            neighborhood_pose_indices = boolean_vector_to_indices_set(neighborhood_vector)
                            pymol_res_pdb_indices = [uaa_pdb_index, nucleophile_pdb_index]
                            for pose_index in neighborhood_pose_indices:
                                info = pose.pdb_info().pose2pdb(int(pose_index)).split(" ")
                                # if info[1] in chains:
                                redesigned_res_pdb_indices.add(info[0])
                                pymol_res_pdb_indices.append(info[0])
                            redesigned_res_pdb_indices -= set(str(i) for i in range(3, 11))
                            redesigned_res_pdb_indices -= set(str(i) for i in range(133, 159))
                            with open(os.path.join(variant_path, "pymol.txt"), "w") as pf:
                                pf.write("show sticks, resi " + "+".join(pymol_res_pdb_indices) + "; hide sticks, /*//*/*/N; hide sticks, /*//*/*/O; hide sticks, /*//*/*/C; hide sticks, h.")
                    elif args.option == "mpnnall":
                        mpnn_path = os.path.join(variant_path, args.output_directory)
                        for pdb_index in range(11, 133):
                            redesigned_res_pdb_indices.add(str(pdb_index))
                    fixed_residues_set = {uaa_pdb_index, nucleophile_pdb_index}.union(set(args.catalytic_residues))
                    redesigned_res_pdb_indices -= fixed_residues_set
                    if uaa_pdb_index == "23" and nucleophile_pdb_index == "136":
                        redesigned_res_pdb_indices.add("134")
                    design_residues = " --redesigned_residues '"
                    symm_residues = " --symmetry_residues '"
                    symm_weight = " --symmetry_weights '"
                    for pdb_index in redesigned_res_pdb_indices:
                        design_residues += " A" + pdb_index + " B" + pdb_index
                        symm_residues += "A" + str(pdb_index) + ",B" + str(pdb_index) + "|"
                        symm_weight += "0.5,0.5|"
                    design_range = design_residues + "'" + symm_residues[:-1] + "'" + symm_weight[:-1] + "'"
                    design_range += " --fixed_residues '"
                    for pdb_index in fixed_residues_set:
                        design_range += " A" + pdb_index + " B" + pdb_index
                    design_range = design_range[:-1] + "'"
                    if os.path.isdir(mpnn_path):
                        shutil.rmtree(mpnn_path)
                    os.mkdir(mpnn_path)
                    with open(os.path.join(mpnn_path, "run.pbs"), "w") as pf:
                        pf.write("#!/bin/bash\n")
                        pf.write("#SBATCH --clusters=amarel\n")
                        pf.write("#SBATCH --partition=p_sdk94_1\n")
                        pf.write("#SBATCH --job-name=" + variant_pdb[:-4] + "\n")
                        pf.write("#SBATCH --nodes=1\n")
                        pf.write("#SBATCH --ntasks=5\n")
                        pf.write("#SBATCH --cpus-per-task=1\n")
                        pf.write("#SBATCH --gres=gpu:1\n")
                        pf.write("#SBATCH --mem=10000\n")
                        pf.write("#SBATCH --time=3-00:00:00\n")
                        pf.write("#SBATCH --output=slurm.%N.%j.log\n")
                        pf.write("#SBATCH --error=slurm.%N.%j.err\n")
                        pf.write("#SBATCH --requeue\n")
                        pf.write("#SBATCH --export=ALL\n")
                        pf.write("#SBATCH --begin=now\n")
                        pf.write("#SBATCH --open-mode=append\n")
                        pf.write("\n")
                        pf.write("srun -n 1 python /projects/f_sdk94_1/Tools/LigandMPNN/run.py" +
                                 " --model_type ligand_mpnn" +
                                 " --checkpoint_path_sc /projects/f_sdk94_1/Tools/LigandMPNN/model_params/ligandmpnn_sc_v_32_002_16.pt" +
                                 " --pdb_path ../" + variant_pdb.split(".")[0] + ".MPNN.pdb" +
                                 " --out_folder ./ --save_stats 1 --pack_side_chains 1 --number_of_packs_per_design 1" +
                                 " --batch_size 20 --number_of_batches 5" + design_range +
                                 " --ligand_mpnn_use_side_chain_context 1\n")
                        pf.write("exit\n")

                elif args.option == "barplot" and lowest_fasc is not None:
                    score_types = ["ddG", "UAA", "UAA_cst"]
                    path = lowest_fasc[:fasc.rfind("/")]
                    wt_fasc = os.path.join(path, pdb_chain + "_apo.sc")
                    wt_total_score = 0
                    if os.path.isfile(wt_fasc):
                        wt_total_score = read_scores_from_fasc(wt_fasc)[0]["total_score"]
                    scores = read_scores_from_fasc(fasc, n=args.read_n_decoys)
                    
                    ddg_values = list()
                    residue_energies = list()
                    residue_cst_scores = list()
                    for i, score_dict in enumerate(scores):
                        # decoy = os.path.join(path, score_dict["decoy"])
                        ddg = score_dict["total_score"] - wt_total_score
                        residue_energy = score_dict["substrates"]
                        residue_cst_score = score_dict["atom_pair_constraint"] + \
                                score_dict["angle_constraint"] + score_dict["dihedral_constraint"]
                        ddg_values.append(ddg)
                        residue_energies.append(residue_energy)
                        residue_cst_scores.append(residue_cst_score)
                        # residue_energy = 0
                        # residue_cst_terms = list()
                        # with open(decoy, "r") as pf:
                        #     for line in pf:
                        #         if line.startswith("pose "):
                        #             info = line[:-1].split(" ")
                        #             ddg_values.append(float(info[-1]) - float(info[-2]) - float(info[-13]) - wt_total_score)
                        #         elif line.startswith("TYZ:MP-OH-connect_") or line.startswith("ASX:MP-OD2-connect_") or \
                        #                 line.startswith("GLX:MP-OE2-connect_") or line.startswith("HIX:MP-NE2-connect_") or \
                        #                 line.startswith("CYX:MP-SG-connect_") or line.startswith("LYX:MP-NZ-connect_") or \
                        #                 line.startswith("SER:MP-SG-connect_") or line.startswith("THR:MP-NZ-connect_") or \
                        #                 line.startswith("TYR:MP-SG-connect_") or line.startswith("MET:MP-NZ-connect_") or \
                        #                 line.startswith("OBY:"):
                        #             info = line[:-1].split(" ")
                        #             residue_energy += (float(info[-1]) - float(info[-2]) - float(info[-13]))
                        #             residue_cst_terms.append(float(info[-11]))
                        #             residue_cst_terms.append(float(info[-12]))
                        #             residue_cst_terms.append(float(info[-14]))
                        # residue_energies.append(residue_energy)
                        # residue_cst.append(max(residue_cst_terms))
                    x = list(range(len(ddg_values)))
                    xticks = list(range(1, len(x) + 1))
                    for i, y_values in enumerate([ddg_values, residue_energies, residue_cst_scores]):
                        plt.bar(x, y_values, width=0.5, color='blue')
                        plt.xticks(x, xticks)
                        # plt.yticks(list(range(-7, 11)), list(range(-7, 11)))
                        plt.xlabel('trajectories')
                        plt.ylabel(score_types[i])
                        plt.title(variant_pdb_numbering)
                        # plt.show()
                        fig = plt.gcf()
                        # fig.set_size_inches(8, 6)
                        fig.savefig(os.path.join(pdb_ligand_path, variant_pdb_numbering + "." + score_types[i] + ".png"), dpi=100)
                        plt.figure().clear()
                        plt.close()
                        plt.cla()
                        plt.clf()
                    images = [Image.open(os.path.join(pdb_ligand_path, variant_pdb_numbering + "." + x + ".png")) for x in score_types]
                    widths, heights = zip(*(i.size for i in images))
                    if args.concatenate == "horizontal":
                        total_width = sum(widths)
                        max_height = max(heights)
                        new_im = Image.new('RGB', (total_width, max_height))
                        x_offset = 0
                        for im in images:
                            new_im.paste(im, (x_offset, 0))
                            x_offset += im.size[0]
                    elif args.concatenate == "vertical":
                        max_width = max(widths)
                        total_height = sum(heights)
                        new_im = Image.new('RGB', (max_width, total_height))
                        y_offset = 0
                        for im in images:
                            new_im.paste(im, (0, y_offset))
                            y_offset += im.size[1]
                    new_im.save(os.path.join(pdb_ligand_path, variant_pdb_numbering + ".png"))
                # break
            if args.option == "csv" or args.option == "mpnn" or args.option == "mpnnall":
                df_list = list(sorted(df_list, key=lambda var_dict: var_dict["ddG"][0]))
                df_list.insert(0, pd.DataFrame({c: pd.Series(dtype=t) for c, t in {"variant": "str", "ddG": "float", "xlink": "float", "xlink_cst": "float", "mutations": "str"}.items()}))
                df = pd.concat(df_list, sort=False)
                df.to_csv(os.path.join(pdb_ligand_path, pdb_chain_ligand + "." + args.step + ".csv"))
            # break
        # if args.option == "csv" or args.option == "mpnn":
        #     df.to_csv(os.path.join(pdb_ligand_path, pdb_chain + ".scores.csv"))
        # if args.option == "point":
        #     df.score_type = df.score_type.astype('category')
        #     df['decoy'] = range(1, len(df) + 1)
        #     sns.barplot(x='decoy', y='score', hue='score_type', hue_order=score_types, data=df, estimator=np.median)
        #     sns.swarmplot(x='decoy', y="score", hue='score_type', hue_order=score_types, data=df, color="0", alpha=.35)
        #     plt.savefig("test.png")
