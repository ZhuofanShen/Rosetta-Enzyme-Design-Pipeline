#!/usr/bin/python3
import argparse
import numpy as np
import os
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueNameSelector, ResidueIndexSelector, ResiduePropertySelector, ChainSelector, \
    NeighborhoodResidueSelector, InterGroupInterfaceByVectorSelector
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
import shutil


def aa_selector(aa_type, npz):
    threshold_phe = 2
    if "CYS" in aa_type:
        threshold_cys = 0.1
    threshold_his = 2
    if "HIS" in aa_type:
        threshold_his = 0.15
    threshold_glu = 2
    if "GLU" in aa_type:
        threshold_glu = 0.15
    threshold_asp = 2
    if "ASP" in aa_type:
        threshold_asp = 0.15
    threshold_lys = 2
    if "LYS" in aa_type:
        threshold_lys = 0.15
    if "PHE" or "TYR" in aa_type:
        threshold_phe = 0.15
    threshold_trp = 2
    if "TRP" in aa_type:
        threshold_trp = 0.15
    threshold_leu = 2
    if "LEU" in aa_type:
        threshold_leu = 0.15
    threshold_met = 2
    if "MET" in aa_type:
        threshold_met = 0.15
    threshold_cys = 2
    info = np.load(npz)
    mask = info["mask"]
    position_aa_pssm = info["log_p"][0][mask == 1]
    p_cys = np.exp(position_aa_pssm[:,1])
    p_his = np.exp(position_aa_pssm[:,6])
    p_glu = np.exp(position_aa_pssm[:,3])
    p_asp = np.exp(position_aa_pssm[:,2])
    p_lys = np.exp(position_aa_pssm[:,8])
    p_phe = np.exp(position_aa_pssm[:,4]) + np.exp(position_aa_pssm[:,19])
    p_trp = np.exp(position_aa_pssm[:,18])
    p_leu = np.exp(position_aa_pssm[:,9])
    p_met = np.exp(position_aa_pssm[:,10])
    p_aa_vector = list(filter(lambda index: index is not None, map(lambda i, j, k, u, v, w, x, y, z: str(w[0] + 1) \
            if i > threshold_cys or j > threshold_his or k > threshold_glu or u > threshold_asp or v > threshold_lys or \
            w[1] > threshold_phe or x > threshold_trp or y > threshold_leu or z > threshold_met else None, \
            p_cys, p_his, p_glu, p_asp, p_lys, enumerate(p_phe), p_trp, p_leu, p_met)))
    p_aa_selector = ResidueIndexSelector(",".join(p_aa_vector))
    aa_selector = ResidueNameSelector(",".join(aa_type))
    return OrResidueSelector(p_aa_selector, aa_selector)

def select_inter_group_interface(group1_selector, group2_selector):
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.nearby_atom_cut(7) # 5.5
    interface_selector.cb_dist_cut(12.0) # 11.0
    interface_selector.vector_dist_cut(10.5) # 9.0
    interface_selector.vector_angle_cut(76.0) # 75.0
    interface_selector.group1_selector(group1_selector)
    interface_selector.group2_selector(group2_selector)
    group1_interface_selector = AndResidueSelector(interface_selector, group1_selector)
    group2_interface_selector = AndResidueSelector(interface_selector, group2_selector)
    return group1_interface_selector, group2_interface_selector

def index_boolean_filter(index, boolean):
    if type(boolean) == np.ndarray:
        for bool in boolean:
            if bool:
                return index
    elif boolean:
        return index
    return None

def boolean_vector_to_indices_set(boolean_vector):
    return set(filter(lambda x: x is not None, map(index_boolean_filter, \
            range(1, len(boolean_vector) + 1), boolean_vector)))

def convert_vector_to_dict(pose, position, res_selector):
    res_vectors = list()
    res_vectors.append(res_selector.apply(pose))
    pose_copy = Pose(pose)
    mutator = MutateResidue()
    mutator.set_selector(ResidueIndexSelector(str(position)))
    mutator.set_res_name("TYR")
    mutator.apply(pose_copy)
    pose_copy.set_chi(1, position, 60)
    for chi2 in [90, -90]:
        pose_copy.set_chi(2, position, chi2)
        res_vectors.append(res_selector.apply(pose_copy))
    pose_copy.set_chi(1, position, 180)
    for chi2 in [75, -105]:
        pose_copy.set_chi(2, position, chi2)
        res_vectors.append(res_selector.apply(pose_copy))
    pose_copy.set_chi(1, position, -60)
    for chi2 in [105, -75]:
        pose_copy.set_chi(2, position, chi2)
        res_vectors.append(res_selector.apply(pose_copy))
    res_vectors = np.stack(res_vectors, axis=1)
    res_dict = dict()
    for res_index in boolean_vector_to_indices_set(res_vectors):
        nucleophile_name1 = pose.residue(res_index).name1()
        nucleophile_positions = res_dict.get(nucleophile_name1)
        if nucleophile_positions:
            nucleophile_positions.append(res_index)
        else:
            res_dict[nucleophile_name1] = [res_index]
    return res_dict

def find_nucleophile_positions(pose, uaa, nucleophile, group1_selector, group2_selector, npz):
    if npz is not None:
        if len(uaa) > 0:
            uaa_selector = aa_selector(uaa, npz)
        if len(nucleophile) > 0:
            nucleophile_selector = aa_selector(nucleophile, npz)
    else:
        if len(uaa) > 0:
            uaa_selector = ResidueNameSelector(",".join(uaa))
        if len(nucleophile) > 0:
            nucleophile_selector = ResidueNameSelector(",".join(nucleophile))
    if len(uaa) > 0:
        uaa_group1_selector = AndResidueSelector(group1_selector, uaa_selector)
    else:
        uaa_group1_selector = group1_selector
    uaa_group1_interface_positions = boolean_vector_to_indices_set(uaa_group1_selector.apply(pose))
    uaa_nucleophile_dict1 = dict()
    for position in uaa_group1_interface_positions:
        _, nucleophile_group2_selector = select_inter_group_interface(ResidueIndexSelector(str(position)), group2_selector)
        if len(nucleophile) > 0:
            nucleophile_group2_selector = AndResidueSelector(nucleophile_group2_selector, nucleophile_selector)
        # print(pose.residue(position).name1() + str(position))
        nucleophile_group2_dict = convert_vector_to_dict(pose, position, nucleophile_group2_selector)
        # print(nucleophile_group2_dict)
        if len(nucleophile_group2_dict) > 0:
            wildtype = pose.residue(position).name1()
            uaa_nucleophile_dict1[wildtype + str(position)] = nucleophile_group2_dict
    return uaa_nucleophile_dict1

def find_interface_uaa_nucleophile_positions(pose, chains, uaa, nucleophile, npz):
    protein_selector = ResiduePropertySelector(ResidueProperty.PROTEIN)
    chain1, chain2 = chains.split("-")
    group1_selector = AndResidueSelector(ChainSelector(chain1), protein_selector)
    group2_selector = AndResidueSelector(ChainSelector(chain2), protein_selector)
    # group1_interface_selector, group2_interface_selector = select_inter_group_interface(group1_selector, group2_selector)
    uaa_nucleophile_dict1 = find_nucleophile_positions(pose, uaa, nucleophile, group1_selector, group2_selector, npz)
    uaa_nucleophile_dict2 = find_nucleophile_positions(pose, uaa, nucleophile, group2_selector, group1_selector, npz)
    return uaa_nucleophile_dict1, uaa_nucleophile_dict2

def find_uaa_nucleophile_terminal_positions(pose, chain, uaa, nucleophile, ter, npz):
    protein_selector = ResiduePropertySelector(ResidueProperty.PROTEIN)
    chain_selector = AndResidueSelector(ChainSelector(chain), protein_selector)
    chain_positions = boolean_vector_to_indices_set(chain_selector.apply(pose))
    n_ter_site = min(chain_positions)
    n_ter_cutoff = n_ter_site + ter
    c_ter_site = max(chain_positions)
    c_ter_cutoff = c_ter_site - ter
    n_ter_selector = ResidueIndexSelector(",".join(str(r) for r in filter(lambda x: x >= n_ter_site and x < n_ter_cutoff, chain_positions)))
    c_ter_selector = ResidueIndexSelector(",".join(str(r) for r in filter(lambda x: x > c_ter_cutoff and x <= c_ter_site, chain_positions)))
    uaa_nucleophile_dict1 = find_nucleophile_positions(pose, uaa, nucleophile, n_ter_selector, c_ter_selector, npz)
    uaa_nucleophile_dict2 = find_nucleophile_positions(pose, uaa, nucleophile, c_ter_selector, n_ter_selector, npz)
    return uaa_nucleophile_dict1, uaa_nucleophile_dict2

def find_uaa_nucleophile_positions(pose, chain, uaa, nucleophile, npz):
    protein_selector = ResiduePropertySelector(ResidueProperty.PROTEIN)
    chain_selector = AndResidueSelector(ChainSelector(chain), protein_selector)
    if npz is not None:
        if len(uaa) > 0:
            uaa_selector = aa_selector(uaa, npz)
        if len(nucleophile) > 0:
            nucleophile_selector = aa_selector(nucleophile, npz)
    else:
        if len(uaa) > 0:
            uaa_selector = ResidueNameSelector(",".join(uaa))
        if len(nucleophile) > 0:
            nucleophile_selector = ResidueNameSelector(",".join(nucleophile))
    uaa_nucleophile_dict = dict()
    if len(uaa) > 0:
        uaa_chain_selector = AndResidueSelector(chain_selector, uaa_selector)
    else:
        uaa_chain_selector = chain_selector
    uaa_chain = uaa_chain_selector.apply(pose)
    for position in filter(lambda x: x is not None, map(lambda x, y: x + 1 if y else None, range(len(uaa_chain)), uaa_chain)):
        uaa_position_selector = ResidueIndexSelector(str(position))
        chain_not_uaa_position_selector = AndResidueSelector(chain_selector, NotResidueSelector(uaa_position_selector))
        _, nucleophile_chain_selector = select_inter_group_interface(uaa_position_selector, chain_not_uaa_position_selector)
        # if position == 5:
        #     print(position)
        #     nucleophile_chain_vector = nucleophile_chain_selector.apply(pose)
        #     print(nucleophile_chain_vector)
        #     s = str()
        #     for i, v in enumerate(nucleophile_chain_vector):
        #         if v:
        #             s += str(i + 1) + "+"
        #     print(s)
        #     print(nucleophile_chain_vector[126])
        #     nucleophile_vector = nucleophile_selector.apply(pose)
        #     print(nucleophile_vector[126])
        # elif position == 36:
        #     print(position)
        #     nucleophile_chain_vector = nucleophile_chain_selector.apply(pose)
        #     print(nucleophile_chain_vector)
        #     s = str()
        #     for i, v in enumerate(nucleophile_chain_vector):
        #         if v:
        #             s += str(i + 1) + "+"
        #     print(s)
        #     print(nucleophile_chain_vector[109])
        #     nucleophile_vector = nucleophile_selector.apply(pose)
        #     print(nucleophile_vector[109])
        if len(nucleophile) > 0:
            nucleophile_chain_selector = AndResidueSelector(nucleophile_chain_selector, nucleophile_selector)
        else:
            nucleophile_chain_selector = nucleophile_chain_selector
        nucleophile_chain_dict = convert_vector_to_dict(pose, position, nucleophile_chain_selector)
        if len(nucleophile_chain_dict) > 0:
            wildtype = pose.residue(position).name1()
            uaa_nucleophile_dict[wildtype + str(position)] = nucleophile_chain_dict
    
    return uaa_nucleophile_dict

def dict_to_position_files(prefix, pdb_chain, uaa_nucleophile_dict):
    if not os.path.isdir(os.path.join(prefix, pdb_chain)):
        os.mkdir(os.path.join(prefix, pdb_chain))
    for wildtype_position, nucleophile_position in uaa_nucleophile_dict.items():
        for nucleophile, position_list in nucleophile_position.items():
            for position in position_list:
                with open(os.path.join(prefix, pdb_chain, pdb_chain + "_" + wildtype_position + "_" + nucleophile + str(position) + '.pos'), 'w') as p_pos:
                    p_pos.write('N_CST 2\n')
                    p_pos.write('1: ' + wildtype_position[1:] + '\n')
                    p_pos.write('2: ' + str(position) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', type=str)
    parser.add_argument('-uaa', type=str, nargs="*", default=["ALL"]) # "PHE", "TYR", "TRP", "MET", "LEU", "HIS"
    parser.add_argument('-np', '--nucleophile', type=str, nargs="*", default=["CYS", "HIS", "ASP", "GLU", "LYS"])
    parser.add_argument('-ter', type=int, default=None)
    parser.add_argument('-pssm', '--use_pssm', action="store_true")
    args = parser.parse_args()

    if "ALL" in args.uaa:
        args.uaa = list()
    
    if "ALL" in args.nucleophile:
        args.nucleophile = list()

    init()

    for pdb_id in os.listdir(args.dir):
        prefix = os.path.join(args.dir, pdb_id)
        old_npz_path = os.path.join(prefix, "unconditional_probs_only")
        for pdb_chain in filter(lambda x: x.endswith(".pdb") and not x.endswith("_relaxed.pdb"), \
                os.listdir(os.path.join(args.dir, pdb_id))):
            try:
                pose = pose_from_pdb(os.path.join(prefix, pdb_chain))
            except:
                with open("ROSETTA_CLASH." + pdb_chain[:-4] + ".log", "w") as p_clash:
                    p_clash.write(pdb_chain + "\n")
                continue
            npz = None
            if args.use_pssm:
                old_npz = os.path.join(old_npz_path, pdb_chain[:-3] + "npz")
                new_npz = os.path.join(prefix, pdb_chain[:-3] + "npz")
                if os.path.isfile(old_npz):
                    os.rename(old_npz, new_npz)
                    npz = new_npz
                elif os.path.isfile(new_npz):
                    npz = new_npz
            i = pdb_chain.find('-')
            j = pdb_chain.rfind('-')
            if i == j:
                chain = pdb_chain[j + 1:-4]
                if args.ter:
                    uaa_nucleophile_dict1, uaa_nucleophile_dict2 = find_uaa_nucleophile_terminal_positions(pose, chain, args.uaa, args.nucleophile, args.ter, npz)
                    dict_to_position_files(prefix, pdb_chain[:-4], uaa_nucleophile_dict1)
                    dict_to_position_files(prefix, pdb_chain[:-4], uaa_nucleophile_dict2)
                else:
                    uaa_nucleophile_dict = find_uaa_nucleophile_positions(pose, chain, args.uaa, args.nucleophile, npz)
                    dict_to_position_files(prefix, pdb_chain[:-4], uaa_nucleophile_dict)
            else:
                chain1 = pdb_chain[i + 1:j]
                chain2 = pdb_chain[j + 1:-4]
                uaa_nucleophile_dict1, uaa_nucleophile_dict2 = find_interface_uaa_nucleophile_positions(pose, chain1 + "-" + chain2, args.uaa, args.nucleophile, npz)
                dict_to_position_files(prefix, pdb_chain[:-4], uaa_nucleophile_dict1)
                dict_to_position_files(prefix, pdb_chain[:-4], uaa_nucleophile_dict2)
        if os.path.isdir(old_npz_path):
            shutil.rmtree(old_npz_path)
        if os.path.isdir(os.path.join(prefix, "seq")):
            shutil.rmtree(os.path.join(prefix, "seq"))
        if os.path.isfile(os.path.join(prefix, "parsed_pdbs.jsonl")):
            os.remove(os.path.join(prefix, "parsed_pdbs.jsonl"))
