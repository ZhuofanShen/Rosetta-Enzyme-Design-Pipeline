import argparse
import json
import numpy as np
import os
import pandas as pd
import torch


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=str)
    parser.add_argument("-o", "--override", action="store_true")
    parser.add_argument("-t", "--temperature", type=float, default=4)
    args = parser.parse_args()
    return args

def get_pseudoWT_sequence(directory, override=False):
    position_pseudoWT_dict = dict()
    if not override:
        for npy in filter(lambda x: x.endswith(".npy"), os.listdir(directory)):
            pseudoWT_position = npy.split(".")[-2]
            position = pseudoWT_position[1:]
            if position.isdigit():
                position_pseudoWT_dict[position] = pseudoWT_position[0]
    sc_folder = os.path.join(directory, directory.split("_")[0])
    if os.path.isdir(sc_folder):
        for sc in filter(lambda x: x.endswith("A.sc"), \
                os.listdir(sc_folder)):
            pseudoWT_position = sc.split(".")[-2][:-1]
            pseudoWT = position_pseudoWT_dict.get(pseudoWT_position[1:])
            if pseudoWT is not None and pseudoWT != pseudoWT_position[0]:
                raise Exception("Pseudo-wildtype of the sc file " + pseudoWT_position[0] + " at position " + \
                        pseudoWT_position[1:] + " is not consistent with the npy file " + pseudoWT + ".")
            position_pseudoWT_dict[pseudoWT_position[1:]] = pseudoWT_position[0]
    all_AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    pseudoWT_position_list = list()
    pseudoWT_one_hot = list()
    for position, pseudoWT in sorted(position_pseudoWT_dict.items(), key=lambda item: int(item[0])):
        pseudoWT_position_list.append(pseudoWT + position)
        pseudoWT_one_hot.append(all_AA.index(pseudoWT))
    return pseudoWT_position_list, np.array(pseudoWT_one_hot)

def get_pseudoWT_dG(protein_path):
    pseudoWT_dG_substrate = None
    with open(protein_path, "r") as pf:
        for line in pf:
            if line.startswith("pose "):
                tmp = line[:-1].split(" ")
                pseudoWT_dG = float(tmp[-1]) - float(tmp[-13])
            elif line.startswith("RRT") or line.startswith("SST") or \
                line.startswith("RST") or line.startswith("SRT"):
                tmp = line[:-1].split(" ")
                pseudoWT_dG_substrate = float(tmp[-1]) - float(tmp[-13])
                break
    return pseudoWT_dG, pseudoWT_dG_substrate

def read_score_file(sc_path):
    with open(sc_path, "r") as pf:
        for line in pf:
            break
    try:
        scores = json.loads(line)
    except:
        print(sc_path)
    dG = scores["total_score"] - scores["coordinate_constraint"]
    dG_substrate = scores.get("substrates")
    return dG, dG_substrate

def create_dG_matrix(directory, pseudoWT_position_list, override=False):
    protein = directory.split("_")[0]
    pseudoWT_npy = os.path.join(directory, protein + ".pseudoWT.npy")
    if os.path.isfile(pseudoWT_npy) and not override:
        state_specific_pseudoWT_dGs = np.load(pseudoWT_npy, allow_pickle=True)
    else:
        state_specific_pseudoWT_dGs = np.empty((0, 2), dtype=np.float32)
        for configuration in ["1R2R", "1S2S", "1R2S", "1S2R"]:
            for carbene in ["rot1", "rot2", "rot3", "rot4"]:
                for ester in ["+", "-"]:
                    state = configuration + "-" + carbene + ester
                    protein_state = protein + "_" + state
                    pseudoWT_dG, pseudoWT_dG_substrate = get_pseudoWT_dG(\
                            os.path.join(directory, "input", protein_state + ".pdb"))
                    state_specific_pseudoWT_dGs = np.concatenate((state_specific_pseudoWT_dGs, \
                            np.asarray([[pseudoWT_dG, pseudoWT_dG_substrate]])), axis=0)
        pseudoWT_dG_fold, _ = get_pseudoWT_dG(os.path.join(directory, "input", protein + ".pdb"))
        state_specific_pseudoWT_dGs = np.concatenate((state_specific_pseudoWT_dGs, \
                np.asarray([[pseudoWT_dG_fold, None]])), axis=0)
        np.save(pseudoWT_npy, state_specific_pseudoWT_dGs)
    state_specific_pseudoWT_dGs = state_specific_pseudoWT_dGs.astype(np.float32)
    state_specific_pseudoWT_dG = torch.from_numpy(state_specific_pseudoWT_dGs[:, 0])
    state_specific_pseudoWT_dG_substrate = torch.from_numpy(state_specific_pseudoWT_dGs[:-1, 1])
    all_AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    n_AA = len(all_AA)
    site_state_aa_specific_dG = np.empty((0, 33, n_AA), dtype=np.float32)
    site_state_aa_specific_dG_substrate = np.empty((0, 32, n_AA), dtype=np.float32)
    for pseudoWT_position in pseudoWT_position_list:
        site_npy = os.path.join(directory, protein + "." + pseudoWT_position + ".npy")
        if os.path.isfile(site_npy) and not override:
            state_aa_specific_dGs = np.load(site_npy, allow_pickle=True)
        else:
            state_aa_specific_dGs = np.empty((0, n_AA, 2), dtype=np.float32)
            for configuration in ["1R2R", "1S2S", "1R2S", "1S2R"]:
                for carbene in ["rot1", "rot2", "rot3", "rot4"]:
                    for ester in ["+", "-"]:
                        state = configuration + "-" + carbene + ester
                        protein_state = protein + "_" + state
                        aa_specific_dGs = np.empty((0, 2), dtype=np.float32)
                        for aa in all_AA:
                            dG, dG_substrate = read_score_file(os.path.join(directory, \
                                    protein_state, protein_state + "." + pseudoWT_position + aa + ".sc"))
                            aa_specific_dGs = np.concatenate((aa_specific_dGs, \
                                    np.asarray([[dG, dG_substrate]])), axis=0)
                        state_aa_specific_dGs = np.concatenate((state_aa_specific_dGs, \
                                aa_specific_dGs[None, :, :]), axis=0)
            aa_specific_dGs = np.empty((0, 2), dtype=np.float32)
            for aa in all_AA:
                dG_fold, _ = read_score_file(os.path.join(directory, protein, \
                        protein + "." + pseudoWT_position + aa + ".sc"))
                aa_specific_dGs = np.concatenate((aa_specific_dGs, np.asarray([[dG_fold, None]])), axis=0)
            state_aa_specific_dGs = np.concatenate((state_aa_specific_dGs, \
                    aa_specific_dGs[None, :, :]), axis=0)
            np.save(site_npy, state_aa_specific_dGs)
        site_state_aa_specific_dG = np.concatenate((site_state_aa_specific_dG, \
                state_aa_specific_dGs[None, :, :, 0]), axis=0)
        site_state_aa_specific_dG_substrate = np.concatenate((site_state_aa_specific_dG_substrate, \
                state_aa_specific_dGs[None, :-1, :, 1]), axis=0)
    state_site_aa_specific_dG = torch.from_numpy(np.transpose(\
            site_state_aa_specific_dG, (1, 0, 2)).astype(np.float32))
    state_site_aa_specific_dG_substrate = torch.from_numpy(np.transpose(\
            site_state_aa_specific_dG_substrate, (1, 0, 2)).astype(np.float32))
    return state_specific_pseudoWT_dG, state_specific_pseudoWT_dG_substrate, \
            state_site_aa_specific_dG, state_site_aa_specific_dG_substrate

def normalize_dG(pseudoWT_one_hot, state_specific_pseudoWT_dG, state_site_aa_specific_dG):
    state_site_specific_pseudoWT_dG = state_site_aa_specific_dG\
            [:, torch.arange(pseudoWT_one_hot.shape[0]), pseudoWT_one_hot]
    state_site_aa_specific_ddG = state_site_aa_specific_dG - \
            state_site_specific_pseudoWT_dG[:, :, None]
    state_site_aa_specific_dG = state_site_aa_specific_ddG + \
            state_specific_pseudoWT_dG[:, None, None]
    return state_site_aa_specific_dG

def multistate_design_fitness(state_specific_pseudoWT_dG, state_site_aa_specific_dG, \
        temperature=4):
    site_aa_specific_dG_fold = state_site_aa_specific_dG[-1, :, :]
    state_site_aa_specific_dG = state_site_aa_specific_dG[:-1, :, :]
    state_site_aa_specific_dG_bind = state_site_aa_specific_dG - site_aa_specific_dG_fold
    state_site_aa_specific_partition = torch.exp(torch.neg(state_site_aa_specific_dG_bind)/temperature)
    pseudoWT_dG_fold = state_specific_pseudoWT_dG[-1]
    state_specific_pseudoWT_dG = state_specific_pseudoWT_dG[:-1]
    state_specific_pseudoWT_dG_bind = state_specific_pseudoWT_dG - pseudoWT_dG_fold
    state_specific_pseudoWT_partition = torch.exp(torch.neg(state_specific_pseudoWT_dG_bind)/temperature)
    site_aa_specific_ddG_fold = site_aa_specific_dG_fold - pseudoWT_dG_fold
    site_aa_specific_ddG_fold_boltzmann_factor = torch.exp(torch.neg(site_aa_specific_ddG_fold)/temperature)
    state_site_aa_specific_ddG_bind = torch.transpose(torch.sub(torch.transpose(state_site_aa_specific_dG_bind, 0, 2), state_specific_pseudoWT_dG_bind), 0, 2)
    state_site_aa_specific_ddG_bind_boltzmann_factor = torch.exp(torch.neg(state_site_aa_specific_ddG_bind)/temperature)
    state_site_aa_specific_ddG = torch.transpose(torch.sub(torch.transpose(state_site_aa_specific_dG, 0, 2), state_specific_pseudoWT_dG), 0, 2)
    state_site_aa_specific_ddG_boltzmann_factor = torch.exp(torch.neg(state_site_aa_specific_ddG)/temperature)
    n_pos = state_site_aa_specific_dG.shape[1]
    n_AA = state_site_aa_specific_dG.shape[2]
    configuration_site_aa_specific_partition_sum = torch.empty(0, n_pos, n_AA)
    configuration_specific_pseudoWT_partition_sum = torch.empty(0)
    configuration_site_aa_specific_ddG_bind_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    for i_state in range(0, state_site_aa_specific_partition.shape[0], 8):
        site_aa_specific_conformer_partition_sum = torch.sum(state_site_aa_specific_partition[i_state: i_state + 8, :, :], dim=0)
        configuration_site_aa_specific_partition_sum = torch.cat(\
            (configuration_site_aa_specific_partition_sum, site_aa_specific_conformer_partition_sum[None, :, :]), \
            dim=0)
        pseudoWT_conformer_partition_sum = torch.sum(state_specific_pseudoWT_partition[i_state: i_state + 8], dim=0)
        configuration_specific_pseudoWT_partition_sum = torch.cat(\
            (configuration_specific_pseudoWT_partition_sum, pseudoWT_conformer_partition_sum[None]), \
            dim=0)
        site_aa_specific_ddG_bind_boltzmann_factor = torch.div(site_aa_specific_conformer_partition_sum, pseudoWT_conformer_partition_sum)
        configuration_site_aa_specific_ddG_bind_boltzmann_factor = torch.cat(\
            (configuration_site_aa_specific_ddG_bind_boltzmann_factor, site_aa_specific_ddG_bind_boltzmann_factor[None, :, :]), \
            dim=0)
    configuration_site_aa_specific_ddG_boltzmann_factor = torch.mul(\
            configuration_site_aa_specific_ddG_bind_boltzmann_factor, \
            site_aa_specific_ddG_fold_boltzmann_factor)
    site_aa_specific_unfavorable_trans_diastereomers_partition_sum = torch.sum(configuration_site_aa_specific_partition_sum[:2], dim=0)
    site_aa_specific_unfavorable_cis_diastereomers_partition_sum = torch.sum(configuration_site_aa_specific_partition_sum[2:], dim=0)
    configuration_site_aa_specific_enantio_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    configuration_site_aa_specific_diastereo_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    configuration_site_aa_specific_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    state_site_aa_specific_enantio_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    state_site_aa_specific_diastereo_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    state_site_aa_specific_fitness_boltzmann_factor = torch.empty(0, n_pos, n_AA)
    for i_configuration in range(configuration_site_aa_specific_ddG_bind_boltzmann_factor.shape[0]):
        trans_configurations = {0, 1}
        cis_configurations = {2, 3}
        if i_configuration in trans_configurations:
            site_aa_specific_unfavorable_enantiomers_partition_sum = configuration_site_aa_specific_partition_sum[next(iter(trans_configurations - {i_configuration}))]
            site_aa_specific_unfavorable_diastereomers_partition_sum = site_aa_specific_unfavorable_cis_diastereomers_partition_sum
        elif i_configuration in cis_configurations:
            site_aa_specific_unfavorable_enantiomers_partition_sum = configuration_site_aa_specific_partition_sum[next(iter(cis_configurations - {i_configuration}))]
            site_aa_specific_unfavorable_diastereomers_partition_sum = site_aa_specific_unfavorable_trans_diastereomers_partition_sum
        site_aa_specific_enantioselectivity_boltzmann_factor = torch.div(configuration_site_aa_specific_partition_sum[i_configuration, :, :], \
                site_aa_specific_unfavorable_enantiomers_partition_sum)
        site_aa_specific_enantio_fitness_boltzmann_factor = torch.mul(\
                configuration_site_aa_specific_ddG_boltzmann_factor[i_configuration], \
                site_aa_specific_enantioselectivity_boltzmann_factor)
        configuration_site_aa_specific_enantio_fitness_boltzmann_factor = torch.cat(\
                (configuration_site_aa_specific_enantio_fitness_boltzmann_factor, \
                site_aa_specific_enantio_fitness_boltzmann_factor[None, :, :]), dim=0)
        site_aa_specific_diastereoselectivity_boltzmann_factor = torch.div(configuration_site_aa_specific_partition_sum[i_configuration, :, :], \
                site_aa_specific_unfavorable_diastereomers_partition_sum)
        site_aa_specific_diastereo_fitness_boltzmann_factor = torch.mul(\
                configuration_site_aa_specific_ddG_boltzmann_factor[i_configuration], \
                site_aa_specific_diastereoselectivity_boltzmann_factor)
        configuration_site_aa_specific_diastereo_fitness_boltzmann_factor = torch.cat(\
                (configuration_site_aa_specific_diastereo_fitness_boltzmann_factor, \
                site_aa_specific_diastereo_fitness_boltzmann_factor[None, :, :]), dim=0)
        site_aa_specific_unfavorable_stereoisomers_partition_sum = site_aa_specific_unfavorable_enantiomers_partition_sum + site_aa_specific_unfavorable_diastereomers_partition_sum
        site_aa_specific_stereoselectivity_boltzmann_factor = torch.div(configuration_site_aa_specific_partition_sum[i_configuration, :, :], \
                site_aa_specific_unfavorable_stereoisomers_partition_sum)
        site_aa_specific_fitness_boltzmann_factor = torch.mul(\
                configuration_site_aa_specific_ddG_boltzmann_factor[i_configuration], \
                site_aa_specific_stereoselectivity_boltzmann_factor)
        configuration_site_aa_specific_fitness_boltzmann_factor = torch.cat(\
                (configuration_site_aa_specific_fitness_boltzmann_factor, \
                site_aa_specific_fitness_boltzmann_factor[None, :, :]), dim=0)
        for i_state in range(8 * i_configuration, 8 * (i_configuration + 1)):
            site_aa_specific_enantioselectivity_boltzmann_factor = torch.div(\
                state_site_aa_specific_partition[i_state, :, :], \
                site_aa_specific_unfavorable_enantiomers_partition_sum)
            site_aa_specific_enantio_fitness_boltzmann_factor = torch.mul(\
                    state_site_aa_specific_ddG_boltzmann_factor[i_state], \
                    site_aa_specific_enantioselectivity_boltzmann_factor)
            state_site_aa_specific_enantio_fitness_boltzmann_factor = torch.cat(\
                    (state_site_aa_specific_enantio_fitness_boltzmann_factor, \
                    site_aa_specific_enantio_fitness_boltzmann_factor[None, :, :]), dim=0)
            site_aa_specific_diastereoselectivity_boltzmann_factor = torch.div(\
                state_site_aa_specific_partition[i_state, :, :], \
                site_aa_specific_unfavorable_diastereomers_partition_sum)
            site_aa_specific_diastereo_fitness_boltzmann_factor = torch.mul(\
                    state_site_aa_specific_ddG_boltzmann_factor[i_state], \
                    site_aa_specific_diastereoselectivity_boltzmann_factor)
            state_site_aa_specific_diastereo_fitness_boltzmann_factor = torch.cat(\
                    (state_site_aa_specific_diastereo_fitness_boltzmann_factor, \
                    site_aa_specific_diastereo_fitness_boltzmann_factor[None, :, :]), dim=0)
            site_aa_specific_stereoselectivity_boltzmann_factor = torch.div(\
                state_site_aa_specific_partition[i_state, :, :], \
                site_aa_specific_unfavorable_stereoisomers_partition_sum)
            site_aa_specific_fitness_boltzmann_factor = torch.mul(\
                    state_site_aa_specific_ddG_boltzmann_factor[i_state], \
                    site_aa_specific_stereoselectivity_boltzmann_factor)
            state_site_aa_specific_fitness_boltzmann_factor = torch.cat(\
                    (state_site_aa_specific_fitness_boltzmann_factor, \
                    site_aa_specific_fitness_boltzmann_factor[None, :, :]), dim=0)
    return site_aa_specific_ddG_fold, site_aa_specific_ddG_fold_boltzmann_factor, \
            state_site_aa_specific_ddG_bind_boltzmann_factor, \
            state_site_aa_specific_ddG_boltzmann_factor, \
            state_site_aa_specific_enantio_fitness_boltzmann_factor, \
            state_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            state_site_aa_specific_fitness_boltzmann_factor, \
            configuration_site_aa_specific_ddG_bind_boltzmann_factor, \
            configuration_site_aa_specific_ddG_boltzmann_factor, \
            configuration_site_aa_specific_enantio_fitness_boltzmann_factor, \
            configuration_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            configuration_site_aa_specific_fitness_boltzmann_factor

def dump_fitness_matrix(matrix, pseudoWT_position_list, prefix, temperature=None):
    all_AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    if temperature is not None:
        matrix = torch.neg(torch.log(matrix) * temperature)
    df = pd.DataFrame(matrix.numpy(), columns = all_AA, index = pseudoWT_position_list)
    df.to_csv(prefix + ".csv")

def write_position_specific_scoring_matrix(directory, pseudoWT_position_list, temperature, \
        site_aa_specific_ddG_fold, site_aa_specific_ddG_fold_boltzmann_factor, \
        state_site_aa_specific_ddG_bind_boltzmann_factor, \
        state_site_aa_specific_ddG_boltzmann_factor, \
        state_site_aa_specific_enantio_fitness_boltzmann_factor, \
        state_site_aa_specific_diastereo_fitness_boltzmann_factor, \
        state_site_aa_specific_fitness_boltzmann_factor, \
        configuration_site_aa_specific_ddG_bind_boltzmann_factor, \
        configuration_site_aa_specific_ddG_boltzmann_factor, \
        configuration_site_aa_specific_enantio_fitness_boltzmann_factor, \
        configuration_site_aa_specific_diastereo_fitness_boltzmann_factor, \
        configuration_site_aa_specific_fitness_boltzmann_factor):
    pdb = directory.split("_")[0]
    sub_dir = os.path.join(directory, pdb)
    if not os.path.isdir(sub_dir):
        os.mkdir(sub_dir)
    dump_fitness_matrix(site_aa_specific_ddG_fold, pseudoWT_position_list, \
            os.path.join(sub_dir, directory + ".stability_score"))
    dump_fitness_matrix(site_aa_specific_ddG_fold_boltzmann_factor, pseudoWT_position_list, \
            os.path.join(sub_dir, directory + ".stability"))
    i_state = 0
    for i_configuration, configuration in enumerate(["1R2R", "1S2S", "1R2S", "1S2R"]):
        sub_dir = os.path.join(directory, pdb + "_" + configuration)
        if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)
        prefix = os.path.join(sub_dir, directory + "_" + configuration)
        dump_fitness_matrix(configuration_site_aa_specific_ddG_bind_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".activity_score", temperature=temperature)
        dump_fitness_matrix(configuration_site_aa_specific_ddG_bind_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".activity")
        dump_fitness_matrix(configuration_site_aa_specific_ddG_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".positive_design_score", temperature=temperature)
        dump_fitness_matrix(configuration_site_aa_specific_ddG_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".positive_design_fitness")
        dump_fitness_matrix(configuration_site_aa_specific_enantio_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".enantio_MSD_score", temperature=temperature)
        dump_fitness_matrix(configuration_site_aa_specific_enantio_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".enantio_MSD_fitness")
        dump_fitness_matrix(configuration_site_aa_specific_diastereo_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".diastereo_MSD_score", temperature=temperature)
        dump_fitness_matrix(configuration_site_aa_specific_diastereo_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".diastereo_MSD_fitness")
        dump_fitness_matrix(configuration_site_aa_specific_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".multi-state_design_score", temperature=temperature)
        dump_fitness_matrix(configuration_site_aa_specific_fitness_boltzmann_factor[i_configuration], \
                pseudoWT_position_list, prefix + ".multi-state_design_fitness")
        for carbene in ["rot1", "rot2", "rot3", "rot4"]:
            for ester in ["+", "-"]:
                state = configuration + "-" + carbene + ester
                sub_dir = os.path.join(directory, pdb + "_" + state)
                if not os.path.isdir(sub_dir):
                    os.mkdir(sub_dir)
                prefix = os.path.join(sub_dir, directory + "_" + state)
                dump_fitness_matrix(state_site_aa_specific_ddG_bind_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".activity_score", temperature=temperature)
                dump_fitness_matrix(state_site_aa_specific_ddG_bind_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".activity")
                dump_fitness_matrix(state_site_aa_specific_ddG_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".positive_design_score", temperature=temperature)
                dump_fitness_matrix(state_site_aa_specific_ddG_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".positive_design_fitness")
                dump_fitness_matrix(state_site_aa_specific_enantio_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".enantio_MSD_score", temperature=temperature)
                dump_fitness_matrix(state_site_aa_specific_enantio_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".enantio_MSD_fitness")
                dump_fitness_matrix(state_site_aa_specific_diastereo_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".diastereo_MSD_score", temperature=temperature)
                dump_fitness_matrix(state_site_aa_specific_diastereo_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".diastereo_MSD_fitness")
                dump_fitness_matrix(state_site_aa_specific_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".multi-state_design_score", temperature=temperature)
                dump_fitness_matrix(state_site_aa_specific_fitness_boltzmann_factor[i_state], \
                        pseudoWT_position_list, prefix + ".multi-state_design_fitness")
                i_state += 1

def main(args):
    directory = args.directory
    override = args.override
    temperature = args.temperature
    pseudoWT_position_list, pseudoWT_one_hot = get_pseudoWT_sequence(directory, override=override)
    state_specific_pseudoWT_dG, state_specific_pseudoWT_dG_substrate, \
            state_site_aa_specific_unnormalized_dG, \
            state_site_aa_specific_unnormalized_dG_substrate = \
            create_dG_matrix(directory, pseudoWT_position_list, override=override)
    state_site_aa_specific_dG = normalize_dG(pseudoWT_one_hot, \
            state_specific_pseudoWT_dG, state_site_aa_specific_unnormalized_dG)
    site_aa_specific_ddG_fold, site_aa_specific_ddG_fold_boltzmann_factor, \
            state_site_aa_specific_ddG_bind_boltzmann_factor, \
            state_site_aa_specific_ddG_boltzmann_factor, \
            state_site_aa_specific_enantio_fitness_boltzmann_factor, \
            state_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            state_site_aa_specific_fitness_boltzmann_factor, \
            configuration_site_aa_specific_ddG_bind_boltzmann_factor, \
            configuration_site_aa_specific_ddG_boltzmann_factor, \
            configuration_site_aa_specific_enantio_fitness_boltzmann_factor, \
            configuration_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            configuration_site_aa_specific_fitness_boltzmann_factor = \
            multistate_design_fitness(state_specific_pseudoWT_dG, \
                    state_site_aa_specific_dG, temperature=temperature)
    write_position_specific_scoring_matrix(directory, pseudoWT_position_list, temperature, \
            site_aa_specific_ddG_fold, site_aa_specific_ddG_fold_boltzmann_factor, \
            state_site_aa_specific_ddG_bind_boltzmann_factor, \
            state_site_aa_specific_ddG_boltzmann_factor, \
            state_site_aa_specific_enantio_fitness_boltzmann_factor, \
            state_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            state_site_aa_specific_fitness_boltzmann_factor, \
            configuration_site_aa_specific_ddG_bind_boltzmann_factor, \
            configuration_site_aa_specific_ddG_boltzmann_factor, \
            configuration_site_aa_specific_enantio_fitness_boltzmann_factor, \
            configuration_site_aa_specific_diastereo_fitness_boltzmann_factor, \
            configuration_site_aa_specific_fitness_boltzmann_factor)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
