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
                pseudoWT_dG_entity = float(tmp[-1]) - float(tmp[-13])
            elif line.startswith("RRT") or line.startswith("SST") or \
                line.startswith("RST") or line.startswith("SRT"):
                tmp = line[:-1].split(" ")
                pseudoWT_dG_substrate = float(tmp[-1]) - float(tmp[-13])
                break
    return pseudoWT_dG_entity, pseudoWT_dG_substrate

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

def multi_state_fitness(pseudoWT_dG_fold, state_site_aa_specific_dG_entity, state_site_aa_specific_dG_substrate):
    site_aa_specific_dG_fold = state_site_aa_specific_dG_entity[-1, :, :]
    site_aa_specific_ddG_fold = site_aa_specific_dG_fold - pseudoWT_dG_fold
    state_site_aa_specific_dG_entity = state_site_aa_specific_dG_entity[:-1, :, :]
    state_site_aa_specific_dG_bind = state_site_aa_specific_dG_entity - site_aa_specific_dG_fold
    state_site_aa_specific_score = state_site_aa_specific_dG_bind
    npos = site_aa_specific_dG_fold.shape[0]
    configuration_site_aa_specific_lowest_score = torch.empty(0, npos, 20)
    configuration_site_aa_specific_dG_entity = torch.empty(0, npos, 20)
    configuration_site_aa_specific_dG_substrate = torch.empty(0, npos, 20)
    for i_state in range(0, state_site_aa_specific_score.shape[0], 8):
        site_aa_specific_lowest_score, site_aa_specific_most_stable_conformation = torch.min(state_site_aa_specific_score[i_state: i_state + 8, :, :], dim=0)
        site_aa_specific_most_favorable_state_dG_entity = torch.squeeze(torch.gather(state_site_aa_specific_dG_entity[i_state: i_state + 8, :, :], 0, \
                torch.unsqueeze(site_aa_specific_most_stable_conformation, 0)), 0)
        site_aa_specific_most_favorable_state_dG_substrate = torch.squeeze(torch.gather(state_site_aa_specific_dG_substrate[i_state: i_state + 8, :, :], 0, \
                torch.unsqueeze(site_aa_specific_most_stable_conformation, 0)), 0)
        configuration_site_aa_specific_lowest_score = torch.cat(\
            (configuration_site_aa_specific_lowest_score, site_aa_specific_lowest_score[None, :, :]), \
            dim=0)
        configuration_site_aa_specific_dG_entity = torch.cat(\
            (configuration_site_aa_specific_dG_entity, site_aa_specific_most_favorable_state_dG_entity[None, :, :]), \
            dim=0)
        configuration_site_aa_specific_dG_substrate = torch.cat(\
            (configuration_site_aa_specific_dG_substrate, site_aa_specific_most_favorable_state_dG_substrate[None, :, :]), \
            dim=0)
    site_aa_specific_most_competing_trans_diastereomer_score, site_aa_specific_most_competing_trans_diastereomer = torch.min(configuration_site_aa_specific_lowest_score[:2], dim=0)
    site_aa_specific_most_competing_trans_diastereomer_dG_entity = torch.squeeze(torch.gather(configuration_site_aa_specific_dG_entity[:2], 0, \
        torch.unsqueeze(site_aa_specific_most_competing_trans_diastereomer, 0)), 0)
    site_aa_specific_most_competing_trans_diastereomer_dG_substrate = torch.squeeze(torch.gather(configuration_site_aa_specific_dG_substrate[:2], 0, \
        torch.unsqueeze(site_aa_specific_most_competing_trans_diastereomer, 0)), 0)
    site_aa_specific_most_competing_cis_diastereomer_score, site_aa_specific_most_competing_cis_diastereomer = torch.min(configuration_site_aa_specific_lowest_score[2:], dim=0)
    site_aa_specific_most_competing_cis_diastereomer_dG_entity = torch.squeeze(torch.gather(configuration_site_aa_specific_dG_entity[2:], 0, \
        torch.unsqueeze(site_aa_specific_most_competing_cis_diastereomer, 0)), 0)
    site_aa_specific_most_competing_cis_diastereomer_dG_substrate = torch.squeeze(torch.gather(configuration_site_aa_specific_dG_substrate[2:], 0, \
        torch.unsqueeze(site_aa_specific_most_competing_cis_diastereomer, 0)), 0)
    state_site_aa_specific_score_stereoselectivity = state_site_aa_specific_score.clone()
    state_site_aa_specific_dG_stereoselectivity = state_site_aa_specific_dG_entity.clone()
    state_site_aa_specific_dG_substrate_stereoselectivity = state_site_aa_specific_dG_substrate.clone()
    state_site_aa_specific_score_enantioselectivity = state_site_aa_specific_score.clone()
    state_site_aa_specific_dG_enantioselectivity = state_site_aa_specific_dG_entity.clone()
    state_site_aa_specific_dG_substrate_enantioselectivity = state_site_aa_specific_dG_substrate.clone()
    state_site_aa_specific_score_diastereoselectivity = state_site_aa_specific_score.clone()
    state_site_aa_specific_dG_diastereoselectivity = state_site_aa_specific_dG_entity.clone()
    state_site_aa_specific_dG_substrate_diastereoselectivity = state_site_aa_specific_dG_substrate.clone()
    for i_configuration in range(configuration_site_aa_specific_lowest_score.shape[0]):
        trans_configurations = {0, 1}
        cis_configurations = {2, 3}
        if i_configuration in trans_configurations:
            enantiomer_configuration = next(iter(trans_configurations - {i_configuration}))
            site_aa_specific_most_competing_diastereomer_score = site_aa_specific_most_competing_cis_diastereomer_score
            site_aa_specific_most_competing_diastereomer_dG_entity = site_aa_specific_most_competing_cis_diastereomer_dG_entity
            site_aa_specific_most_competing_diastereomer_dG_substrate = site_aa_specific_most_competing_cis_diastereomer_dG_substrate
        elif i_configuration in cis_configurations:
            enantiomer_configuration = next(iter(cis_configurations - {i_configuration}))
            site_aa_specific_most_competing_diastereomer_score = site_aa_specific_most_competing_trans_diastereomer_score
            site_aa_specific_most_competing_diastereomer_dG_entity = site_aa_specific_most_competing_trans_diastereomer_dG_entity
            site_aa_specific_most_competing_diastereomer_dG_substrate = site_aa_specific_most_competing_trans_diastereomer_dG_substrate
        site_aa_specific_most_competing_enantiomer_score = configuration_site_aa_specific_lowest_score[enantiomer_configuration]
        site_aa_specific_most_competing_enantiomer_dG_entity = configuration_site_aa_specific_dG_entity[enantiomer_configuration]
        site_aa_specific_most_competing_enantiomer_dG_substrate = configuration_site_aa_specific_dG_substrate[enantiomer_configuration]
        site_aa_specific_most_competing_stereoisomer_score, site_aa_specific_most_competing_stereoisomer = torch.min( \
                torch.stack((site_aa_specific_most_competing_enantiomer_score, site_aa_specific_most_competing_diastereomer_score)), dim=0)
        site_aa_specific_most_competing_stereoisomer_dG_entity = torch.squeeze( \
            torch.gather(torch.stack((site_aa_specific_most_competing_enantiomer_dG_entity, site_aa_specific_most_competing_diastereomer_dG_entity)), 0, \
            torch.unsqueeze(site_aa_specific_most_competing_stereoisomer, 0)), 0)
        site_aa_specific_most_competing_stereoisomer_dG_substrate = torch.squeeze( \
            torch.gather(torch.stack((site_aa_specific_most_competing_enantiomer_dG_substrate, site_aa_specific_most_competing_diastereomer_dG_substrate)), 0, \
            torch.unsqueeze(site_aa_specific_most_competing_stereoisomer, 0)), 0)
        for i_state in range(8 * i_configuration, 8 * (i_configuration + 1)):
            state_site_aa_specific_score_stereoselectivity[i_state] = state_site_aa_specific_score_stereoselectivity[i_state] \
                    - site_aa_specific_most_competing_stereoisomer_score
            state_site_aa_specific_dG_stereoselectivity[i_state] = state_site_aa_specific_dG_stereoselectivity[i_state] \
                    - site_aa_specific_most_competing_stereoisomer_dG_entity
            state_site_aa_specific_dG_substrate_stereoselectivity[i_state] = state_site_aa_specific_dG_substrate_stereoselectivity[i_state] \
                    - site_aa_specific_most_competing_stereoisomer_dG_substrate
            state_site_aa_specific_score_enantioselectivity[i_state] = state_site_aa_specific_score_enantioselectivity[i_state] \
                    - site_aa_specific_most_competing_enantiomer_score
            state_site_aa_specific_dG_enantioselectivity[i_state] = state_site_aa_specific_dG_enantioselectivity[i_state] \
                    - site_aa_specific_most_competing_enantiomer_dG_entity
            state_site_aa_specific_dG_substrate_enantioselectivity[i_state] = state_site_aa_specific_dG_substrate_enantioselectivity[i_state] \
                    - site_aa_specific_most_competing_enantiomer_dG_substrate
            state_site_aa_specific_score_diastereoselectivity[i_state] = state_site_aa_specific_score_diastereoselectivity[i_state] \
                    - site_aa_specific_most_competing_diastereomer_score
            state_site_aa_specific_dG_diastereoselectivity[i_state] = state_site_aa_specific_dG_diastereoselectivity[i_state] \
                    - site_aa_specific_most_competing_diastereomer_dG_entity
            state_site_aa_specific_dG_substrate_diastereoselectivity[i_state] = state_site_aa_specific_dG_substrate_diastereoselectivity[i_state] \
                    - site_aa_specific_most_competing_diastereomer_dG_substrate
    state_site_aa_specific_fitness = state_site_aa_specific_dG_entity - pseudoWT_dG_fold + state_site_aa_specific_dG_stereoselectivity
    return site_aa_specific_ddG_fold, state_site_aa_specific_fitness, \
            state_site_aa_specific_dG_entity, state_site_aa_specific_dG_bind, \
            state_site_aa_specific_dG_stereoselectivity, \
            state_site_aa_specific_dG_substrate_stereoselectivity, \
            state_site_aa_specific_dG_enantioselectivity, \
            state_site_aa_specific_dG_substrate_enantioselectivity, \
            state_site_aa_specific_dG_diastereoselectivity, \
            state_site_aa_specific_dG_substrate_diastereoselectivity

def dump_fitness_matrix(fitness_score_matrix, position_aa_list, file_name, temperature):
    all_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    df = pd.DataFrame(fitness_score_matrix.numpy(), columns = all_aa, index = position_aa_list)
    df.to_csv(file_name + "_score.csv")
    if temperature:
        probability_matrix = torch.exp(torch.neg(fitness_score_matrix)/temperature)
        df = pd.DataFrame(probability_matrix.numpy(), columns = all_aa, index = position_aa_list)
        df.to_csv(file_name + ".csv")

def write_position_specific_scoring_matrix(directory, position_aa_list, temperature, \
            site_aa_specific_ddG_fold, state_site_aa_specific_fitness, state_site_aa_specific_dG_entity, \
            state_site_aa_specific_dG_bind, state_site_aa_specific_dG_stereoselectivity, \
            state_site_aa_specific_dG_enantioselectivity, state_site_aa_specific_dG_diastereoselectivity):
    protein = directory.split("_")[0]
    dump_fitness_matrix(site_aa_specific_ddG_fold, position_aa_list, \
            directory + "/" + protein + "/" + directory + ".stability", temperature)
    i_state = 0
    for configuration in ["1R2R", "1S2S", "1R2S", "1S2R"]:
        for carbene in ["rot1", "rot2", "rot3", "rot4"]:
            for ester in ["+", "-"]:
                state = configuration + "-" + carbene + ester
                protein_state = protein + "_" + state
                prefix = directory + "/" + protein_state + "/" + directory + "_" + state
                dump_fitness_matrix(state_site_aa_specific_fitness[i_state], position_aa_list, prefix + ".multi-state_fitness", temperature)
                dump_fitness_matrix(state_site_aa_specific_dG_entity[i_state], position_aa_list, prefix + ".single-state_fitness", temperature)
                dump_fitness_matrix(state_site_aa_specific_dG_bind[i_state], position_aa_list, prefix + ".activity", temperature)
                dump_fitness_matrix(state_site_aa_specific_dG_stereoselectivity[i_state], position_aa_list, prefix + ".stereoselectivity", temperature)
                dump_fitness_matrix(state_site_aa_specific_dG_enantioselectivity[i_state], position_aa_list, prefix + ".enantioselectivity", temperature)
                dump_fitness_matrix(state_site_aa_specific_dG_diastereoselectivity[i_state], position_aa_list, prefix + ".diastereoselectivity", temperature)
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
    site_aa_specific_ddG_fold, state_site_aa_specific_fitness, \
            state_site_aa_specific_dG_entity, state_site_aa_specific_dG_bind, \
            state_site_aa_specific_dG_stereoselectivity, \
            state_site_aa_specific_dG_substrate_stereoselectivity, \
            state_site_aa_specific_dG_enantioselectivity, \
            state_site_aa_specific_dG_substrate_enantioselectivity, \
            state_site_aa_specific_dG_diastereoselectivity, \
            state_site_aa_specific_dG_substrate_diastereoselectivity = \
            multi_state_fitness(state_specific_pseudoWT_dG[-1], \
                    state_site_aa_specific_dG, state_site_aa_specific_unnormalized_dG_substrate)
    write_position_specific_scoring_matrix(directory, pseudoWT_position_list, temperature, \
            site_aa_specific_ddG_fold, state_site_aa_specific_fitness, state_site_aa_specific_dG_entity, \
            state_site_aa_specific_dG_bind, state_site_aa_specific_dG_stereoselectivity, \
            state_site_aa_specific_dG_enantioselectivity, state_site_aa_specific_dG_diastereoselectivity)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
