#!/usr/bin/python3
import argparse
import math
import os
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueIndexSelector, ResiduePropertySelector, ChainSelector, \
    NeighborhoodResidueSelector, InterGroupInterfaceByVectorSelector

def detect_inter_group_interface(pose, chains):
    protein_selector = ResiduePropertySelector(ResidueProperty.PROTEIN)
    group1_selector = AndResidueSelector(ChainSelector(chains[0]), protein_selector)
    group2_selector = AndResidueSelector(ChainSelector(chains[1]), protein_selector)
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.group1_selector(group1_selector)
    interface_selector.group2_selector(group2_selector)
    group1_interface_selector = AndResidueSelector(interface_selector, group1_selector)
    group2_interface_selector = AndResidueSelector(interface_selector, group2_selector)
    group1_interface_vector = group1_interface_selector.apply(pose)
    group2_interface_vector = group2_interface_selector.apply(pose)
    return group1_interface_vector, group2_interface_vector

def convert_vector_to_list(res_vector):
    res_list = list()
    for index, res_bool in enumerate(res_vector):
        if res_bool:
            res_list.append(index + 1)
    return res_list

def write_position_files(prefix, res_list1, res_list2, workload):
    total_files = math.ceil(len(res_list1) * len(res_list2) / workload)
    pos_2_len = math.ceil(len(res_list2) / total_files)
    for i in range(0, total_files):
        with open(prefix + '_' + str(i + 1) + '.pos', 'w+') as p_pos:
            p_pos.write('N_CST 2\n')
            p_pos.write('1:')
            for res in res_list1:
                p_pos.write(' ' + str(res))
            p_pos.write('\n')
            p_pos.write('2:')
            if i < total_files - 1:
                for res in res_list2[i * pos_2_len: (i + 1) * pos_2_len]:
                    p_pos.write(' ' + str(res))
            else:
                for res in res_list2[i * pos_2_len:]:
                    p_pos.write(' ' + str(res))
            p_pos.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('-chain', type=str, help='Make staples within one chain.')
    parser.add_argument('-chains', type=str, nargs=2, help='Make staples across the interface of two chains.')
    parser.add_argument('-wl', '--workload', type=int, default=1000)
    args = parser.parse_args()

    params = list()
    for f in os.listdir():
        if f.endswith('.params'):
            params.append(f)
    if len(params) > 0:
        init('-extra_res_fa ' + ' '.join(params))
    else:
        init()
    pose = pose_from_pdb(args.pdb)
    prefix = args.pdb[args.pdb.rfind('/') + 1:-4]
    i = prefix.find('_')
    if i != -1:
        prefix = args.pdb[:i]
    if args.chain:
        # Residues having the same pdb chain id do not necessarily share the same pose chain id.
        # Convert pdb chain index to pose chain index.
        chain_selector = AndResidueSelector(ChainSelector(args.chain), ResiduePropertySelector(ResidueProperty.PROTEIN))
        chain_vector = chain_selector.apply(pose)
        chain_begin = None
        chain_end = None
        for res_index, res in enumerate(chain_vector):
            if res:
                if not chain_begin:
                    chain_begin = res_index + 1
                # chain_index = pose.chain(res_index + 1)
                # break
            else:
                if chain_begin:
                    chain_end = res_index
                    break
        if not chain_end:
            chain_end = res_index
        # Get the begining and the end position of the chain.
        # chain_begin = pose.chain_begin(chain_index)
        # chain_end = pose.chain_end(chain_index)
        chain_begin_res_list = list(range(chain_begin, chain_begin + 50))
        chain_end_res_list = list(range(chain_end - 49, chain_end + 1))
        prefix += '-' + args.chain
        os.mkdir(prefix)
        write_position_files(prefix + '/' + prefix + '_1', chain_begin_res_list, chain_end_res_list, args.workload)
        write_position_files(prefix + '/' + prefix + '_2', chain_end_res_list, chain_begin_res_list, args.workload)
    elif args.chains:
        group1_interface_vector, group2_interface_vector = detect_inter_group_interface(pose, args.chains)
        group1_interface_list = convert_vector_to_list(group1_interface_vector)
        group2_interface_list = convert_vector_to_list(group2_interface_vector)
        prefix += '-' + args.chains[0] + args.chains[1]
        os.mkdir(prefix)
        write_position_files(prefix + '/' + prefix, group1_interface_list, group2_interface_list, args.workload)
