#!/usr/bin/python3
'''
score function
symmetry
coord cst
task factory
    read tyz cyx from remark lines
move map
'''
import argparse
import os
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    InitializeFromCommandline, IncludeCurrent, ExtraRotamers, \
    OperateOnResidueSubset, RestrictToRepackingRLT, \
    RestrictAbsentCanonicalAASRLT, PreventRepackingRLT, \
    RestrictToRepacking
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueIndexSelector, NeighborhoodResidueSelector, \
    InterGroupInterfaceByVectorSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    RMSDMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, \
    AddOrRemoveMatchCsts, DetectProteinLigandInterface
from pyrosetta.rosetta.protocols.protein_interface_design \
    import FavorNativeResidue
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('-params', type=str, nargs='*')
    parser.add_argument('-ft', '--fold_tree', type=str, nargs='*')
    parser.add_argument('-chis', type=str, nargs='*')
    parser.add_argument('-sf', '--score_function', type=str)
    parser.add_argument('-symm', '--symmetry', type=str)
    # parser.add_argument('-dup', '--duplicated_chains', type=str, nargs='*')
    parser.add_argument('-enzdescst', type=str)
    parser.add_argument('-cst', type=str)
    parser.add_argument('-nbh', type=float)
    parser.add_argument('-nataa', type=bool)
    parser.add_argument('--enzdes', action='store_true')
    parser.add_argument('-muts', type=str, nargs='*')
    parser.add_argument('-des', type=str, nargs='*')
    parser.add_argument('-decoys', type=int)
    parser.add_argument('-rmsd', type=str)
    args = parser.parse_args()
    return args

def init_pyrosetta_with_opts(args):
    opts = '-ex1 -ex2 -ignore_zero_occupancy false \
        -use_input_sc' # -nblist_autoupdate
    if args.score_function.startswith('beta_nov16'):
        opts += ' -corrections::beta_nov16'
    if args.params:
        opts += " -extra_res_fa"
        for params_file in args.params:
            opts = opts + ' ' + params_file
    if args.enzdescst:
        opts += ' -enzdes:cstfile {} -run:preserve_header'.format(args.enzdescst)
    if args.cst:
        opts += ' -constraints:cst_fa_file {}'.format(args.cst)
    init(opts)

def create_fold_tree(edge_list):
    fold_tree = FoldTree()
    for edge_str in edge_list:
        edge_num = edge_str.split(',')
        if len(edge_num) is 3:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), int(edge_num[2]))
        elif len(edge_num) is 4:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), edge_num[2], edge_num[3])
        else:
            raise Exception("The number of arguments in add_edge function should be 3 or 4")
    return fold_tree

def create_score_fc(sfxn='ref2015', weights=None):
    score_function = create_score_function(sfxn)
    if sfxn.startswith('ref2015'):
        score_function.set_weight(ScoreType.fa_intra_rep_nonprotein, 0.545)
    return score_function

def get_match_substrate_pose_indexes(pdb):
    match_substrate_pose_indexes = set()
    info = pose.pdb_info()
    flag = False
    with open(pdb, 'r') as pdb:
        for line in pdb:
            if line.startswith('REMARK 666 MATCH TEMPLATE'):
                ligand_chain_id = line[26]
                ligand_res_id = int(line[32:36])
                ligand_pose_id = info.pdb2pose(ligand_chain_id, ligand_res_id)
                match_substrate_pose_indexes.add(ligand_pose_id)
                # if duplicated_chains:
                #     for duplicated_chain in filter(lambda x: x != ligand_chain_id, duplicated_chains):
                #         ligand_pose_id = info.pdb2pose(duplicated_chain, ligand_res_id)
                #         match_substrate_pose_indexes.add(ligand_pose_id)
                flag = True
            elif flag == True:
                break
    return match_substrate_pose_indexes

def create_coord_cst(pose, ligand_res_indexes):
    coord_cst_gen = CoordinateConstraintGenerator()
    if ligand_res_indexes:
        ligand_res_index_str = ','.join(str(ligand_res_index) for ligand_res_index in ligand_res_indexes)
        ligand_selector = ResidueIndexSelector(ligand_res_index_str)
        coord_cst_gen.set_residue_selector(NotResidueSelector(ligand_selector))
    return coord_cst_gen

def create_enzdes_cst():
# Add enzdes constraints
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

def create_relax_task_factory(point_mutations, nbh=None):
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    if point_mutations and point_mutations[0] != 'None':
        mutated_selector = OrResidueSelector()
        for point_mutation in point_mutations:
            mutation_info = point_mutation.split(',')
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(mutation_info[1])
            point_mutation_selector = ResidueIndexSelector(mutation_info[0])
            task_factory.push_back(OperateOnResidueSubset(restriction, point_mutation_selector))
            mutated_selector.add_residue_selector(point_mutation_selector)
        if nbh:
            repacking_selector = NeighborhoodResidueSelector()
            repacking_selector.set_focus_selector(mutated_selector)
            repacking_selector.set_distance(args.nbh)
            repacking_selector.set_include_focus_in_subset(False)
            repacking = RestrictToRepackingRLT()
            task_factory.push_back(OperateOnResidueSubset(repacking, repacking_selector))
            mutated_and_repacking_selector = OrResidueSelector(mutated_selector, repacking_selector)
            prevent = PreventRepackingRLT()
            task_factory.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_selector, True))
        else:
            repacking = RestrictToRepackingRLT()
            task_factory.push_back(OperateOnResidueSubset(repacking, mutated_selector, True))
    else:
        task_factory.push_back(RestrictToRepacking())
    return task_factory

def create_enzdes_task_factory(pdb, pose, nbh=None):
    # Mutatable
    info = pose.pdb_info()
    matching_res_indexes = set()
    flag = False
    with open(pdb, 'r') as pdb:
        for line in pdb:
            if line.startswith('REMARK 666 MATCH TEMPLATE'):
                # the linker is repackable
                # ligand_chain_id = line[26]
                ligand_res_id = int(line[32:36])
                ligand_pose_id = info.pdb2pose('A', ligand_res_id)
                matching_res_indexes.add(ligand_pose_id)
                # matching residues are repackable
                # motif_chain_id = line[49]
                motif_res_id = int(line[55:59])
                motif_pose_id = info.pdb2pose('A', motif_res_id)
                matching_res_indexes.add(motif_pose_id)
                flag = True
            elif flag == True:
                break
    matching_res_indexes_str = ','.join(str(idx) for idx in matching_res_indexes)
    matching_res_selector = ResidueIndexSelector(matching_res_indexes_str)
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.group1_selector(matching_res_selector)
    interface_selector.group2_selector(NotResidueSelector(matching_res_selector))
    mutatable_selector = AndResidueSelector(interface_selector, NotResidueSelector(matching_res_selector))
    # Repacking
    if nbh:
        repacking_selector = NeighborhoodResidueSelector()
        repacking_selector.set_focus_selector(mutatable_selector)
        repacking_selector.set_distance(nbh)
        repacking_selector.set_include_focus_in_subset(False)
    # Task Factory
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    restriction = RestrictAbsentCanonicalAASRLT()
    restriction.aas_to_keep('GASTVLIMPFYWDEHQNKR') # No cystine
    task_factory.push_back(OperateOnResidueSubset(restriction, mutatable_selector))
    repack = RestrictToRepackingRLT()
    if nbh:
        task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
        prevent = PreventRepackingRLT()
        task_factory.push_back(OperateOnResidueSubset(prevent, \
            OrResidueSelector(mutatable_selector, repacking_selector), True))
    else:
        task_factory.push_back(OperateOnResidueSubset(repack, mutatable_selector, True))
    return task_factory

def create_design_task_factory(mutatable_res_indexes, nbh=None):
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    mutatable_res_idx_str = ','.join(mutatable_res_indexes)
    mutatable_selector = ResidueIndexSelector(mutatable_res_idx_str)
    restriction = RestrictAbsentCanonicalAASRLT()
    restriction.aas_to_keep('GASTVLIMPFYWDEHQNKR') # no cystine
    task_factory.push_back(OperateOnResidueSubset(restriction, mutatable_selector))
    repack = RestrictToRepackingRLT()
    if nbh:
        repacking_selector = NeighborhoodResidueSelector()
        repacking_selector.set_focus_selector(mutatable_selector)
        repacking_selector.set_distance(nbh)
        repacking_selector.set_include_focus_in_subset(False)
        task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
        prevent = PreventRepackingRLT()
        task_factory.push_back(OperateOnResidueSubset(prevent, \
            OrResidueSelector(mutatable_selector, repacking_selector), True))
    else:
        task_factory.push_back(OperateOnResidueSubset(repack, mutatable_selector, True))
    return task_factory

def create_move_map(pose, ligand_res_indexes):
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    if ligand_res_indexes:
        for ligand_pose_id in ligand_res_indexes:
            edge = pose.fold_tree().get_residue_edge(ligand_pose_id)
            if not edge.is_jump():
                raise Exception('Edge of the linker is not a jump edge.')
            mm.set_jump(edge.label(), True)
    return mm

def create_fast_relax_mover(score_function, task_factory, move_map):
    # Make FastRelax mover
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)
    return fast_relax

def run_job_distributor(pose, name, decoys, sfxn, *movers):
    if decoys:
        job_distributor = PyJobDistributor(name, decoys, sfxn)
    else:
        job_distributor = PyJobDistributor(name, 5, sfxn)
    while not job_distributor.job_complete:
        pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        job_distributor.output_decoy(pose_copy)

if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    pose = pose_from_pdb(args.pdb)
    if args.fold_tree:
        fold_tree = create_fold_tree(args.fold_tree)
        pose.fold_tree(fold_tree)
    if args.chis:
        for chi in args.chis:
            chi_info = chi.split(',')
            pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
    # Applying symmetry if specified
    if args.symmetry:
        sfxn = SymmetricScoreFunction()
        sfxn.add_weights_from_file(args.score_function)
        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.set_keep_pdb_info_labels(True)
        sfsm.apply(pose)
    else:
        sfxn = create_score_fc(args.score_function)
    # get pose index of the linkers
    if args.enzdes:
        match_substrate_pose_indexes = get_match_substrate_pose_indexes(args.pdb)
    else:
        match_substrate_pose_indexes = None
    # coordinate constraint
    coord_cst_gen = create_coord_cst(pose, match_substrate_pose_indexes)
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    add_csts.apply(pose)
    # constraint
    if args.cst:
        add_fa_constraints_from_cmdline(pose, sfxn)
    # enzyme design constraint
    if args.enzdescst:
        enzdes_cst = create_enzdes_cst()
        enzdes_cst.apply(pose)
    # declare movers
    movers = list()
    if args.rmsd:
        if args.rmsd == 'True':
            rmsdm = RMSDMetric(pose)
        else:
            rmsdm = RMSDMetric(pose, ResidueIndexSelector(args.rmsd))
        movers.append(rmsdm)
    if args.des or args.enzdes:
        if args.nataa:
            favor_nataa = FavorNativeResidue(pose, 1.5)
        if args.des:
            tf = create_design_task_factory(args.des, nbh=args.nbh)
        elif args.enzdes:
            tf = create_enzdes_task_factory(args.pdb, pose, nbh=args.nbh)
    else:
        tf = create_relax_task_factory(args.muts, nbh=args.nbh)
    mm = create_move_map(pose, match_substrate_pose_indexes)
    fr = create_fast_relax_mover(sfxn, tf, mm)
    movers.append(fr)
    name = args.pdb.split('/')[-1][:-4]
    run_job_distributor(pose, name, args.decoys, sfxn, *movers)
