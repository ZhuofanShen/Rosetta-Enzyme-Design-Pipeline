#!/usr/bin/python3
import argparse
import os
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    InitializeFromCommandline, IncludeCurrent, ExtraRotamers, \
    OperateOnResidueSubset, RestrictToRepackingRLT, \
    RestrictAbsentCanonicalAASRLT, PreventRepackingRLT, \
    RestrictToRepacking, PreventRepacking
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
    AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.protein_interface_design \
    import FavorNativeResidue
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('-ref', '--reference_pose', type=str)
    parser.add_argument('-params', '--parameters_files', type=str, nargs='*')
    parser.add_argument('-ft', '--fold_tree', type=str, nargs='*')
    parser.add_argument('-chis', '--chi_dihedrals', type=str, nargs='*', default=list())
    parser.add_argument('-sf', '--score_function', type=str, default='ref2015_cst')
    parser.add_argument('--score_terms', type=str, nargs='*', default=list())
    parser.add_argument('-symm', '--symmetry', type=str)
    parser.add_argument('-cst', '--constraints', type=str)
    parser.add_argument('-nbh', '--neighborhood', type=float, default=0)
    parser.add_argument('-relax', '--fast_relax', action='store_true')
    parser.add_argument('-muts', '--mutations', type=str, nargs='*', default=list())
    parser.add_argument('-nataa', '--favor_native_residue', type=float)
    parser.add_argument('-des', '--designable_sites', type=int, nargs='*')
    parser.add_argument('-enzdes_cst', '--enzyme_design_constraints', type=str)
    parser.add_argument('-subs', '--substrates', type=int, nargs='*')
    parser.add_argument('-no_cys', '--cystine', action='store_false')
    parser.add_argument('-ncaa', '--noncanonical_amino_acids', type=str, nargs='*', default=list(), help='name3')
    parser.add_argument('-xform', '--substrate_rigid_body_transformations', type=bool, default=False)
    parser.add_argument('-n', '--decoys', type=int, default=50)
    parser.add_argument('-rmsd', type=int, nargs='*', default=None)
    parser.add_argument('--annotated_name', action='store_true')
    parser.add_argument('-debug', '--debug_mode', action='store_true')
    args = parser.parse_args()
    return args

def init_pyrosetta_with_opts(args):
    opts = '-ex1 -ex2 -ignore_zero_occupancy false -use_input_sc -no_optH false -flip_HNQ'
    if args.score_function.startswith('beta_nov16'):
        opts += ' -corrections::beta_nov16'
    if args.parameters_files:
        opts += " -extra_res_fa " + ' '.join(args.parameters_files)
    if args.enzyme_design_constraints:
        opts += ' -enzdes:cstfile {} -run:preserve_header'.format(args.enzyme_design_constraints)
    if args.constraints:
        opts += ' -constraints:cst_fa_file {}'.format(args.constraints)
    init(opts)

def create_fold_tree(edge_list):
    fold_tree = FoldTree()
    for edge_str in edge_list:
        edge_num = edge_str.split(',')
        if len(edge_num) == 3:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), int(edge_num[2]))
        elif len(edge_num) == 4:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), edge_num[2], edge_num[3])
        else:
            raise Exception("The number of arguments in add_edge function should be 3 or 4")
    return fold_tree

def get_match_pose_indexes(info, pdb, symmetry):
    match_substrate_pose_indexes = set()
    match_res_pose_indexes = set()
    flag = False
    main_chain = None
    with open(pdb, 'r') as pdb:
        for line in pdb:
            if line.startswith('REMARK 666 MATCH TEMPLATE'):
                # substrate
                substrate_chain_id = line[26]
                if not main_chain:
                    main_chain = substrate_chain_id
                if main_chain == substrate_chain_id or not symmetry:
                    substrate_res_id = int(line[32:36])
                    substrate_pose_id = info.pdb2pose(substrate_chain_id, substrate_res_id)
                    match_substrate_pose_indexes.add(substrate_pose_id)
                # match residues
                motif_chain_id = line[49]
                if main_chain == motif_chain_id or not symmetry:
                    motif_res_id = int(line[55:59])
                    motif_pose_id = info.pdb2pose(motif_chain_id, motif_res_id)
                    match_res_pose_indexes.add(motif_pose_id)
                flag = True
            elif flag == True:
                break
    return match_substrate_pose_indexes, match_res_pose_indexes

def create_coord_cst(ref_pose=None, no_coord_cst=None):
    coord_cst_gen = CoordinateConstraintGenerator()
    if ref_pose:
        coord_cst_gen.set_reference_pose(ref_pose)
    if no_coord_cst:
        no_coord_cst_selector = ResidueIndexSelector(','.join(str(no_coord_cst_res) for no_coord_cst_res in no_coord_cst))
        coord_cst_gen.set_residue_selector(NotResidueSelector(no_coord_cst_selector))
    return coord_cst_gen

def create_enzdes_cst():
# Add enzdes constraints
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

def create_relax_task_factory(point_mutations: list = list(), active_site_positions: list = list(), neighborhood: float = 0):
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    repack = RestrictToRepackingRLT()
    prevent = PreventRepackingRLT()
    if len(point_mutations) > 0:
        # Mutate
        mutated_selector = OrResidueSelector()
        for point_mutation in point_mutations:
            mutation_info = point_mutation.split(',')
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(mutation_info[1])
            point_mutation_selector = ResidueIndexSelector(mutation_info[0])
            task_factory.push_back(OperateOnResidueSubset(restriction, point_mutation_selector))
            mutated_selector.add_residue_selector(point_mutation_selector)
        # Repack and static
        if neighborhood > 0:
            if len(active_site_positions) > 0:
                # Repack
                enzyme_core_selector = ResidueIndexSelector(','.join(str(active_site_position) for active_site_position in active_site_positions))
                designable_repacking_selector = NeighborhoodResidueSelector()
                designable_repacking_selector.set_focus_selector(OrResidueSelector(mutated_selector, enzyme_core_selector))
                designable_repacking_selector.set_distance(neighborhood)
                designable_repacking_selector.set_include_focus_in_subset(True)
                repacking_selector = AndResidueSelector(designable_repacking_selector, NotResidueSelector(mutated_selector))
                task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
                # Static
                task_factory.push_back(OperateOnResidueSubset(prevent, designable_repacking_selector, True))
            else:
                # Repack
                repacking_selector = NeighborhoodResidueSelector()
                repacking_selector.set_focus_selector(mutated_selector)
                repacking_selector.set_distance(args.neighborhood)
                repacking_selector.set_include_focus_in_subset(False)
                task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
                # Static
                mutated_and_repacking_selector = OrResidueSelector(mutated_selector, repacking_selector)
                task_factory.push_back(OperateOnResidueSubset(prevent, mutated_and_repacking_selector, True))
        else:
            # Repack
            task_factory.push_back(OperateOnResidueSubset(repack, mutated_selector, True))
    else:
        if neighborhood > 0:
            if len(active_site_positions) > 0:
                # Repack
                enzyme_core_selector = ResidueIndexSelector(','.join(str(active_site_position) for active_site_position in active_site_positions))
                repacking_selector = NeighborhoodResidueSelector()
                repacking_selector.set_focus_selector(enzyme_core_selector)
                repacking_selector.set_distance(neighborhood)
                repacking_selector.set_include_focus_in_subset(True)
                task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
                # Static
                task_factory.push_back(OperateOnResidueSubset(prevent, repacking_selector, True))
            else:
                # Static
                task_factory.push_back(PreventRepacking())
        else:
            # Repack
            task_factory.push_back(RestrictToRepacking())
    return task_factory

def create_design_task_factory(designable_res_indexes: list, active_site_positions: list = list(), \
        neighborhood: float = 0, cystine: bool = True, noncanonical_amino_acids: list = list()):
    task_factory = TaskFactory()
    if len(noncanonical_amino_acids) > 0:
        ncaa_palette = CustomBaseTypePackerPalette()
        for ncaa in noncanonical_amino_acids:
            ncaa_palette.add_type(ncaa)
        task_factory.set_packer_palette(ncaa_palette)
    task_factory.push_back(IncludeCurrent())
    # Mutatable
    designable_selector = ResidueIndexSelector(','.join(str(designable_res_index) for designable_res_index in designable_res_indexes))
    if not cystine:
        restriction = RestrictAbsentCanonicalAASRLT()
        restriction.aas_to_keep('AGILPVFWYDERHKSTCMNQ')
        task_factory.push_back(OperateOnResidueSubset(restriction, designable_selector))
    # Repack
    repack = RestrictToRepackingRLT()
    prevent = PreventRepackingRLT()
    if neighborhood > 0:
        if len(active_site_positions) > 0:
            enzdes_cst_selector = ResidueIndexSelector(','.join(str(enzdes_cst_index) for enzdes_cst_index in active_site_positions))
            designable_repacking_selector = NeighborhoodResidueSelector()
            designable_repacking_selector.set_focus_selector(OrResidueSelector(designable_selector, enzdes_cst_selector))
            designable_repacking_selector.set_distance(neighborhood)
            designable_repacking_selector.set_include_focus_in_subset(True)
            repacking_selector = AndResidueSelector(designable_repacking_selector, NotResidueSelector(designable_selector))
            task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
            task_factory.push_back(OperateOnResidueSubset(prevent, designable_repacking_selector, True))
        else:
            repacking_selector = NeighborhoodResidueSelector()
            repacking_selector.set_focus_selector(designable_selector)
            repacking_selector.set_distance(neighborhood)
            repacking_selector.set_include_focus_in_subset(False)
            task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
            task_factory.push_back(OperateOnResidueSubset(prevent, \
                OrResidueSelector(designable_selector, repacking_selector), True))
    else:
        task_factory.push_back(OperateOnResidueSubset(repack, designable_selector, True))
    return task_factory

def create_enzyme_design_task_factory(active_site_positions: list, neighborhood: float = 0, \
        cystine: bool = True, noncanonical_amino_acids: list = list()):
    task_factory = TaskFactory()
    if len(noncanonical_amino_acids) > 0:
        ncaa_palette = CustomBaseTypePackerPalette()
        for ncaa in noncanonical_amino_acids:
            ncaa_palette.add_type(ncaa)
        task_factory.set_packer_palette(ncaa_palette)
    task_factory.push_back(IncludeCurrent())
    # Mutatable
    enzyme_core_selector = ResidueIndexSelector(','.join(str(active_site_position) for active_site_position in active_site_positions))
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.group1_selector(enzyme_core_selector)
    interface_selector.group2_selector(NotResidueSelector(enzyme_core_selector))
    designable_selector = AndResidueSelector(interface_selector, NotResidueSelector(enzyme_core_selector))
    if not cystine:
        restriction = RestrictAbsentCanonicalAASRLT()
        restriction.aas_to_keep('AGILPVFWYDERHKSTCMNQ')
        task_factory.push_back(OperateOnResidueSubset(restriction, designable_selector))
    # Repack
    repack = RestrictToRepackingRLT()
    if neighborhood > 0:
        repacking_wo_enzyme_core_selector = NeighborhoodResidueSelector()
        repacking_wo_enzyme_core_selector.set_focus_selector(interface_selector)
        repacking_wo_enzyme_core_selector.set_distance(neighborhood)
        repacking_wo_enzyme_core_selector.set_include_focus_in_subset(False)
        repacking_selector = OrResidueSelector(repacking_wo_enzyme_core_selector, enzyme_core_selector)
        task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
        prevent = PreventRepackingRLT()
        task_factory.push_back(OperateOnResidueSubset(prevent, OrResidueSelector(designable_selector, repacking_selector), True))
    else:
        task_factory.push_back(OperateOnResidueSubset(repack, designable_selector, True))
    return task_factory

def create_move_map(pose, ligand_res_indexes: list = list()):
    mm = MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    if len(ligand_res_indexes) > 0:
        for ligand_pose_id in ligand_res_indexes:
            edge = pose.fold_tree().get_residue_edge(ligand_pose_id)
            if not edge.is_jump():
                raise Exception('Edge of the ligand is not a jump edge.')
            mm.set_jump(edge.label(), True)
    return mm

def create_fast_relax_mover(score_function, task_factory, move_map=None):
    # Make FastRelax mover
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    if move_map:
        fast_relax.set_movemap(move_map)
    return fast_relax

def run_job_distributor(pose, sfxn, decoys, name, annotated_name, *movers):
    if decoys:
        job_distributor = PyJobDistributor(name, decoys, sfxn)
    else:
        job_distributor = PyJobDistributor(name, 5, sfxn)
    pdb_name = job_distributor.pdb_name
    while not job_distributor.job_complete:
        current_id = job_distributor.current_id
        pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        job_distributor.output_decoy(pose_copy)
        if annotated_name:
            mutations_str = str()
            for index, wt_res in pose.sequence():
                res = pose_copy.sequence()[index]
                if res != wt_res:
                    mutations_str += '_' + wt_res + str(index) + res
            os.rename(pdb_name + '_' + str(current_id) + '.pdb', pdb_name + '_' + str(current_id) + mutations_str + '.pdb')

if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    pose = pose_from_pdb(args.pdb)
    if args.reference_pose:
        ref_pose = pose_from_pdb(args.reference_pose)
    if args.fold_tree:
        fold_tree = create_fold_tree(args.fold_tree)
        pose.fold_tree(fold_tree)
        if args.reference_pose:
            ref_pose.fold_tree(fold_tree)
    for chi in args.chi_dihedrals:
        chi_info = chi.split(',')
        pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
        if args.reference_pose:
            ref_pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
    # Applying symmetry if specified
    if args.symmetry:
        sfxn = SymmetricScoreFunction()
        sfxn.add_weights_from_file(args.score_function)
        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.set_keep_pdb_info_labels(True)
        sfsm.apply(pose)
        if args.reference_pose:
            sfsm.apply(ref_pose)
    else:
        sfxn = create_score_function(args.score_function)
    # Add score terms
    for score_term in args.score_terms:
        term_weight = score_term.split(':')
        exec("sfxn.set_weight(ScoreType.{}, {})".format(term_weight[0], term_weight[1]))
    # get pose index of the substrate and catalytic residues
    active_site_positions = set()
    if args.enzyme_design_constraints:
        pdb_info = pose.pdb_info()
        match_substrate_pose_indexes, match_res_pose_indexes = get_match_pose_indexes(pdb_info, args.pdb, args.symmetry)
        match_indexes = match_substrate_pose_indexes.union(match_res_pose_indexes)
        active_site_positions.update(match_index for match_index in match_indexes)
    if args.substrates:
        active_site_positions.update(args.substrates)
    # coordinate constraint
    if args.reference_pose:
        coord_cst_gen = create_coord_cst(ref_pose = ref_pose)
    else:
        coord_cst_gen = create_coord_cst()
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    add_csts.apply(pose)
    # constraint
    if args.constraints:
        add_fa_constraints_from_cmdline(pose, sfxn)
    # enzyme design constraint
    if args.enzyme_design_constraints:
        enzdes_cst = create_enzdes_cst()
        enzdes_cst.apply(pose)
    # declare movers
    movers = list()
    if not args.rmsd:
        rmsdm = RMSDMetric(pose)
    else:
        rmsdm = RMSDMetric(pose, ResidueIndexSelector(','.join(str(res) for res in args.rmsd)))
    movers.append(rmsdm)
    if args.fast_relax:
        tf = create_relax_task_factory(point_mutations = args.mutations, active_site_positions = active_site_positions, \
                neighborhood = args.neighborhood)
    elif args.designable_sites or args.substrates or args.enzyme_design_constraints:
        if args.favor_native_residue:
            favor_nataa = FavorNativeResidue(pose, args.favor_native_residue)
        if args.designable_sites:
            tf = create_design_task_factory(args.designable_sites, active_site_positions = active_site_positions, \
                    neighborhood = args.neighborhood, cystine = args.cystine, noncanonical_amino_acids = args.noncanonical_amino_acids)
        else:
            tf = create_enzyme_design_task_factory(active_site_positions, neighborhood = args.neighborhood, \
                    cystine = args.cystine, noncanonical_amino_acids = args.noncanonical_amino_acids)
    else:
        tf = create_relax_task_factory(point_mutations = args.mutations, neighborhood = args.neighborhood)
    # If the substrate is allowed to perform rigid body transformations (under development)
    if args.substrate_rigid_body_transformations:
        if args.substrates:
            mm = create_move_map(pose, ligand_res_indexes = args.substrates)
        elif args.enzyme_design_constraints:
            mm = create_move_map(pose, ligand_res_indexes = match_substrate_pose_indexes)
        fr = create_fast_relax_mover(sfxn, tf, move_map = mm)
    else:
        fr = create_fast_relax_mover(sfxn, tf)
    movers.append(fr)
    name = args.pdb.split('/')[-1][:-4]
    if args.debug_mode:
        print(pose.fold_tree())
        if args.enzyme_design_constraints:
            print(match_substrate_pose_indexes)
            print(match_res_pose_indexes)
        print(tf.create_task_and_apply_taskoperations(pose))
        print(mm)
    else:
        run_job_distributor(pose, sfxn, args.decoys, name, args.annotated_name, *movers)
