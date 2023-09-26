#!/usr/bin/python3
import argparse
import math
import os
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    InitializeFromCommandline, IncludeCurrent, ExtraRotamers, \
    OperateOnResidueSubset, RestrictToRepackingRLT, \
    RestrictAbsentCanonicalAASRLT, PreventRepackingRLT, \
    RestrictToRepacking, PreventRepacking
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import \
    add_fa_constraints_from_cmdline, ConstraintSet, AmbiguousConstraint, \
    AtomPairConstraint, AngleConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import CircularHarmonicFunc, \
    FlatHarmonicFunc, HarmonicFunc
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueIndexSelector, InterGroupInterfaceByVectorSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    RMSDMetric, TotalEnergyMetric
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
    parser.add_argument("pdb", type=str)
    parser.add_argument("-ref", "--coordinate_reference_pdb", type=str)
    parser.add_argument("-params", "--parameters_files", type=str, nargs="*")
    parser.add_argument("-ft", "--fold_tree", type=str, nargs="*")
    parser.add_argument("-chis", "--chi_dihedrals", type=str, nargs="*", default=list())
    parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst")
    parser.add_argument("--score_terms", type=str, nargs="*", default=list())
    parser.add_argument("-symm", "--symmetry", type=str)
    parser.add_argument("-coord_boundary", "--coordinate_constraint_bounded_width", type=float)
    parser.add_argument("-no_coord_cst", "--no_coordinate_constraint_residues", type=str, nargs="*", default=list())
    parser.add_argument("-cst", "--constraints", type=str)
    parser.add_argument("-muts", "--mutations", type=str, nargs="*", default=list())
    parser.add_argument("-nataa", "--favor_native_residue", type=float)
    parser.add_argument("-des", "--design_positions", type=str, nargs="*", default=list())
    parser.add_argument("-no_opt_coord_cst", "--no_optimization_region_coordinate_constraint", action="store_true")
    parser.add_argument("-enzdes_cst", "--enzyme_design_constraints", type=str)
    parser.add_argument("--design_active_site", action="store_true")
    parser.add_argument("-subs", "--substrates", type=str, nargs="*")
    parser.add_argument("-rpk_intf", "--repack_interface_only", action="store_true")
    parser.add_argument("-min_intf", "--minimize_interface_only", action="store_true")
    parser.add_argument("-ddG_ref", "--ddG_reference_pdb", type=str)
    parser.add_argument("-xform", "--substrate_rigid_body_transformations", action="store_true")
    parser.add_argument("-premin", "--pre_minimization", action="store_true")
    parser.add_argument("-no_cys", "--no_cystine", action="store_true")
    parser.add_argument("-ncaa", "--noncanonical_amino_acids", type=str, nargs="*", default=list(), help="name3")
    parser.add_argument("-optH", "--optimize_protonation_state", action="store_true")
    parser.add_argument("-no_rmsd", "--no_rmsd_residues", type=str, nargs="*", default=list())
    parser.add_argument("-n", "--decoys", type=int, default=50)
    parser.add_argument("--annotated_name", action="store_true")
    parser.add_argument("-debug", "--debug_mode", action="store_true")
    args = parser.parse_args()
    return args

def init_pyrosetta_with_opts(args):
    if args.optimize_protonation_state:
        no_optH = "false -flip_HNQ"
    else:
        no_optH = "true"
    opts = "-ex1 -ex2 -ignore_zero_occupancy false -use_input_sc -no_optH " + no_optH
    if args.score_function.startswith("beta_nov16"):
        opts += " -corrections::beta_nov16"
    if args.parameters_files:
        opts += " -extra_res_fa " + " ".join(args.parameters_files)
    if args.enzyme_design_constraints:
        opts += " -enzdes:cstfile {} -run:preserve_header".format(args.enzyme_design_constraints)
    if args.constraints:
        opts += " -constraints:cst_fa_file {}".format(args.constraints)
    init(opts)

def set_symmetry(symmetry, *pose_list):
    if symmetry:
        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.set_keep_pdb_info_labels(True)
        for pose in pose_list:
            if pose:
                sfsm.apply(pose)

def set_score_function(score_function_weight_file, symmetry=False, score_terms=list()):
    if symmetry:
        score_function = SymmetricScoreFunction()
        score_function.add_weights_from_file(score_function_weight_file)
    else:
        score_function = create_score_function(score_function_weight_file)
    # Add score terms
    for score_term in score_terms:
        term_weight = score_term.split(":")
        exec("score_function.set_weight(ScoreType.{}, {})".format(term_weight[0], term_weight[1]))
    return score_function

def customize_score_function(rep_type="hard", cartesian=False, symmetry=False, membrane=False, constraint_weight=1.0, score_terms=list()):
    """
    Custimize an appropriate score function on top of a L-J potential 
    that is either hard (ref2015) or soft (ref2015_soft), and having  
    symmetry and/or membrane and/or constraints settings.
    """
    assert rep_type in ["hard", "soft"]

    if symmetry: # Declare the symmetry score function
        score_function = SymmetricScoreFunction()
        if rep_type == "hard":
            if membrane:
                score_function.add_weights_from_file("franklin2019")
            else:
                score_function.add_weights_from_file("ref2015")
        elif rep_type == "soft":
            if membrane: # Set up a soft-rep version of franklin2019
                score_function.add_weights_from_file("ref2015_soft")
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function.add_weights_from_file("ref2015_soft")
    else: # Declare the ordinary score function
        if rep_type == "hard":
            if membrane:
                score_function = create_score_function("franklin2019")
            else:
                score_function = create_score_function("ref2015")
        elif rep_type == "soft":
            if membrane: # Set up a soft-rep version of franklin2019
                score_function = create_score_function("ref2015_soft")
                score_function.set_weight(ScoreType.fa_water_to_bilayer, 1.0)
            else:
                score_function = create_score_function("ref2015_soft")
    # The score functions do not have constraint weights.
    # If requisited, the constraint weights are added.
    if constraint_weight is not None:
        score_function.set_weight(ScoreType.atom_pair_constraint, constraint_weight)
        score_function.set_weight(ScoreType.coordinate_constraint, constraint_weight)
        score_function.set_weight(ScoreType.angle_constraint, constraint_weight)
        score_function.set_weight(ScoreType.dihedral_constraint, constraint_weight)
        score_function.set_weight(ScoreType.metalbinding_constraint, constraint_weight)
        score_function.set_weight(ScoreType.chainbreak, constraint_weight)
        score_function.set_weight(ScoreType.res_type_constraint, constraint_weight)
    # Set the Cartesian score term.
    if cartesian:
        score_function.set_weight(ScoreType.cart_bonded, 0.5)
        score_function.set_weight(ScoreType.pro_close, 0)
        score_function.set_weight(ScoreType.metalbinding_constraint, 0)
        score_function.set_weight(ScoreType.chainbreak, 0)
    # Add score terms
    for score_term in score_terms:
        term_weight = score_term.split(":")
        exec("score_function.set_weight(ScoreType.{}, {})".format(term_weight[0], term_weight[1]))

    return score_function

def create_fold_tree(edge_list):
    fold_tree = FoldTree()
    for edge_str in edge_list:
        edge_num = edge_str.split(",")
        if len(edge_num) == 3:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), int(edge_num[2]))
        elif len(edge_num) == 4:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), edge_num[2], edge_num[3])
        else:
            raise Exception("The number of arguments in add_edge function should be 3 or 4")
    return fold_tree

def get_match_pose_indices(info, pdb, symmetry):
    match_substrate_pose_indices = set()
    match_res_pose_indices = set()
    flag = False
    main_chain = None
    with open(pdb, "r") as pdb:
        for line in pdb:
            if line.startswith("REMARK 666 MATCH TEMPLATE"):
                # substrate
                substrate_chain_id = line[26]
                if not main_chain:
                    main_chain = substrate_chain_id
                if main_chain == substrate_chain_id or not symmetry:
                    substrate_res_id = int(line[32:36])
                    substrate_pose_id = info.pdb2pose(substrate_chain_id, substrate_res_id)
                    match_substrate_pose_indices.add(substrate_pose_id)
                # match residues
                motif_chain_id = line[49]
                if main_chain == motif_chain_id or not symmetry:
                    motif_res_id = int(line[55:59])
                    motif_pose_id = info.pdb2pose(motif_chain_id, motif_res_id)
                    match_res_pose_indices.add(motif_pose_id)
                flag = True
            elif flag == True:
                break
    return match_substrate_pose_indices, match_res_pose_indices

def create_coord_cst(coord_ref_pose=None, coord_boundary = None, no_coord_cst_selection = False, sidechain=False):
    coord_cst_gen = CoordinateConstraintGenerator()
    if coord_boundary:
        coord_cst_gen.set_bounded(True)
        coord_cst_gen.set_bounded_width(coord_boundary)
    if coord_ref_pose:
        coord_cst_gen.set_reference_pose(coord_ref_pose)
    if no_coord_cst_selection is not None:
        if no_coord_cst_selection is not False:
            coord_cst_gen.set_residue_selector(NotResidueSelector(no_coord_cst_selection))
    else:
        coord_cst_gen.set_residue_selector(OrResidueSelector())
    if sidechain:
        coord_cst_gen.set_sidechain(True)
    return coord_cst_gen

def apply_dihedral_constraint(pose, atom_name_list, residue_id_list, degree, standard_deviation):
    atom_list = list()
    for atom_name, residue_id in zip(atom_name_list, residue_id_list):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(degree) == list:
        circular_harmonic_cst = AmbiguousConstraint()
        for deg, sd in zip(degree, standard_deviation):
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * deg, sd)
            cst = DihedralConstraint(atom_list[0], atom_list[1], atom_list[2], atom_list[3], circular_harmonic_fc)
            circular_harmonic_cst.add_individual_constraint(cst)
    else:
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * degree, standard_deviation)
        circular_harmonic_cst = DihedralConstraint(atom_list[0], atom_list[1], atom_list[2], atom_list[3], circular_harmonic_fc)
    pose.add_constraint(circular_harmonic_cst)

def create_enzdes_cst():
# Add enzdes constraints
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

def pre_minimization(pose, substrate_pose_indices):
    move_map = MoveMap()
    move_map.set_bb(False)
    substrate_catalytic_sidechain_selector = ResidueIndexSelector(",".join(str(substrate_pose_index) for substrate_pose_index in substrate_pose_indices))
    move_map.set_chi(substrate_catalytic_sidechain_selector.apply(pose))
    pre_minimizer = MinMover()
    pre_minimizer.min_type("lbfgs_armijo_nonmonotone")
    pre_minimizer.movemap(move_map)
    return pre_minimizer

def pdb_to_pose_numbering(pose, chain_id_pdb_indices):
    pose_indices = set()
    for chain_id_pdb_index in chain_id_pdb_indices:
        if "," in chain_id_pdb_index:
            chain_id_pdb_index, final_aa_type = chain_id_pdb_index.split(",")
            chain_id, pdb_index = chain_id_pdb_index[0], int(chain_id_pdb_index[1:])
            pose_index = pose.pdb_info().pdb2pose(chain_id, pdb_index)
            if pose_index > 0:
                pose_indices.add(str(pose_index) + "," + final_aa_type)
        else:
            if "-" in chain_id_pdb_index:
                chain_id = chain_id_pdb_index[0]
                start_pdb_index, end_pdb_index = chain_id_pdb_index[1:].split("-")
                pdb_indices = range(int(start_pdb_index), int(end_pdb_index) + 1)
            else:
                chain_id, pdb_index = chain_id_pdb_index[0], int(chain_id_pdb_index[1:])
                pdb_indices = [pdb_index]
            for pdb_index in pdb_indices:
                pose_index = pose.pdb_info().pdb2pose(chain_id, pdb_index)
                if pose_index > 0:
                    pose_indices.add(pose_index)
    return pose_indices

def create_task_factory(point_mutations: set = set(), design_positions: set = set(), \
        theozyme_positions: set = set(), design_active_site: bool = False, repacking_range: bool = False, \
        ddG_ref_pose = None, no_cystine: bool = False, noncanonical_amino_acids: list = None):
    """
    Priority: point mutations > specified design positions > protein-ligand interface design
    """
    # Create task factory
    task_factory = TaskFactory()
    if noncanonical_amino_acids:
        ncaa_palette = CustomBaseTypePackerPalette()
        for ncaa in noncanonical_amino_acids:
            ncaa_palette.add_type(ncaa)
        task_factory.set_packer_palette(ncaa_palette)
    
    # Every positions are designable by defalut other than specification.
    task_factory.push_back(IncludeCurrent())

    # Specify point mutations
    if len(point_mutations) > 0:
        mutated_position_list = set()
        for point_mutation in point_mutations:
            position_aa = point_mutation.split(",")
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(position_aa[1])
            task_factory.push_back(OperateOnResidueSubset(restriction, ResidueIndexSelector(position_aa[0])))
            mutated_position_list.add(position_aa[0])
        mutation_selection = ResidueIndexSelector(",".join(mutated_position_list))
    else:
        mutation_selection = OrResidueSelector()

    # Specify design positions
    design_selection = OrResidueSelector()
    # Site directed design
    if len(design_positions) > 0:
        site_spec_design_selection = ResidueIndexSelector(",".join(str(design_position) for design_position in design_positions))
        if len(point_mutations) > 0:
            site_spec_design_selection = AndResidueSelector(site_spec_design_selection, NotResidueSelector(mutation_selection))
        design_selection.add_residue_selector(site_spec_design_selection)
    # Theozyme protein interface design
    if len(theozyme_positions) > 0:
        theozyme_selection = ResidueIndexSelector(",".join(str(theozyme_position) for theozyme_position in theozyme_positions))
        # Find protein-theozyme interface design positions
        if design_active_site:
            theozyme_protein_interface_selection = InterGroupInterfaceByVectorSelector()
            # theozyme_protein_interface_selection.nearby_atom_cut(7.0) # 5.5
            # theozyme_protein_interface_selection.cb_dist_cut(12.0) # 11.0
            # theozyme_protein_interface_selection.vector_dist_cut(10.5) # 9.0
            # theozyme_protein_interface_selection.vector_angle_cut(76.0) # 75.0
            theozyme_protein_interface_selection.group1_selector(theozyme_selection)
            theozyme_protein_interface_selection.group2_selector(NotResidueSelector(theozyme_selection))
            protein_interface_selection = AndResidueSelector(theozyme_protein_interface_selection, \
                    NotResidueSelector(theozyme_selection))
            if len(point_mutations) > 0:
                protein_interface_selection = AndResidueSelector(protein_interface_selection, NotResidueSelector(mutation_selection))
            design_selection.add_residue_selector(protein_interface_selection)
    else:
        theozyme_selection = OrResidueSelector()

    if no_cystine:
        restriction = RestrictAbsentCanonicalAASRLT()
        restriction.aas_to_keep("AGILPVFWYDERHKSTMNQ")
        task_factory.push_back(OperateOnResidueSubset(restriction, design_selection))

    # Repack and static
    if len(point_mutations) > 0 or len(design_positions) > 0 or len(theozyme_positions) > 0:
        substitution_selection = OrResidueSelector(mutation_selection, design_selection)
        substitution_theozyme_selection = OrResidueSelector(substitution_selection, theozyme_selection)
        if repacking_range:
            # Repack
            substitution_repacking_selection = InterGroupInterfaceByVectorSelector()
            # substitution_repacking_selection.nearby_atom_cut(7.0) # 5.5
            # substitution_repacking_selection.cb_dist_cut(12.0) # 11.0
            # substitution_repacking_selection.vector_dist_cut(10.5) # 9.0
            # substitution_repacking_selection.vector_angle_cut(76.0) # 75.0
            substitution_repacking_selection.group1_selector(substitution_theozyme_selection)
            substitution_repacking_selection.group2_selector(NotResidueSelector(substitution_theozyme_selection))
            substitution_repacking_selection = OrResidueSelector(substitution_repacking_selection, \
                    theozyme_selection) # ensure to include theozymes
            repacking_selection = AndResidueSelector(substitution_repacking_selection, \
                    NotResidueSelector(substitution_selection)) # repack w/o mutations
            if ddG_ref_pose:
                repacking_vector = repacking_selection.apply(ddG_ref_pose)
                repacking_selection = ResidueIndexSelector(",".join(filter(lambda x: x is not None, \
                        map(lambda x, y: str(x) if y == 1 else None, range(1, len(repacking_vector) + 1), repacking_vector))))
            substitution_repacking_selection = OrResidueSelector(repacking_selection, \
                    substitution_selection) # ensure to include mutations
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), repacking_selection))
            # Static
            task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), OrResidueSelector(substitution_selection, \
                    repacking_selection), True))
        else:
            # Repack
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), substitution_selection, True))
            substitution_repacking_selection = None # all residue
    else:
        if repacking_range:
            # Static
            task_factory.push_back(PreventRepacking())
            substitution_repacking_selection = False # no residue
        else:
            # Repack
            task_factory.push_back(RestrictToRepacking())
            substitution_repacking_selection = None # all residue

    return task_factory, substitution_repacking_selection

def create_move_map(pose, substitution_repacking_selection = None, ddG_ref_pose = None, ligand_res_indices: set = set()):
    move_map = MoveMap()
    if substitution_repacking_selection:
        interface_selector = InterGroupInterfaceByVectorSelector()
        interface_selector.nearby_atom_cut(10.0) # 5.5
        interface_selector.cb_dist_cut(15.0) # 11.0
        interface_selector.vector_dist_cut(12.0) # 9.0
        interface_selector.vector_angle_cut(80.0) # 75.0
        interface_selector.group1_selector(substitution_repacking_selection)
        interface_selector.group2_selector(NotResidueSelector(substitution_repacking_selection))
        minimization_selection = OrResidueSelector(substitution_repacking_selection, interface_selector)
        if ddG_ref_pose: # passing residue indices vectors to a movemap
            minimization_vector = minimization_selection.apply(ddG_ref_pose)
            move_map.set_bb(minimization_vector)
            move_map.set_chi(minimization_vector)
        else: # passing residue selectors to a movemap factory
            minimization_vector = minimization_selection.apply(pose)
            move_map.set_bb(minimization_vector)
            move_map.set_chi(minimization_vector)
    elif substitution_repacking_selection is None:
        move_map.set_bb(True)
        move_map.set_chi(True)
    elif substitution_repacking_selection is False:
        move_map.set_bb(False)
        move_map.set_chi(False)
    if len(ligand_res_indices) > 0:
        for ligand_pose_id in ligand_res_indices:
            edge = pose.fold_tree().get_residue_edge(ligand_pose_id)
            if not edge.is_jump():
                raise Exception("Edge of the ligand is not a jump edge.")
            move_map.set_jump(edge.label(), True)
    return move_map

def create_fast_relax_mover(score_function, task_factory, move_map=None):
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    if move_map:
        fast_relax.set_movemap(move_map)
    return fast_relax

def calculate_energy(score_function, pose, selector=None, score_type=None):
    # Create the metric
    metric = TotalEnergyMetric()
    metric.set_scorefunction(score_function)
    if score_type:
        exec("metric.set_scoretype(ScoreType.{})".format(score_type))
    # Add the selector
    if selector:
        metric.set_residue_selector(selector)
    return metric.calculate(pose)

def run_job_distributor(pose, score_function, decoys, name, annotated_name, *movers):
    if decoys:
        job_distributor = PyJobDistributor(name, decoys, score_function)
    else:
        job_distributor = PyJobDistributor(name, 5, score_function)
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
                    mutations_str += "_" + wt_res + str(index) + res
            os.rename(pdb_name + "_" + str(current_id) + ".pdb", pdb_name + "_" + str(current_id) + mutations_str + ".pdb")

def main(args):
    pose = pose_from_pdb(args.pdb)
    coord_ref_pose = None
    if args.coordinate_reference_pdb:
        coord_ref_pose = pose_from_pdb(args.coordinate_reference_pdb)
    ddG_ref_pose = None
    if args.ddG_reference_pdb:
        ddG_ref_pose = pose_from_pdb(args.ddG_reference_pdb)
    if args.fold_tree:
        fold_tree = create_fold_tree(args.fold_tree)
        pose.fold_tree(fold_tree)
        if args.coordinate_reference_pdb:
            coord_ref_pose.fold_tree(fold_tree)
        if args.ddG_coordinate_reference_pdb:
            ddG_ref_pose.fold_tree(fold_tree)
    for chi in args.chi_dihedrals:
        chi_info = chi.split(",")
        pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
        # if args.coordinate_reference_pdb:
        #     coord_ref_pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
        # if args.ddG_reference_pdb:
        #     ddG_ref_pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
    # Applying symmetry if specified
    set_symmetry(args.symmetry, [pose, coord_ref_pose, ddG_ref_pose])
    # create score function
    score_function = set_score_function(args.score_function, symmetry=args.symmetry, score_terms=args.score_terms)
    # get pose index of the substrate and catalytic residues
    theozyme_positions = set() # including substrates and catalytic residues
    if args.enzyme_design_constraints:
        pdb_info = pose.pdb_info()
        match_substrate_pose_indices, match_res_pose_indices = get_match_pose_indices(pdb_info, args.pdb, args.symmetry)
        match_indices = set(match_substrate_pose_indices + match_res_pose_indices)
        theozyme_positions.update(match_index for match_index in match_indices)
    if args.substrates:
        substrate_pose_indices = pdb_to_pose_numbering(pose, args.substrates)
        theozyme_positions.update(substrate_pose_indices)
    # constraint
    if args.constraints:
        add_fa_constraints_from_cmdline(pose, score_function)
    # enzyme design constraint
    if args.enzyme_design_constraints:
        enzdes_cst = create_enzdes_cst()
        enzdes_cst.apply(pose)
    # define rmsd mover
    if args.no_rmsd_residues:
        no_rmsd_residues = pdb_to_pose_numbering(pose, args.no_rmsd_residues)
        rmsd_metric = RMSDMetric(pose, NotResidueSelector(ResidueIndexSelector(",".join(str(res) for res in no_rmsd_residues))))
    else:
        rmsd_metric = RMSDMetric(pose)
    # native amino acid type
    if args.favor_native_residue and (args.design_positions or args.design_active_site):
        favor_nataa = FavorNativeResidue(pose, args.favor_native_residue)
        favor_nataa.apply(pose)
    # pre minimization
    if args.pre_minimization and args.substrate_rigid_body_transformations and len(theozyme_positions) > 0:
        pre_minimizer = pre_minimization(pose, theozyme_positions)
        pre_minimizer.apply(pose)
    # convert pdb numbering to pose numbering
    mutation_pose_indices = pdb_to_pose_numbering(pose, args.mutations)
    design_pose_indices = pdb_to_pose_numbering(pose, args.design_positions)
    # create task factory
    task_factory, substitution_repacking_selection = create_task_factory(point_mutations = mutation_pose_indices, design_positions = design_pose_indices, \
            theozyme_positions = theozyme_positions, design_active_site = args.design_active_site, \
            repacking_range = args.repack_interface_only, ddG_ref_pose = ddG_ref_pose, no_cystine = args.no_cystine, \
            noncanonical_amino_acids = args.noncanonical_amino_acids)
    # coordinate constraint
    if len(args.no_coordinate_constraint_residues) > 0:
        no_coord_cst_residues = pdb_to_pose_numbering(pose, args.no_coordinate_constraint_residues)
        no_coord_cst_selection = ResidueIndexSelector(",".join(str(no_coord_cst_res) for no_coord_cst_res in no_coord_cst_residues))
    else:
        no_coord_cst_selection = False # All residues are restrained.
    if args.no_optimization_region_coordinate_constraint:
        if substitution_repacking_selection:
            if no_coord_cst_selection:
                no_coord_cst_selection = OrResidueSelector(no_coord_cst_selection, substitution_repacking_selection)
            else:
                no_coord_cst_selection = substitution_repacking_selection
        elif substitution_repacking_selection is None:
            no_coord_cst_selection = None # Nothing is restrained.
    coord_cst_gen = create_coord_cst(coord_ref_pose = coord_ref_pose, \
            coord_boundary = args.coordinate_constraint_bounded_width, no_coord_cst_selection = no_coord_cst_selection)
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    add_csts.apply(pose)
    # create move map
    if not args.minimize_interface_only:
        substitution_repacking_selection = None
    all_substrate_pose_indices = set()
    if args.substrate_rigid_body_transformations and (args.substrates or args.enzyme_design_constraints):
        if args.substrates:
            all_substrate_pose_indices.update(substrate_pose_indices)
        if args.enzyme_design_constraints:
            all_substrate_pose_indices.update(match_substrate_pose_indices)
    move_map = create_move_map(pose, substitution_repacking_selection = substitution_repacking_selection, ddG_ref_pose = ddG_ref_pose, ligand_res_indices = all_substrate_pose_indices)
    # create fast relax mover
    fr = create_fast_relax_mover(score_function, task_factory, move_map = move_map)
    name = args.pdb.split("/")[-1][:-4]
    if args.debug_mode:
        print(pose.fold_tree())
        if args.enzyme_design_constraints:
            print(match_substrate_pose_indices)
            print(match_res_pose_indices)
        print(task_factory.create_task_and_apply_taskoperations(pose))
        print(move_map)
    else:
        run_job_distributor(pose, score_function, args.decoys, name, args.annotated_name, fr, rmsd_metric)


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    main(args)
