#!/usr/bin/python3
import argparse
import math
import os
from pyrosetta import *
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictToRepackingRLT, RestrictAbsentCanonicalAASRLT, \
    PreventRepackingRLT, PreventRepacking
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import \
    add_fa_constraints_from_cmdline, ConstraintSet, AmbiguousConstraint, \
    AtomPairConstraint, AngleConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc, \
    FlatHarmonicFunc, CircularHarmonicFunc
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueIndexSelector, InterGroupInterfaceByVectorSelector, \
    NeighborhoodResidueSelector
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
    parser.add_argument("-ddG_ref", "--ddG_reference_pdb", type=str)
    parser.add_argument("-params", "--parameters_files", type=str, nargs="*")
    parser.add_argument("-optH", "--optimize_protonation_state", action="store_true")
    parser.add_argument("-symm", "--symmetry", type=str)
    parser.add_argument("-ft", "--fold_tree", type=str, nargs="*")
    parser.add_argument("-chis", "--chi_dihedrals", type=str, nargs="*", default=list())
    parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst")
    parser.add_argument("--score_terms", type=str, nargs="*", default=list())
    parser.add_argument("-coord_boundary", "--coordinate_constraint_bounded_width", type=float)
    parser.add_argument("-no_coord_cst", "--no_coordinate_constraint_residues", \
            type=str, nargs="*", default=list())
    parser.add_argument("-no_opt_coord_cst", "--no_coordinate_constraint_on_optimization_region", \
            action="store_true")
    parser.add_argument("-cst", "--constraints", type=str)
    parser.add_argument("-enzdes_cst", "--enzyme_design_constraints", type=str)
    parser.add_argument("-static", "--static_positions", type=str, nargs="*", default=list())
    parser.add_argument("-muts", "--mutations", type=str, nargs="*", default=list())
    parser.add_argument("-des", "--design_positions", type=str, nargs="*", default=list())
    parser.add_argument("-subs", "--substrates", type=str, nargs="*")
    parser.add_argument("-cat", "--catalytic_residues", type=str, nargs="*")
    parser.add_argument("--design_binding_site", action="store_true")
    parser.add_argument("-nataa", "--favor_native_residue", type=float)
    parser.add_argument("-noaa", "--excluded_amino_acid_types", type=str, \
            help="String of one-letter AA codes.")
    parser.add_argument("-ncaa", "--noncanonical_amino_acids", type=str, \
            nargs="*", default=list(), help="name3")
    parser.add_argument("-premin", "--pre_minimization", action="store_true")
    parser.add_argument("-tform", "--substrate_rigid_body_transformations", action="store_true")
    parser.add_argument("-rpk_nbh", "--repack_neighborhood_only", action="store_true")
    parser.add_argument("-min_nbh", "--minimize_neighborhood_only", action="store_true")
    parser.add_argument("-no_rmsd", "--no_rmsd_residues", type=str, nargs="*", default=list())
    parser.add_argument("-n", "--n_decoys", type=int, default=50)
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

def customize_score_function(rep_type="hard", cartesian=False, symmetry=False, \
        membrane=False, constraint_weight=1.0, score_terms=list()):
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
            elif flag is True:
                break
    return match_substrate_pose_indices, match_res_pose_indices

def create_coord_cst(coord_ref_pose=None, coord_boundary = None, \
            no_coord_cst_selection = False, sidechain=False):
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

def apply_atom_pair_constraint(pose, atom_name_list, residue_id_list, distance, \
        standard_deviation, bound=None):
    atom_list = list()
    for atom_name, residue_id in zip(atom_name_list, residue_id_list):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(distance) is list:
        harmonic_cst = AmbiguousConstraint()
        if bound:
            for i, dis_sd in enumerate(zip(distance, standard_deviation)):
                if bound:
                    harmonic_fc = FlatHarmonicFunc(dis_sd[0], dis_sd[1], bound[i])
                else:
                    harmonic_fc = HarmonicFunc(dis_sd[0], dis_sd[1])
                cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
                harmonic_cst.add_individual_constraint(cst)
    else:
        if bound:
            harmonic_fc = FlatHarmonicFunc(distance, standard_deviation, bound)
        else:
            harmonic_fc = HarmonicFunc(distance, standard_deviation)
        harmonic_cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
    pose.add_constraint(harmonic_cst)

def apply_angle_constraint(pose, atom_name_list, residue_id_list, degree, standard_deviation):
    atom_list = list()
    for atom_name, residue_id in zip(atom_name_list, residue_id_list):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(degree) is list:
        circular_harmonic_cst = AmbiguousConstraint()
        for deg, sd in zip(degree, standard_deviation):
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * deg, sd)
            cst = AngleConstraint(atom_list[0], atom_list[1], atom_list[2], circular_harmonic_fc)
            circular_harmonic_cst.add_individual_constraint(cst)
    else:
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * degree, standard_deviation)
        circular_harmonic_cst = AngleConstraint(atom_list[0], atom_list[1], \
                atom_list[2], circular_harmonic_fc)
    pose.add_constraint(circular_harmonic_cst)

def apply_dihedral_constraint(pose, atom_name_list, residue_id_list, degree, standard_deviation):
    atom_list = list()
    for atom_name, residue_id in zip(atom_name_list, residue_id_list):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(degree) is list:
        circular_harmonic_cst = AmbiguousConstraint()
        for deg, sd in zip(degree, standard_deviation):
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * deg, sd)
            cst = DihedralConstraint(atom_list[0], atom_list[1], atom_list[2], atom_list[3], circular_harmonic_fc)
            circular_harmonic_cst.add_individual_constraint(cst)
    else:
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * degree, standard_deviation)
        circular_harmonic_cst = DihedralConstraint(atom_list[0], atom_list[1], \
                atom_list[2], atom_list[3], circular_harmonic_fc)
    pose.add_constraint(circular_harmonic_cst)

def create_enzdes_cst():
# Add enzdes constraints
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

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

def get_jump_edges(fold_tree, rigid_body_tform_pose_indices):
    jump_edges = set()
    for rigid_body_tform_pose_index in rigid_body_tform_pose_indices:
        jump_edge = fold_tree.get_residue_edge(rigid_body_tform_pose_index)
        if not jump_edge.is_jump():
            raise Exception("Edge of the ligand is not a jump edge.")
        jump_edges.add(jump_edge.label())
    return jump_edges

def select_neighborhood_region(focus_selection, include_focus: bool, method="vector"):
    if method == "vector":
        neighborhood_selection = InterGroupInterfaceByVectorSelector()
        neighborhood_selection.group1_selector(focus_selection)
        neighborhood_selection.group2_selector(NotResidueSelector(focus_selection))
        if include_focus:
            neighborhood_selection = OrResidueSelector(neighborhood_selection, \
                    focus_selection)
        else:
            neighborhood_selection = AndResidueSelector(neighborhood_selection, \
                    NotResidueSelector(focus_selection))
    elif method == "distance":
        neighborhood_selection = NeighborhoodResidueSelector()
        neighborhood_selection.set_focus_selector(focus_selection)
        neighborhood_selection.set_distance(8.0)
        neighborhood_selection.set_include_focus_in_subset(include_focus)
    return neighborhood_selection

def fix_residue_selection(ddG_ref_pose, selector):
    selection_vector = selector.apply(ddG_ref_pose)
    selection = ResidueIndexSelector(",".join(filter(lambda x: x is not None, \
            map(lambda x, y: str(x) if y == 1 else None, \
            range(1, len(selection_vector) + 1), selection_vector))))
    return selection

def pre_minimization(pose, theozyme_positions, jump_edges: set = set()):
    pre_minimization_selector = ResidueIndexSelector(",".join(str(theozyme_position) \
            for theozyme_position in theozyme_positions))
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(pre_minimization_selector.apply(pose))
    for jump_edge in jump_edges:
        move_map.set_jump(jump_edge, True)
    pre_minimizer = MinMover()
    pre_minimizer.min_type("lbfgs_armijo_nonmonotone")
    pre_minimizer.movemap(move_map)
    return pre_minimizer

def create_task_factory(static_positions: set = set(), point_mutations: set = set(), \
        design_positions: set = set(), theozyme_positions: set = set(), \
        design_binding_site: bool = False, only_repack_neighborhood: bool = False, \
        ddG_ref_pose = None, excluded_amino_acid_types: str = None, \
        noncanonical_amino_acids: list = None):
    """
    Priority:
        1. Freeze static residues.
        2. Introduce point mutations.
        3. Perform site-directed design.
        4. Design protein-ligand interfaces.
    """
    # Create a task factory.
    task_factory = TaskFactory()
    if noncanonical_amino_acids:
        ncaa_palette = CustomBaseTypePackerPalette()
        for ncaa in noncanonical_amino_acids:
            ncaa_palette.add_type(ncaa)
        task_factory.set_packer_palette(ncaa_palette)

    # Every position is designable by defalut other than explicit specification.
    task_factory.push_back(IncludeCurrent())

    # Specify the point mutations.
    mutation_positions = set()
    if len(point_mutations) > 0:
        for point_mutation in point_mutations:
            mutating_position, target_aa = point_mutation.split(",")
            if not mutating_position in static_selection:
                restriction = RestrictAbsentCanonicalAASRLT()
                restriction.aas_to_keep(target_aa)
                task_factory.push_back(OperateOnResidueSubset(restriction, \
                        ResidueIndexSelector(mutating_position)))
                mutation_positions.add(mutating_position)
    mutation_selection = ResidueIndexSelector(",".join(mutation_positions) + ",")
    interface_design_excluding_positions = static_positions + mutation_positions

    # Select the design positions.
    design_selection = OrResidueSelector()
    # Specify the site-directed design positions.
    if len(static_positions) > 0:
        design_positions = design_positions - static_positions
    if len(mutation_positions) > 0:
        design_positions = design_positions - mutation_positions
    if len(design_positions) > 0:
        site_directed_design_selection = ResidueIndexSelector(\
            ",".join(str(design_position) for design_position in design_positions))
        design_selection.add_residue_selector(site_directed_design_selection)
    interface_design_excluding_positions = interface_design_excluding_positions + design_positions
    # Specify the theozyme positions.
    theozyme_selection = ResidueIndexSelector(",".join(str(theozyme_position)\
            for theozyme_position in theozyme_positions) + ",")
    if len(theozyme_positions) > 0:
        # Design the theozyme-protein interface.
        if design_binding_site:
            protein_binding_site_selection = select_neighborhood_region(theozyme_selection, False)
            if len(interface_design_excluding_positions) > 0:
                protein_binding_site_selection = AndResidueSelector(protein_binding_site_selection, \
                        NotResidueSelector(ResidueIndexSelector(",".join(interface_design_excluding_positions))))
            if ddG_ref_pose: # Fix the design positions.
                protein_binding_site_selection = fix_residue_selection(protein_binding_site_selection, ddG_ref_pose)
            design_selection.add_residue_selector(protein_binding_site_selection)

    # Exclude some AA types if specified.
    if excluded_amino_acid_types:
        all_AAs = set("AGILPVFWYDERHKSTMNQ")
        excluded_AAs = set(excluded_amino_acid_types)
        restriction = RestrictAbsentCanonicalAASRLT("".join(\
                available_AA for available_AA in all_AAs - excluded_AAs))
        restriction.aas_to_keep()
        task_factory.push_back(OperateOnResidueSubset(restriction, design_selection))

    # Specify the static positions.
    static_selection = ResidueIndexSelector(",".join(str(static_position) \
            for static_position in static_positions) + ",")

    # Identify the repack and static regions.
    if len(mutation_positions) > 0 or len(design_positions) > 0 or len(theozyme_positions) > 0:
        # Site-directed AA substitutions positions.
        substitution_selection = OrResidueSelector(mutation_selection, design_selection)
        if only_repack_neighborhood:
            # Repack the neighborhood region of the site-directed AA substitutions and theozyme positions.
            substitution_theozyme_selection = OrResidueSelector(substitution_selection, theozyme_selection)
            substitution_repacking_selection = select_neighborhood_region(substitution_theozyme_selection, True)
            # Exclude the site-directed AA substitutions from the repacking region.
            repacking_selection = AndResidueSelector(substitution_repacking_selection, \
                    NotResidueSelector(substitution_selection))
            if ddG_ref_pose: # Fix the repacking positions.
                repacking_selection = fix_residue_selection(repacking_selection, ddG_ref_pose)
                substitution_repacking_selection = OrResidueSelector(substitution_selection, \
                        repacking_selection)
            # Exclude the static positions from the repacking region.
            repacking_selection = AndResidueSelector(substitution_repacking_selection, \
                    NotResidueSelector(static_selection))
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), repacking_selection))
            # Static region.
            task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), \
                    OrResidueSelector(substitution_selection, repacking_selection), True))
        else:
            # Exclude the site-directed AA substitutions and static positions from repacking.
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), \
                    OrResidueSelector(substitution_selection, static_selection), True))
            substitution_repacking_selection = None
    else:
        if only_repack_neighborhood:
            # Do not repack any part of the pose.
            task_factory.push_back(PreventRepacking())
            substitution_repacking_selection = False
        else:
            # Exclude the static positions from the repacking region.
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), static_selection, True))
            substitution_repacking_selection = None
    # When setting up the movemap , "substitution_repacking_selection" is used as the 
    # focus selector so that its neighborhood region is also subject to minimization.
    # The "substitution_repacking_selection" selector may includes static positions.
    # All these static positions will be excluded from the minimization region.
    return task_factory, substitution_repacking_selection

def create_move_map(pose, focus_selection = None, only_minimize_neighborhood: bool = False, \
        static_positions: set = set(), ddG_ref_pose = None, jump_edges: set = set()):
    '''
    When "only_minimize_neighborhood" is True, only the "focus_selection" along with its 
    neighborhood region are subject to minimization except the static positions.
    Freeze a residue does not mean not to repack and minimize the residues around it, 
    especially when the freezed residue is included in the theozyme positions.
    '''
    move_map = MoveMap()
    static_selection = ResidueIndexSelector(",".join(str(static_position) \
            for static_position in static_positions) + ",")
    if only_minimize_neighborhood:
        if focus_selection:
            minimization_selection = select_neighborhood_region(focus_selection, True)
            if len(static_positions) > 0:
                minimization_selection = AndResidueSelector(minimization_selection, \
                        NotResidueSelector(static_selection))
            if ddG_ref_pose: # Pass residue indices vectors to the movemap.
                minimization_vector = minimization_selection.apply(ddG_ref_pose)
                move_map.set_bb(minimization_vector)
                move_map.set_chi(minimization_vector)
            else: # Pass residue selectors to a movemap factory. Under development.
                minimization_vector = minimization_selection.apply(pose)
                move_map.set_bb(minimization_vector)
                move_map.set_chi(minimization_vector)
        elif focus_selection is None:
            raise Exception("Using -min_nbh is not allowed without using -rpk_nbh at the same time.")
        elif focus_selection is False:
            move_map.set_bb(False)
            move_map.set_chi(False)
    else:
        if len(static_positions) > 0:
            minimization_selection = NotResidueSelector(static_selection)
            move_map.set_bb(minimization_selection.apply(pose))
            move_map.set_chi(minimization_selection.apply(pose))
        else:
            move_map.set_bb(True)
            move_map.set_chi(True)
    for jump_edge in jump_edges > 0:
        move_map.set_jump(jump_edge, True)
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

def run_job_distributor(pose, score_function, name, annotated_name, *movers, n_decoys=50):
    job_distributor = PyJobDistributor(name, n_decoys, score_function)
    pdb_name = job_distributor.pdb_name
    while not job_distributor.job_complete:
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
            current_id = job_distributor.current_id
            os.rename(pdb_name + "_" + str(current_id) + ".pdb", \
                    pdb_name + "_" + str(current_id) + mutations_str + ".pdb")

def main(args):
    pose = pose_from_pdb(args.pdb)
    coord_ref_pose = None
    if args.coordinate_reference_pdb:
        coord_ref_pose = pose_from_pdb(args.coordinate_reference_pdb)
    ddG_ref_pose = None
    if args.ddG_reference_pdb:
        if args.ddG_reference_pdb is "True":
            ddG_ref_pose = Pose(pose)
        else:
            ddG_ref_pose = pose_from_pdb(args.ddG_reference_pdb)
    if args.fold_tree:
        fold_tree = create_fold_tree(args.fold_tree)
        pose.fold_tree(fold_tree)
    for chi in args.chi_dihedrals:
        chi_info = chi.split(",")
        pose.set_chi(int(chi_info[0]), int(chi_info[1]), float(chi_info[2]))
    # Applying the symmetry if specified.
    set_symmetry(args.symmetry, [pose, coord_ref_pose, ddG_ref_pose])
    # Create the score function.
    score_function = set_score_function(args.score_function, \
            symmetry=args.symmetry, score_terms=args.score_terms)
    # Add constraint files from the command line.
    if args.constraints:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Add coordinate constraints.
    if len(args.no_coordinate_constraint_residues) > 0:
        no_coord_cst_residues = pdb_to_pose_numbering(pose, args.no_coordinate_constraint_residues)
        no_coord_cst_selection = ResidueIndexSelector(",".join(str(no_coord_cst_res) \
                for no_coord_cst_res in no_coord_cst_residues))
    else:
        no_coord_cst_selection = False
    coord_cst_gen = create_coord_cst(coord_ref_pose = coord_ref_pose, \
            coord_boundary = args.coordinate_constraint_bounded_width, \
            no_coord_cst_selection = no_coord_cst_selection)
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    add_csts.apply(pose)
    # Add enzyme design constraints.
    if args.enzyme_design_constraints:
        enzdes_cst = create_enzdes_cst()
        enzdes_cst.apply(pose)
    # Create the RMSD metric.
    if args.no_rmsd_residues:
        no_rmsd_residues = pdb_to_pose_numbering(pose, args.no_rmsd_residues)
        rmsd_metric = RMSDMetric(pose, NotResidueSelector(ResidueIndexSelector(\
                ",".join(str(res) for res in no_rmsd_residues))))
    else:
        rmsd_metric = RMSDMetric(pose)
    # Favor native AA types.
    if args.favor_native_residue and (args.design_positions or args.design_binding_site):
        favor_nataa = FavorNativeResidue(pose, args.favor_native_residue)
        favor_nataa.apply(pose)
    # Get pose indices of substrates and catalytic residues.
    theozyme_pose_indices = set()
    if args.enzyme_design_constraints:
        pdb_info = pose.pdb_info()
        match_substrate_pose_indices, match_res_pose_indices = \
                get_match_pose_indices(pdb_info, args.pdb, args.symmetry)
        match_indices = set(match_substrate_pose_indices + match_res_pose_indices)
        theozyme_pose_indices.update(match_index for match_index in match_indices)
    if args.substrates:
        substrate_pose_indices = pdb_to_pose_numbering(pose, args.substrates)
        theozyme_pose_indices.update(substrate_pose_indices)
    if args.catalytic_residues:
        catalytic_residue_pose_indices = pdb_to_pose_numbering(pose, args.catalytic_residues)
        theozyme_pose_indices.update(catalytic_residue_pose_indices)
    # Convert pdb numberings to pose numberings.
    static_pose_indices = pdb_to_pose_numbering(pose, args.static_positions)
    mutation_pose_indices = pdb_to_pose_numbering(pose, args.mutations)
    design_pose_indices = pdb_to_pose_numbering(pose, args.design_positions)
    # Rigid body transformations.
    rigid_body_tform_pose_indices = set()
    if args.substrate_rigid_body_transformations:
        if args.substrates:
            rigid_body_tform_pose_indices.update(substrate_pose_indices)
        if args.enzyme_design_constraints:
            rigid_body_tform_pose_indices.update(match_substrate_pose_indices)
    rigid_body_tform_pose_indices = rigid_body_tform_pose_indices - static_pose_indices
    free_jump_edges = get_jump_edges(pose.fold_tree(), rigid_body_tform_pose_indices)
    # Perform pre-minimization.
    free_theozyme_pose_indices = theozyme_pose_indices - static_pose_indices
    if args.pre_minimization and len(free_theozyme_pose_indices) > 0:
        pre_minimizer = pre_minimization(pose, free_theozyme_pose_indices, jump_edges=free_jump_edges)
        pre_minimizer.apply(pose)
    # Create the task factory.
    task_factory, substitution_repacking_selection = create_task_factory(\
            static_positions = static_pose_indices, point_mutations = mutation_pose_indices, \
            design_positions = design_pose_indices, theozyme_positions = theozyme_pose_indices, \
            design_binding_site = args.design_binding_site, only_repack_neighborhood = args.repack_neighborhood_only, \
            ddG_ref_pose = ddG_ref_pose, excluded_amino_acid_types = args.excluded_amino_acid_types, \
            noncanonical_amino_acids = args.noncanonical_amino_acids)
    # Create the move map.
    move_map = create_move_map(pose, focus_selection = substitution_repacking_selection, \
            only_minimize_neighborhood = args.minimize_neighborhood_only, \
            static_positions = static_pose_indices, ddG_ref_pose = ddG_ref_pose, \
            jump_edges = free_jump_edges)
    # Create the fast relax mover.
    fr = create_fast_relax_mover(score_function, task_factory, move_map = move_map)
    # Update coordinate constraints. Under development.
    if args.no_coordinate_constraint_on_optimization_region:
        if substitution_repacking_selection:
            if no_coord_cst_selection:
                no_coord_cst_selection = OrResidueSelector(\
                        no_coord_cst_selection, substitution_repacking_selection)
            else:
                no_coord_cst_selection = substitution_repacking_selection
        elif substitution_repacking_selection is None:
            no_coord_cst_selection = None
    coord_cst_gen = create_coord_cst(coord_ref_pose = coord_ref_pose, \
            coord_boundary = args.coordinate_constraint_bounded_width, \
            no_coord_cst_selection = no_coord_cst_selection)
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    add_csts.apply(pose)
    # Run jobs.
    name = args.pdb.split("/")[-1][:-4]
    if args.debug_mode:
        print(pose.fold_tree())
        if args.enzyme_design_constraints:
            print(match_substrate_pose_indices)
            print(match_res_pose_indices)
        print(task_factory.create_task_and_apply_taskoperations(pose))
        print(move_map)
    else:
        run_job_distributor(pose, score_function, name, args.annotated_name, \
                fr, rmsd_metric, n_decoys=args.n_decoys)


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    main(args)
