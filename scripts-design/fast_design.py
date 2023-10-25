#!/usr/bin/python3
import argparse
import math
import numpy as np
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
    ResidueIndexSelector, ResidueNameSelector, \
    InterGroupInterfaceByVectorSelector, NeighborhoodResidueSelector
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
    parser.add_argument("-ft", "--fold_tree", type=str, nargs="*")
    parser.add_argument("-edges", "--alter_jump_edges", type=str, nargs="*", help=\
            "[$edge,$atom1,$atom2] or [$edge,$upstream_edge] or [$edge,$upstream_edge,$atom1,$atom2] * n")
    parser.add_argument("-chis", "--chi_dihedrals", type=str, nargs="*", default=list(), \
            help="$chain$residue_index,$chi,$degree or $residue_name3,$chi,$degree")
    parser.add_argument("-symm", "--symmetry", type=str)
    parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst")
    parser.add_argument("--score_terms", type=str, nargs="*", default=list())
    parser.add_argument("-coord_boundary", "--coordinate_constraint_bounded_width", type=float)
    parser.add_argument("-no_coord_cst", "--no_coordinate_constraint_residues", \
            type=str, nargs="*", default=list())
    parser.add_argument("-no_opt_coord_cst", "--no_coordinate_constraint_on_optimization_region", \
            action="store_true")
    parser.add_argument("-cst", "--constraints", type=str)
    parser.add_argument("-dis_atoms", "--distance_constraint_atoms", type=str, nargs="*", \
            help="$chain$residue_index,$atom_name or $residue_name3,$atom_name * 2n")
    parser.add_argument("-dis_params", "--distance_constraint_parameters", type=str, nargs="*", \
            help="$distance,$standard_deviation or $dis,$sd,$dis,$sd,... * n")
    parser.add_argument("-dis_bounds", "--distance_constraint_boundaries", type=str, nargs="*", \
            help="$boundary or $boundary,$boundary,... or placeholder * n")
    parser.add_argument("-angle_atoms", "--angle_constraint_atoms", type=str, nargs="*", \
            help="$chain$residue_index,$atom_name or $residue_name3,$atom_name * 3n")
    parser.add_argument("-angle_params", "--angle_constraint_parameters", type=str, nargs="*", \
            help="$degree,$standard_deviation or $deg,$sd,$deg,$sd,... * n")
    parser.add_argument("-dihe_atoms", "--dihedral_constraint_atoms", type=str, nargs="*", \
            help="$chain$residue_index,$atom_name or $residue_name3,$atom_name * 4n")
    parser.add_argument("-dihe_params", "--dihedral_constraint_parameters", type=str, nargs="*", \
            help="$degree,$standard_deviation or $deg,$sd,$deg,$sd,... * n")
    parser.add_argument("-enzdes_cst", "--enzyme_design_constraints", type=str)
    parser.add_argument("-static", "--static_residues", type=str, nargs="*", default=list())
    parser.add_argument("-static_ids", "--static_residue_identities", type=str, nargs="*", help="name3")
    parser.add_argument("-cat", "--catalytic_residues", type=str, nargs="*")
    parser.add_argument("-subs", "--substrates", type=str, nargs="*")
    parser.add_argument("-sub_ids", "--substrate_identities", type=str, nargs="*", help="name3")
    parser.add_argument("-premin", "--pre_minimization", action="store_true")
    parser.add_argument("-muts", "--mutations", type=str, nargs="*", default=list())
    parser.add_argument("-des", "--design_residues", type=str, nargs="*", default=list())
    parser.add_argument("-des_bs", "--design_binding_site", action="store_true")
    parser.add_argument("-nataa", "--favor_native_residue", type=float)
    parser.add_argument("-noaa", "--excluded_amino_acid_types", type=str, \
            help="String of one-letter AA codes.")
    parser.add_argument("-ncaa", "--noncanonical_amino_acids", type=str, nargs="*", \
            default=list(), help="name3")
    parser.add_argument("-rpk_nbh", "--repack_neighborhood_only", action="store_true")
    parser.add_argument("-rpk_bs", "--repack_binding_site", action="store_true")
    parser.add_argument("-min_nbh", "--minimize_neighborhood_only", action="store_true")
    parser.add_argument("-tform", "--substrate_rigid_body_transformations", action="store_true")
    parser.add_argument("-tform_enzdes_subs", "--enzyme_design_constraints_substrate_transformations", \
            action="store_true")
    parser.add_argument("-no_rmsd", "--no_rmsd_residues", type=str, nargs="*", default=list())
    parser.add_argument("-n", "--n_decoys", type=int, default=50)
    parser.add_argument("-o", "--output_filename_prefix", type=str)
    parser.add_argument("--annotated_name", action="store_true")
    parser.add_argument("-nosave", "--no_save_decoys", action="store_true")
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

def load_pdb_as_pose(pdb, coordinate_reference_pdb, ddG_reference_pdb):
    pose = pose_from_pdb(pdb)
    sequence_length = len(pose.sequence())
    coord_ref_pose = None
    if coordinate_reference_pdb:
        coord_ref_pose = pose_from_pdb(coordinate_reference_pdb)
    ddG_ref_pose = None
    if ddG_reference_pdb:
        if ddG_reference_pdb == "True":
            ddG_ref_pose = Pose(pose)
        else:
            ddG_ref_pose = pose_from_pdb(ddG_reference_pdb)
        ddG_reference_sequence_length = len(ddG_ref_pose.sequence())
        truncate_sequence_end = ddG_reference_sequence_length - sequence_length
        assert truncate_sequence_end >= 0
    else:
        truncate_sequence_end = 0
    return pose, coord_ref_pose, ddG_ref_pose, sequence_length, truncate_sequence_end

def pdb_to_pose_numbering(pose, chain_id_pdb_indices):
    '''
    Return a set of residue position strings.
    '''
    pose_indices = set()
    jump_pose_indices = set()
    for chain_id_pdb_index in chain_id_pdb_indices:
        chain_id = chain_id_pdb_index[0]
        if "," in chain_id_pdb_index:
            pdb_index, suffix = chain_id_pdb_index[1:].split(",")
            pose_index = pose.pdb_info().pdb2pose(chain_id, int(pdb_index))
            if pose_index > 0:
                pose_indices.add(str(pose_index) + "," + suffix)
        else:
            if "-" in chain_id_pdb_index:
                start_pdb_index, end_pdb_index = chain_id_pdb_index[1:].split("-")
                pdb_indices = range(int(start_pdb_index), int(end_pdb_index) + 1)
            else:
                pdb_indices = [int(chain_id_pdb_index[1:])]
            for i, pdb_index in enumerate(pdb_indices):
                pose_index = pose.pdb_info().pdb2pose(chain_id, pdb_index)
                if pose_index > 0:
                    pose_indices.add(str(pose_index))
                    if i == 0:
                        jump_pose_indices.add(str(pose_index))
    return pose_indices, jump_pose_indices

def residue_name3_selector(pose, name3_list):
    pose_indices = set()
    for pose_index in range(1, len(pose.sequence()) + 1):
        if pose.residue(pose_index).name3() in name3_list:
            pose_indices.add(str(pose_index))
    return pose_indices

def create_fold_tree(edges):
    fold_tree = FoldTree()
    for edge_str in edges:
        edge = edge_str.split(",")
        if len(edge) == 3:
            fold_tree.add_edge(int(edge[0]), int(edge[1]), int(edge[2]))
        elif len(edge) == 4:
            fold_tree.add_edge(int(edge[0]), int(edge[1]), edge[2], edge[3])
        else:
            raise Exception("An edge infomation list should contain 3 or 4 elements.")
    return fold_tree

def alter_fold_tree_jump_edges(fold_tree, alter_jump_edges):
    jump_edge_labels = list(range(1, fold_tree.num_jump() + 1))
    alter_jump_edges_dict = dict()
    for alter_jump_edge in alter_jump_edges:
        alter_jump_edge = alter_jump_edge.strip("[").strip("]").split(",")
        jump_edge_label = int(alter_jump_edge[0])
        if jump_edge_label < 0:
            jump_edge_label = jump_edge_labels[jump_edge_label]
        if len(alter_jump_edge) == 2:
            upstream_edge_label = alter_jump_edge[1:]
            alter_jump_edges_dict[jump_edge_label] = (int(upstream_edge_label))
        elif len(alter_jump_edge) == 3:
            start_atom, stop_atom = alter_jump_edge[1:]
            alter_jump_edges_dict[jump_edge_label] = (0, start_atom, stop_atom)
        elif len(alter_jump_edge) == 4:
            upstream_edge_label, start_atom, stop_atom = alter_jump_edge[1:]
            alter_jump_edges_dict[jump_edge_label] = (int(upstream_edge_label), start_atom, stop_atom)
        else:
            raise Exception("An alter edge list should contain 3 or 4 elements.")
    edge_strings = list()
    for edge_str in fold_tree.to_string().split("EDGE")[1:]:
        edge = tuple(filter(lambda x: x != "", edge_str.strip(" ").split(" ")))
        downstream_residue = int(edge[1])
        edge_label = fold_tree.get_residue_edge(downstream_residue).label()
        alter_jump_edge = alter_jump_edges_dict.get(edge_label)
        if alter_jump_edge:
            upstream_edge_label = alter_jump_edge[0]
            if upstream_edge_label < 0:
                upstream_edge_label = jump_edge_labels[upstream_edge_label]
            if upstream_edge_label == 0:
                upstream_residue = int(edge[0])
            else:
                upstream_residue = fold_tree.jump_edge(upstream_edge_label).stop()
            if len(alter_jump_edge) == 1:
                edge = (str(upstream_residue), str(downstream_residue), str(edge_label))
            elif len(alter_jump_edge) == 3:
                edge = (str(upstream_residue), str(downstream_residue), \
                        alter_jump_edge[1], alter_jump_edge[2])
        edge_strings.append(",".join(edge))
    return create_fold_tree(edge_strings)

def set_chi_dihedral(pose, chi_dihedrals):
    for chi in chi_dihedrals:
        residue_info, chi_index, dihedral_value = chi.split(",")
        try:
            chi_residue_pose_indices, _ = pdb_to_pose_numbering(pose, [residue_info])
        except:
            chi_residue_pose_indices = residue_name3_selector(pose, [residue_info])
        for chi_residue_pose_index in chi_residue_pose_indices:
            pose.set_chi(int(chi_index), int(chi_residue_pose_index), float(dihedral_value))

def create_coord_cst(coord_ref_pose=None, coord_boundary=None, no_coord_cst_selection=False, \
        sidechain=False):
    coord_cst_gen = CoordinateConstraintGenerator()
    if coord_boundary:
        coord_cst_gen.set_bounded(True)
        coord_cst_gen.set_bounded_width(coord_boundary)
    if coord_ref_pose is not None:
        coord_cst_gen.set_reference_pose(coord_ref_pose)
    if no_coord_cst_selection is not None:
        if no_coord_cst_selection is not False:
            coord_cst_gen.set_residue_selector(NotResidueSelector(no_coord_cst_selection))
    else:
        coord_cst_gen.set_residue_selector(OrResidueSelector())
    if sidechain:
        coord_cst_gen.set_sidechain(True)
    return coord_cst_gen

def create_harmonic_constraint(pose, atom_names, pose_indices, distance, standard_deviation, \
        boundary=None):
    atom_list = list()
    for atom_name, residue_id in zip(atom_names, pose_indices):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(distance) is list or type(distance) is np.ndarray:
        assert len(distance) == len(standard_deviation)
        harmonic_cst = AmbiguousConstraint()
        for i, dis_sd in enumerate(zip(distance, standard_deviation)):
            if boundary:
                harmonic_fc = FlatHarmonicFunc(dis_sd[0], dis_sd[1], boundary[i])
            else:
                harmonic_fc = HarmonicFunc(dis_sd[0], dis_sd[1])
            cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
            harmonic_cst.add_individual_constraint(cst)
    else:
        if boundary:
            harmonic_fc = FlatHarmonicFunc(distance, standard_deviation, boundary)
        else:
            harmonic_fc = HarmonicFunc(distance, standard_deviation)
        harmonic_cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
    return harmonic_cst

def create_circular_harmonic_constraint(pose, atom_names, pose_indices, degree, standard_deviation):
    atom_list = list()
    for atom_name, residue_id in zip(atom_names, pose_indices):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(degree) is list or type(degree) is np.ndarray:
        assert len(degree) == len(standard_deviation)
        circular_harmonic_cst = AmbiguousConstraint()
        for deg, sd in zip(degree, standard_deviation):
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * deg, sd)
            if len(atom_list) == 3:
                cst = AngleConstraint(atom_list[0], atom_list[1], atom_list[2], \
                        circular_harmonic_fc)
            elif len(atom_list) == 4:
                cst = DihedralConstraint(atom_list[0], atom_list[1], atom_list[2], \
                        atom_list[3], circular_harmonic_fc)
            else:
                raise Exception("More than four atoms are specified in the circular harmonic constraint.")
            circular_harmonic_cst.add_individual_constraint(cst)
    else:
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * degree, standard_deviation)
        if len(atom_list) == 3:
            circular_harmonic_cst = AngleConstraint(atom_list[0], atom_list[1], \
                    atom_list[2], circular_harmonic_fc)
        elif len(atom_list) == 4:
            circular_harmonic_cst = DihedralConstraint(atom_list[0], atom_list[1], \
                    atom_list[2], atom_list[3], circular_harmonic_fc)
    return circular_harmonic_cst

def create_constraints(pose, constraint_atoms, constraint_parameters, boundaries=None):
    assert len(constraint_atoms) % len(constraint_parameters) == 0
    constraints = list()
    length = int(len(constraint_atoms) / len(constraint_parameters))
    if boundaries:
        assert len(boundaries) == len(constraint_parameters)
    # Enumerate multiple constraint specifications.
    for cst_index, cst_first_atom in enumerate(range(0, len(constraint_atoms), length)):
        # Get value(s) and standard deviation(s).
        current_cst_parameters = constraint_parameters[cst_index].split(",")
        n_parameters = len(current_cst_parameters)
        if n_parameters == 2:
            value, standard_deviation = float(current_cst_parameters[0]), float(current_cst_parameters[1])
        elif n_parameters > 2:
            assert n_parameters % 2 == 0
            current_cst_parameters = np.array(list(float(param) for param in current_cst_parameters))\
                    .reshape(-1,2).transpose()
            value, standard_deviation = current_cst_parameters[0], current_cst_parameters[1]
        # Get boundary(s).
        if boundaries:
            boundary = boundaries[cst_index]
            if "," in boundary:
                boundary = list()
                for bound in boundary.split(","):
                    try:
                        bound = float(bound)
                    except:
                        bound = None
                boundary.append(bound)
            try:
                boundary = float(boundary)
            except:
                boundary = None
        else:
            boundary = None
        # Get atoms.
        constraints_pose_indices = list()
        current_cst_atoms = constraint_atoms[cst_first_atom:cst_first_atom + length]
        atom_names = list()
        try:
            pose_index_atom_name_set, _ = pdb_to_pose_numbering(pose, current_cst_atoms)
            pose_indices = list()
            for pose_index_atom_name in pose_index_atom_name_set:
                pose_index, atom_name = pose_index_atom_name.split(",")
                pose_indices.append(int(pose_index))
                atom_names.append(atom_name)
            constraints_pose_indices.append(pose_indices)
        except: # Select possibly multiple residues by residue identity.
            cst_residue_name3 = None
            for cst_atom in current_cst_atoms:
                residue_name3, atom_name = cst_atom.split(",")
                if not cst_residue_name3:
                    cst_residue_name3 = residue_name3
                elif cst_residue_name3 != residue_name3:
                    raise Exception("The residue identity assigned to the constraint is not consistent.")
                atom_names.append(atom_name)
            for pose_index in residue_name3_selector(pose, [cst_residue_name3]):
                constraints_pose_indices.append([int(pose_index)] * len(atom_names))
        finally:
            for pose_indices in constraints_pose_indices:
                if len(atom_names) == 2:
                    constraint = create_harmonic_constraint(pose, atom_names, pose_indices, value, \
                            standard_deviation, boundary=boundary)
                else:
                    constraint = create_circular_harmonic_constraint(pose, atom_names, pose_indices, \
                            value, standard_deviation)
                constraints.append(constraint)
    return constraints

def create_enzdes_cst():
# Add enzdes constraints
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

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
                    match_substrate_pose_indices.add(str(substrate_pose_id))
                # match residues
                motif_chain_id = line[49]
                if main_chain == motif_chain_id or not symmetry:
                    motif_res_id = int(line[55:59])
                    motif_pose_id = info.pdb2pose(motif_chain_id, motif_res_id)
                    match_res_pose_indices.add(str(motif_pose_id))
                flag = True
            elif flag is True:
                break
    return match_substrate_pose_indices, match_res_pose_indices

def pose_indices_to_jump_edges(fold_tree, rigid_body_tform_pose_indices):
    jump_edges = set()
    for rigid_body_tform_pose_index in rigid_body_tform_pose_indices:
        jump_edge = fold_tree.get_residue_edge(int(rigid_body_tform_pose_index))
        if not jump_edge.is_jump():
            raise Exception("The specified edge is not a jump edge.")
        jump_edges.add(jump_edge.label())
    return jump_edges

def index_boolean_filter(index, boolean):
    for bool in boolean:
        if bool:
            return index
    return None

def boolean_vector_to_indices_set(boolean_vector, n_monomers=1, truncate_sequence_end=0):
    boolean_vector = np.array(boolean_vector)
    sequence_length = len(boolean_vector) // n_monomers
    boolean_vector = boolean_vector[:n_monomers*sequence_length]\
            .reshape(n_monomers, sequence_length).transpose()
    if truncate_sequence_end > 0:
        boolean_vector = boolean_vector[:-truncate_sequence_end,:]
    indices_iterator = filter(lambda x: x is not None, map(index_boolean_filter, \
            range(1, len(boolean_vector) + 1), boolean_vector))
    return set(str(index) for index in indices_iterator)

def set_symmetry(symmetry, sequence_length, *pose_list):
    n_monomer = 1
    if symmetry:
        sfsm = SetupForSymmetryMover(args.symmetry)
        sfsm.set_keep_pdb_info_labels(True)
        for pose in filter(lambda p: p is not None, pose_list):
            sfsm.apply(pose)
            if not n_monomer:
                n_monomer = round(len(pose.sequence().strip("X")) / sequence_length)
    return n_monomer

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

def pre_minimization(pose, pre_min_positions, jump_edges:set=set()):
    pre_minimization_selector = ResidueIndexSelector(",".join(pre_min_positions))
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(pre_minimization_selector.apply(pose))
    for jump_edge in jump_edges:
        move_map.set_jump(jump_edge, True)
    pre_minimizer = MinMover()
    pre_minimizer.min_type("lbfgs_armijo_nonmonotone")
    pre_minimizer.movemap(move_map)
    return pre_minimizer

def create_task_factory(static_positions:set=set(), point_mutations:set=set(), \
        design_positions:set=set(), theozyme_positions:set=set(), \
        design_binding_site:bool=False, repack_neighborhood_only:bool=False, \
        repack_binding_site:bool=False, ddG_ref_pose=None, n_monomers:int=1, \
        truncate_sequence_end:int=0, excluded_amino_acid_types:str=None, \
        noncanonical_amino_acids:list=None):
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
    task_factory.push_back(ExtraRotamers(0, 1, 1))
    task_factory.push_back(ExtraRotamers(0, 2, 1))

    # Site-directed AA substitutions positions.
    site_directed_substitution_positions = set()

    # Specify the point mutations.
    mutation_positions = set()
    for point_mutation in point_mutations:
        mutating_position, target_aa = point_mutation.split(",")
        if not mutating_position in static_positions:
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(target_aa)
            task_factory.push_back(OperateOnResidueSubset(restriction, \
                    ResidueIndexSelector(mutating_position)))
            mutation_positions.add(mutating_position)
    site_directed_substitution_positions.update(mutation_positions)

    # Select the design positions.
    design_selection = OrResidueSelector()
    # Specify the site-directed design positions.
    if len(static_positions) > 0:
        design_positions = design_positions - static_positions
    if len(mutation_positions) > 0:
        design_positions = design_positions - mutation_positions
    if len(design_positions) > 0:
        site_directed_design_selection = ResidueIndexSelector(",".join(design_positions))
        design_selection.add_residue_selector(site_directed_design_selection)
    site_directed_substitution_positions.update(design_positions)

    # Site-directed AA substitutions positions.
    substitution_selection = ResidueIndexSelector(",".join(site_directed_substitution_positions) + ",")
    # Exclude the AA substitutions and static positions from repacking.
    substitution_static_selection = ResidueIndexSelector(",".join(\
            site_directed_substitution_positions.union(static_positions)) + ",")

    # Design the theozyme-protein interface.
    theozyme_selection = ResidueIndexSelector(",".join(theozyme_positions) + ",")
    if design_binding_site and len(theozyme_positions) > 0:
        binding_site_selection = select_neighborhood_region(theozyme_selection, False)
        binding_site_selection = AndResidueSelector(binding_site_selection, \
                NotResidueSelector(substitution_static_selection))
        if ddG_ref_pose is not None: # Fix the design positions.
            binding_site_positions = boolean_vector_to_indices_set(\
                    binding_site_selection.apply(ddG_ref_pose), \
                    n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)
            binding_site_selection = ResidueIndexSelector(",".join(binding_site_positions) + ",")
        design_selection.add_residue_selector(binding_site_selection)
        substitution_selection = OrResidueSelector(substitution_selection, binding_site_selection)
        substitution_static_selection = OrResidueSelector(substitution_static_selection, \
                binding_site_selection)

    # Exclude some AA types if specified.
    if excluded_amino_acid_types:
        all_AAs = set("AGILPVFWYDERHKSTMNQ")
        excluded_AAs = set(excluded_amino_acid_types)
        restriction = RestrictAbsentCanonicalAASRLT(",".join(\
                available_AA for available_AA in all_AAs - excluded_AAs))
        restriction.aas_to_keep()
        task_factory.push_back(OperateOnResidueSubset(restriction, design_selection))

    # Identify the repack and static regions.
    if len(mutation_positions) > 0 or len(design_positions) > 0 or len(theozyme_positions) > 0:
        if repack_neighborhood_only:
            # Repack the neighborhood region of the AA substitutions and theozyme positions.
            if repack_binding_site:
                substitution_repacking_selection = select_neighborhood_region(\
                        OrResidueSelector(substitution_selection, theozyme_selection), True)
            else:
                substitution_repacking_selection = select_neighborhood_region(substitution_selection, True)
                substitution_repacking_selection = OrResidueSelector(\
                        substitution_repacking_selection, theozyme_selection)
            if ddG_ref_pose is not None: # Fix the repacking positions.
                substitution_repacking_positions = boolean_vector_to_indices_set(\
                        substitution_repacking_selection.apply(ddG_ref_pose), \
                        n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)
                substitution_repacking_selection = ResidueIndexSelector(\
                        ",".join(substitution_repacking_positions))
            # Exclude the AA substitutions and static positions from repacking.
            repacking_selection = AndResidueSelector(substitution_repacking_selection, \
                    NotResidueSelector(substitution_static_selection))
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), repacking_selection))
            # No repacking region.
            task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), \
                    OrResidueSelector(substitution_selection, repacking_selection), True))
        else:
            # Exclude the AA substitutions and static positions from repacking.
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), \
                    substitution_static_selection, True))
            substitution_repacking_selection = None
    else:
        if repack_neighborhood_only:
            # Do not repack any part of the pose.
            task_factory.push_back(PreventRepacking())
            substitution_repacking_selection = False
        else:
            # Exclude the static positions from the repacking region.
            static_selection = ResidueIndexSelector(",".join(static_positions) + ",")
            task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), static_selection, True))
            substitution_repacking_selection = None
    return task_factory, substitution_repacking_selection

def create_move_map(pose, focus_selection=None, minimize_neighborhood_only:bool=False, \
        static_positions:set=set(), ddG_ref_pose=None, n_monomers=1, \
        truncate_sequence_end=0, jump_edges:set=set()):
    '''
    When "minimize_neighborhood_only" is True, only the "focus_selection" along with 
    its neighborhood region are subject to minimization except the static positions.
    Set a residue to static does not mean not to repack and minimize the residues 
    around it, especially when the static residue is included in the theozyme positions.
    '''
    move_map = MoveMap()
    static_selection = ResidueIndexSelector(",".join(static_positions) + ",")
    if minimize_neighborhood_only:
        if focus_selection:
            minimization_selection = select_neighborhood_region(focus_selection, True)
            if len(static_positions) > 0:
                minimization_selection = AndResidueSelector(minimization_selection, \
                        NotResidueSelector(static_selection))
            if ddG_ref_pose is not None: # Fix the minimization positions.
                minimization_positions = boolean_vector_to_indices_set(\
                        minimization_selection.apply(ddG_ref_pose), \
                        n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)
                minimization_selection = ResidueIndexSelector(",".join(minimization_positions) + ",")
                minimization_vector = minimization_selection.apply(pose)
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
            if ddG_ref_pose is not None: # Fix the minimization positions.
                minimization_positions = boolean_vector_to_indices_set(\
                        minimization_selection.apply(ddG_ref_pose), \
                        n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)
                minimization_selection = ResidueIndexSelector(",".join(minimization_positions) + ",")
                minimization_vector = minimization_selection.apply(pose)
                move_map.set_bb(minimization_vector)
                move_map.set_chi(minimization_vector)
            else: # Pass residue selectors to a movemap factory. Under development.
                minimization_vector = minimization_selection.apply(pose)
                move_map.set_bb(minimization_vector)
                move_map.set_chi(minimization_vector)
        else:
            move_map.set_bb(True)
            move_map.set_chi(True)
    for jump_edge in jump_edges:
        move_map.set_jump(jump_edge, True)
    return move_map

def create_fast_relax_mover(score_function, task_factory, move_map=None):
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    if move_map:
        fast_relax.set_movemap(move_map)
    return fast_relax

def calculate_energy(score_function, pose, selector=None, score_type:str=None):
    # Create the metric
    metric = TotalEnergyMetric()
    metric.set_scorefunction(score_function)
    if score_type:
        exec("metric.set_scoretype(ScoreType.{})".format(score_type))
    # Add the selector
    if selector:
        metric.set_residue_selector(selector)
    return metric.calculate(pose)

def calculate_energy(score_function, pose, selector=None, score_type:str=None):
    # Create the metric
    metric = TotalEnergyMetric()
    metric.set_scorefunction(score_function)
    if score_type:
        exec("metric.set_scoretype(ScoreType.{})".format(score_type))
    # Add the selector
    if selector:
        metric.set_residue_selector(selector)
    return metric.calculate(pose)

def read_energies_from_pdb(pdb_path, substrate_identities):
    dG_substrate = None
    with open(pdb_path, "r") as pf:
        for line in pf:
            if line.startswith("label "):
                score_terms = line[:-1].split(" ")
            elif line.startswith("pose "):
                scores = line[:-1].split(" ")
                dG_total = float(scores[score_terms.index("total")]) - \
                        float(scores[score_terms.index("coordinate_constraint")])
            elif line[:3] in substrate_identities:
                dG_substrate += float(line[:-1].split(" ")[-1])
    return dG_total, dG_substrate

def run_jobs(pose, score_function, *movers, n_decoys=5, theozyme_positions:set=set(), \
        output_filename=None):
    if not output_filename:
        output_filename = pose.pdb_info().name().split("/")[-1][:-4]
    finished_decoys = 0
    already_finished = False
    if os.path.isfile(output_filename + ".pdb"):
        already_finished = True
        finished_decoys = n_decoys
        dG_total, dG_substrate = read_energies_from_pdb(output_filename + ".pdb")
    else:
        for decoy in range(1, n_decoys + 1):
            checkpoint = output_filename + "." + str(decoy) + ".in_progress.pdb"
            if os.path.isfile(checkpoint):
                finished_decoys = decoy
                dG_total, dG_substrate = read_energies_from_pdb(checkpoint)
    for decoy in range(finished_decoys + 1, n_decoys + 1):
        pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        if decoy == 1:
            best_decoy = pose_copy
            best_score = score_function(pose_copy)
            pose_copy.dump_pdb(output_filename + ".1.in_progress.pdb")
            dG_total = None
        else:
            current_score = score_function(pose_copy)
            if current_score < best_score:
                best_decoy = pose_copy
                best_score = current_score
                pose_copy.dump_pdb(output_filename + "." + str(decoy) + ".in_progress.pdb")
                os.remove(output_filename + "." + str(decoy - 1) + ".in_progress.pdb")
                dG_total = None
            else:
                os.rename(output_filename + "." + str(decoy - 1) + ".in_progress.pdb", \
                        output_filename + "." + str(decoy) + ".in_progress.pdb")
        if dG_total == None:
            dG_total = calculate_energy(score_function, best_decoy) - \
                    calculate_energy(score_function, best_decoy, score_type="coordinate_constraint")
            if len(theozyme_positions) > 0:
                substrate_selector = ResidueIndexSelector(",".join(theozyme_positions))
                dG_substrate = calculate_energy(score_function, best_decoy, selector=substrate_selector) - \
                        calculate_energy(score_function, best_decoy, selector=substrate_selector, score_type="coordinate_constraint")
    if not already_finished:
        os.rename(output_filename + "." + str(decoy) + ".in_progress.pdb", output_filename + ".pdb")
    if len(theozyme_positions) > 0:
        with open(output_filename + ".dat", "w") as pf:
            pf.write(str(dG_total) + "\n" + str(dG_substrate) + "\n")
    else:
        with open(output_filename + ".dat", "w") as pf:
            pf.write(str(dG_total) + "\n")

def run_job_distributor(pose, score_function, *movers, n_decoys=50, \
        output_filename_prefix:str=None, annotated_name:str=None):
    if not output_filename_prefix:
        output_filename_prefix = pose.pdb_info().name().split("/")[-1][:-4]
    job_distributor = PyJobDistributor(output_filename_prefix, n_decoys, score_function)
    while not job_distributor.job_complete:
        pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        job_distributor.output_decoy(pose_copy)
        if annotated_name:
            mutations_str = str()
            for pose_index, wt_aa in pose.sequence():
                aa = pose_copy.sequence()[pose_index]
                if aa != wt_aa:
                    mutations_str += "_" + wt_aa + str(pose_index) + aa
            current_id = job_distributor.current_id
            os.rename(job_distributor.pdb_name + "_" + str(current_id) + ".pdb", \
                    job_distributor.pdb_name + mutations_str + "_" + str(current_id) + ".pdb")

def main(args):
    # Create the score function.
    score_function = set_score_function(args.score_function, symmetry=args.symmetry, \
            score_terms=args.score_terms)
    # Load pdb as pose.
    pose, coord_ref_pose, ddG_ref_pose, sequence_length, truncate_sequence_end = \
            load_pdb_as_pose(args.pdb, args.coordinate_reference_pdb, args.ddG_reference_pdb)
    # Group fragments in the fold tree.
    if args.fold_tree:
        pose.fold_tree(create_fold_tree(args.fold_tree))
    elif args.alter_jump_edges:
        pose.fold_tree(alter_fold_tree_jump_edges(pose.fold_tree(), args.alter_jump_edges))
    # Set chi dihedrals.
    set_chi_dihedral(pose, args.chi_dihedrals)
    # Read constraint files from the command line.
    if args.constraints:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Add constraints on-the-fly.
    constraints = list()
    if args.distance_constraint_atoms:
        constraints.extend(create_constraints(pose, args.distance_constraint_atoms, \
                args.distance_constraint_parameters, boundaries=args.distance_constraint_boundaries))
    if args.angle_constraint_atoms:
        constraints.extend(create_constraints(pose, args.angle_constraint_atoms, \
                args.angle_constraint_parameters))
    if args.dihedral_constraint_atoms:
        constraints.extend(create_constraints(pose, args.dihedral_constraint_atoms, \
                args.dihedral_constraint_parameters))
    for constraint in constraints:
        pose.add_constraint(constraint)
    # Add coordinate constraints.
    if len(args.no_coordinate_constraint_residues) > 0:
        no_coord_cst_residues, _ = pdb_to_pose_numbering(pose, args.no_coordinate_constraint_residues)
        no_coord_cst_selection = ResidueIndexSelector(",".join(no_coord_cst_residues))
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
    # Favor native AA types.
    if args.favor_native_residue and (args.design_residues or args.design_binding_site):
        favor_nataa = FavorNativeResidue(pose, args.favor_native_residue)
        favor_nataa.apply(pose)
    # Create the RMSD metric.
    no_rmsd_residues, _ = pdb_to_pose_numbering(pose, args.no_rmsd_residues)
    rmsd_metric = RMSDMetric(pose, NotResidueSelector(ResidueIndexSelector(\
                ",".join(no_rmsd_residues) + ",")))
    # Convert pdb numberings to pose numberings.
    static_pose_indices, _ = pdb_to_pose_numbering(pose, args.static_residues)
    if args.static_residue_identities:
        static_pose_indices.update(residue_name3_selector(pose, args.static_residue_identities))
    mutation_pose_indices, _ = pdb_to_pose_numbering(pose, args.mutations)
    design_pose_indices, _ = pdb_to_pose_numbering(pose, args.design_residues)
    # Get pose indices of substrates and catalytic residues.
    theozyme_pose_indices = set()
    rigid_body_tform_pose_indices = set()
    if args.catalytic_residues:
        catalytic_residue_pose_indices, _ = pdb_to_pose_numbering(pose, args.catalytic_residues)
        theozyme_pose_indices.update(catalytic_residue_pose_indices)
    if args.substrates:
        substrate_pose_indices, jump_pose_indices = pdb_to_pose_numbering(pose, args.substrates)
        theozyme_pose_indices.update(substrate_pose_indices)
        rigid_body_tform_pose_indices.update(jump_pose_indices)
    if args.enzyme_design_constraints:
        pdb_info = pose.pdb_info()
        match_substrate_pose_indices, match_res_pose_indices = \
                get_match_pose_indices(pdb_info, args.pdb, args.symmetry)
        theozyme_pose_indices.update(match_substrate_pose_indices.union(match_res_pose_indices))
        if args.enzyme_design_constraints_substrate_transformations:
            rigid_body_tform_pose_indices.update(match_substrate_pose_indices)
    if args.substrate_identities:
        substrate_by_id_pose_indices = residue_name3_selector(pose, args.substrate_identities)
        theozyme_pose_indices.update(substrate_by_id_pose_indices)
        nonredundant_theozyme_pose_indices = theozyme_pose_indices.copy()
        rigid_body_tform_pose_indices.update(substrate_by_id_pose_indices)
        if args.ddG_reference_pdb:
            theozyme_pose_indices.update(residue_name3_selector(ddG_ref_pose, args.substrate_identities))
    # Verify the redundant length of ddG_ref_pose.
    assert len(theozyme_pose_indices) - len(nonredundant_theozyme_pose_indices) == truncate_sequence_end
    # Rigid body transformations.
    rigid_body_tform_jump_edges = set()
    if args.substrate_rigid_body_transformations:
        rigid_body_tform_pose_indices = rigid_body_tform_pose_indices - static_pose_indices
        rigid_body_tform_jump_edges = pose_indices_to_jump_edges(pose.fold_tree(), \
                rigid_body_tform_pose_indices)
    # Applying the symmetry if specified.
    n_monomers = set_symmetry(args.symmetry, sequence_length, pose, coord_ref_pose, ddG_ref_pose)
    # Perform pre-minimization.
    minimization_theozyme_pose_indices = nonredundant_theozyme_pose_indices - static_pose_indices
    if args.pre_minimization and len(minimization_theozyme_pose_indices) > 0:
        pre_minimizer = pre_minimization(pose, minimization_theozyme_pose_indices, \
                jump_edges=rigid_body_tform_jump_edges)
        pre_minimizer.apply(pose)
    # Create the task factory.
    task_factory, substitution_repacking_selection = create_task_factory(\
            static_positions=static_pose_indices, point_mutations=mutation_pose_indices, \
            design_positions=design_pose_indices, theozyme_positions=theozyme_pose_indices, \
            design_binding_site=args.design_binding_site, \
            repack_neighborhood_only=args.repack_neighborhood_only, \
            repack_binding_site=args.repack_binding_site, \
            ddG_ref_pose=ddG_ref_pose, n_monomers=n_monomers, \
            truncate_sequence_end=truncate_sequence_end, \
            excluded_amino_acid_types=args.excluded_amino_acid_types, \
            noncanonical_amino_acids=args.noncanonical_amino_acids)
    # Create the move map.
    move_map = create_move_map(pose, focus_selection=substitution_repacking_selection, \
            minimize_neighborhood_only=args.minimize_neighborhood_only, \
            static_positions=static_pose_indices, ddG_ref_pose=ddG_ref_pose, \
            n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end, \
            jump_edges=rigid_body_tform_jump_edges)
    # Create the fast relax mover.
    fr = create_fast_relax_mover(score_function, task_factory, move_map=move_map)
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
        # Unload old coordinate constraints and add new coordinate constraints.
        coord_cst_gen = create_coord_cst(coord_ref_pose=coord_ref_pose, \
                coord_boundary=args.coordinate_constraint_bounded_width, \
                no_coord_cst_selection=no_coord_cst_selection)
        add_csts = AddConstraints()
        add_csts.add_generator(coord_cst_gen)
        add_csts.apply(pose)
    # Run jobs.
    if args.debug_mode:
        print(pose.fold_tree())
        if args.enzyme_design_constraints:
            print(match_substrate_pose_indices)
            print(match_res_pose_indices)
        print(task_factory.create_task_and_apply_taskoperations(pose))
        print(move_map)
        print(score_function.show(pose))
    elif args.no_save_decoys:
        run_jobs(pose, score_function, fr, rmsd_metric, n_decoys=5, \
                theozyme_positions=theozyme_pose_indices, \
                output_filename=args.output_filename_prefix)
    else:
        run_job_distributor(pose, score_function, fr, rmsd_metric, \
                n_decoys=args.n_decoys, output_filename_prefix=args.output_filename_prefix, \
                annotated_name=args.annotated_name)


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    main(args)
