#!/usr/bin/python3
import argparse
import json
import numpy as np
import os
from pyrosetta import init, create_score_function, FoldTree, \
    MoveMap, pose_from_pdb, Pose, PyJobDistributor
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.palette import CustomBaseTypePackerPalette
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictToRepackingRLT, RestrictAbsentCanonicalAASRLT, \
    PreventRepackingRLT
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import \
    add_fa_constraints_from_cmdline, ConstraintSet, AmbiguousConstraint, \
    AtomPairConstraint, AngleConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc, \
    FlatHarmonicFunc, CircularHarmonicFunc
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.jump_selector import JumpIndexSelector
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, \
    move_map_action
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    TrueResidueSelector, ResidueIndexSelector, \
    InterGroupInterfaceByVectorSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    RMSDMetric, TotalEnergyMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, \
    AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.protein_interface_design import \
    FavorNativeResidue
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.utility import vector1_bool


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=str)
    parser.add_argument("-cloud", "--cloud_pdb", action="store_true")
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
    parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst", \
            help="$base_$soft_$cart_$cst")
    parser.add_argument("--score_terms", type=str, nargs="*", default=list())
    parser.add_argument("-coord_cst_sd", "--coordinate_constraints_standard_deviation", type=float, \
            default=0.5)
    parser.add_argument("-bounded_coord_cst", "--bounded_coordinate_constraints", type=float)
    parser.add_argument("-ca_coord_cst", "--only_CA_coordinate_constraints", action="store_true")
    parser.add_argument("-all_coord_cst", "--all_atom_coordinate_constraints_positions", \
            type=str, nargs="*", default=list())
    parser.add_argument("-no_coord_cst", "--no_backbone_coordinate_constraints_residues", \
            type=str, nargs="*", default=list())
    parser.add_argument("-no_pack_coord_cst", "--no_coordinate_constraints_on_packing_region", \
            action="store_true")
    parser.add_argument("-cst", "--constraint_file", type=str)
    parser.add_argument("-dis_atoms", "--distance_constraint_atoms", type=str, nargs="*", \
            help="$chain$residue_index,$atom_name or $residue_name3,$atom_name * 2n")
    parser.add_argument("-dis_params", "--distance_constraint_parameters", type=str, nargs="*", \
            help="$distance,$standard_deviation or $dis,$sd,$dis,$sd,... * n")
    parser.add_argument("-bounded_dis", "--bounded_distances", type=str, nargs="*", \
            help="$distance or $distance,$distance,... or placeholder * n")
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
    parser.add_argument("-cat_ids", "--catalytic_residue_identities", type=str, nargs="*", help="name3")
    parser.add_argument("-subs", "--substrates", type=str, nargs="*")
    parser.add_argument("-sub_ids", "--substrate_identities", type=str, nargs="*", help="name3")
    parser.add_argument("-premin", "--pre_minimization", action="store_true")
    parser.add_argument("-ddg_wt", "--ddG_wildtype", action="store_true", \
            help="Keep all -muts, -des, -des_bs and -des_enzdes positions as wildtype.")
    parser.add_argument("-muts", "--mutations", type=str, nargs="*", default=list(), 
            help="Site-directed AA substitution. Highests priority.")
    parser.add_argument("-des", "--design_residues", type=str, nargs="*", default=list(), \
            help="Site-directed AA design. Lower priority than -muts.")
    parser.add_argument("-des_bs", "--design_binding_site", action="store_true", \
            help="Substrate-binding site AA design. Lower priority than -muts and -des.")
    parser.add_argument("-des_enzdes", "--design_enzdes_shell", action="store_true", \
            help="Enzdes shell AA design. Lower priority than -muts and -des.")
    parser.add_argument("-nataa", "--favor_native_residue", type=float)
    parser.add_argument("-noaa", "--excluded_amino_acid_types", type=str, \
            help="Concatenated string of excluded AA one-letter codes.")
    parser.add_argument("-ncaa", "--noncanonical_amino_acids", type=str, nargs="*", \
            default=list(), help="$name1,$name3")
    parser.add_argument("-ncaa_des", "--allow_ncaa_in_design", action="store_true")
    parser.add_argument("-rpk_nbh", "--repack_neighborhood_only", action="store_true")
    parser.add_argument("-no_rpk_bs", "--no_repack_binding_site", action="store_true")
    parser.add_argument("-no_rpk_enzdes", "--no_repack_enzdes_shell", action="store_true")
    parser.add_argument("-min_nbh", "--minimize_neighborhood_only", action="store_true")
    parser.add_argument("-tform", "--substrate_rigid_body_transformations", action="store_true")
    parser.add_argument("-tform_enzdes", "--enzdes_substrates_transformations", \
            action="store_true")
    parser.add_argument("-no_rmsd", "--no_rmsd_residues", type=str, nargs="*", default=list())
    parser.add_argument("-n", "--n_decoys", type=int, default=50)
    parser.add_argument("-prefix", "--output_filename_prefix", type=str)
    parser.add_argument("-suffix", "--output_filename_mutations_suffix", action="store_true")
    parser.add_argument("-save_n", "--save_n_decoys", type=int)
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
    if args.constraint_file:
        opts += " -constraints:cst_fa_file {}".format(args.constraint_file)
    init(opts)

def set_score_function(score_function_name:str, symmetry:bool=False, score_terms=list()):
    """
    Custimize an appropriate score function on top of a L-J potential that is 
    either hard or soft, and having symmetry, cartesian and/or constraints settings.
    """
    membrane = False
    if score_function_name.startswith("franklin2019") and ("_soft" in score_function_name or \
            "_cart" in score_function_name or "_cst" in score_function_name):
        score_function_name.replace("franklin2019", "ref2015")
        membrane = True
    cartesian = False
    constraints = False
    if "_soft" in score_function_name:
        if "_cart" in score_function_name:
            score_function_name.replace("_cart", str())
            cartesian = True
        if "_cst" in score_function_name:
            score_function_name.replace("_cst", str())
            constraints = True
    elif score_function_name.startswith("franklin2019") and \
            "_cart" in score_function_name and "_cst" in score_function_name:
        score_function_name.replace("_cst", str())
        constraints = True
    if symmetry:
        score_function = SymmetricScoreFunction()
        score_function.add_weights_from_file(score_function_name)
    else:
        score_function = create_score_function(score_function_name)
    if membrane:
        score_function.set_weight(ScoreType.fa_water_to_bilayer, 1)
    if cartesian:
        score_function.set_weight(ScoreType.cart_bonded, 0.5)
        # score_function.set_weight(ScoreType.pro_close, 0)
        # score_function.set_weight(ScoreType.metalbinding_constraint, 0)
        # score_function.set_weight(ScoreType.chainbreak, 0)
    if constraints:
        if not cartesian:
            score_function.set_weight(ScoreType.pro_close, 1.25)
            score_function.set_weight(ScoreType.metalbinding_constraint, 1)
            score_function.set_weight(ScoreType.chainbreak, 1)
        score_function.set_weight(ScoreType.atom_pair_constraint, 1)
        score_function.set_weight(ScoreType.coordinate_constraint, 1)
        score_function.set_weight(ScoreType.angle_constraint, 1)
        score_function.set_weight(ScoreType.dihedral_constraint, 1)
        score_function.set_weight(ScoreType.res_type_constraint, 1)
    for score_term in score_terms:
        term_weight = score_term.split(":")
        exec("score_function.set_weight(ScoreType.{}, {})".format(term_weight[0], term_weight[1]))
    return score_function

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

def residue_name3_selector(pose, name3_list, sequence_length:int=None):
    pose_indices = set()
    for pose_index in range(1, len(pose.sequence()) + 1):
        if pose.residue(pose_index).name3() in name3_list:
            if sequence_length and pose_index > sequence_length:
                continue
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
    new_jump_edge_label = 0
    for edge_str in fold_tree.to_string().split("EDGE")[1:]:
        edge_info = tuple(filter(lambda x: x != "", edge_str.strip(" ").split(" ")))
        downstream_residue = int(edge_info[1])
        edge = fold_tree.get_residue_edge(downstream_residue)
        if edge.is_jump():
            upstream_residue = int(edge_info[0])
            edge_label = edge.label()
            assert edge_label == int(edge_info[2])
            new_jump_edge_label += 1
            alter_jump_edge = alter_jump_edges_dict.get(edge_label)
            if alter_jump_edge:
                upstream_edge_label = alter_jump_edge[0]
                if upstream_edge_label < 0:
                    upstream_edge_label = jump_edge_labels[upstream_edge_label]
                if upstream_edge_label > 0:
                    upstream_residue = fold_tree.jump_edge(upstream_edge_label).stop()
                if len(alter_jump_edge) == 1:
                    edge_info = (str(upstream_residue), str(downstream_residue), str(new_jump_edge_label))
                elif len(alter_jump_edge) == 3:
                    edge_info = (str(upstream_residue), str(downstream_residue), \
                            alter_jump_edge[1], alter_jump_edge[2])
                    new_jump_edge_label -= 1
            elif new_jump_edge_label < edge_label:
                edge_info = (str(upstream_residue), str(downstream_residue), str(new_jump_edge_label))
        edge_strings.append(",".join(edge_info))
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

def create_coordinate_constraints(reference_pose=None, selection=None, \
        standard_deviation:float=0.5, bounded:float=None, ca_only:bool=False, \
        side_chain:bool=False):
    coord_cst_gen = CoordinateConstraintGenerator()
    if reference_pose is not None:
        coord_cst_gen.set_reference_pose(reference_pose)
    if selection is not None:
        coord_cst_gen.set_residue_selector(selection)
    coord_cst_gen.set_sd(standard_deviation)
    if bounded:
        coord_cst_gen.set_bounded(True)
        coord_cst_gen.set_bounded_width(bounded)
    coord_cst_gen.set_ca_only(ca_only)
    coord_cst_gen.set_sidechain(side_chain)
    return coord_cst_gen

def create_harmonic_constraint(pose, atom_names, pose_indices, distance, standard_deviation, \
        bounded_distance:float=None):
    atom_list = list()
    for atom_name, residue_id in zip(atom_names, pose_indices):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    if type(distance) is list or type(distance) is np.ndarray:
        assert len(distance) == len(standard_deviation)
        harmonic_cst = AmbiguousConstraint()
        for i_nested, dis_sd in enumerate(zip(distance, standard_deviation)):
            if bounded_distance:
                harmonic_fc = FlatHarmonicFunc(dis_sd[0], dis_sd[1], bounded_distance[i_nested])
            else:
                harmonic_fc = HarmonicFunc(dis_sd[0], dis_sd[1])
            cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
            harmonic_cst.add_individual_constraint(cst)
    else:
        if bounded_distance:
            harmonic_fc = FlatHarmonicFunc(distance, standard_deviation, bounded_distance)
        else:
            harmonic_fc = HarmonicFunc(distance, standard_deviation)
        harmonic_cst = AtomPairConstraint(atom_list[0], atom_list[1], harmonic_fc)
    return harmonic_cst

def create_circular_harmonic_constraint(pose, atom_names, pose_indices, degrees, standard_deviation):
    atom_list = list()
    for atom_name, residue_id in zip(atom_names, pose_indices):
        atom_list.append(AtomID(pose.residue(residue_id).atom_index(atom_name), residue_id))
    radians = np.radians(degrees)
    standard_deviation = np.radians(standard_deviation)
    if type(radians) is np.ndarray:
        assert len(radians) == len(standard_deviation)
        circular_harmonic_cst = AmbiguousConstraint()
        for rad, sd in zip(radians, standard_deviation):
            circular_harmonic_fc = CircularHarmonicFunc(rad, sd)
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
        circular_harmonic_fc = CircularHarmonicFunc(radians, standard_deviation)
        if len(atom_list) == 3:
            circular_harmonic_cst = AngleConstraint(atom_list[0], atom_list[1], \
                    atom_list[2], circular_harmonic_fc)
        elif len(atom_list) == 4:
            circular_harmonic_cst = DihedralConstraint(atom_list[0], atom_list[1], \
                    atom_list[2], atom_list[3], circular_harmonic_fc)
    return circular_harmonic_cst

def create_constraints(pose, constraint_atoms, constraint_parameters, bounded_distances=None):
    assert len(constraint_atoms) % len(constraint_parameters) == 0
    constraints = list()
    length = int(len(constraint_atoms) / len(constraint_parameters))
    if bounded_distances:
        assert len(bounded_distances) == len(constraint_parameters)
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
        # Get bounded distances.
        if bounded_distances:
            bounded_distance = bounded_distances[cst_index]
            if "," in bounded_distance:
                bounded_distance = list()
                for bounded_dis_i in bounded_distance.split(","):
                    try:
                        bounded_dis_i = float(bounded_dis_i)
                    except:
                        bounded_dis_i = None
                bounded_distance.append(bounded_dis_i)
            try:
                bounded_distance = float(bounded_distance)
            except:
                bounded_distance = None
        else:
            bounded_distance = None
        # Get atoms.
        current_cst_atoms = constraint_atoms[cst_first_atom:cst_first_atom + length]
        atom_names = list()
        try:
            pose_index_atom_name_set, _ = pdb_to_pose_numbering(pose, current_cst_atoms)
            pose_indices = list()
            for pose_index_atom_name in pose_index_atom_name_set:
                pose_index, atom_name = pose_index_atom_name.split(",")
                pose_indices.append(int(pose_index))
                atom_names.append(atom_name)
            pose_indices_list = [pose_indices]
        except: # Select possibly multiple residues by residue identity.
            pose_indices_list = list()
            for cst_atom in current_cst_atoms:
                residue_name3, atom_name = cst_atom.split(",")
                atom_names.append(atom_name)
                new_pose_indices = residue_name3_selector(pose, [residue_name3])
                if len(pose_indices_list) == 0:
                    pose_indices_list = list([int(new_pose_index)] for new_pose_index in new_pose_indices)
                else:
                    new_pose_indices_list = list()
                    for pose_indices in pose_indices_list:
                        for new_pose_index in new_pose_indices:
                            new_pose_indices_list.append(pose_indices + [int(new_pose_index)])
                    pose_indices_list = new_pose_indices_list
        finally:
            for pose_indices in pose_indices_list:
                if len(atom_names) == 2:
                    constraint = create_harmonic_constraint(pose, atom_names, pose_indices, value, \
                            standard_deviation, bounded_distance=bounded_distance)
                else:
                    constraint = create_circular_harmonic_constraint(pose, atom_names, pose_indices, \
                            value, standard_deviation)
                constraints.append(constraint)
    return constraints

def create_enzdes_constraints():
    enz_cst = AddOrRemoveMatchCsts()
    enz_cst.set_cst_action(ADD_NEW)
    return enz_cst

def get_enzdes_pose_indices(pdb_info, pdb, symmetry):
    '''
    Get pose indices of enzdes positions by reading the pdb file REMARK 666 header.
    '''
    enzdes_substrate_pose_indices = set()
    enzdes_res_pose_indices = set()
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
                    substrate_pose_id = pdb_info.pdb2pose(substrate_chain_id, substrate_res_id)
                    enzdes_substrate_pose_indices.add(str(substrate_pose_id))
                # enzdes residues
                motif_chain_id = line[49]
                if main_chain == motif_chain_id or not symmetry:
                    motif_res_id = int(line[55:59])
                    motif_pose_id = pdb_info.pdb2pose(motif_chain_id, motif_res_id)
                    enzdes_res_pose_indices.add(str(motif_pose_id))
                flag = True
            elif flag:
                break
    return enzdes_substrate_pose_indices, enzdes_res_pose_indices

def pose_indices_to_jump_edges(fold_tree, rigid_body_tform_pose_indices):
    jump_edges = set()
    for rigid_body_tform_pose_index in rigid_body_tform_pose_indices:
        jump_edge = fold_tree.get_residue_edge(int(rigid_body_tform_pose_index))
        if not jump_edge.is_jump():
            raise Exception("The specified edge is not a jump edge.")
        jump_edges.add(jump_edge.label())
    return jump_edges

def index_boolean_filter(index, boolean):
    if type(boolean) == np.ndarray:
        for bool in boolean:
            if bool:
                return str(index)
    elif boolean:
        return str(index)
    return False

def boolean_vector_to_indices_set(boolean_vector, n_monomers:int=1, \
        truncate_sequence_end:int=0, output_vector:bool=False):
    boolean_vector = np.array(boolean_vector)
    sequence_length = len(boolean_vector) // n_monomers
    boolean_matrix = boolean_vector[:n_monomers*sequence_length]\
            .reshape(n_monomers, sequence_length).transpose()
    if truncate_sequence_end > 0:
        boolean_matrix = boolean_matrix[:-truncate_sequence_end,:]
    indices = set(filter(lambda x: x, map(index_boolean_filter, \
            range(1, len(boolean_matrix) + 1), boolean_matrix)))
    if output_vector:
        boolean_vector = vector1_bool(len(boolean_matrix))
        for index in indices:
            boolean_vector[int(index)] = True
        return boolean_vector
    else:
        return indices

def set_symmetry(symmetry:str, *pose_list, sequence_length:int=0):
    if symmetry:
        sfsm = SetupForSymmetryMover(symmetry)
        sfsm.set_keep_pdb_info_labels(True)
        n_monomers = None
        for pose in filter(lambda p: p is not None, pose_list):
            sfsm.apply(pose)
            if not n_monomers and sequence_length > 0:
                n_monomers = round(len(pose.sequence().strip("X")) / sequence_length)
    else:
        n_monomers = 1
    return n_monomers

def select_neighborhood_region(focus_selection, include_focus: bool, method:str="vector"):
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

def create_pre_minimizer(score_function, assembly_length, pre_minimization_pose_indices, \
        jump_edges:set=set()):
    pre_minimization_vector = vector1_bool(assembly_length)
    for pose_index in pre_minimization_pose_indices:
        pre_minimization_vector[int(pose_index)] = True
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(pre_minimization_vector)
    for jump_edge in jump_edges:
        move_map.set_jump(jump_edge, True)
    pre_minimizer = MinMover()
    pre_minimizer.score_function(score_function)
    pre_minimizer.min_type("lbfgs_armijo_nonmonotone")
    if score_function.get_weight(ScoreType.cart_bonded) > 0:
        pre_minimizer.cartesian(True)
    pre_minimizer.movemap(move_map)
    return pre_minimizer

def create_task_factory(nopack_positions:set=set(), point_mutations:set=set(), \
        design_positions:set=set(), theozyme_positions:set=set(), enzdes_positions:set=set(), \
        design_binding_site:bool=False, design_enzdes_shell:bool=False, \
        repack_neighborhood_only:bool=False, repack_binding_site:bool=True, \
        repack_enzdes_shell:bool=True, ddG_ref_pose=None, ddG_wildtype:bool=False, \
        n_monomers:int=1, truncate_sequence_end:int=0, excluded_amino_acid_types:str=None, \
        noncanonical_amino_acids:list=list(), allow_ncaa_in_design:bool=False):
    """
    Priority:
        1. No pack residues.
        2. Introduce point mutations.
        3. Perform site-directed design.
        4. Design protein-ligand interfaces.
    """

    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    task_factory.push_back(ExtraRotamers(0, 1, 1))
    task_factory.push_back(ExtraRotamers(0, 2, 1))

    # Compatible amino acid types.
    AA_1to3_dict = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", \
            "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", \
            "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", \
            "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    all_AAs = set(AA_1to3_dict.keys())
    # Include noncanonical amino acids.
    ncaa_palette = CustomBaseTypePackerPalette()
    for ncaa in noncanonical_amino_acids:
        ncaa_name1, ncaa_name3 = ncaa.split(",")
        AA_1to3_dict[ncaa_name1] = ncaa_name3
        ncaa_palette.add_type(ncaa_name3)
    if allow_ncaa_in_design:
        all_AAs = set(AA_1to3_dict.keys())
        task_factory.set_packer_palette(ncaa_palette)

    # Site-directed AA substitutions positions.
    # Repack its surrounding shell w/o including the focus region itself.
    site_directed_substitution_positions = set()

    # Repack its surrounding shell along with the focus region itself if not static.
    repack_shell_focus_positions = set()

    # Specify the point mutations.
    mutators = list()
    mutation_positions = set()
    site_directed_ncaa_positions = set()
    canonical_AAs = set("ACDEFGHIKLMNPQRSTVWY")
    for point_mutation in point_mutations:
        mutating_position, target_AA = point_mutation.split(",")
        if target_AA not in AA_1to3_dict.keys():
            raise Exception("Specified AA type " + target_AA + " is not declared in -ncaa.")
        point_mutation_selection = ResidueIndexSelector(mutating_position)
        if mutating_position in nopack_positions or not target_AA in canonical_AAs:
            if not ddG_wildtype:
                mutator = MutateResidue()
                mutator.set_selector(point_mutation_selection)
                target_AA_name3 = AA_1to3_dict[target_AA]
                mutator.set_res_name(target_AA_name3)
                mutators.append(mutator)
            repack_shell_focus_positions.add(mutating_position)
            if not mutating_position in nopack_positions:
            # Repack the introduced noncanonical AA.
                site_directed_ncaa_positions.add(mutating_position)
        elif not ddG_wildtype:
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(target_AA)
            task_factory.push_back(OperateOnResidueSubset(restriction, \
                    point_mutation_selection))
        mutation_positions.add(mutating_position)
    site_directed_substitution_positions.update(mutation_positions)

    # Select the design positions.
    design_selection = OrResidueSelector()
    # Specify the site-directed design positions.
    design_positions = design_positions - mutation_positions
    if len(design_positions.intersection(nopack_positions)) > 0:
        raise Exception("Site-directed design positions cannot contain nopack positions.")
    site_directed_design_selection = ResidueIndexSelector(",".join(design_positions) + ",")
    design_selection.add_residue_selector(site_directed_design_selection)
    site_directed_substitution_positions.update(design_positions)

    # Site-directed AA substitutions positions. May include nopack positions or NCAA mutations.
    substitution_selection = ResidueIndexSelector(",".join(site_directed_substitution_positions) + ",")

    # Design its surrounding shell and repack the focus region itself.
    design_shell_focus_positions = set()
    if design_binding_site:
        design_shell_focus_positions.update(theozyme_positions)
    if design_enzdes_shell:
        design_shell_focus_positions.update(enzdes_positions)
    # Repack its surrounding shell along with the focus region itself.
    if repack_binding_site:
        repack_shell_focus_positions.update(theozyme_positions)
    if repack_enzdes_shell:
        repack_shell_focus_positions.update(enzdes_positions)
    # Repack this region in addition to the surrounding shell repack region.
    additional_repacking_positions = theozyme_positions.union(enzdes_positions)
    # The following 3 selections may contain site-directed substitution positions 
    # or nopack positions or redundant ddG_ref_pose positions. Use with caution.
    design_shell_focus_selection = ResidueIndexSelector(",".join(\
            design_shell_focus_positions) + ",")
    repack_shell_focus_selection = ResidueIndexSelector(",".join(\
            repack_shell_focus_positions) + ",")
    additional_repacking_selection = ResidueIndexSelector(",".join(\
            additional_repacking_positions) + ",")

    # Design the theozyme-protein interface.
    design_shell_selection = select_neighborhood_region(design_shell_focus_selection, False)
    # Exclude the site-directed AA substitutions and nopack positions from the design shell.
    design_shell_selection = AndResidueSelector(design_shell_selection, \
            NotResidueSelector(ResidueIndexSelector(",".join(\
            site_directed_substitution_positions.union(nopack_positions)) + ",")))
    if ddG_ref_pose is not None:
        # Fix the design positions and truncate redundant ddG_ref_pose positions.
        binding_site_positions = boolean_vector_to_indices_set(\
                design_shell_selection.apply(ddG_ref_pose), \
                n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)
        design_shell_selection = ResidueIndexSelector(",".join(binding_site_positions) + ",")
    design_selection.add_residue_selector(design_shell_selection)
    substitution_selection = OrResidueSelector(substitution_selection, design_shell_selection)

    # Every position is designable by defalut in the task factory other than specification.
    if excluded_amino_acid_types and not ddG_wildtype:
        # Exclude some AA types if specified.
        excluded_AAs = set(excluded_amino_acid_types)
        restriction = RestrictAbsentCanonicalAASRLT(",".join(all_AAs - excluded_AAs))
        restriction.aas_to_keep()
        task_factory.push_back(OperateOnResidueSubset(restriction, design_selection))

    # Exclude nopack positions from repacking.
    nopack_selection = ResidueIndexSelector(",".join(nopack_positions) + ",")
    # Need to exclude NCAA substitutions from the substitution selection when repacking.
    # In other words, explicitly include NCAA substitutions in repacking.
    caa_substitution_selection = AndResidueSelector(substitution_selection, NotResidueSelector(\
            ResidueIndexSelector(",".join(site_directed_ncaa_positions) + ",")))
    if repack_neighborhood_only:
        # Repack the neighborhood region of the AA substitutions and theozyme positions.
        # min_shell_focus_selection may contain nopack positions or redundant ddG_ref_pose positions.
        min_shell_focus_selection = select_neighborhood_region(OrResidueSelector(\
                substitution_selection, repack_shell_focus_selection), True)
        min_shell_focus_selection = OrResidueSelector(min_shell_focus_selection, \
                additional_repacking_selection)
        if ddG_ref_pose is not None:
            # Fix the repacking positions. # type(min_shell_focus_selection) == set
            min_shell_focus_selection = boolean_vector_to_indices_set(\
                    min_shell_focus_selection.apply(ddG_ref_pose), n_monomers=n_monomers)
            # Filter out redundant ddG_ref_pose positions.
            pack_positions = boolean_vector_to_indices_set(ResidueIndexSelector(\
                    ",".join(min_shell_focus_selection) + ",").apply(ddG_ref_pose), \
                    n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end)\
                    - nopack_positions
            pack_selection = ResidueIndexSelector(",".join(pack_positions) + ",")
        else:
            pack_selection = AndResidueSelector(min_shell_focus_selection, \
                    NotResidueSelector(nopack_selection))
        if not ddG_wildtype:
            # Exclude the canonical AA substitutions and nopack positions from repacking.
            repacking_selection = AndResidueSelector(pack_selection, NotResidueSelector(\
                    caa_substitution_selection))
        else:
            repacking_selection = pack_selection
        task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), repacking_selection))
        # No repacking region.
        task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), pack_selection, True))
    else:
        if not ddG_wildtype:
            # Exclude the canonical AA substitutions and nopack positions from repacking.
            not_repacking_selection = OrResidueSelector(caa_substitution_selection, nopack_selection)
        else:
            not_repacking_selection = nopack_selection
        task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), \
                not_repacking_selection, True))
        task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), nopack_selection))
        min_shell_focus_selection = None

    return task_factory, mutators, min_shell_focus_selection

def create_move_map(focus_selection, minimize_neighborhood_only:bool=False, \
        static_positions:set=set(), ddG_ref_pose=None, n_monomers=1, \
        truncate_sequence_end=0, jump_edges:set=set()):
    '''
    When "minimize_neighborhood_only" is True, only the "focus_selection" along with 
    its neighborhood region are subject to minimization except the static positions.
    Set a residue to static does not mean not to repack and minimize the residues 
    around it, especially when the static residue is included in the theozyme positions.
    '''
    static_selection = ResidueIndexSelector(",".join(static_positions) + ",")
    if minimize_neighborhood_only:
        if focus_selection == None:
            raise Exception("Using -min_nbh is not allowed without using -rpk_nbh at the same time.")
        if type(focus_selection) == set:
            assert ddG_ref_pose is not None
            focus_selection = ResidueIndexSelector(",".join(focus_selection) + ",")
        minimization_selection = select_neighborhood_region(focus_selection, True)
        if len(static_positions) > 0:
            minimization_selection = AndResidueSelector(minimization_selection, \
                    NotResidueSelector(static_selection))
        if ddG_ref_pose is not None:
            # Fix the minimization positions and truncate redundant ddG_ref_pose positions.
            minimization_vector = boolean_vector_to_indices_set(\
                    minimization_selection.apply(ddG_ref_pose), n_monomers=n_monomers, \
                    truncate_sequence_end=truncate_sequence_end, output_vector=True)
            move_map = MoveMap()
            move_map.set_bb(minimization_vector)
            move_map.set_chi(minimization_vector)
        else: # Pass residue selectors to a movemap factory.
            move_map = MoveMapFactory()
            move_map.add_bb_action(move_map_action.mm_enable, minimization_selection)
            move_map.add_chi_action(move_map_action.mm_enable, minimization_selection)
    else:
        minimization_vector = vector1_bool(len(focus_selection))
        for pose_index in range(1, len(focus_selection) + 1):
            minimization_vector[pose_index] = True
        for static_position in static_positions:
            minimization_vector[str(static_position)] = False
        move_map = MoveMap()
        move_map.set_bb(minimization_vector)
        move_map.set_chi(minimization_vector)
    if type(move_map) == MoveMap:
        for jump_edge in jump_edges:
            move_map.set_jump(jump_edge, True)
    else:
        for jump_edge in jump_edges:
            move_map.add_jump_action(move_map_action.mm_enable, JumpIndexSelector(jump_edge))
    return move_map

def create_fast_relax_mover(score_function, task_factory, move_map=None):
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    if type(move_map) == MoveMap:
        fast_relax.set_movemap(move_map)
    elif type(move_map) == MoveMapFactory:
        fast_relax.set_movemap_factory(move_map)
    if score_function.get_weight(ScoreType.cart_bonded) > 0:
        fast_relax.cartesian(True)
    return fast_relax

def load_pdb_as_pose(score_function, pdb:str, fold_tree, chi_dihedrals:list, \
        symmetry:str, constraint_file:str, geometry_constraints, constraints):
    # Load pdb as pose.
    pose = pose_from_pdb(pdb)
    if fold_tree:
        pose.fold_tree(fold_tree)
    set_chi_dihedral(pose, chi_dihedrals)
    # Apply symmetry if specified.
    set_symmetry(symmetry, pose)
    # Read constraint files from the command line and apply to pose.
    if constraint_file:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Apply geometry constraints.
    for geometry_constraint in geometry_constraints:
        pose.add_constraint(geometry_constraint)
    # Apply coordinate constraints, EnzDes constraints and AA type constraints.
    for constraint in constraints:
        constraint.apply(pose)
    return pose

def parse_cloud_pdb(cloud_pdb):
    cloud_pdb_lines = list()
    pdb_lines = list()
    with open(cloud_pdb, "r") as pf:
        substrate_name3 = None
        read_rotamer = None
        for line in pf:
            if line.startswith("REMARK 666 MATCH TEMPLATE"):
                if not substrate_name3:
                    chain_id = line[26]
                    substrate_name3 = line[28:31]
                    pdb_index = line[32:36]
                pdb_lines.append(line)
                read_rotamer = False
            elif line.startswith("ATOM  ") or line.startswith("HETATM"):
                if read_rotamer:
                    if line[17:20] == substrate_name3:
                        rotamer_lines.append(line[:21] + chain_id + pdb_index + line[26:])
                else:
                    if line[17:20] == substrate_name3 and line[21] == chain_id and \
                            line[22:26] == pdb_index:
                        pass
                    else:
                        pdb_lines.append(line)
            elif line.startswith("MODEL ") and read_rotamer == False:
                read_rotamer = True
                rotamer_lines = list()
            elif line.startswith("ENDMDL"):
                read_rotamer = False
                cloud_pdb_lines.append(rotamer_lines)
    cloud_pdb_lines.insert(0, pdb_lines)
    return cloud_pdb_lines

def read_scores_from_pdb(pdb_path, theozyme_positions:set=set(), \
        score_terms:set={"total_score", "coordinate_constraint", "atom_pair_constraint", \
        "angle_constraint", "dihedral_constraint"}):
    scores = dict()
    dG_substrate = None
    with open(pdb_path, "r") as pf:
        for line in pf:
            if line.startswith("label "):
                all_score_terms = line[:-1].split(" ")
            elif line.startswith("pose "):
                all_scores = line[:-1].split(" ")
                for score_term in score_terms:
                    key_word = score_term
                    if score_term == "total_score":
                        key_word = "total"
                    scores[score_term] = float(all_scores[all_score_terms.index(key_word)])
            elif line.split(" ")[0].split("_")[-1] in theozyme_positions:
                if dG_substrate == None:
                    dG_substrate = 0
                dG_substrate += float(line[:-1].split(" ")[-1])
    if dG_substrate is not None:
        scores["substrates"] = dG_substrate
    return scores

def energy_metric(score_function, pose, selection=None, score_type:str=None):
    # Create the metric
    metric = TotalEnergyMetric()
    metric.set_scorefunction(score_function)
    if score_type:
        exec("metric.set_scoretype(ScoreType.{})".format(score_type))
    # Add the selector
    if selection:
        metric.set_residue_selector(selection)
    return metric.calculate(pose)

def calculate_pose_scores(pose, score_function, theozyme_positions:set=set(), \
        score_terms:set={"coordinate_constraint", "atom_pair_constraint", \
        "angle_constraint", "dihedral_constraint"}):
    scores = dict()
    for score_term in score_terms:
        scores[score_term] = energy_metric(score_function, pose, score_type=score_term)
    if len(theozyme_positions) > 0:
        substrate_selection = ResidueIndexSelector(",".join(theozyme_positions))
        scores["substrates"] = energy_metric(score_function, pose, selection=substrate_selection)
    return scores

def run_job(score_function, pose, fold_tree, chi_dihedrals:list, constraint_file:str, \
        geometry_constraints, constraints, symmetry:str, movers, cloud_pdb_lines:list=None, \
        n_decoys:int=50, save_n_decoys:int=1, theozyme_positions:set=set(), \
        output_filename_prefix:str=None, wildtype_sequence:str=str()):
    if not output_filename_prefix:
        output_filename_prefix = pose.pdb_info().name().split("/")[-1][:-4]
    n_finished_decoys = 0
    decoy_scores_list = [None] * save_n_decoys
    decoy_filenames_list = [None] * save_n_decoys
    checkpoint_filename_scores = dict()
    for filename in filter(lambda f: f.startswith(output_filename_prefix) and \
            f.endswith(".pdb"), os.listdir()):
        scores = read_scores_from_pdb(filename, theozyme_positions=theozyme_positions)
        filename_split = filename.split(".")
        if filename.endswith("_checkpoint.pdb"):
            checkpoint_filename_scores[filename] = read_scores_from_pdb(filename)
        elif len(filename_split) >= 3 and filename_split[-2].isdigit():
            ith_ranked_decoy = int(filename_split[-2])
            now_saved_n_decoys = len(decoy_scores_list)
            if ith_ranked_decoy <= now_saved_n_decoys:
                decoy_scores_list[ith_ranked_decoy - 1] = scores
                decoy_filenames_list[ith_ranked_decoy - 1] = filename
            else:
                decoy_scores_list = decoy_scores_list + [None] * (ith_ranked_decoy - now_saved_n_decoys \
                        - 1) + [scores]
                decoy_filenames_list = decoy_filenames_list + [None] * (ith_ranked_decoy - \
                        now_saved_n_decoys - 1) + [filename]
        else:
            if filename.endswith("_tmp.pdb"):
                os.remove(filename)
            else:
                checkpoint_filename_scores[filename] = read_scores_from_pdb(filename)
    first_checkpoint = True
    for filename, scores in sorted(checkpoint_filename_scores.items(), key=lambda x: x[1]["total_score"]):
        if not first_checkpoint:
            os.remove(filename)
        if filename.endswith("_checkpoint.pdb"):
            n_finished_decoys = min(int(filename[:-4].split(".")[-1].split("_")[0]), n_decoys)
            output_filename = filename[:filename[:-4].rfind(".")] + ".pdb"
            if filename != output_filename:
                os.rename(filename, output_filename)
                filename = output_filename
        else:
            n_finished_decoys = n_decoys
        decoy_scores_list[0] = scores
        decoy_filenames_list[0] = filename
        first_checkpoint = False
    checkpoint_saved = True
    if decoy_filenames_list[0] == None:
        checkpoint_saved = False
    decoy_scores_list = list(filter(lambda sc: sc, decoy_scores_list))
    decoy_filenames_list = list(filter(lambda fn: fn, decoy_filenames_list))
    if len(decoy_filenames_list) == 0:
        checkpoint_saved = True
    else:
        with open(output_filename_prefix + ".sc", "w") as pf:
            for filename, scores in zip(decoy_filenames_list, decoy_scores_list):
                scores["decoy"] = filename
                pf.write(json.dumps(scores) + "\n")
    if not checkpoint_saved:
        n_finished_decoys = len(decoy_filenames_list)
        ranked_decoy_filename = decoy_filenames_list[0]
        checkpoint_filename = ranked_decoy_filename[:ranked_decoy_filename[:-4].rfind(".")] + \
                "." + str(n_finished_decoys) + "_checkpoint.pdb"
        os.rename(ranked_decoy_filename, checkpoint_filename)
        decoy_filenames_list[0] = checkpoint_filename
    for load_index, ranked_decoy_filename in enumerate(decoy_filenames_list[1:]):
        new_ranked_decoy_filename = ranked_decoy_filename[:ranked_decoy_filename[:-4].rfind(".")] + \
                "." + str(load_index + 2) + ".pdb"
        if load_index + 2 > save_n_decoys:
            os.remove(ranked_decoy_filename)
        elif new_ranked_decoy_filename != ranked_decoy_filename:
            os.rename(ranked_decoy_filename, new_ranked_decoy_filename)
            decoy_filenames_list[load_index + 1] = new_ranked_decoy_filename
    if len(decoy_filenames_list) > save_n_decoys:
        decoy_scores_list = decoy_scores_list[:save_n_decoys]
        decoy_filenames_list = decoy_filenames_list[:save_n_decoys]
    for i_decoy in range(n_finished_decoys + 1, n_decoys + 1):
        if cloud_pdb_lines:
            if len(cloud_pdb_lines) > 2:
                rotamer_index = np.random.randint(1, len(cloud_pdb_lines) - 1)
            else:
                rotamer_index = 1
            rotamer_lines = cloud_pdb_lines[rotamer_index]
            del cloud_pdb_lines[rotamer_index]
            tmp_pdb = output_filename_prefix + "." + str(i_decoy) + "_tmp.pdb"
            with open(tmp_pdb, "w") as p_pdb:
                p_pdb.writelines(cloud_pdb_lines[0])
                p_pdb.writelines(rotamer_lines)
            pose_copy = load_pdb_as_pose(score_function, tmp_pdb, fold_tree, chi_dihedrals, \
                    symmetry, constraint_file, geometry_constraints, constraints)
            os.remove(tmp_pdb)
        else:
            pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        mutated_sequence = pose_copy.sequence()
        pdb_info = pose_copy.pdb_info()
        insert_index = 0
        current_score = energy_metric(score_function, pose_copy)
        for index in range(len(decoy_scores_list)):
            if current_score >= decoy_scores_list[index]["total_score"]:
                insert_index = index + 1
            else:
                break
        current_output_filename = output_filename_prefix
        for index, wt_aa in enumerate(wildtype_sequence):
            aa = mutated_sequence[index]
            if aa != wt_aa:
                current_output_filename += "_" + wt_aa + pdb_info.pose2pdb(index + 1).split(" ")[0] + aa
        if i_decoy < n_decoys:
            checkpoint_file_suffix = "." + str(i_decoy) + "_checkpoint.pdb"
        else:
            checkpoint_file_suffix = ".pdb"
        if insert_index == 0:
            if i_decoy > 1:
                for move_index in range(len(decoy_filenames_list), 1, -1):
                    ranked_decoy_filename = decoy_filenames_list[move_index - 1]
                    new_ranked_decoy_filename = ranked_decoy_filename\
                            [:ranked_decoy_filename[:-4].rfind(".")] + "." + str(move_index + 1) + ".pdb"
                    os.rename(ranked_decoy_filename, new_ranked_decoy_filename)
                    decoy_filenames_list[move_index - 1] = new_ranked_decoy_filename
                old_checkpoint_filename = decoy_filenames_list[0]
                ranked_decoy_filename = old_checkpoint_filename[:old_checkpoint_filename[:-4].rfind(".")]\
                        + ".2.pdb"
                os.rename(old_checkpoint_filename, ranked_decoy_filename)
                decoy_filenames_list[0] = ranked_decoy_filename
            current_output_filename += checkpoint_file_suffix
        else:
            old_checkpoint_filename = decoy_filenames_list[0]
            checkpoint_filename = old_checkpoint_filename[:old_checkpoint_filename[:-4].rfind(".")] + \
                    checkpoint_file_suffix
            os.rename(old_checkpoint_filename, checkpoint_filename)
            decoy_filenames_list[0] = checkpoint_filename
            if insert_index < save_n_decoys:
                decoy_filenames_list[0] = checkpoint_filename
                for move_index in range(len(decoy_filenames_list), insert_index, -1):
                    ranked_decoy_filename = decoy_filenames_list[move_index - 1]
                    new_ranked_decoy_filename = ranked_decoy_filename\
                            [:ranked_decoy_filename[:-4].rfind(".")] + "." + str(move_index + 1) + ".pdb"
                    os.rename(ranked_decoy_filename, new_ranked_decoy_filename)
                    decoy_filenames_list[move_index - 1] = new_ranked_decoy_filename
                current_output_filename += "." + str(insert_index + 1) + ".pdb"
            else:
                continue
        scores = calculate_pose_scores(pose_copy, score_function, theozyme_positions=theozyme_positions)
        scores["total_score"] = current_score
        decoy_scores_list = decoy_scores_list[:insert_index] + [scores] + decoy_scores_list[insert_index:]
        pose_copy.dump_pdb(current_output_filename)
        decoy_filenames_list = decoy_filenames_list[:insert_index] + [current_output_filename] + \
                decoy_filenames_list[insert_index:]
        if len(decoy_filenames_list) > save_n_decoys:
            decoy_scores_list = decoy_scores_list[:save_n_decoys]
            os.remove(decoy_filenames_list[save_n_decoys])
            decoy_filenames_list = decoy_filenames_list[:save_n_decoys]
        with open(output_filename_prefix + ".sc", "w") as pf:
            for filename, scores in zip(decoy_filenames_list, decoy_scores_list):
                scores["decoy"] = filename
                pf.write(json.dumps(scores) + "\n")

def run_job_distributor(score_function, pose, fold_tree, chi_dihedrals:list, \
        constraint_file:str, geometry_constraints, constraints, symmetry:str, movers, \
        cloud_pdb_lines:list=None, n_decoys:int=5, output_filename_prefix:str=None, \
        wildtype_sequence:str=str()):
    if not output_filename_prefix:
        output_filename_prefix = pose.pdb_info().name().split("/")[-1][:-4]
    job_distributor = PyJobDistributor(output_filename_prefix, n_decoys, score_function)
    while not job_distributor.job_complete:
        i_decoy = job_distributor.current_id
        if cloud_pdb_lines:
            rotamer_index = np.random.randint(1, len(cloud_pdb_lines) - 1)
            rotamer_lines = cloud_pdb_lines[rotamer_index]
            del cloud_pdb_lines[rotamer_index]
            tmp_pdb = output_filename_prefix + "." + str(i_decoy) + "_tmp.pdb"
            with open(tmp_pdb, "w") as p_pdb:
                p_pdb.writelines(cloud_pdb_lines[0])
                p_pdb.writelines(rotamer_lines)
            pose_copy = load_pdb_as_pose(score_function, tmp_pdb, fold_tree, chi_dihedrals, \
                    symmetry, constraint_file, geometry_constraints, constraints)
            os.remove(tmp_pdb)
        else:
            pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        job_distributor.output_decoy(pose_copy)
        mutated_sequence = pose_copy.sequence()
        pdb_info = pose_copy.pdb_info()
        suffix = str()
        for index, wt_aa in enumerate(wildtype_sequence):
            aa = mutated_sequence[index]
            if aa != wt_aa:
                suffix += "_" + wt_aa + pdb_info.pose2pdb(index + 1).split(" ")[0] + aa
        os.rename(job_distributor.pdb_name + "_" + str(i_decoy) + ".pdb", \
                job_distributor.pdb_name + suffix + "_" + str(i_decoy) + ".pdb")

def main(args):
    # Create the score function.
    score_function = set_score_function(args.score_function, symmetry=args.symmetry, \
            score_terms=args.score_terms)
    # Load pdb files as poses.
    pose = pose_from_pdb(args.pdb)
    coord_ref_pose = None
    if args.coordinate_reference_pdb:
        coord_ref_pose = pose_from_pdb(args.coordinate_reference_pdb)
    ddG_ref_pose = None
    if args.ddG_reference_pdb:
        if args.ddG_reference_pdb == "True":
            ddG_ref_pose = pose
        else:
            ddG_ref_pose = pose_from_pdb(args.ddG_reference_pdb)
    # Parse cloud pdb.
    if args.cloud_pdb:
        cloud_pdb_lines = parse_cloud_pdb(args.pdb)
        args.n_decoys = min(args.n_decoys, len(cloud_pdb_lines[1:]))
        with open(args.pdb[:-4] + ".tmp.pdb", "w") as pf:
            pf.writelines(cloud_pdb_lines[0])
            pf.writelines(cloud_pdb_lines[1])
        pose = pose_from_pdb(args.pdb[:-4] + ".tmp.pdb")
        os.remove(args.pdb[:-4] + ".tmp.pdb")
    else:
        cloud_pdb_lines = None
    # Get sequences information.
    wildtype_sequence = pose.sequence()
    sequence_length = len(wildtype_sequence)
    if not args.output_filename_mutations_suffix:
        wildtype_sequence = str()
    if args.ddG_reference_pdb:
        ddG_reference_sequence_length = len(ddG_ref_pose.sequence())
        truncate_sequence_end = ddG_reference_sequence_length - sequence_length
        assert truncate_sequence_end >= 0
    else:
        truncate_sequence_end = 0
    # Convert pdb numberings to pose numberings.
    # Any ddG_ref_pose redundant positions in the site-directed residue actions, 
    # i.e., -static, -mut and -des, will be ignored.
    static_pose_indices, _ = pdb_to_pose_numbering(pose, args.static_residues)
    if args.static_residue_identities:
        static_pose_indices.update(residue_name3_selector(pose, args.static_residue_identities, \
                sequence_length=sequence_length))
    mutation_pose_indices, _ = pdb_to_pose_numbering(pose, args.mutations)
    design_pose_indices, _ = pdb_to_pose_numbering(pose, args.design_residues)
    # Get pose indices of substrates and catalytic residues.
    theozyme_pose_indices = set()
    rigid_body_tform_pose_indices = set()
    if args.catalytic_residues:
        catalytic_residue_pose_indices, _ = pdb_to_pose_numbering(pose, args.catalytic_residues)
        theozyme_pose_indices.update(catalytic_residue_pose_indices)
    if args.catalytic_residue_identities:
        catalytic_residue_by_id_pose_indices = residue_name3_selector(pose, args.catalytic_residue_identities)
        theozyme_pose_indices.update(catalytic_residue_by_id_pose_indices)
    if args.substrates:
        substrate_pose_indices, jump_pose_indices = pdb_to_pose_numbering(pose, args.substrates)
        theozyme_pose_indices.update(substrate_pose_indices)
        if args.substrate_rigid_body_transformations:
            rigid_body_tform_pose_indices.update(jump_pose_indices)
    if args.substrate_identities:
        substrate_by_id_pose_indices = residue_name3_selector(pose, args.substrate_identities)
        theozyme_pose_indices.update(substrate_by_id_pose_indices)
        if args.substrate_rigid_body_transformations:
            rigid_body_tform_pose_indices.update(substrate_by_id_pose_indices)
        if args.ddG_reference_pdb:
            theozyme_pose_indices.update(residue_name3_selector(ddG_ref_pose, args.substrate_identities))
    # Get pose indices of enzdes positions.
    enzdes_pose_indices = set()
    # if args.enzyme_design_constraints:
    if args.ddG_reference_pdb and args.ddG_reference_pdb != "True":
        enzdes_substrate_pose_indices, enzdes_res_pose_indices = get_enzdes_pose_indices(\
                ddG_ref_pose.pdb_info(), args.ddG_reference_pdb, args.symmetry)
    else:
        enzdes_substrate_pose_indices, enzdes_res_pose_indices = get_enzdes_pose_indices(\
                pose.pdb_info(), args.pdb, args.symmetry)
    enzdes_pose_indices = enzdes_substrate_pose_indices.union(enzdes_res_pose_indices)
    if args.enzdes_substrates_transformations:
        rigid_body_tform_pose_indices.update(enzdes_substrate_pose_indices)
    # Get all theozyme positions.
    nonredundant_theozyme_pose_indices = set(filter(lambda index: int(index) \
            <= sequence_length, theozyme_pose_indices.union(enzdes_pose_indices)))
    # Set fold tree.
    if not args.fold_tree and not args.alter_jump_edges:
        fold_tree = pose.fold_tree()
    else:
        if args.fold_tree:
            fold_tree = create_fold_tree(args.fold_tree)
        elif args.alter_jump_edges:
            fold_tree = alter_fold_tree_jump_edges(pose.fold_tree(), args.alter_jump_edges)
        pose.fold_tree(fold_tree)
    # Rigid body transformations.
    rigid_body_tform_pose_indices = rigid_body_tform_pose_indices - static_pose_indices
    rigid_body_tform_jump_edges = set()
    if len(rigid_body_tform_pose_indices) > 0:
        rigid_body_tform_jump_edges = pose_indices_to_jump_edges(fold_tree, \
                rigid_body_tform_pose_indices)
    # Create geometry constraints on-the-fly.
    geometry_constraints = list()
    if args.distance_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, args.distance_constraint_atoms, \
                args.distance_constraint_parameters, \
                bounded_distances=args.bounded_distances))
    if args.angle_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, args.angle_constraint_atoms, \
                args.angle_constraint_parameters))
    if args.dihedral_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, args.dihedral_constraint_atoms, \
                args.dihedral_constraint_parameters))
    # Set chi dihedrals.
    set_chi_dihedral(pose, args.chi_dihedrals)
    # Apply symmetry if specified.
    n_monomers = set_symmetry(args.symmetry, pose, coord_ref_pose, ddG_ref_pose, \
            sequence_length=sequence_length)
    # Create the task factory.
    task_factory, mutators, min_shell_focus_selection = create_task_factory(\
            nopack_positions=static_pose_indices, \
            point_mutations=mutation_pose_indices, \
            design_positions=design_pose_indices, \
            theozyme_positions=theozyme_pose_indices, \
            enzdes_positions=enzdes_pose_indices, \
            design_binding_site=args.design_binding_site, \
            design_enzdes_shell=args.design_enzdes_shell, \
            repack_neighborhood_only=args.repack_neighborhood_only, \
            repack_binding_site=not args.no_repack_binding_site, \
            repack_enzdes_shell=not args.no_repack_enzdes_shell, \
            ddG_ref_pose=ddG_ref_pose, ddG_wildtype=args.ddG_wildtype, \
            n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end, \
            excluded_amino_acid_types=args.excluded_amino_acid_types, \
            noncanonical_amino_acids=args.noncanonical_amino_acids, \
            allow_ncaa_in_design=args.allow_ncaa_in_design)
    # List of coordinate constraints, EnzDes constraints and AA type constraints.
    constraints = list()
    # Add coordinate constraints.
    no_coord_cst_selection = OrResidueSelector()
    no_coord_cst_residues, _ = pdb_to_pose_numbering(pose, \
            args.no_backbone_coordinate_constraints_residues)
    no_coord_cst_selection.add_residue_selector(ResidueIndexSelector(\
            ",".join(no_coord_cst_residues) + ","))
    if args.no_coordinate_constraints_on_packing_region:
        if not min_shell_focus_selection:
            no_coord_cst_selection = TrueResidueSelector()
        elif type(min_shell_focus_selection) == set:
            pack_positions = set(filter(lambda index: int(index) <= sequence_length, \
                    min_shell_focus_selection)) - static_pose_indices
            no_coord_cst_selection.add_residue_selector(ResidueIndexSelector(\
                    ",".join(pack_positions) + ","))
        else:
            no_coord_cst_selection.add_residue_selector(min_shell_focus_selection)
    add_csts = AddConstraints()
    bb_coord_cst_gen = create_coordinate_constraints(reference_pose=coord_ref_pose, \
            selection=NotResidueSelector(no_coord_cst_selection), \
            standard_deviation=args.coordinate_constraints_standard_deviation, \
            bounded=args.bounded_coordinate_constraints, \
            ca_only=args.only_CA_coordinate_constraints)
    add_csts.add_generator(bb_coord_cst_gen)
    all_atom_coord_cst_selection = ResidueIndexSelector(",".join(\
            args.all_atom_coordinate_constraints_positions) + ",")
    all_atom_coord_cst_gen = create_coordinate_constraints(reference_pose=coord_ref_pose, \
            selection=all_atom_coord_cst_selection, \
            standard_deviation=args.coordinate_constraints_standard_deviation, \
            bounded=args.bounded_coordinate_constraints, \
            ca_only=args.only_CA_coordinate_constraints, side_chain=True)
    add_csts.add_generator(all_atom_coord_cst_gen)
    add_csts.apply(pose)
    constraints.append(add_csts)
    # Read constraint files from the command line and apply to pose.
    if args.constraint_file:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Apply geometry constraints.
    for geometry_constraint in geometry_constraints:
        pose.add_constraint(geometry_constraint)
    # Apply enzyme design constraints.
    if args.enzyme_design_constraints:
        enzdes_cst = create_enzdes_constraints()
        enzdes_cst.apply(pose)
        constraints.append(enzdes_cst)
    # Favor native AA types.
    if args.favor_native_residue and (args.design_residues or args.design_binding_site):
        favor_nataa = FavorNativeResidue(pose, args.favor_native_residue)
        favor_nataa.apply(pose)
        constraints.append(favor_nataa)
    # List of movers.
    movers = list()
    # Create the RMSD metric.
    no_rmsd_residues, _ = pdb_to_pose_numbering(pose, args.no_rmsd_residues)
    rmsd_metric = RMSDMetric(pose, NotResidueSelector(ResidueIndexSelector(\
            ",".join(no_rmsd_residues) + ",")))
    # Make some point mutations if not args.ddG_wildtype.
    movers.extend(mutators)
    # Perform pre-minimization.
    assembly_length = len(pose.sequence())
    pre_minimization_pose_indices = nonredundant_theozyme_pose_indices - static_pose_indices
    if args.pre_minimization and len(pre_minimization_pose_indices) > 0:
        pre_minimizer = create_pre_minimizer(score_function, assembly_length, \
                pre_minimization_pose_indices, jump_edges=rigid_body_tform_jump_edges)
        movers.append(pre_minimizer)
    # Create the move map.
    if not args.minimize_neighborhood_only:
        min_shell_focus_selection = assembly_length
    move_map = create_move_map(min_shell_focus_selection, \
            minimize_neighborhood_only=args.minimize_neighborhood_only, \
            static_positions=static_pose_indices, ddG_ref_pose=ddG_ref_pose, \
            n_monomers=n_monomers, truncate_sequence_end=truncate_sequence_end, \
            jump_edges=rigid_body_tform_jump_edges)
    # Create the fast relax mover.
    fast_relax = create_fast_relax_mover(score_function, task_factory, move_map=move_map)
    movers.append(fast_relax)
    # Add the RMSD metric.
    movers.append(rmsd_metric)
    # Run jobs.
    if args.debug_mode:
        print(pose.fold_tree())
        if args.enzyme_design_constraints:
            print(enzdes_substrate_pose_indices)
            print(enzdes_res_pose_indices)
        packer_task = task_factory.create_task_and_apply_taskoperations(pose)
        packer_task.show()
        if type(move_map) == MoveMapFactory:
            move_map = move_map.create_movemap_from_pose(pose)
        move_map.show()
        print(score_function.show(pose))
    elif args.save_n_decoys:
        run_job(score_function, pose, fold_tree, args.chi_dihedrals, \
                args.constraint_file, geometry_constraints, constraints, \
                args.symmetry, movers, cloud_pdb_lines=cloud_pdb_lines, \
                n_decoys=args.n_decoys, save_n_decoys=args.save_n_decoys, \
                theozyme_positions=nonredundant_theozyme_pose_indices, \
                output_filename_prefix=args.output_filename_prefix, \
                wildtype_sequence=wildtype_sequence)
    else:
        run_job_distributor(score_function, pose, fold_tree, args.chi_dihedrals, \
                args.constraint_file, geometry_constraints, constraints, \
                args.symmetry, movers, cloud_pdb_lines=cloud_pdb_lines, \
                n_decoys=args.n_decoys, output_filename_prefix=args.output_filename_prefix, \
                wildtype_sequence=wildtype_sequence)


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    main(args)
