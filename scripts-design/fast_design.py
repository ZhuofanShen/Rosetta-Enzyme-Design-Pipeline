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
    add_fa_constraints_from_cmdline, AmbiguousConstraint, \
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
from pyrosetta.rosetta.protocols.membrane import \
    AddMembraneMover
from pyrosetta.rosetta.protocols.membrane.symmetry import \
    SymmetricAddMembraneMover
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.protein_interface_design import \
    FavorNativeResidue
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
    parser.add_argument("-symm", "--symmetry", type=str)
    parser.add_argument("-span", "--membrane_span_file", type=str)
    parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst", \
            help="$base_$soft_$cart_$cst")
    parser.add_argument("--score_terms", type=str, nargs="*", default=list())
    parser.add_argument("-ft", "--fold_tree", type=str, nargs="*")
    parser.add_argument("-ddg_wt", "--ddG_wildtype", action="store_true", \
            help="Only repack instead of -premuts, -muts, -des, -des_bs and -des_enzdes.")
    parser.add_argument("-premuts", "--pre_mutations", type=str, nargs="*", default=list(), 
            help="Site-directed AA substitution. Applied before everything.")
    parser.add_argument("-edges", "--alter_jump_edges", type=str, nargs="*", help=\
            "[$edge,$atom1,$atom2] or [$edge,$upstream_edge] or [$edge,$upstream_edge,$atom1,$atom2] * n")
    parser.add_argument("-chis", "--chi_dihedrals", type=str, nargs="*", default=list(), \
            help="$chain$residue_index,$chi,$degree or $residue_name3,$chi,$degree")
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
    elif constraints:
        score_function.set_weight(ScoreType.pro_close, 1.25)
        score_function.set_weight(ScoreType.metalbinding_constraint, 1)
        score_function.set_weight(ScoreType.chainbreak, 1)
    if constraints:
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
    name3_chain_id_dict = dict()
    for name3 in name3_list:
        if "-" in name3:
            chain_id, name3 = name3.split("-")
            name3_chain_id_dict[name3] = set(chain_id)
        else:
            name3_chain_id_dict[name3] = set()
    pdb_info = pose.pdb_info()
    for pose_index in list(range(1, len(pose.sequence()) + 1))[:sequence_length]:
        name3 = pose.residue(pose_index).name3().strip(" ")
        chain_id_set = name3_chain_id_dict.get(name3)
        if chain_id_set is not None:
            chain_id = pdb_info.pose2pdb(pose_index).split(" ")[1]
            if len(chain_id_set) == 0 or chain_id in chain_id_set:
                pose_indices.add(str(pose_index))
    return pose_indices

def create_residue_mutators(mutations, noncanonical_amino_acids:list=list(), \
        ddG_wildtype:bool=False):
    # Compatible amino acid alphabet.
    AA_1to3_dict = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", \
            "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", \
            "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", \
            "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    # Include noncanonical amino acids.
    for ncaa in noncanonical_amino_acids:
        ncaa_name1, ncaa_name3 = ncaa.split(",")
        AA_1to3_dict[ncaa_name1] = ncaa_name3
    # Create residue mutators.
    mutators = list()
    mutation_pose_indices = set()
    for mutation in mutations:
        mutating_position, target_AA = mutation.split(",")
        mutation_pose_indices.add(mutating_position)
        if not ddG_wildtype:
            mutator = MutateResidue()
            mutator.set_selector(ResidueIndexSelector(mutating_position))
            target_AA_name3 = AA_1to3_dict.get(target_AA)
            if target_AA_name3 is None:
                raise Exception("Specified AA type " + target_AA + " is not declared.")
            mutator.set_res_name(target_AA_name3)
            mutators.append(mutator)
    return mutators, mutation_pose_indices

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
    if type(distance) is np.ndarray:
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

def get_enzdes_pose_indices(pdb_info, pdb, symmetry):
    '''
    Get pose indices of enzdes positions by reading the pdb file REMARK 666 header.
    '''
    enzdes_substrate_pose_indices = set()
    enzdes_res_pose_indices = set()
    read_remark666_lines = False
    main_chain = None
    with open(pdb, "r") as fpdb:
        for line in fpdb:
            if line.startswith("REMARK 666 MATCH TEMPLATE"):
                line_split = list(filter(lambda x: x, line.split(" ")))
                # substrate
                substrate_chain_id = line_split[4]
                if not main_chain:
                    main_chain = substrate_chain_id
                if main_chain == substrate_chain_id or not symmetry:
                    substrate_res_id = int(line_split[6])
                    substrate_pose_id = pdb_info.pdb2pose(substrate_chain_id, substrate_res_id)
                    enzdes_substrate_pose_indices.add(str(substrate_pose_id))
                # enzdes residues
                motif_chain_id = line_split[9]
                if main_chain == motif_chain_id or not symmetry:
                    motif_res_id = int(line_split[11])
                    motif_pose_id = pdb_info.pdb2pose(motif_chain_id, motif_res_id)
                    enzdes_res_pose_indices.add(str(motif_pose_id))
                read_remark666_lines = True
            elif read_remark666_lines:
                break
    return enzdes_substrate_pose_indices, enzdes_res_pose_indices

def pose_indices_to_jump_edges(fold_tree, rigid_body_tform_pose_indices):
    jump_edges = set()
    for rigid_body_tform_pose_index in rigid_body_tform_pose_indices:
        jump_edge = fold_tree.get_residue_edge(int(rigid_body_tform_pose_index))
        if jump_edge.is_jump():
            jump_edges.add(jump_edge.label())
    return jump_edges

def apply_symmetry_membrane(symmetry:str, membrane_span_file:str, *pose_list, sequence_length:int=0):
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
    if membrane_span_file:
        if symmetry:
            add_membrane = SymmetricAddMembraneMover(membrane_span_file)
        else:
            add_membrane = AddMembraneMover(membrane_span_file)
        for pose in filter(lambda p: p is not None, pose_list):
            add_membrane.apply(pose)
    return n_monomers

def index_boolean_filter(index, boolean):
    if type(boolean) == np.ndarray:
        for bool in boolean:
            if bool:
                return str(index)
    elif boolean:
        return str(index)
    return False

def boolean_vector_to_indices_set(boolean_vector, n_monomers:int=1):
    boolean_vector = np.array(boolean_vector)
    sequence_length = len(boolean_vector) // n_monomers
    boolean_matrix = boolean_vector[:n_monomers*sequence_length]\
            .reshape(n_monomers, sequence_length).transpose()
    indices = set(filter(lambda x: x, map(index_boolean_filter, \
            range(1, len(boolean_matrix) + 1), boolean_matrix)))
    return indices

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

def create_pre_minimization_move_map(assembly_length, pre_minimization_pose_indices, \
        jump_edges:set=set()):
    pre_minimization_vector = vector1_bool(assembly_length)
    for pose_index in pre_minimization_pose_indices:
        pre_minimization_vector[int(pose_index)] = True
    move_map = MoveMap()
    move_map.set_bb(False)
    move_map.set_chi(pre_minimization_vector)
    for jump_edge in jump_edges:
        move_map.set_jump(jump_edge, True)
    return move_map

def create_task_factory(specified_static_positions:set=set(), pre_mutation_positions:set=set(), \
        point_mutations:set=set(), design_positions:set=set(), theozyme_positions:set=set(), \
        enzdes_positions:set=set(), design_binding_site:bool=False, design_enzdes_shell:bool=False, \
        repack_neighborhood_only:bool=False, repack_binding_site:bool=True, \
        repack_enzdes_shell:bool=True, ddG_ref_pose=None, n_monomers:int=1, \
        sequence_length:int=0, ddG_wildtype:bool=False, excluded_amino_acid_types:str=None, \
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
    all_AAs = set("ACDEFGHIKLMNPQRSTVWY")
    # Include noncanonical amino acids.
    if allow_ncaa_in_design:
        ncaa_palette = CustomBaseTypePackerPalette()
        for ncaa in noncanonical_amino_acids:
            ncaa_name1, ncaa_name3 = ncaa.split(",")
            ncaa_palette.add_type(ncaa_name3)
            all_AAs.add(ncaa_name1)
        task_factory.set_packer_palette(ncaa_palette)

    # Repack its surrounding shell w/o including the focus region itself.
    # Specify the point mutations. Not include pre mutation positions.
    packer_mutation_positions = set()
    for point_mutation in point_mutations:
        mutating_position, target_AA = point_mutation.split(",")
        packer_mutation_positions.add(mutating_position)
        if not ddG_wildtype:
            restriction = RestrictAbsentCanonicalAASRLT()
            restriction.aas_to_keep(target_AA)
            task_factory.push_back(OperateOnResidueSubset(restriction, \
                    ResidueIndexSelector(mutating_position)))
    mutation_positions = packer_mutation_positions.union(pre_mutation_positions)

    # Specify the site-directed design positions.
    design_positions = design_positions - packer_mutation_positions
    if len(design_positions.intersection(specified_static_positions)) > 0:
        raise Exception("Site-directed design positions cannot contain specified static positions.")
    design_selection = ResidueIndexSelector(",".join(design_positions) + ",")
    packer_substitution_selection = ResidueIndexSelector(",".join(\
            packer_mutation_positions.union(design_positions)) + ",")

    # Design the surrounding shell of the focus region.
    design_shell_focus_positions = set()
    if design_binding_site:
        design_shell_focus_positions.update(theozyme_positions)
    if design_enzdes_shell:
        design_shell_focus_positions.update(enzdes_positions)
    # Repack the surrounding shell of the focus region.
    repack_shell_focus_positions = pre_mutation_positions.copy()
    if repack_binding_site:
        repack_shell_focus_positions.update(theozyme_positions)
    if repack_enzdes_shell:
        repack_shell_focus_positions.update(enzdes_positions)
    # Always repack the focus region itself.
    additional_repacking_positions = theozyme_positions.union(enzdes_positions)
    # The following 3 selections may contain pre-mutation positions 
    # or specified static positions or redundant ddG_ref_pose positions.
    design_shell_focus_selection = ResidueIndexSelector(",".join(\
            design_shell_focus_positions) + ",")
    repack_shell_focus_selection = ResidueIndexSelector(",".join(\
            repack_shell_focus_positions) + ",")
    additional_repacking_selection = ResidueIndexSelector(",".join(\
            additional_repacking_positions) + ",")

    # Design the theozyme-protein interface.
    design_shell_selection = select_neighborhood_region(design_shell_focus_selection, False)
    # Fix the design positions (or not).
    # Exclude the site-directed substitution positions and static positions.
    if ddG_ref_pose is not None:
        assert sequence_length > 0
        binding_site_positions = set(filter(lambda x: int(x) <= sequence_length, \
                boolean_vector_to_indices_set(design_shell_selection.apply(ddG_ref_pose), \
                n_monomers=n_monomers))) - mutation_positions - specified_static_positions
        design_shell_selection = ResidueIndexSelector(",".join(binding_site_positions) + ",")
    else:
        design_shell_selection = AndResidueSelector(design_shell_selection, \
                NotResidueSelector(ResidueIndexSelector(",".join(\
                mutation_positions.union(specified_static_positions)) + ",")))
    packer_substitution_selection = OrResidueSelector(packer_substitution_selection, design_shell_selection)

    # Every position is designable by defalut in the task factory other than specification.
    if excluded_amino_acid_types and not ddG_wildtype:
        # Exclude some AA types if specified.
        excluded_AAs = set(excluded_amino_acid_types)
        restriction = RestrictAbsentCanonicalAASRLT(",".join(all_AAs - excluded_AAs))
        restriction.aas_to_keep()
        task_factory.push_back(OperateOnResidueSubset(restriction, \
                OrResidueSelector(design_selection, design_shell_selection)))

    # Repacking region w/o AA type change.
    if repack_neighborhood_only:
        # Repack the neighborhood region of the AA substitutions and theozyme positions.
        # May contain specified static positions or redundant ddG_ref_pose positions.
        min_shell_focus_selection = select_neighborhood_region(OrResidueSelector(\
                packer_substitution_selection, repack_shell_focus_selection), True)
        min_shell_focus_selection = OrResidueSelector(min_shell_focus_selection, \
                additional_repacking_selection)
    else:
        min_shell_focus_selection = TrueResidueSelector()
    # Fix the repacking positions (or not).
    # Exclude the specified static positions.
    if ddG_ref_pose is not None:
        min_shell_focus_positions = boolean_vector_to_indices_set(\
                min_shell_focus_selection.apply(ddG_ref_pose), n_monomers=n_monomers)
        min_shell_focus_selection = ResidueIndexSelector(",".join(min_shell_focus_positions) + ",")
        packer_sampling_positions = set(filter(lambda x: int(x) <= sequence_length, \
                min_shell_focus_positions)) - specified_static_positions
        packer_sampling_selection = ResidueIndexSelector(",".join(packer_sampling_positions) + ",")
    else:
        packer_sampling_selection = AndResidueSelector(min_shell_focus_selection, \
                NotResidueSelector(ResidueIndexSelector(",".join(specified_static_positions) + ",")))

    # Set the repacking and static region in the task factory.
    if ddG_wildtype:
        # Repack instead of mutate or design.
        repacking_selection = packer_sampling_selection
    else:
        # Exclude the packer substitution positions from repacking.
        repacking_selection = AndResidueSelector(packer_sampling_selection, \
                NotResidueSelector(packer_substitution_selection))
    task_factory.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), repacking_selection))
    task_factory.push_back(OperateOnResidueSubset(PreventRepackingRLT(), packer_sampling_selection, True))

    return task_factory, min_shell_focus_selection, packer_sampling_selection

def create_move_map(focus_selection, static_positions:set=set(), ddG_ref_pose=None, \
        n_monomers:int=1, sequence_length:int=0, assembly_length:int=0, jump_edges:set=set()):
    # Only the focus selection along with its neighborhood shell wll be subject to 
    # backbone and sidechain minimization except the static positions.
    minimization_selection = select_neighborhood_region(focus_selection, True)
    if ddG_ref_pose is not None: # Create a move map with fixed minimization positions.
        assert sequence_length > 0
        assert assembly_length > 0
        minimization_positions = set(filter(lambda x: int(x) <= sequence_length, \
                boolean_vector_to_indices_set(minimization_selection.apply(ddG_ref_pose), \
                n_monomers=n_monomers))) - static_positions
        minimization_vector = vector1_bool(assembly_length)
        for minimization_position in minimization_positions:
            minimization_vector[int(minimization_position)] = True
        move_map = MoveMap()
        move_map.set_bb(minimization_vector)
        move_map.set_chi(minimization_vector)
        for jump_edge in jump_edges:
            move_map.set_jump(jump_edge, True)
    else: # Create a movemap factory using residue and jump selectors.
        minimization_selection = AndResidueSelector(minimization_selection, \
                NotResidueSelector(ResidueIndexSelector(",".join(static_positions) + ",")))
        move_map = MoveMapFactory()
        move_map.add_bb_action(move_map_action.mm_enable, minimization_selection)
        move_map.add_chi_action(move_map_action.mm_enable, minimization_selection)
        for jump_edge in jump_edges:
            move_map.add_jump_action(move_map_action.mm_enable, JumpIndexSelector(jump_edge))
    return move_map

def create_packer(score_function, task_factory):
    packer = PackRotamersMover(score_function)
    packer.task_factory(task_factory)
    return packer

def create_minimizer(score_function, move_map=None):
    minimizer = MinMover()
    minimizer.min_type("lbfgs_armijo_nonmonotone")
    minimizer.score_function(score_function)
    if score_function.get_weight(ScoreType.cart_bonded) > 0:
        minimizer.cartesian(True)
    if type(move_map) == MoveMap:
        minimizer.set_movemap(move_map)
    elif type(move_map) == MoveMapFactory:
        minimizer.set_movemap_factory(move_map)
    return minimizer

def create_fast_relax_mover(score_function, task_factory, move_map=None):
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    if score_function.get_weight(ScoreType.cart_bonded) > 0:
        fast_relax.cartesian(True)
    if type(move_map) == MoveMap:
        fast_relax.set_movemap(move_map)
    elif type(move_map) == MoveMapFactory:
        fast_relax.set_movemap_factory(move_map)
    return fast_relax

def load_pdb_as_pose(score_function, pdb:str, pre_mutators, fold_tree, chi_dihedrals:list, \
        symmetry:str, membrane_span_file:str, constraint_file:str, geometry_constraints, \
        constraints, favor_native_residue):
    # Load pdb as pose.
    pose = pose_from_pdb(pdb)
    for pre_mutator in pre_mutators:
        pre_mutator.apply(pose)
    if fold_tree:
        pose.fold_tree(fold_tree)
    set_chi_dihedral(pose, chi_dihedrals)
    # Apply symmetry if specified.
    apply_symmetry_membrane(symmetry, membrane_span_file, pose)
    # Read constraint files from the command line and apply to pose.
    if constraint_file:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Apply geometry constraints.
    for geometry_constraint in geometry_constraints:
        pose.add_constraint(geometry_constraint)
    # Apply coordinate constraints and EnzDes constraints.
    for constraint in constraints:
        constraint.apply(pose)
    # Apply AA type constraints.
    if favor_native_residue:
        FavorNativeResidue(pose, favor_native_residue)
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
                    line_split = list(filter(lambda x: x, line.split(" ")))
                    chain_id = line_split[4]
                    substrate_name3 = line_split[5]
                    pdb_index = line_split[6]
                    space = " " * (4 - len(pdb_index))
                pdb_lines.append(line)
                read_rotamer = False
            elif line.startswith("ATOM  ") or line.startswith("HETATM"):
                res_name3 = line[17:20].strip(" ")
                if read_rotamer:
                    if res_name3 == substrate_name3:
                        rotamer_lines.append(line[:21] + chain_id + space + pdb_index + line[26:])
                else:
                    if res_name3 == substrate_name3 and line[21] == chain_id and \
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
    if selection is not None:
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

def run_job(score_function, pose, pre_mutators, fold_tree, chi_dihedrals:list, \
        symmetry:str, membrane_span_file:str, constraint_file:str, geometry_constraints, \
        constraints, favor_native_residue, movers, cloud_pdb_lines:list=None, \
        n_decoys:int=50, save_n_decoys:int=1, theozyme_positions:set=set(), \
        output_filename_prefix:str=None, wildtype_sequence:str=str()):
    if not output_filename_prefix:
        output_filename_prefix = os.path.basename(pose.pdb_info().name())\
                .rstrip(".pdb").rstrip(".tmp") + "_relaxed"
    if os.path.abspath(pose.pdb_info().name().rstrip(".pdb").rstrip(".tmp")) == \
            os.path.abspath("") + "/" + output_filename_prefix:
        raise Exception("Please assign an output filename prefix different than the input pdb filename \
                or store the input pdb file in a folder other than the current directory.")
    if os.path.isfile(output_filename_prefix + ".sc"):
        decoy_scores_list = list()
        with open(output_filename_prefix + ".sc", "r") as pf:
            for line in pf:
                scores = json.loads(line)
                scores.pop("decoy", None)
                decoy_scores_list.append(scores)
    else:
        decoy_scores_list = [None] * save_n_decoys
    decoy_filenames_list = [None] * save_n_decoys
    checkpoint_filename_scores = dict()
    for filename in filter(lambda f: f.startswith(output_filename_prefix) and \
            f.endswith(".pdb"), os.listdir()):
        scores = read_scores_from_pdb(filename, theozyme_positions=theozyme_positions)
        filename_split = filename.split(".")
        if filename.endswith("_checkpoint.pdb"):
            checkpoint_filename_scores[filename] = scores
        elif len(filename_split) >= 3 and filename_split[-2].isdigit():
            ith_ranked_decoy = int(filename_split[-2])
            if ith_ranked_decoy <= len(decoy_scores_list):
                existing_scores = decoy_scores_list[ith_ranked_decoy - 1]
                if existing_scores:
                    assert existing_scores["total_score"] - scores["total_score"] < 0.1
                else:
                    decoy_scores_list[ith_ranked_decoy - 1] = scores
            else:
                decoy_scores_list = decoy_scores_list + [None] * (ith_ranked_decoy - \
                        len(decoy_scores_list) - 1) + [scores]
            if ith_ranked_decoy <= len(decoy_filenames_list):
                decoy_filenames_list[ith_ranked_decoy - 1] = filename
            else:
                decoy_filenames_list = decoy_filenames_list + [None] * (ith_ranked_decoy - \
                        len(decoy_filenames_list) - 1) + [filename]
        else:
            if filename.endswith("_tmp.pdb"):
                os.remove(filename)
            else:
                checkpoint_filename_scores[filename] = scores
    n_finished_decoys = 0
    checkpoint_saved = False
    for filename, scores in sorted(checkpoint_filename_scores.items(), key=lambda x: x[1]["total_score"]):
        if checkpoint_saved:
            raise Exception("More than one checkpoint file found in the current directory!")
        if filename.endswith("_checkpoint.pdb"):
            n_finished_decoys = min(int(filename[:-4].split(".")[-1].split("_")[0]), n_decoys)
        decoy_scores_list[0] = scores
        decoy_filenames_list[0] = filename
        checkpoint_saved = True
    decoy_scores_list_1 = list(filter(lambda s: s,map(lambda s, f: s if f else None, \
            decoy_scores_list[:len(decoy_filenames_list)], decoy_filenames_list)))
    decoy_scores_list_2 = decoy_scores_list[len(decoy_filenames_list):]
    decoy_filenames_list = list(filter(lambda fn: fn, decoy_filenames_list))
    if len(decoy_filenames_list) >= save_n_decoys:
        decoy_scores_list = decoy_scores_list_1 + decoy_scores_list_2
        for decoy_filename in decoy_filenames_list[save_n_decoys:]:
            os.remove(decoy_filename)
        decoy_filenames_list = decoy_filenames_list[:save_n_decoys]
    else:
        decoy_scores_list = decoy_scores_list_1
    if not n_finished_decoys:
        n_finished_decoys = len(decoy_scores_list_1) + len(decoy_scores_list_2)
    now_unsaved_n_decoys = save_n_decoys - len(decoy_filenames_list)
    n_unfinished_decoys = n_decoys - n_finished_decoys
    if now_unsaved_n_decoys > n_unfinished_decoys:
        n_decoys += now_unsaved_n_decoys - n_unfinished_decoys
    if len(decoy_filenames_list) > 0:
        checkpoint_filename = decoy_filenames_list[0]
        new_checkpoint_filename = None
        i = checkpoint_filename[:-4].rfind(".")
        if i == -1:
            i = -4
        if n_finished_decoys < n_decoys and not checkpoint_filename.endswith("_checkpoint.pdb"):
            new_checkpoint_filename = checkpoint_filename[:i] + "." + \
                    str(n_finished_decoys) + "_checkpoint.pdb"
        elif n_finished_decoys == n_decoys and i != -4:
            new_checkpoint_filename = checkpoint_filename[:i] + ".pdb"
        if new_checkpoint_filename:
            os.rename(checkpoint_filename, new_checkpoint_filename)
            decoy_filenames_list[0] = new_checkpoint_filename
    for load_index, ranked_decoy_filename in enumerate(decoy_filenames_list[1:]):
        new_ranked_decoy_filename = ranked_decoy_filename[:ranked_decoy_filename[:-4].rfind(".")] + \
                "." + str(load_index + 2) + ".pdb"
        if new_ranked_decoy_filename != ranked_decoy_filename:
            os.rename(ranked_decoy_filename, new_ranked_decoy_filename)
            decoy_filenames_list[load_index + 1] = new_ranked_decoy_filename
    if len(decoy_scores_list) > 0:
        with open(output_filename_prefix + ".sc", "w") as pf:
            for ith_decoy, scores in enumerate(decoy_scores_list):
                if ith_decoy < len(decoy_filenames_list):
                    filename = decoy_filenames_list[ith_decoy]
                    scores = scores.copy()
                    scores["decoy"] = filename
                pf.write(json.dumps(scores) + "\n")
    for i_decoy in range(n_finished_decoys + 1, n_decoys + 1):
        if cloud_pdb_lines:
            if len(cloud_pdb_lines) > 2:
                rotamer_index = np.random.randint(1, len(cloud_pdb_lines) - 1)
            elif len(cloud_pdb_lines) == 2:
                rotamer_index = 1
            else:
                i_decoy = n_decoys
                break
            rotamer_lines = cloud_pdb_lines[rotamer_index]
            del cloud_pdb_lines[rotamer_index]
            tmp_pdb = output_filename_prefix + "." + str(i_decoy) + "_tmp.pdb"
            with open(tmp_pdb, "w") as p_pdb:
                p_pdb.writelines(cloud_pdb_lines[0])
                p_pdb.writelines(rotamer_lines)
            pose_copy = load_pdb_as_pose(score_function, tmp_pdb, pre_mutators, \
                    fold_tree, chi_dihedrals, symmetry, membrane_span_file, \
                    constraint_file, geometry_constraints, constraints, favor_native_residue)
            os.remove(tmp_pdb)
        else:
            pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        mutated_sequence = pose_copy.sequence()
        insert_index = 0
        current_score = energy_metric(score_function, pose_copy)
        for insert_index in range(len(decoy_scores_list)):
            if current_score < decoy_scores_list[insert_index]["total_score"]:
                break
            insert_index += 1
        scores = calculate_pose_scores(pose_copy, score_function, theozyme_positions=theozyme_positions)
        scores["total_score"] = current_score
        decoy_scores_list = decoy_scores_list[:insert_index] + [scores] + decoy_scores_list[insert_index:]
        current_output_filename = output_filename_prefix
        for index, wt_aa in enumerate(wildtype_sequence):
            aa = mutated_sequence[index]
            if aa != wt_aa:
                current_output_filename += "_" + wt_aa + \
                        pose_copy.pdb_info().pose2pdb(index + 1).split(" ")[0] + aa
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
                for move_index in range(len(decoy_filenames_list), insert_index, -1):
                    ranked_decoy_filename = decoy_filenames_list[move_index - 1]
                    new_ranked_decoy_filename = ranked_decoy_filename\
                            [:ranked_decoy_filename[:-4].rfind(".")] + "." + str(move_index + 1) + ".pdb"
                    os.rename(ranked_decoy_filename, new_ranked_decoy_filename)
                    decoy_filenames_list[move_index - 1] = new_ranked_decoy_filename
                current_output_filename += "." + str(insert_index + 1) + ".pdb"
            else:
                current_output_filename = None
        if current_output_filename:
            pose_copy.dump_pdb(current_output_filename)
            decoy_filenames_list = decoy_filenames_list[:insert_index] + [current_output_filename] + \
                    decoy_filenames_list[insert_index:]
            if len(decoy_filenames_list) > save_n_decoys:
                os.remove(decoy_filenames_list[save_n_decoys])
                decoy_filenames_list = decoy_filenames_list[:save_n_decoys]
        with open(output_filename_prefix + ".sc", "w") as pf:
            for ith_decoy, scores in enumerate(decoy_scores_list):
                if ith_decoy < len(decoy_filenames_list):
                    filename = decoy_filenames_list[ith_decoy]
                    scores = scores.copy()
                    scores["decoy"] = filename
                pf.write(json.dumps(scores) + "\n")

def run_job_distributor(score_function, pose, pre_mutators, fold_tree, chi_dihedrals:list, \
        symmetry:str, membrane_span_file:str, constraint_file:str, geometry_constraints, \
        constraints, favor_native_residue, movers, cloud_pdb_lines:list=None, \
        n_decoys:int=5, output_filename_prefix:str=None, wildtype_sequence:str=str()):
    if not output_filename_prefix:
        output_filename_prefix = os.path.basename(pose.pdb_info().name())\
                .rstrip(".pdb").rstrip(".tmp") + "_relaxed"
    job_distributor = PyJobDistributor(output_filename_prefix, n_decoys, score_function)
    while not job_distributor.job_complete:
        i_decoy = job_distributor.current_id
        if cloud_pdb_lines:
            if len(cloud_pdb_lines) > 2:
                rotamer_index = np.random.randint(1, len(cloud_pdb_lines) - 1)
            elif len(cloud_pdb_lines) == 2:
                rotamer_index = 1
            else:
                break
            rotamer_lines = cloud_pdb_lines[rotamer_index]
            del cloud_pdb_lines[rotamer_index]
            tmp_pdb = output_filename_prefix + "." + str(i_decoy) + "_tmp.pdb"
            with open(tmp_pdb, "w") as p_pdb:
                p_pdb.writelines(cloud_pdb_lines[0])
                p_pdb.writelines(rotamer_lines)
            pose_copy = load_pdb_as_pose(score_function, tmp_pdb, pre_mutators, \
                    fold_tree, chi_dihedrals, symmetry, membrane_span_file, \
                    constraint_file, geometry_constraints, constraints, favor_native_residue)
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
    if args.cloud_pdb: # Parse cloud pdb.
        cloud_pdb_lines = parse_cloud_pdb(args.pdb)
        args.n_decoys = min(args.n_decoys, len(cloud_pdb_lines[1:]))
        tmp_pdb = args.pdb.rstrip(".pdb") + ".tmp.pdb"
        with open(tmp_pdb, "w") as pf:
            pf.writelines(cloud_pdb_lines[0])
            pf.writelines(cloud_pdb_lines[1])
        pose = pose_from_pdb(tmp_pdb)
        os.remove(tmp_pdb)
    else:
        cloud_pdb_lines = None
        pose = pose_from_pdb(args.pdb)
    coord_ref_pose = None
    if args.coordinate_reference_pdb:
        coord_ref_pose = pose_from_pdb(args.coordinate_reference_pdb)
    ddG_ref_pose = None
    if args.ddG_reference_pdb:
        ddG_ref_pose = pose_from_pdb(args.ddG_reference_pdb)
    # Get sequences information.
    wildtype_sequence = pose.sequence()
    sequence_length = len(wildtype_sequence)
    if not args.output_filename_mutations_suffix:
        wildtype_sequence = str()
    # Convert pdb numberings to pose numberings.
    # Any ddG_ref_pose redundant positions in -static, -mut and -des will be ignored.
    static_pose_indices, _ = pdb_to_pose_numbering(pose, args.static_residues)
    if args.static_residue_identities:
        static_pose_indices.update(residue_name3_selector(pose, args.static_residue_identities, \
                sequence_length=sequence_length))
    pre_mutations, _ = pdb_to_pose_numbering(pose, args.pre_mutations)
    mutations, _ = pdb_to_pose_numbering(pose, args.mutations)
    mutations = mutations - pre_mutations
    design_pose_indices, _ = pdb_to_pose_numbering(pose, args.design_residues)
    # Classify NCAA mutations and static mutations into pre-mutations.
    static_ncaa_mutations = set()
    canonical_AAs = set("ACDEFGHIKLMNPQRSTVWY")
    for mutation in mutations:
        mutating_position, target_AA = mutation.split(",")
        if mutating_position in static_pose_indices or not target_AA in canonical_AAs:
            static_ncaa_mutations.add(mutation)
    mutations = mutations - static_ncaa_mutations
    pre_mutations.update(static_ncaa_mutations)
    # Apply pre-mutations.
    pre_mutators, pre_mutation_pose_indices = create_residue_mutators(pre_mutations, \
            noncanonical_amino_acids=args.noncanonical_amino_acids, ddG_wildtype=args.ddG_wildtype)
    for pre_mutator in pre_mutators:
        pre_mutator.apply(pose)
    # Set fold tree.
    if not args.fold_tree and not args.alter_jump_edges:
        fold_tree = pose.fold_tree()
    else:
        if args.fold_tree:
            fold_tree = create_fold_tree(args.fold_tree)
        elif args.alter_jump_edges:
            fold_tree = alter_fold_tree_jump_edges(pose.fold_tree(), args.alter_jump_edges)
        pose.fold_tree(fold_tree)
    # Set chi dihedrals.
    set_chi_dihedral(pose, args.chi_dihedrals)
    # Select residue by name3 reference pose.
    if args.ddG_reference_pdb:
        res_name3_ref_pose = ddG_ref_pose
        enzdes_ref_pdb = args.ddG_reference_pdb
    else:
        res_name3_ref_pose = pose
        enzdes_ref_pdb = args.pdb
    rigid_body_tform_pose_indices = set()
    # Get substrates and catalytic residues pose indices.
    # Theozyme_pose_indices will include ddG_ref_pose redundant positions.
    theozyme_pose_indices = set()
    if args.catalytic_residues:
        catalytic_residue_pose_indices, _ = pdb_to_pose_numbering(pose, args.catalytic_residues)
        theozyme_pose_indices.update(catalytic_residue_pose_indices)
    if args.catalytic_residue_identities:
        cat_id_pose_indices = residue_name3_selector(res_name3_ref_pose, args.catalytic_residue_identities)
        theozyme_pose_indices.update(cat_id_pose_indices)
    if args.substrates:
        substrate_pose_indices, jump_pose_indices = pdb_to_pose_numbering(pose, args.substrates)
        theozyme_pose_indices.update(substrate_pose_indices)
        if args.substrate_rigid_body_transformations:
            rigid_body_tform_pose_indices.update(jump_pose_indices)
    if args.substrate_identities:
        substrate_id_pose_indices = residue_name3_selector(res_name3_ref_pose, args.substrate_identities)
        theozyme_pose_indices.update(substrate_id_pose_indices)
        if args.substrate_rigid_body_transformations:
            rigid_body_tform_pose_indices.update(substrate_id_pose_indices)
    # Get EnzDes positions pose indices.
    # Enzdes_pose_indices will include ddG_ref_pose redundant positions.
    enzdes_pose_indices = set()
    enzdes_substrate_pose_indices, enzdes_res_pose_indices = get_enzdes_pose_indices(\
            res_name3_ref_pose.pdb_info(), enzdes_ref_pdb, args.symmetry)
    enzdes_pose_indices = enzdes_substrate_pose_indices.union(enzdes_res_pose_indices)
    if args.enzdes_substrates_transformations:
        rigid_body_tform_pose_indices.update(enzdes_substrate_pose_indices)
    # Exclude ddG_ref_pose redundant positions and static positions.
    pre_minimization_pose_indices = set(filter(lambda index: int(index) \
            <= sequence_length, theozyme_pose_indices.union(enzdes_pose_indices)))
    pre_minimization_pose_indices = pre_minimization_pose_indices - static_pose_indices
    rigid_body_tform_pose_indices = set(filter(lambda index: int(index) \
            <= sequence_length, rigid_body_tform_pose_indices))
    rigid_body_tform_pose_indices = rigid_body_tform_pose_indices - static_pose_indices
    # Rigid body transformations.
    rigid_body_tform_jump_edges = set()
    if len(rigid_body_tform_pose_indices) > 0:
        rigid_body_tform_jump_edges = pose_indices_to_jump_edges(fold_tree, \
                rigid_body_tform_pose_indices)
    # Apply symmetry if specified.
    n_monomers = apply_symmetry_membrane(args.symmetry, args.membrane_span_file, \
            pose, coord_ref_pose, ddG_ref_pose, sequence_length=sequence_length)
    # Read constraint files from the command line and apply to pose.
    if args.constraint_file:
        add_fa_constraints_from_cmdline(pose, score_function)
    # Add geometry constraints on-the-fly.
    geometry_constraints = list()
    if args.distance_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, \
                args.distance_constraint_atoms, args.distance_constraint_parameters, \
                bounded_distances=args.bounded_distances))
    if args.angle_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, \
                args.angle_constraint_atoms, args.angle_constraint_parameters))
    if args.dihedral_constraint_atoms:
        geometry_constraints.extend(create_constraints(pose, \
                args.dihedral_constraint_atoms, args.dihedral_constraint_parameters))
    for geometry_constraint in geometry_constraints:
        pose.add_constraint(geometry_constraint)
    # Create the task factory.
    task_factory, min_shell_focus_selection, packer_sampling_selection = \
            create_task_factory(specified_static_positions=static_pose_indices, \
            pre_mutation_positions=pre_mutation_pose_indices, \
            point_mutations=mutations, design_positions=design_pose_indices, \
            theozyme_positions=theozyme_pose_indices, \
            enzdes_positions=enzdes_pose_indices, \
            design_binding_site=args.design_binding_site, \
            design_enzdes_shell=args.design_enzdes_shell, \
            repack_neighborhood_only=args.repack_neighborhood_only, \
            repack_binding_site=not args.no_repack_binding_site, \
            repack_enzdes_shell=not args.no_repack_enzdes_shell, \
            ddG_ref_pose=ddG_ref_pose, n_monomers=n_monomers, \
            sequence_length=sequence_length, ddG_wildtype=args.ddG_wildtype, \
            excluded_amino_acid_types=args.excluded_amino_acid_types, \
            noncanonical_amino_acids=args.noncanonical_amino_acids, \
            allow_ncaa_in_design=args.allow_ncaa_in_design)
    # List of coordinate constraints and EnzDes constraints.
    constraints = list()
    # Add coordinate constraints.
    no_coord_cst_selection = OrResidueSelector()
    no_coord_cst_residues, _ = pdb_to_pose_numbering(pose, \
            args.no_backbone_coordinate_constraints_residues)
    no_coord_cst_selection.add_residue_selector(ResidueIndexSelector(\
            ",".join(no_coord_cst_residues) + ","))
    if args.no_coordinate_constraints_on_packing_region:
        no_coord_cst_selection.add_residue_selector(packer_sampling_selection)
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
    # Add enzyme design constraints.
    if args.enzyme_design_constraints:
        enzdes_cst = AddOrRemoveMatchCsts()
        enzdes_cst.set_cst_action(ADD_NEW)
        enzdes_cst.apply(pose)
        constraints.append(enzdes_cst)
    # Add AA type constraints.
    if args.favor_native_residue and (args.design_residues or args.design_binding_site):
        FavorNativeResidue(pose, args.favor_native_residue)
    # List of movers.
    movers = list()
    # Create the RMSD metric.
    no_rmsd_residues, _ = pdb_to_pose_numbering(pose, args.no_rmsd_residues)
    rmsd_metric = RMSDMetric(pose, NotResidueSelector(ResidueIndexSelector(\
            ",".join(no_rmsd_residues) + ",")))
    # Perform pre-minimization.
    assembly_length = len(pose.sequence())
    if args.pre_minimization and len(pre_minimization_pose_indices) > 0:
        pre_minimization_move_map = create_pre_minimization_move_map(assembly_length, \
                pre_minimization_pose_indices, jump_edges=rigid_body_tform_jump_edges)
        pre_minimizer = create_minimizer(score_function, move_map=pre_minimization_move_map)
        movers.append(pre_minimizer)
    # Create the move map.
    if args.minimize_neighborhood_only:
        if not args.repack_neighborhood_only:
            raise Exception("Using -min_nbh is not allowed without using -rpk_nbh at the same time.")
    else: # Minimize the whole pose.
        min_shell_focus_selection = TrueResidueSelector()
    move_map = create_move_map(min_shell_focus_selection, static_positions=static_pose_indices, \
            ddG_ref_pose=ddG_ref_pose, n_monomers=n_monomers, sequence_length=sequence_length, \
            assembly_length=assembly_length, jump_edges=rigid_body_tform_jump_edges)
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
        run_job(score_function, pose, pre_mutators, fold_tree, args.chi_dihedrals, args.symmetry, \
                args.membrane_span_file, args.constraint_file, geometry_constraints, \
                constraints, args.favor_native_residue, movers, cloud_pdb_lines=cloud_pdb_lines, \
                n_decoys=args.n_decoys, save_n_decoys=args.save_n_decoys, \
                theozyme_positions=pre_minimization_pose_indices, \
                output_filename_prefix=args.output_filename_prefix, \
                wildtype_sequence=wildtype_sequence)
    else:
        run_job_distributor(score_function, pose, pre_mutators, fold_tree, args.chi_dihedrals, \
                args.symmetry, args.membrane_span_file, args.constraint_file, geometry_constraints, \
                constraints, args.favor_native_residue, movers, cloud_pdb_lines=cloud_pdb_lines, \
                n_decoys=args.n_decoys, output_filename_prefix=args.output_filename_prefix, \
                wildtype_sequence=wildtype_sequence)


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    main(args)
