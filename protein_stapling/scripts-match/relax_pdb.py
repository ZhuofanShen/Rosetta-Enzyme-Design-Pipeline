#!/usr/bin/python3
import argparse
from fast_design import parse_arguments, init_pyrosetta_with_opts, create_score_function, create_coord_cst
import os
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
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
    add_fa_constraints_from_cmdline, ConstraintSet, \
    AmbiguousConstraint, AngleConstraint, AtomPairConstraint, DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import CircularHarmonicFunc, FlatHarmonicFunc, HarmonicFunc
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, ChainSelector, \
    ResidueIndexSelector, ResidueNameSelector, ResiduePropertySelector, \
    NeighborhoodResidueSelector, InterGroupInterfaceByVectorSelector
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


def run_job(pose, sfxn, decoys, pdb_file_name, *movers):
    best_score = None
    for _ in range(decoys):
        pose_copy = Pose(pose)
        for mover in movers:
            mover.apply(pose_copy)
        current_score = sfxn(pose_copy)
        if not best_score or current_score < best_score:
            best_score = current_score
            pose_copy.dump_pdb(pdb_file_name[:-3] + "relaxed.pdb")


if __name__ == "__main__":
    args = parse_arguments()
    init_pyrosetta_with_opts(args)
    sfxn = create_score_function(args.score_function)
    # define movers
    movers = list()
    # coordinate constraint
    coord_cst_gen = create_coord_cst()
    add_csts = AddConstraints()
    add_csts.add_generator(coord_cst_gen)
    movers.append(add_csts)
    # rmsd
    pose = pose_from_pdb(args.pdb)
    rmsdm = RMSDMetric(pose)
    movers.append(rmsdm)
    # create task factory
    tf = standard_task_factory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(IncludeCurrent())
    tf.push_back(ExtraRotamers(0, 1, 1))
    tf.push_back(ExtraRotamers(0, 2, 1))
    fr = FastRelax()
    fr.set_scorefxn(sfxn)
    fr.set_task_factory(tf)
    movers.append(fr)
    # print(tf.create_task_and_apply_taskoperations(pose))
    run_job(pose, sfxn, args.decoys, args.pdb, *movers)
