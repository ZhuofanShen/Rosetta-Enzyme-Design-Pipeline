# Rosetta Enzyme Design Pipeline

A Rosetta/PyRosetta–based workflow for designing and optimizing enzyme active sites.

The core of the workflow is a PyRosetta module/script **`enzdes_utils.py`**. All operations—including introducing point mutations and optimizing protein–substrate complexes—are performed by calling this script.  
Example projects are included in the repository (e.g., `hemoprotein_biocatalysts_design/` and `protein_stapling/`). 

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Repository Layout](#repository-layout)
- [Requirements](#requirements)
- [Installation](#installation)
- [Prepare Your Inputs](#prepare-your-inputs)
- [Quickstart](#quickstart)
- [Typical Workflows](#typical-workflows)
- [Inputs & Outputs](#inputs--outputs)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Overview

This pipeline wraps commonly used Rosetta enzyme‑design steps into a practical, reproducible workflow. It focuses on two tasks:

1. **Design around a bound substrate/ligand** using enzyme‑design constraints (the “EnzDes” application family in Rosetta). 
2. **Systematic mutational exploration** (single or combinatorial) with local packing/minimization and optional sequence design.

Where helpful, it also accommodates the classic **Dock → Design** pattern: use **RosettaMatch** or other docking practices to place catalytic geometries (the “theozyme”) on a scaffold, then redesign around the placed geometry with **FastDesign** or **LigandMPNN**.

---

## Key Features

- **Single or batched point mutations** at specified positions, with local repacking/relax and scoring.
- **Ligand‑aware design**: read ligand `.params`, apply enzyme‑design geometric constraints (`.cst`), and optimize the active site.
- **Design cycles**: alternate packing/minimization and sequence optimization (FastDesign‑style) around the active site (where applicable).
- **Results you can filter**: write out structures and tabular scores so you can rank designs by stability, interface energy, constraint satisfaction, etc.
- **Scriptable with PyRosetta** for custom pipelines; designed to be hackable.

---

## Repository Layout
Rosetta-Enzyme-Design-Pipeline/

├─ enzdes_utils/ # Core PyRosetta utilities, exposes the main entry point (enzdes_utils.py)

├─ hemoprotein_biocatalysts_design/ # Example project: heme protein design (inputs, scripts)

├─ protein_stapling/ # Example project: staple/design utilities (inputs, scripts)

└─ README.md


---

## Requirements

- **Python 3.8+** (tested with modern 3.11)
- **PyRosetta 4** (download & install from the official site; see below).
- **Rosetta command‑line tools (optional)** if you plan to run **Match** externally or mix in RosettaScripts XML flows.

> **Note on PyRosetta**: PyRosetta is distributed as licensed wheels/tarballs. Install from the official PyRosetta website following their instructions (Windows/macOS/Linux guidance provided).

---

## Installation

1. **Clone this repo**
   ```bash
   git clone https://github.com/ZhuofanShen/Rosetta-Enzyme-Design-Pipeline.git
   cd Rosetta-Enzyme-Design-Pipeline
2. **Create a clean Python environment (example with conda)**
   ```bash
   conda create -n rosetta-enzdes python=3.10 -y
   conda activate rosetta-enzdes

## Prepare Your Inputs

1. **Protein-substrate complex structures**

Dock the substrate/TS into the protein active site using selected docking software or RosettaMatch, or simply overlay the TS onto the cofactor. The `enzdes_utils.py` script can explicitly consider substrate/TS isomers and conformers by loading a cloud PDB file following a specific format.

2. **Ligand params & conformers**

Generate a ligand .params file (and optional rotamer library PDB) with molfile_to_params.py from Rosetta/PyRosetta tooling. Use --parameters_files to load the params at run time. 

3. **Enzyme‑design constraints**

Recommended: Use --distance/angle/dihedral_constraint_atoms along with --distance/angle/dihedral_constraint_parameters to define desired catalytic geometries, i.e., distances/angles/dihedrals.
Alternative: Prepare an enzyme‑design .cst (Match/EnzDes geometric constraint) file that encodes desired catalytic geometries.

Rosetta’s official docs provide short, practical tutorials on preparing ligands and constraints.

## Quickstart

The exact CLI flags for enzdes_utils.py may differ—run python enzdes_utils/enzdes_utils.py -h to list all supported options in your version.

Example: Design around a bound ligand with catalytic constraints
From the repository root:
python enzdes_utils/enzdes_utils.py \
  pdb                                  inputs/scaffold_with_ligand.pdb \
  or
  --cloud_pdb                          inputs/scaffold_with_ligand_conformers.pdb \
  --index_reference_pdb                inputs/scaffold_with_ligand_stereoisomers_and_conformers.pdb \
  --parameters_files                   ligands/LIG.params \
  --score_function                     ref2015_cst \
  --score_terms                        fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1 \
  --fold_tree                          1,153,-1 1,154,1 154,155,FE,N1 \
  --alter_jump_edges                   [-1,-2,FE,N1] \
  --dihedral_constraint_atoms          LIG,N1 LIG,FE LIG,CARB LIG,CEST \
  --dihedral_constraint_parameters     0,180,34.95 \
  --mutations                          A42,D A85,G \
  --catalytic_residues                 A42 A67-69 A85 \
  or
  --design_residues                    A42 A67-69 A85 \
  or
  --design_binding_site \
  --substrate_identities               LIG \
  --repack_neighborhood_only \
  --minimize_neighborhood_only \
  --substrate_rigid_body_transformations \
  --n_decoys                           20 \
  --save_n_decoys                      5 \
  --output_filename_prefix             results/run_01
What this does (typical behavior):

Loads the input complex, ligand params, and enzyme‑design constraints.

Applies the requested point mutations (single or multiple).

Runs a pack/minimize/design cycle around the active site.

Writes designed models and a score table into `results/` with a filename prefix "run_01".

## Typical Workflows
1) De novo active‑site placement → design

Match: Place catalytic residues and the ligand geometry on a scaffold using Rosetta Match with your .cst file. 
Rosetta 

Design: Use this pipeline (or Rosetta FastDesign/RosettaScripts) to optimize sequences and side chains around the placed theozyme. 

2) Mutational scanning around a known active site

Provide a list of single mutations (e.g., A42D,A85G) or a small Cartesian product.

For each variant, the script repacks/minimizes locally and scores; output is a table plus designed structures.

3) Substrate/transition‑state analog optimization

Load ligand .params and a .cst describing catalytic geometry.

Run design cycles to improve binding and geometry compliance (e.g., constraint scores). 

## Inputs & Outputs

Inputs

PDB: protein (optionally with ligand already placed).

Ligand .params and (optionally) conformer PDB library. 

.cst file: enzyme‑design / matcher constraints. 

Mutation list and/or resfile (if your build supports it).

Outputs

Designed structures (PDB or silent files) for each job.

Tabular scores (CSV/TSV, e.g., total score, InterfaceDelta, constraint energies).

Logs describing run settings and Rosetta/PyRosetta versions.

## Troubleshooting
run `enzdes_utils.py` with an additional flag --debug_mode.

## Citations

Rosetta Enzyme Design (EnzDes) docs — overview and usage. 

Rosetta Match — place catalytic geometries (theozyme) on scaffolds. 

Constraint (.cst) file format — match/enzdes geometric constraints. 

Ligand preparation — generating .params with molfile_to_params.py. 

PyRosetta installation — official guidance & platform notes. 

Background reading: Richter et al. “De Novo Enzyme Design Using Rosetta3.” Protein Sci. (2011). 

## License

The code is released under the MIT License.

## Acknowledgements

Built on Rosetta and PyRosetta; thanks to RosettaCommons for software and documentation. 

Repository authors and contributors listed in GitHub insights.
