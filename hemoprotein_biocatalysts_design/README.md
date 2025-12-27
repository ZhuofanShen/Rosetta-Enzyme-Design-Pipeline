# Rosetta Enzyme Design Pipeline

This repository implements a **Rosetta / PyRosetta–based enzyme design workflow** developed for **heme-dependent biocatalysts**, with a focus on:

- Raw .pdb files preprocessing
- Transition-state (TS)–guided Rosetta **FastDesign** and **FastRelax** design
- **Multi-State Design (MSD)** for stereoselectivity optimization  

The pipeline was used to generate the data reported in our *Nature Communications* study and is designed as a **step-wise, restartable workflow** suitable for large-scale **Slurm** execution.

---

## Repository Structure

After cloning the repository, place all pipeline scripts in a single `scripts/` folder and raw PDB files in `pdb/`:

```text
$HOME/Rosetta-Enzyme-Design-Pipeline/
├── enzdes_utils/
│   └── fast_design.py
├── hemoprotein_biocatalysts_design/
│   └── pdb/
│       └── ${PDB_ID}.pdb
└── scripts/


enzdes_utils/
Core PyRosetta utilities (FastDesign, relax, scoring).

scripts/
Pipeline driver scripts and Slurm submission templates.

hemoprotein_biocatalysts_design/
Project-specific working directory (inputs, intermediates, outputs).
```

## Step 1 — Preprocess Raw PDBs and WT Relaxation
**1.1 Clean and preprocess raw PDB files**

From `hemoprotein_biocatalysts_design/`, run:

`python clean_hemoprotein_pdb.py -i pdb -o processed_pdb -link`

`python classify_hemoprotein_folds.py`

This step:

Cleans raw PDB files

Standardizes atom naming
Classifies hemoproteins into fold families

**1.2 Resulting directory structure**
```text
hemoprotein_biocatalysts_design/
├── pdb/
├── processed_pdb/
│   └── ${PDB_ID}.pdb
├── HEA_monomer/
├── HEB_monomer/
├── HEC_monomer/
├── HEM_monomer/
│   └── ${PDB_ID}/
│       ├── ${PDB_ID}_clean.pdb
│       └── DIOXYGENASE.fold
├── HEM_homo-oligomer/
├── HEM_hetero-oligomer_mono-heme/
└── ...
```
**1.3 WT pre-relaxation**

Submit Slurm jobs:
`sbatch scripts/relax_raw_WT_pdb.pbs`

After all jobs finish, select the best relaxed WT structure:
`python scripts/get_best_relaxed_decoy.py HEM_monomer -s 1`

-s 1 indicates step-1 WT relaxation, not variant–TS relaxation.

Final WT structure layout:
```text
HEM_monomer/
└── ${PDB_ID}/
    ├── ${PDB_ID}_clean.pdb
    ├── ${PDB_ID}_relaxed.pdb
    └── DIOXYGENASE.fold
```

## Step 2 — Superimpose TS Stereoisomers and Conformers

Align TS stereoisomers and conformers onto the heme cofactor of the pre-relaxed WT scaffold:
```bash
python scripts/generate_hemoprotein_substrate_complexes.py \
    HEM \
    -f DIOXYGENASE \
    -sub substrates/cyclopropanation_styrene_EDA
```

`-f/--fold` is optional.

Output:
```text
HEM/
└── ${PDB_ID}/
    ├── ${PDB_ID}_clean.pdb
    ├── ${PDB_ID}_relaxed.pdb
    ├── DIOXYGENASE.fold
    └── cyclopropanation_styrene_EDA_distal/
        └── complexes/
            └── ${PDB_ID}_{stereoisomer}_{conformer}.pdb
```

These structures serve as WT–TS complex starting models.

## Step 3 — FastDesign of WT Scaffolds Toward TS Binding

Generate FastDesign Slurm scripts:
```bash
python scripts/generate_design_slurm_scripts.py \
    HEM \
    -f DIOXYGENASE \
    -stereo 1R2S \
    -sub substrates/cyclopropanation_styrene_EDA \
    -n 5
```

If `-f/--fold` is omitted, the default script name is default.sh.

Execute all generated scripts:

```bash
for path in `ls HEM/*/cyclopropanation_styrene_EDA_*/*_FastDesign/*rot*/DIOXYGENASE.sh`; do
    cd ${path%/*}
    bash ${path##*/}
    sleep 0.01
    cd ../../../../..
done
```

After all jobs finish, collect design scores:
```python
python scripts/generate_design_scores_table.py \
    HEM \
    -fold DIOXYGENASE \
    -stereo 1R2S
```

Output:

`HEM_DIOXYGENASE_1R2S.csv`

Corresponds to **Supplementary Data 2** in the paper

## Step 4 — FastRelax of Designed Variants in TS Complexes

Generate relax jobs:
```bash
python scripts/fast_design_to_relax.py \
    HEM_DIOXYGENASE_1R2S.csv \
    -d HEM \
    -sub substrates/cyclopropanation_styrene_EDA \
    -stereo 2 \
    -th 0
```

Run all relax scripts:
```bash
for path in `ls HEM/*/cyclopropanation_styrene_EDA_*/FastDesign_Relax/run_generate_relax_slurm_scripts.sh`; do
    cd ${path%/*}
    bash ${path##*/}
    sleep 0.01
    cd ../../../..
done

sbatch fast_design_relax.pbs
```

After all relax jobs complete:
```python
python scripts/generate_relax_scores_table_multi_scaf.py \
    HEM \
    -stereo 2 \
    -sub cyclopropanation_styrene_EDA \
    -rank s
```

Output:

`HEM_stereoselectivity.csv`

Corresponds to **Supplementary Data 3**

Extract best relaxed decoys:
```bash
cd HEM/${PDB_ID}/cyclopropanation_styrene_EDA_distal
python ../../../../scripts/get_best_relaxed_decoys.py FastDesign_Relax -s 2
```

-s 2 indicates variant–TS relaxation.

## Step 5 — Multi-State Design (MSD)

Prepare MSD inputs:
```bash
mkdir Multi-State_Design
mkdir Multi-State_Design/${PDB_ID}_${mutations}/input
mv ${PDB_ID}_{stereoisomer}_{conformer}.pdb \
   Multi-State_Design/${PDB_ID}_${mutations}/input
```

Edit mutation sites:

`vi ../../../../scripts/site-saturation_mutagenesis.pbs`


Submit MSD jobs:

`sbatch ../../../../scripts/site-saturation_mutagenesis.pbs ${PDB_ID}_${mutations}`


After completion:

`python ../../../../scripts/multistate_design_fitness.py ${PDB_ID}_${mutations}`


Outputs:
```text
${PDB_ID}_${mutations}/${PDB_ID}_${stereoisomer}/
└── ${PDB_ID}_${mutations}_${stereoisomer}.{fitness_type}.csv
```

Generate sequence logos:

`python logo_plot.py ${PDB_ID}_${mutations}_${stereoisomer}.{fitness_type}.csv`


Use the logo plots to propose new variants
(e.g. 6F0A_C129A_S167V_R231L).

## Step 6 — Iterative FastRelax ↔ MSD Cycles

For a new variant:
```bash
cd HEM/${PDB_ID}/cyclopropanation_styrene_EDA_distal
mkdir FastRelax
cd FastRelax
```

Generate relax scripts:
```bash
python ../../../../scripts/generate_relax_slurm_scripts.py \
    6F0A_C129A_S167V \
    -muts A129,A A167,V \
    -script relax \
    -n 30
```

Submit jobs:

sbatch fast_relax_6M8F.pbs


Analyze results:
```bash
python ../../../../scripts/generate_relax_scores_table.py \
    -baseline 6F0A_WT \
    -stereo 2
```

Restart MSD from a selected variant:
```python
python ../../../../scripts/get_best_relaxed_decoys.py \
    6F0A_C129A_S167V_R231L \
    -s 2
```

Repeat Step 5 for further refinement.
