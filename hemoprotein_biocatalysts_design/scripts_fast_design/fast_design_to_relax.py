import argparse
import os
import pandas as pd


def pdb_2_seq_dict(pdb_path, restype_3to1):
    position_aa_dict = dict()
    with open(pdb_path, "r") as pf:
        for line in pf:
            current_resi = None
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                resi = int(line[22:26])
                aa_type = restype_3to1.get(line[17:20])
                if resi != current_resi and aa_type:
                    position_aa_dict[line[21] + str(resi)] = aa_type
                    current_resi = resi
    return position_aa_dict

restype_3to1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', type=str)
    parser.add_argument('-d', '--directory', type=str, default="HEM")
    parser.add_argument("-sub", "--substrate", type=str, default="substrates/cyclopropanation_styrene_EDA")
    parser.add_argument('-stereo', '--preferred_stereoisomer', type=int, choices=[0, 1, 2, 3], default=2)
    parser.add_argument('-th', '--threshold', type=float, default=0)
    args = parser.parse_args()
    
    args.substrate = args.substrate.strip("/")
    substrate = args.substrate.split("/")[-1]

    columns = ["PDB_ID", "WT"]
    stereoisomer = ["RR", "SS", "RS", "SR"][args.preferred_stereoisomer]
    for rot in range(1, 5):
        for ester in ["p", "m"]:
            prefix = stereoisomer + str(rot) + ester + "_"
            columns.append(prefix + "TS")
            columns.append(prefix + "ddG")
    df = pd.read_csv(args.csv, usecols=columns)
    for _, row in df.iterrows():
        pdb = row.PDB_ID
        pdb_path = args.directory + "/" + pdb + "/" + pdb + "_relaxed.pdb"
        if os.path.isfile(pdb_path):
            wt_seq_dict = pdb_2_seq_dict(pdb_path, restype_3to1)
        else:
            with open(args.directory + ".err", "a") as pf:
                pf.write(pdb + "\n")
            continue
        for i in range(2, 18, 2):
            exec("ts_sc = row." + columns[i])
            exec("ddg = row." + columns[i+1])
            if not ts_sc < args.threshold or not ddg < args.threshold:
                continue
            isomer = columns[i].split("_")[0]
            if isomer[3] == "p":
                ester = "+"
            elif isomer[3] == "m":
                ester = "-"
            isomer = "1" + isomer[0] + "2" + isomer[1] + "-rot" + isomer[2] + ester
            design_traj_path = args.directory + "/" + pdb + "/" + substrate + "_proximal/FastDesign/" + pdb + "_" + isomer
            if_proximal = False
            if os.path.isdir(design_traj_path):
                for design_traj_pdb in os.listdir(design_traj_path):
                    if design_traj_pdb.startswith(pdb + "_" + isomer) and design_traj_pdb.endswith(".pdb"):
                        if_proximal = True
                        break
                if if_proximal:
                    seq_dict = pdb_2_seq_dict(design_traj_path + "/" + design_traj_pdb, restype_3to1)
                    mutations = list()
                    for chain_pos, aa_type in wt_seq_dict.items():
                        mut_aa_type = seq_dict.get(chain_pos)
                        if not mut_aa_type or aa_type == mut_aa_type:
                            continue
                        mutations.append(chain_pos + "," + mut_aa_type)
                    if len(mutations) > 0:
                        mutations = " -muts " + " ".join(mutations)
                    else:
                        mutations = ""
                    relax_path = args.directory + "/" + pdb + "/" + substrate + "_proximal/FastDesign_Relax"
                    if not os.path.isdir(relax_path):
                        os.mkdir(relax_path)
                    with open(relax_path + "/run_generate_relax_slurm_scripts.sh", "a") as pf:
                        pf.write("python ../../../../scripts/generate_relax_slurm_scripts.py " + pdb + "_" + isomer + " -sub ../../../../" + args.substrate + mutations + " -script screening -n 5 -t 8\n")
            design_traj_path = args.directory + "/" + pdb + "/" + substrate + "_distal/FastDesign/" + pdb + "_" + isomer
            if_distal = False
            if os.path.isdir(design_traj_path):
                for design_traj_pdb in os.listdir(design_traj_path):
                    if design_traj_pdb.startswith(pdb + "_" + isomer) and design_traj_pdb.endswith(".pdb"):
                        if_distal = True
                        break
                if if_distal:
                    seq_dict = pdb_2_seq_dict(design_traj_path + "/" + design_traj_pdb, restype_3to1)
                    mutations = list()
                    for chain_pos, aa_type in wt_seq_dict.items():
                        mut_aa_type = seq_dict.get(chain_pos)
                        if not mut_aa_type or aa_type == mut_aa_type:
                            continue
                        mutations.append(chain_pos + "," + mut_aa_type)
                    if len(mutations) > 0:
                        mutations = " -muts " + " ".join(mutations)
                    else:
                        mutations = ""
                    relax_path = args.directory + "/" + pdb + "/" + substrate + "_distal/FastDesign_Relax"
                    if not os.path.isdir(relax_path):
                        os.mkdir(relax_path)
                    with open(relax_path + "/run_generate_relax_slurm_scripts.sh", "a") as pf:
                        pf.write("python ../../../../scripts/generate_relax_slurm_scripts.py " + pdb + "_" + isomer + " -sub ../../../../" + args.substrate + mutations + " -script screening -n 5 -t 8\n")
            if if_proximal and if_distal:
                with open(args.directory + ".log", "a") as pf:
                    pf.write(pdb + "_" + isomer + "\n")
