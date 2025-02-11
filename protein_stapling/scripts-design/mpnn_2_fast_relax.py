import argparse
import json
import os
import shutil


xlink_name3_dict = {"O2beY": "OBY"}
uaa_name3_dict = {"O2beY": "TYZ"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str)
    parser.add_argument('-lig', "--ligand", type=str, default="O2beY-CYS")
    parser.add_argument('--process_variants', type=str, nargs="*", default=list()) # e.g., 1ysb-A_S6_E126
    parser.add_argument('-s', "--step", type=int, choices=[1, 2], required=True)
    parser.add_argument('-od', '--output_directory', type=str, default='mpnn')
    parser.add_argument('-esm_n', '--best_n_esm_scored_sequences', type=int)
    parser.add_argument('-n', "--n_decoys", type=int, default=20)
    args = parser.parse_args()

    uaa_type, nucleophile_type = args.ligand.split("-")

    # improper dihedral 1
    dihe_atoms="-dihe_atoms OBY,O1 OBY,C1 OBY,C2 TYZ,OH"
    dihe_params="-dihe_params 0,1"

    if nucleophile_type == "ASP":
        nucleophile_name3="ASX"
        # improper dihedral 2
        improper2=" OBY,S1 OBY,C2 OBY,C1 ASX,OD2"
        
        # dihedral 2 AB
        dihe_atoms+=" OBY,C1 OBY,C2 ASX,OD2 ASX,CG"
        dihe_params+=" 180,59,80,35,280,35"
        # dihedral 2 B
        dihe_atoms+=" OBY,C2 ASX,OD2 ASX,CG ASX,CB"
        dihe_params+=" 180,20"
    elif nucleophile_type == "GLU":
        nucleophile_name3="GLX"
        # improper dihedral 2
        improper2=" OBY,S1 OBY,C2 OBY,C1 GLX,OE2"
        # dihedral 2 AB
        dihe_atoms+=" OBY,C1 OBY,C2 GLX,OE2 GLX,CD"
        dihe_params+=" 180,59,80,35,280,35"
        # dihedral 2 B
        dihe_atoms+=" OBY,C2 GLX,OE2 GLX,CD GLX,CG"
        dihe_params+=" 180,20"
    elif nucleophile_type == "CYS":
        nucleophile_name3="CYX"
        # improper dihedral 2
        improper2=" OBY,S1 OBY,C2 OBY,C1 CYX,SG"
        # dihedral 2 AB
        dihe_atoms+=" OBY,C1 OBY,C2 CYX,SG CYX,CB"
        dihe_params+=" 180,42,70,33,290,33"
        # dihedral 2 B
        dihe_atoms+=" OBY,C2 CYX,SG CYX,CB CYX,CA"
        dihe_params+=" 180,42,70,33,290,33"
    elif nucleophile_type == "HIS":
        nucleophile_name3="HIX"
        # improper dihedral 2
        improper2=" OBY,S1 OBY,C2 OBY,C1 HIX,NE2"
        # # dihedral 2 AB
        # dihe_atoms+=" OBY,C1 OBY,C2 HIX,NE2 HIX,CD2"
        # dihe_params+=" 90,58,270,58"
        # dihedral 2 B
        dihe_atoms+=" OBY,C2 HIX,NE2 HIX,CD2 HIX,CG"
        dihe_params+=" 180,1"
    elif nucleophile_type == "LYS":
        nucleophile_name3="LYX"
        # improper dihedral 2
        improper2=" OBY,S1 OBY,C2 OBY,C1 LYX,NZ"
        # improper dihedral lysine
        dihe_atoms+=" OBY,C2 LYX,CE LYX,NZ LYX,3HZ"
        dihe_params+=" 0,1"
        # dihedral 2 AB
        dihe_atoms+=" OBY,C1 OBY,C2 LYX,NZ LYX,CE"
        dihe_params+=" 180,36,60,33,300,33"
        # dihedral 2 B
        dihe_atoms+=" OBY,C2 LYX,NZ LYX,CE LYX,CD"
        dihe_params+=" 180,36,60,33,300,33"

    # improper dihedral 2
    dihe_atoms+=improper2
    dihe_params+=" 0,1"
    # dihedral 1/2 A
    dihe_atoms+=" OBY,O1 OBY,C1 OBY,C2 OBY,S1"
    dihe_params+=" 180,37,70,33,290,33"
    # dihedral 1 B
    dihe_atoms+=" OBY,C1 TYZ,OH TYZ,CZ TYZ,CE1"
    dihe_params+=" 0,50,180,50"
    # dihedral 1 AB
    dihe_atoms+=" OBY,C2 OBY,C1 TYZ,OH TYZ,CZ"
    dihe_params+=" 180,36,70,33,290,33"

    uaa_name3 = uaa_name3_dict[uaa_type]
    xlink_name3 = xlink_name3_dict[uaa_type]

    params_files_list = ["ligands/" + args.ligand + "/" + uaa_name3 + ".params", 
                         "ligands/" + args.ligand + "/" + xlink_name3 + ".params", 
                         "ligands/" + args.ligand + "/" + nucleophile_name3 + ".params"]

    if args.step == 1:
        from pyrosetta import *
        from pymol import *
        init("-ex1 -ex2 -ignore_zero_occupancy false -use_input_sc -no_optH true " +
            "-extra_res_fa " + " ".join(params_files_list) + " " +
            " -enzdes:cstfile ligands/" + args.ligand + "/" + args.ligand + "_Relax.cst " +
            "-run:preserve_header")

    for pdb in os.listdir(args.directory):
        pdb_path = os.path.join(args.directory, pdb)
        for pdb_chain_ligand in filter(lambda x: x.startswith(pdb) and os.path.isdir(\
                os.path.join(pdb_path, x)), os.listdir(pdb_path)):
            pdb_chain, uaa_nucleophile = pdb_chain_ligand.split("_")
            pdb_ligand_path = os.path.join(pdb_path, pdb_chain_ligand)
            for variant in filter(lambda x: os.path.isdir(os.path.join(pdb_ligand_path, x)) and x.startswith(pdb_chain), os.listdir(pdb_ligand_path)):
                if len(args.process_variants) > 0 and variant not in args.process_variants:
                    continue
                variant_path = os.path.join(pdb_ligand_path, variant)
                mutations = variant.split("_")[1:]
                uaa_pose_index = int(mutations[0][1:])
                nucleophile_pose_index = int(mutations[1][1:])
                for variant_pdb in os.listdir(variant_path):
                    if variant_pdb.startswith("Z") and variant_pdb.endswith(".pdb"):
                        pdb_positions = variant_pdb.split(".")[0][1:].split("X")
                        variant_pdb_numbering = pdb + "_" + mutations[0][0] + pdb_positions[0] + "UAA_" + mutations[1][0] + pdb_positions[1]
                        break
                mpnn_path = os.path.join(variant_path, args.output_directory)
                mpnn_trajs_path = os.path.join(mpnn_path, "packed")
                if not os.path.isdir(mpnn_trajs_path):
                    continue
                for variant_pdb in os.listdir(variant_path):
                    if variant_pdb.startswith("Z") and variant_pdb.endswith(".pdb") and not variant_pdb.endswith(".MPNN.pdb"):
                        variant_pdb_prefix = variant_pdb.split(".")[0]
                        pdb_positions = variant_pdb_prefix[1:].split("X")
                        if args.step == 1:
                            pseudo_wt_pose = pose_from_pdb(os.path.join(variant_path, variant_pdb))
                            pseudo_wt_sequence = pseudo_wt_pose.sequence()
                        break
                if args.step == 1:
                    mpnn_trajs = list()
                    unique_seq_set = set()
                    with open(os.path.join(mpnn_path, "sequences.json"), "r") as pf:
                        seqs_list = json.load(pf)
                    if args.best_n_esm_scored_sequences:
                        for mpnn_traj in filter(lambda x: x.endswith("_1.pdb"), os.listdir(mpnn_trajs_path)):
                            prefix = mpnn_traj.split(".")[0]
                            break
                        with open(os.path.join(mpnn_path, "esm_scores.json"), "r") as pf:
                            total_logp_list = json.load(pf)
                        for mpnn_variant in sorted(zip(total_logp_list[1:], seqs_list[1:], list(range(1, len(total_logp_list) + 1))), key=lambda x: x[0], reverse=True):
                            if mpnn_variant[1] not in unique_seq_set:
                                mpnn_trajs.append(prefix + ".MPNN_packed_" + str(mpnn_variant[2]) + "_1.pdb")
                                unique_seq_set.add(mpnn_variant[1])
                                if len(unique_seq_set) == args.best_n_esm_scored_sequences:
                                    break
                    else:
                        for i_traj, mpnn_traj in enumerate(sorted(filter(lambda x: x.endswith("_1.pdb"), os.listdir(mpnn_trajs_path)), \
                                    key=lambda x: int(x.split("_")[-2]))):
                            if seqs_list[i_traj] not in unique_seq_set:
                                mpnn_trajs.append(mpnn_traj)
                                unique_seq_set.add(seqs_list[i_traj])
                    cmd.load(os.path.join("pdb" + args.directory[3:], pdb, pdb_chain + ".pdb"))
                    traj_mutations_list = list()
                    mutations_pdb_index_set = set()
                    mutations_chain_pdb_index_set = set()
                    for mpnn_traj in mpnn_trajs:
                        mpnn_traj_file = os.path.join(mpnn_trajs_path, mpnn_traj)
                        if nucleophile_name3 == "CYX":
                            with open(mpnn_traj_file, "r") as pf:
                                pdb_lines = pf.readlines()
                            with open(mpnn_traj_file, "w") as pf:
                                for pdb_line in pdb_lines:
                                    pf.write(pdb_line.replace("HIX", "CYX"))
                        cmd.load(mpnn_traj_file, mpnn_traj[:-4])
                        cmd.save(mpnn_traj_file, mpnn_traj[:-4])
                        cmd.delete(mpnn_traj[:-4])
                        mpnn_traj_pose = pose_from_pdb(mpnn_traj_file)
                        mpnn_traj_sequence = mpnn_traj_pose.sequence()
                        mpnn_traj_pose.dump_pdb(mpnn_traj_file)
                        mutations_pdb_index_aa = list()
                        for i in range(len(mpnn_traj_sequence)):
                            if mpnn_traj_sequence[i] != pseudo_wt_sequence[i]:
                                info = pseudo_wt_pose.pdb_info().pose2pdb(i + 1).split(" ")
                                mutations_pdb_index_set.add(int(info[0]))
                                if i+1 != nucleophile_pose_index and i+1 != uaa_pose_index and \
                                        mpnn_traj_sequence[i] != "X" and mpnn_traj_sequence[i] != "Z":
                                    mutations_pdb_index_aa.append(info[1] + info[0] + "," + mpnn_traj_sequence[i])
                                    mutations_chain_pdb_index_set.add(info[1] + info[0])
                        mutations_str = str()
                        if len(mutations_pdb_index_aa) > 0:
                            mutations_str = "-muts " + " ".join(mutations_pdb_index_aa) + " "
                        traj_mutations_list.append(mutations_str + "\n")
                        cmd.load(mpnn_traj_file)
                    with open(os.path.join(mpnn_path, "mutations.txt"), "w") as pf:
                        pf.writelines(traj_mutations_list)
                    if not args.best_n_esm_scored_sequences:
                        with open(os.path.join(mpnn_path, "pymol.json"), "w") as pf:
                            pf.write(json.dumps(list(mutations_pdb_index_set)))
                        cmd.show(representation="sticks", selection="resi " + '+'.join(str(mutation_pdb_index) for mutation_pdb_index in mutations_pdb_index_set))
                        cmd.save(os.path.join(variant_path, variant_pdb_numbering + ".pse"), format='pse')
                        cmd.delete("*")
                        with open(os.path.join(mpnn_path, "mpnn_design_sites.json"), "w") as pf:
                            pf.write(json.dumps(list(mutations_chain_pdb_index_set)))
                elif args.step == 2 and os.path.isfile(os.path.join(mpnn_path, "mutations.txt")):
                    with open(os.path.join(mpnn_path, "mutations.txt"), "r") as pf:
                        traj_mutations_list = pf.readlines()
                    for ith_traj, line in enumerate(filter(lambda x: x.startswith("-muts "), traj_mutations_list)):
                        mutations_pdb_index = list(filter(lambda x: x != "", line.strip("\n")[6:].split(" ")))
                        mutations_str = str()
                        if len(mutations_pdb_index) > 0:
                            mutations_str = "-muts " + " ".join(mutations_pdb_index) + " "
                        relax_path = os.path.join(variant_path, args.output_directory + "_fast_relax_" + str(ith_traj + 1))
                        if os.path.isdir(relax_path):
                            shutil.rmtree(relax_path)
                        os.mkdir(relax_path)
                        with open(os.path.join(relax_path, "relax.sh"), "w") as pf:
                            pf.write('slurmit.py --job ' + variant_pdb_prefix + ' --mem 3000 --time 3-00:00:00 ')
                            pf.write('--command "python ../../../../../../enzdes_utils/fast_design.py ../' + variant_pdb + ' -cloud ')
                            pf.write('-ref ../../../../../pdb' + args.directory[3:] + '/' + pdb + '/' + pdb_chain + '.pdb ')
                            pf.write('-index_ref ../' + variant_pdb + ' -params ../../../../../' + ' ../../../../../'.join(params_files_list) + ' ')
                            pf.write('-sf ref2015_cst -enzdes_cst ../../../../../ligands/' + args.ligand + '/' + args.ligand + '_Relax.cst ')
                            pf.write('-cat_ids ' + uaa_name3 + ' ' + nucleophile_name3 + ' ')
                            pf.write('-sub_ids ' + xlink_name3 + ' ' + dihe_atoms + ' ' + dihe_params + ' ' + mutations_str)
                            pf.write('-premin -rpk_nbh -min_nbh -tform_enzdes -n ' + str(args.n_decoys) + ' -save_n 5 -prefix ' + variant + ' ')
                            pf.write('--score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1;"\n')
                        with open(os.path.join(relax_path, "relax_WT.sh"), "w") as pf:
                            pf.write('slurmit.py --job ' + variant_pdb_prefix + ' --mem 3000 --time 3-00:00:00 ')
                            pf.write('--command "python ../../../../../../enzdes_utils/fast_design.py ')
                            pf.write('../../../../../pdb' + args.directory[3:] + '/' + pdb + '/' + pdb_chain + '_relaxed.pdb ')
                            pf.write('-ref ../../../../../pdb' + args.directory[3:] + '/' + pdb + '/' + pdb_chain + '.pdb ')
                            pf.write('-index_ref ../' + variant_pdb + ' -params ../../../../../' + ' ../../../../../'.join(params_files_list) + ' ')
                            pf.write('-sf ref2015_cst -cat_ids ' + uaa_name3 + ' ' + nucleophile_name3 + ' ')
                            pf.write('-sub_ids ' + xlink_name3 + ' ' + mutations_str + '-ddG_wt -rpk_nbh -min_nbh ')
                            pf.write('-n 3 -save_n 1 -prefix ' + pdb_chain + '_apo ')
                            pf.write('--score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1;"\n')
