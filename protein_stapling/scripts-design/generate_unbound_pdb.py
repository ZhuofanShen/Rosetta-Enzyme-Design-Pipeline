import os


suffix = "antibody"
for pdb in os.listdir("out" + suffix):
    pdb_path = os.path.join("out" + suffix, pdb)
    for pdb_chain_ligand in os.listdir(pdb_path):
        pdb_chain = pdb_chain_ligand.split("_")[0]
        with open(os.path.join("pdb" + suffix, pdb, pdb_chain + "-unbound_relaxed.pdb"), "r") as pf:
            unbound_wt_relaxed_pdb_lines = pf.readlines()
        pdb_chain_ligand_path = os.path.join(pdb_path, pdb_chain_ligand)
        if not os.path.isdir(pdb_chain_ligand_path):
            continue
        for variant in os.listdir(pdb_chain_ligand_path):
            variant_path = os.path.join(pdb_chain_ligand_path, variant)
            if not os.path.isfile(os.path.join(variant_path, variant + ".sc")):
                continue
            for f in os.listdir(variant_path):
                if f.startswith("Z") and f.endswith(".pdb"):
                    uaa_lines = list()
                    nucleophile_lines = list()
                    f_path = os.path.join(variant_path, f)
                    with open(f_path, "r") as pf:
                        for line in pf:
                            if line[17:20] == "TYZ":
                                uaa_chain_res_idx = line[21:26]
                                uaa_lines.append(line)
                            elif line[17:20] in ["CYX", "HIX", "LYX", "ASX", "GLX"]:
                                nucleophile_name3 = line[17:20]
                                nucleophile_chain_res_idx = line[21:26]
                                nucleophile_lines.append(line)

                    with open(f_path[:-4] + "-unbound.pdb", "w") as pf:
                        write_uaa = False
                        write_nucleophile = False
                        for line in filter(lambda x: x.startswith("ATOM  "), unbound_wt_relaxed_pdb_lines):
                            if line[21:26] == uaa_chain_res_idx:
                                if not write_uaa:
                                    for uaa_line in uaa_lines:
                                        if uaa_line[21] == "A":
                                            new_x = str(round(float(uaa_line[31:38]) + 122, 3))
                                            if len(new_x) == 5:
                                                new_x += "0"
                                            elif len(new_x) == 4:
                                                new_x += "00"
                                            elif len(new_x) == 3:
                                                new_x += "000"
                                            elif len(new_x) == 2:
                                                new_x += ".000"
                                            pf.write(uaa_line[:31] + " " + new_x + uaa_line[38:])
                                        else:
                                            pf.write(uaa_line)
                                    write_uaa = True
                            elif line[21:26] == nucleophile_chain_res_idx:
                                if not write_nucleophile:
                                    for nucleophile_line in nucleophile_lines:
                                        if nucleophile_line[21] == "A":
                                            new_x = str(round(float(nucleophile_line[31:38]) + 122, 3))
                                            if len(new_x) == 5:
                                                new_x += "0"
                                            elif len(new_x) == 4:
                                                new_x += "00"
                                            elif len(new_x) == 3:
                                                new_x += "000"
                                            elif len(new_x) == 2:
                                                new_x += ".000"
                                            pf.write(nucleophile_line[:31] + " " + new_x + nucleophile_line[38:])
                                        else:
                                            pf.write(nucleophile_line)
                                    write_nucleophile = True
                            else:
                                pf.write(line)
