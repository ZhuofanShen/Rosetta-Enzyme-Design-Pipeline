import argparse
import os
from pymol import *
import traceback

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=str)
parser.add_argument("-f", "--fold", type=str)
parser.add_argument("-sub", "--substrate_path", type=str, default="substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-iso", "--isomers", type=str, nargs="*", default=["RRT", "SST", "RST", "SRT"])
args = parser.parse_args()

dir = args.directory
fold = args.fold
substrate_path = args.substrate_path
substrate = os.path.basename(substrate_path)

distal_side_nitrogen_order = ["N1", "N2", "N3", "N4"]
proximal_side_nitrogen_order = ["N1", "N4", "N3", "N2"]
for pdb in filter(lambda x: os.path.isfile(os.path.join(dir, x, x + "_relaxed.pdb")) \
        and (not fold or os.path.isfile(os.path.join(dir, x, fold + ".fold"))) \
        , os.listdir(dir)): # and os.path.isfile(os.path.join(dir, x, "main"))
    try:
        link_lines = list()
        ssbond_lines = list()
        remark666_lines = list()
        with open(os.path.join(dir, pdb, pdb + "_relaxed.pdb"), "r") as ppdb:
            for line in ppdb:
                if line.startswith("REMARK 666 MATCH "):
                    remark666_lines.append(line)
                elif line.startswith("LINK"):
                    link_lines.append(line)
                elif line.startswith("SSBOND"):
                    ssbond_lines.append(line)
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    break
        with open(os.path.join(dir, pdb, pdb + ".out"), "r") as pf:
            lines = pf.readlines()
        heme = lines[2].strip("\n").split(" ")
        if heme[0].startswith("HE"):
            heme_nitrogen_atoms = ["NA", "NB", "NC", "ND"]
        elif heme[0] == "WUP":
            heme_nitrogen_atoms = ["N1", "N2", "N3", "N4"]
        else:
            raise Exception("Unknown porphyrin type.")
        if len(lines) == 5:
            line_index = [3, 4]
        elif len(lines) == 4:
            line_index = [3]
        elif len(lines) == 3:
            line_index = []
            with open("no_axial_residue.log", "a") as pf:
                pf.write(pdb + "\n")
        for i in line_index:
            if lines[i][:3] not in ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", \
                        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]:
                continue
            proximal_res = lines[i].strip("\n").split(" ")
            if proximal_res[4] == "PROXIMAL":
                open_coordination_site = "distal"
                nitrogen_order = distal_side_nitrogen_order
            elif proximal_res[4] == "DISTAL":
                open_coordination_site = "proximal"
                nitrogen_order = proximal_side_nitrogen_order
            pdb_sub_path = os.path.join(dir, pdb, substrate + "_" + open_coordination_site)
            if not os.path.isdir(pdb_sub_path):
                os.mkdir(pdb_sub_path)
            if not os.path.isdir(os.path.join(pdb_sub_path, "complexes")):
                os.mkdir(os.path.join(pdb_sub_path, "complexes"))
            another_proximal_res = None
            if len(lines) > 5:
                if i == 3:
                    j = 4
                else:
                    j = 3
                another_proximal_res = lines[j].strip("\n").split(" ")
                if not another_proximal_res[0] in ["ALA", "CYS", "ASP", "GLU", "PHE", \
                        "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", \
                        "ARG", "SER", "THR", "VAL", "TRP", "TYR"]:
                    another_proximal_res = None
                else:
                    with open(os.path.join(pdb_sub_path, another_proximal_res[0]), "w") as _:
                        pass
            cmd.load(os.path.join(dir, pdb, pdb + "_relaxed.pdb"), pdb)
            for stereo in args.isomers:
                cmd.load(os.path.join(substrate_path, stereo + ".pdb"))
            ddG_ref_pdb_string = pdb
            pdb_index = 0
            for stereo in args.isomers:
                for rot in range(4):
                    pdb_index += 1
                    obj = stereo + str(pdb_index)
                    cmd.create(obj, stereo)
                    cmd.alter(obj, "resi='" + str(pdb_index) + "'")
                    cmd.pair_fit("/" + obj + "//X/" + str(pdb_index) + "/" + nitrogen_order[0-rot], \
                                "/" + pdb + "//" + heme[1] + "/" + heme[2] + "/" + heme_nitrogen_atoms[0], \
                            "/" + obj + "//X/" + str(pdb_index) + "/" + nitrogen_order[1-rot], \
                            "/" + pdb + "//" + heme[1] + "/" + heme[2] + "/" + heme_nitrogen_atoms[1], \
                            "/" + obj + "//X/" + str(pdb_index) + "/" + nitrogen_order[2-rot], \
                            "/" + pdb + "//" + heme[1] + "/" + heme[2] + "/" + heme_nitrogen_atoms[2], \
                            "/" + obj + "//X/" + str(pdb_index) + "/" + nitrogen_order[3-rot], \
                            "/" + pdb + "//" + heme[1] + "/" + heme[2] + "/" + heme_nitrogen_atoms[3])
                    ddG_ref_pdb_string += " | " + obj
                    pdb_index += 1
                    obj2 = stereo + str(pdb_index)
                    cmd.create(obj2, obj)
                    cmd.alter(obj2, "resi='" + str(pdb_index) + "'")
                    dihedral = cmd.get_dihedral("/" + obj2 + "//X/" + str(pdb_index) + "/Fe", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/CARB", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/CEST", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/OCAR") + 180.
                    cmd.set_dihedral("/" + obj2 + "//X/" + str(pdb_index) + "/Fe", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/CARB", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/CEST", \
                            "/" + obj2 + "//X/" + str(pdb_index) + "/OCAR", dihedral)
                    ddG_ref_pdb_string += " | " + obj2
            ddG_ref_pdb_name = os.path.join(pdb_sub_path, "complexes", pdb + "_1X2X-rotn.pdb")
            cmd.save(ddG_ref_pdb_name, ddG_ref_pdb_string)
            if len(link_lines) > 0 or len(ssbond_lines) > 0:
                time.sleep(5)
                with open(ddG_ref_pdb_name, "r+") as ppdb:
                    pdb_lines = ppdb.readlines()
                    ppdb.seek(0)
                    if len(ssbond_lines) > 0:
                        ppdb.writelines(ssbond_lines)
                    if len(link_lines) > 0:
                        ppdb.writelines(link_lines)
                    ppdb.writelines(pdb_lines)
                    ppdb.truncate()
            pdb_index = -1
            if len(heme[2]) == 3:
                heme_resi_str = heme[2]
            elif len(heme[2]) == 2:
                heme_resi_str = " " + heme[2]
            elif len(heme[2]) == 1:
                heme_resi_str = "  " + heme[2]
            for stereo in args.isomers:
                for rot in range(4):
                    pdb_index += 2
                    obj = stereo + str(pdb_index)
                    cmd.alter(obj, "resi='1'")
                    pdb_name = os.path.join(pdb_sub_path, "complexes", pdb + \
                            "_1" + stereo[0] + "2" + stereo[1] + "-rot" + str(rot + 1) + ".pdb")
                    cmd.save(pdb_name, pdb + " | " + obj)
                    if len(remark666_lines) > 0 or len(link_lines) > 0 or len(ssbond_lines) > 0:
                        time.sleep(5)
                        with open(pdb_name, "r+") as ppdb:
                            pdb_lines = ppdb.readlines()
                            ppdb.seek(0)
                            if len(remark666_lines) > 0:
                                ppdb.writelines(remark666_lines)
                            ppdb.write("REMARK 666 MATCH TEMPLATE X " + stereo + "    1 MATCH MOTIF " + heme[1] + " " + heme[0] + "  " + heme_resi_str + "  " + str(len(remark666_lines) + 1) + "  1               \n")
                            if len(ssbond_lines) > 0:
                                ppdb.writelines(ssbond_lines)
                            if len(link_lines) > 0:
                                ppdb.writelines(link_lines)
                            ppdb.writelines(pdb_lines)
                            ppdb.truncate()
            cmd.delete("*")
    except Exception:
        with open(substrate + ".err", "a") as pf:
            pf.write(pdb + "\n")
            pf.write(traceback.format_exc())
