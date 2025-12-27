import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("directory", type=str)
parser.add_argument("-f", "--fold", type=str)
parser.add_argument("-sub", "--substrate", type=str, default="substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-stereo", "--stereoisomers", type=str, nargs="*", default=["1R2R", "1S2S", "1R2S", "1S2R"])
parser.add_argument("-n", "--n_decoy", type=int, default=5)
parser.add_argument('-t', '--time', type=str, default='3-00')
args = parser.parse_args()

dir = args.directory
fold = args.fold
substrate = args.substrate.strip("/").split("/")[-1]
substrate_path = os.path.join("../../../../../", args.substrate)
n_decoy = str(args.n_decoy)
time = args.time

config_chi2_dict = {"1R2R+": "80", "1R2R-": "260", "1S2S+": "100", "1S2S-": "280", \
        "1R2S+": "74", "1R2S-": "256", "1S2R+": "104", "1S2R-": "286"}

styrenyl_dihe_params = str()
if not args.substrate.endswith("cyclopropanation_benzofuran_EDA"):
    styrenyl_dihe_params = ' 0,42.97'

variants = 0
jobs = 0
if fold:
    bash_name = fold.replace(" ", "_")
else:
    bash_name = "default"
for pdb in filter(lambda x: os.path.isfile(os.path.join(dir, x, "main")) \
    and (not fold or os.path.isfile(os.path.join(dir, x, fold + ".fold"))), os.listdir(dir)):
    variants += 1
    for open_coordination_site in ["distal", "proximal"]:
        pdb_sub_path = os.path.join(dir, pdb, substrate + "_" + open_coordination_site)
        if not os.path.isdir(pdb_sub_path):
            continue
        design_root_path = os.path.join(pdb_sub_path, "FastDesign")
        jobs += 1
        if not os.path.isdir(design_root_path):
            os.mkdir(design_root_path)
        for stereo in args.stereoisomers:
            for rot in range(1, 5):
                for ester in ["+", "-"]:
                    design_path = os.path.join(design_root_path, pdb + "_" + stereo + "-rot" + str(rot) + ester)
                    if os.path.isdir(design_path):
                        continue
                    os.mkdir(design_path)
                    stereo_abbrv = stereo[1] + stereo[3]
                    chi2 = config_chi2_dict[stereo + ester]
                    styrenyl_dihe_atoms = str()
                    if not args.substrate.endswith("cyclopropanation_benzofuran_EDA"):
                        styrenyl_dihe_atoms = ' ' + stereo_abbrv + 'T,CST1 ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,CAR1 ' + stereo_abbrv + 'T,CAR2'
                    with open(design_path + "/" + bash_name + ".sh", "w") as pf:
                        pf.write('slurmit.py --job ' + pdb + '_' + stereo + '-rot' + str(rot) + ester + ' --partition main --time ' + time + ':00:00 --command "python ../../../../../../enzdes_utils/fast_design.py ../../complexes/' + pdb + '_' + stereo + '-rot' + str(rot) + '.pdb -ref ../../../' + pdb + '_clean.pdb -index_ref ../../complexes/' + pdb + '_1X2X-rotn.pdb -params ' + os.path.join(substrate_path, 'RRT.params') + ' ' + os.path.join(substrate_path, 'SST.params') + ' ' + os.path.join(substrate_path, 'RST.params') + ' ' + os.path.join(substrate_path, 'SRT.params') + ' -edges [-1,-2,FE,N1] -chis ' + stereo_abbrv + 'T,2,' + chi2 + ' -dihe_atoms ' + stereo_abbrv + 'T,N1 ' + stereo_abbrv + 'T,Fe ' + stereo_abbrv + 'T,CARB ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OCAR ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OEST ' + stereo_abbrv + 'T,CET1' + styrenyl_dihe_atoms + ' -dihe_params 0,34.95 0,21.77' + styrenyl_dihe_params + ' -enzdes_cst ' + os.path.join(substrate_path, 'HEM-TS-' + stereo + '.cst') + ' -no_min_sc_ids HEM -no_min_je_ids HEM RRT SST SRT RST -sub_ids RRT SST SRT RST -des_bs -nataa 2' + ' -rpk_nbh -no_rpk_enzdes -min_nbh -n ' + n_decoy + ' -prefix ' + pdb + '_' + stereo + '-rot' + str(rot) + ester + ' --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1"')
print("total variants: " + str(variants))
print("total jobs: " + str(jobs) + " * 8")
