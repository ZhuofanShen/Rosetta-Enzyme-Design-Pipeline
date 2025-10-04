import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("variant_path", type=str)
parser.add_argument("-muts", "--mutations", type=str, nargs="*")
parser.add_argument("-sub", "--substrate", type=str, default="../../../../substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-n", "--n_decoy", type=int, default=30)
parser.add_argument("-script", type=str, default="relax")
parser.add_argument("-t", "--time", type=str, default="15")
parser.add_argument("-sf", "--score_function", type=str, default="ref2015_cst")
args = parser.parse_args()

variant_path = args.variant_path
variant = [s for s in variant_path.split("/")][-1]
pdb = variant.split("_")[0]
mutations = ""
if args.mutations and len(args.mutations) > 0:
    mutations = " -muts " + " ".join(args.mutations)
substrate_path = os.path.join("../..", args.substrate)
n_decoy = str(args.n_decoy)
script = args.script
time = args.time
additional_score_fn_terms = ''
if args.score_function.startswith("ref2015"):
    additional_score_fn_terms = ' --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1'

config_chi2_dict = {"1R2R+": "80", "1R2R-": "260", "1S2S+": "100", "1S2S-": "280", \
        "1R2S+": "74", "1R2S-": "256", "1S2R+": "104", "1S2R-": "286"}

additional_dihe_params = str()
if args.substrate.endswith("styrene_EDA"):
    additional_dihe_params = ' 0,42.97'
elif args.substrate.endswith("cyclopropanation_17_EDA"):
    additional_dihe_params = ' 180,28.00 180,21.77'
elif args.substrate.endswith("cyclopropanation_18_EDA"):
    additional_dihe_params = ' 180,30.00 180,21.77 0,30.00'

if not os.path.isdir(variant_path):
    os.mkdir(variant_path)
for stereo in ["1R2R", "1S2S", "1R2S", "1S2R"]:
    for rot in range(1, 5):
        for ester in ["+", "-"]:
            complex_sub_path = variant_path + "/" + pdb + "_" + stereo + "-rot" + str(rot) + ester
            if not os.path.isdir(complex_sub_path):
                os.mkdir(complex_sub_path)
                stereo_abbrv = stereo[1] + stereo[3]
                chi2 = " " + stereo_abbrv + "T,2," + config_chi2_dict[stereo + ester]
                chi6 = str()
                if substrate_path.endswith("cyclopropanation_p-CF3-styrene_EDA"):
                    chi6 = " " + stereo_abbrv + "T,6,90,270"
                elif substrate_path.endswith("cyclopropanation_17_EDA"):
                    if stereo_abbrv in ["RR", "SS"]:
                        chi6 = " " + stereo_abbrv + "T,7,57,303"
                    elif stereo_abbrv in ["RS", "SR"]:
                        chi6 = " " + stereo_abbrv + "T,5,0,180"
                elif substrate_path.endswith("cyclopropanation_23_EDA"):
                    if stereo_abbrv == "RR":
                        chi6 = " " + stereo_abbrv + "T,5,133,253"
                    elif stereo_abbrv == "SS":
                        chi6 = " " + stereo_abbrv + "T,5,107,227"
                    elif stereo_abbrv == "RS":
                        chi6 = " " + stereo_abbrv + "T,5,104,224"
                    elif stereo_abbrv == "SR":
                        chi6 = " " + stereo_abbrv + "T,5,136,256"
                additional_dihe_atoms = str()
                if args.substrate.endswith("styrene_EDA"):
                    additional_dihe_atoms = ' ' + stereo_abbrv + 'T,CST1 ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,CAR1 ' + stereo_abbrv + 'T,CAR2'
                elif args.substrate.endswith("cyclopropanation_17_EDA"):
                    additional_dihe_atoms = ' ' + stereo_abbrv + 'T,CST1 ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,CES ' + stereo_abbrv + 'T,OES ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,CES ' + stereo_abbrv + 'T,OES ' + stereo_abbrv + 'T,CAR1'
                elif args.substrate.endswith("cyclopropanation_18_EDA"):
                    additional_dihe_atoms = ' ' + stereo_abbrv + 'T,CST1 ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,OES ' + stereo_abbrv + 'T,CES ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,OES ' + stereo_abbrv + 'T,CES ' + stereo_abbrv + 'T,CAR1 ' + stereo_abbrv + 'T,OES ' + stereo_abbrv + 'T,CES ' + stereo_abbrv + 'T,CAR1 ' + stereo_abbrv + 'T,CAR2'
                cofactor_settings = ' -enzdes_cst ' + os.path.join(substrate_path, 'HEM-TS-' + stereo + '.cst') + ' -no_min_sc_ids HEM -no_min_je_ids HEM RRT SST SRT RST'
                if args.substrate.endswith("_DPP"):
                    cofactor_settings = ' ' + os.path.join(substrate_path, 'WUP.params') + ' -enzdes_cst ' + os.path.join(substrate_path, 'WUP-TS-' + stereo + '.cst') + ' -no_min_je_ids WUP RRT SST SRT RST'
                with open(complex_sub_path + "/" + script + ".sh", "w") as pf:
                    pf.write('slurmit.py --job ' + variant + '_' + stereo + '-rot' + str(rot) + ester + ' --partition main --time ' + time + ':00:00 --command "python ../../../../../../../enzdes_utils/fast_design.py ../../../complexes/' + pdb + '_' + stereo + '-rot' + str(rot) + '.pdb -ref ../../../../' + pdb + '_clean.pdb -index_ref ../../../complexes/' + pdb + '_1X2X-rotn.pdb -params ' + os.path.join(substrate_path, 'RRT.params') + ' ' + os.path.join(substrate_path, 'SST.params') + ' ' + os.path.join(substrate_path, 'RST.params') + ' ' + os.path.join(substrate_path, 'SRT.params') + cofactor_settings + ' -edges [-1,-2,FE,N1] -chis' + chi2 + chi6 + ' -dihe_atoms ' + stereo_abbrv + 'T,N1 ' + stereo_abbrv + 'T,Fe ' + stereo_abbrv + 'T,CARB ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OCAR ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OEST ' + stereo_abbrv + 'T,CET1' + additional_dihe_atoms + ' -dihe_params 0,34.95 0,21.77' + additional_dihe_params + ' -sub_ids RRT SST SRT RST' + mutations + ' -rpk_nbh -no_rpk_enzdes -min_nbh -n ' + n_decoy + ' -prefix ' + pdb + '_' + stereo + '-rot' + str(rot) + ester + ' -sf ' + args.score_function + additional_score_fn_terms + '"')
sub_path = variant_path + "/" + pdb
if not os.path.isdir(sub_path):
    os.mkdir(sub_path)
    cofactor_settings = ' -no_min_sc_ids HEM -no_min_je_ids HEM'
    if args.substrate.endswith("_DPP"):
        cofactor_settings = ' ' + os.path.join(substrate_path, 'WUP.params') + ' -enzdes_cst ' + os.path.join(substrate_path, 'WUP.cst') + ' -no_min_je_ids WUP'
    with open(sub_path + "/" + script + ".sh", "w") as pf:
        pf.write('slurmit.py --job ' + variant + ' --partition main --time ' + time + ':00:00 --command "python ../../../../../../../enzdes_utils/fast_design.py ../../../../' + pdb + '_relaxed.pdb -ref ../../../../' + pdb + '_clean.pdb -index_ref ../../../complexes/' + pdb + '_1X2X-rotn.pdb -params ' + os.path.join(substrate_path, 'RRT.params') + ' ' + os.path.join(substrate_path, 'SST.params') + ' ' + os.path.join(substrate_path, 'RST.params') + ' ' + os.path.join(substrate_path, 'SRT.params') + cofactor_settings + ' -sub_ids RRT SST SRT RST' + mutations + ' -rpk_nbh -no_rpk_enzdes -min_nbh -n ' + n_decoy + ' -prefix ' + pdb + ' -sf ' + args.score_function + additional_score_fn_terms + '"')
