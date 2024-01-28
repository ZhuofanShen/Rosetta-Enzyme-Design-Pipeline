import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("variant_path", type=str)
parser.add_argument("-muts", "--mutations", type=str, nargs="*")
parser.add_argument("-sub", "--substrate", type=str, default="../../../../substrates/cyclopropanation_styrene_EDA")
parser.add_argument("-n", "--n_decoy", type=int, default=30)
parser.add_argument("-script", type=str, default="relax")
parser.add_argument("-t", "--time", type=str, default="15")
args = parser.parse_args()

variant_path = args.variant_path
variant = [s for s in variant_path.split("/")][-1]
pdb = variant.split("_")[0]
mutations = " ".join(args.mutations)
substrate_path = os.path.join("../..", args.substrate)
n_decoy = str(args.n_decoy)
script = args.script
time = args.time

config_chi2_dict = {"1R2R+": "80", "1R2R-": "260", "1S2S+": "100", "1S2S-": "280", \
        "1R2S+": "74", "1R2S-": "256", "1S2R+": "104", "1S2R-": "286"}

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
                if substrate_path.endswith("cyclopropanation_p-CF3-styrene_EDA"):
                    chi6 = " " + stereo_abbrv + "T,6,90,270"
                else:
                    chi6 = str()
                with open(complex_sub_path + "/" + script + ".sh", "w") as pf:
                    pf.write('slurmit.py --job ' + variant + '_' + stereo + '-rot' + str(rot) + ester + ' --partition main --time ' + time + ':00:00 --command "python ../../../../../../../enzdes_utils/fast_design.py ../../../complexes/' + pdb + '_' + stereo + '-rot' + str(rot) + '.pdb -ref ../../../../' + pdb + '_clean.pdb -index_ref ../../../complexes/' + pdb + '_1X2X-rotn.pdb -params ' + os.path.join(substrate_path, 'RRT.params') + ' ' + os.path.join(substrate_path, 'SST.params') + ' ' + os.path.join(substrate_path, 'RST.params') + ' ' + os.path.join(substrate_path, 'SRT.params') + ' -edges [-1,-2,FE,N1] -chis' + chi2 + chi6 + ' -dihe_atoms ' + stereo_abbrv + 'T,N1 ' + stereo_abbrv + 'T,Fe ' + stereo_abbrv + 'T,CARB ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OCAR ' + stereo_abbrv + 'T,CEST ' + stereo_abbrv + 'T,OEST ' + stereo_abbrv + 'T,CET1 ' + stereo_abbrv + 'T,CST1 ' + stereo_abbrv + 'T,CST2 ' + stereo_abbrv + 'T,CAR1 ' + stereo_abbrv + 'T,CAR2 -dihe_params 0,34.95 0,21.77 0,42.97 -enzdes_cst ' + os.path.join(substrate_path, 'HEM-TS-' + stereo + '.cst') + ' -no_pack_ids HEM -no_min_sc_ids HEM -no_min_je_ids HEM RRT SST SRT RST -sub_ids RRT SST SRT RST -muts ' + mutations + ' -rpk_nbh -no_rpk_enzdes -min_nbh -n ' + n_decoy + ' -prefix ' + pdb + '_' + stereo + '-rot' + str(rot) + ester + ' --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1"')
sub_path = variant_path + "/" + pdb
if not os.path.isdir(sub_path):
    os.mkdir(sub_path)
    with open(sub_path + "/" + script + ".sh", "w") as pf:
        pf.write('slurmit.py --job ' + variant + ' --partition main --time ' + time + ':00:00 --command "python ../../../../../../../enzdes_utils/fast_design.py ../../../../' + pdb + '_relaxed.pdb -ref ../../../../' + pdb + '_clean.pdb -index_ref ../../../complexes/' + pdb + '_1X2X-rotn.pdb -params ' + os.path.join(substrate_path, 'RRT.params') + ' ' + os.path.join(substrate_path, 'SST.params') + ' ' + os.path.join(substrate_path, 'RST.params') + ' ' + os.path.join(substrate_path, 'SRT.params') + ' -no_pack_ids HEM -no_min_sc_ids HEM -no_min_je_ids HEM -sub_ids RRT SST SRT RST -muts ' + mutations + ' -rpk_nbh -no_rpk_enzdes -min_nbh -n ' + n_decoy + ' -prefix ' + pdb + ' --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1"')
