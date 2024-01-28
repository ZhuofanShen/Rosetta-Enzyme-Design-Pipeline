import argparse
import json
import math
import os
import shutil
import xlwt


def read_scores_from_fasc(var):
    decoy_scores_str = None
    file_names = os.listdir(var)
    for file_name in file_names:
        if file_name.endswith(".fasc"):
            decoy_scores_str = "["
            with open(var + "/" + file_name) as fasc:
                for line in fasc:
                    decoy_scores_str += line[:-1]
                    decoy_scores_str += ","
            decoy_scores_str = decoy_scores_str[:-1] + "]"
    if decoy_scores_str:
        return json.loads(decoy_scores_str)
    else:
        return False

def extract_n_decoys(scores, n=1, is_reversed=False):
    scores = sorted(scores, key=lambda entry: entry["total_score"], reverse=is_reversed)
    lowest_scores = dict()
    for score in scores[0:n]:
        lowest_scores.update({str(score["decoy"]): score["total_score"]})
    return lowest_scores

def read_variant_scores(enz_var, row, sum_sheets, scores_sheets, ts_sheets, \
        green_style, orange_style, preferred_stereoisomer, \
        baseline_dGs=None, baseline_partition_function=None, temperature=4, \
        remove_redundant_decoys=False):
    pdb = baseline_variant.split("_")[0]
    for sum_sheet in sum_sheets:
        sum_sheet.write(row, 0, enz_var)
    for scores_sheet in scores_sheets:
        scores_sheet.write(row, 0, enz_var)
    for ts_sheet in ts_sheets:
        ts_sheet.write(row, 0, enz_var)

    enz_scores_apo = read_scores_from_fasc(enz_var + "/" + pdb)
    if enz_scores_apo:
        lowest_enz_scores = extract_n_decoys(enz_scores_apo)
        for enz_name, enz_dG_fold in lowest_enz_scores.items():
            enz_no = int(enz_name.split("_")[-1][:-4])
            break
        for best_decoy in os.listdir(enz_var + "/" + pdb):
            if best_decoy.endswith("_" + str(enz_no) + ".pdb"):
                enz_name = best_decoy
                break
        if remove_redundant_decoys:
            for decoy in os.listdir(enz_var + "/" + pdb):
                if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                        not decoy.endswith('.fasc'):
                    os.remove(enz_var + "/" + pdb + '/' + decoy)
    elif remove_redundant_decoys:
        shutil.rmtree(enz_var + "/" + pdb)
        enz_dG_fold = 0
    sum_sheets[0].write(row, 1, round(enz_dG_fold, 2))

    diastereo_negative_state_lowest_scores = list()
    diastereo_negative_state_partition_function = 0

    scores_sheet = scores_sheets[0]
    ts_sheet = ts_sheets[0]
    lowest_sc_1r2r = 10000
    partition_function = 0
    for col_1r2r in range(1, 9):
        if col_1r2r % 2 == 1:
            rot_name = str(int((col_1r2r + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1r2r / 2)) + '-'
        enz_scores_1r2r = read_scores_from_fasc(enz_var + "/" + pdb + "_1R2R-rot" + rot_name)
        if enz_scores_1r2r:
            lowest_enz_scores_1r2r = extract_n_decoys(enz_scores_1r2r)
            for enz_name_1r2r, enz_score_1r2r in lowest_enz_scores_1r2r.items():
                enz_no_1r2r = int(enz_name_1r2r.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_var + "/" + pdb + "_1R2R-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1r2r) + ".pdb"):
                    enz_name_1r2r = best_decoy
                    break
            scores_sheet.write(row, col_1r2r, round(enz_score_1r2r, 2))
            partition_function += math.exp(enz_dG_fold - enz_score_1r2r/temperature)
            if enz_score_1r2r < lowest_sc_1r2r:
                lowest_col_1r2r = col_1r2r
                lowest_sc_1r2r = enz_score_1r2r
                # lowest_no_1r2r = enz_no_1r2r
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_var + "/" + pdb + "_1R2R-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_var + "/" + pdb + "_1R2R-rot" + rot_name + '/' + decoy)
        elif remove_redundant_decoys:
            shutil.rmtree(enz_var + "/" + pdb + "_1R2R-rot" + rot_name)
        try:
            with open(enz_var + "/" + pdb + "_1R2R-rot" + rot_name + '/' + enz_name_1r2r, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('RRT'):
                        scores = line[:-1].split(' ')
                        ts_score_1r2r = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                ts_sheet.write(row, col_1r2r, ts_score_1r2r)
        except:
            pass
    scores_sheet.write(row, lowest_col_1r2r, round(lowest_sc_1r2r, 2), green_style)
    try:
        if preferred_stereoisomer == 0:
            positive_state_lowest_score = lowest_sc_1r2r
            positive_state_partition_function = partition_function
        elif preferred_stereoisomer == 1:
            enantio_negative_state_lowest_score = lowest_sc_1r2r
            enantio_negative_state_partition_function = partition_function
        else:
            diastereo_negative_state_lowest_scores.append(lowest_sc_1r2r)
            diastereo_negative_state_partition_function += partition_function
    except:
        pass
    
    scores_sheet = scores_sheets[1]
    ts_sheet = ts_sheets[1]
    lowest_sc_1s2s = 10000
    partition_function = 0
    for col_1s2s in range(1, 9):
        if col_1s2s % 2 == 1:
            rot_name = str(int((col_1s2s + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1s2s / 2)) + '-'
        enz_scores_1s2s = read_scores_from_fasc(enz_var + "/" + pdb + "_1S2S-rot" + rot_name)
        if enz_scores_1s2s:
            lowest_enz_scores_1s2s = extract_n_decoys(enz_scores_1s2s)
            for enz_name_1s2s, enz_score_1s2s in lowest_enz_scores_1s2s.items():
                enz_no_1s2s = int(enz_name_1s2s.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_var + "/" + pdb + "_1S2S-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1s2s) + ".pdb"):
                    enz_name_1s2s = best_decoy
                    break
            scores_sheet.write(row, col_1s2s, round(enz_score_1s2s, 2))
            partition_function += math.exp(enz_dG_fold - enz_score_1s2s/temperature)
            if enz_score_1s2s < lowest_sc_1s2s:
                lowest_col_1s2s = col_1s2s
                lowest_sc_1s2s = enz_score_1s2s
                # lowest_no_1s2s = enz_no_1s2s
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_var + "/" + pdb + "_1S2S-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_var + "/" + pdb + "_1S2S-rot" + rot_name + '/' + decoy)
        elif remove_redundant_decoys:
            shutil.rmtree(enz_var + "/" + pdb + "_1S2S-rot" + rot_name)
        try:
            with open(enz_var + "/" + pdb + "_1S2S-rot" + rot_name + '/' + enz_name_1s2s, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('SST'):
                        scores = line[:-1].split(' ')
                        ts_score_1s2s = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                ts_sheet.write(row, col_1s2s, ts_score_1s2s)
        except:
            pass
    scores_sheet.write(row, lowest_col_1s2s, round(lowest_sc_1s2s, 2), green_style)
    try:
        if preferred_stereoisomer == 1:
            positive_state_lowest_score = lowest_sc_1s2s
            positive_state_partition_function = partition_function
        elif preferred_stereoisomer == 0:
            enantio_negative_state_lowest_score = lowest_sc_1s2s
            enantio_negative_state_partition_function = partition_function
        else:
            diastereo_negative_state_lowest_scores.append(lowest_sc_1s2s)
            diastereo_negative_state_partition_function += partition_function
    except:
        pass
    
    scores_sheet = scores_sheets[2]
    ts_sheet = ts_sheets[2]
    lowest_sc_1r2s = 10000
    partition_function = 0
    for col_1r2s in range(1, 9):
        if col_1r2s % 2 == 1:
            rot_name = str(int((col_1r2s + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1r2s / 2)) + '-'
        enz_scores_1r2s = read_scores_from_fasc(enz_var + "/" + pdb + "_1R2S-rot" + rot_name)
        if enz_scores_1r2s:
            lowest_enz_scores_1r2s = extract_n_decoys(enz_scores_1r2s)
            for enz_name_1r2s, enz_score_1r2s in lowest_enz_scores_1r2s.items():
                enz_no_1r2s = int(enz_name_1r2s.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_var + "/" + pdb + "_1R2S-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1r2s) + ".pdb"):
                    enz_name_1r2s = best_decoy
                    break
            scores_sheet.write(row, col_1r2s, round(enz_score_1r2s, 2))
            partition_function += math.exp(enz_dG_fold - enz_score_1r2s/temperature)
            if enz_score_1r2s < lowest_sc_1r2s:
                lowest_col_1r2s = col_1r2s
                lowest_sc_1r2s = enz_score_1r2s
                # lowest_no_1r2s = enz_no_1r2s
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_var + "/" + pdb + "_1R2S-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_var + "/" + pdb + "_1R2S-rot" + rot_name + '/' + decoy)
        elif remove_redundant_decoys:
            shutil.rmtree(enz_var + "/" + pdb + "_1R2S-rot" + rot_name)
        try:
            with open(enz_var + "/" + pdb + "_1R2S-rot" + rot_name + '/' + enz_name_1r2s, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('RST'):
                        scores = line[:-1].split(' ')
                        ts_score_1r2s = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                ts_sheet.write(row, col_1r2s, ts_score_1r2s)
        except:
            pass
    scores_sheet.write(row, lowest_col_1r2s, round(lowest_sc_1r2s, 2), green_style)
    try:
        if preferred_stereoisomer == 2:
            positive_state_lowest_score = lowest_sc_1r2s
            positive_state_partition_function = partition_function
        elif preferred_stereoisomer == 3:
            enantio_negative_state_lowest_score = lowest_sc_1r2s
            enantio_negative_state_partition_function = partition_function
        else:
            diastereo_negative_state_lowest_scores.append(lowest_sc_1r2s)
            diastereo_negative_state_partition_function += partition_function
    except:
        pass
    
    scores_sheet = scores_sheets[3]
    ts_sheet = ts_sheets[3]
    lowest_sc_1s2r = 10000
    partition_function = 0
    for col_1s2r in range(1, 9):
        if col_1s2r % 2 == 1:
            rot_name = str(int((col_1s2r + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1s2r / 2)) + '-'
        enz_scores_1s2r = read_scores_from_fasc(enz_var + "/" + pdb + "_1S2R-rot" + rot_name)
        if enz_scores_1s2r:
            lowest_enz_scores_1s2r = extract_n_decoys(enz_scores_1s2r)
            for enz_name_1s2r, enz_score_1s2r in lowest_enz_scores_1s2r.items():
                enz_no_1s2r = int(enz_name_1s2r.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_var + "/" + pdb + "_1S2R-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1s2r) + ".pdb"):
                    enz_name_1s2r = best_decoy
                    break
            scores_sheet.write(row, col_1s2r, round(enz_score_1s2r, 2))
            partition_function += math.exp(enz_dG_fold - enz_score_1s2r/temperature)
            if enz_score_1s2r < lowest_sc_1s2r:
                lowest_col_1s2r = col_1s2r
                lowest_sc_1s2r = enz_score_1s2r
                # lowest_no_1s2r = enz_no_1s2r
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_var + "/" + pdb + "_1S2R-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_var + "/" + pdb + "_1S2R-rot" + rot_name + '/' + decoy)
        elif remove_redundant_decoys:
            shutil.rmtree(enz_var + "/" + pdb + "_1S2R-rot" + rot_name)
        try:
            with open(enz_var + "/" + pdb + "_1S2R-rot" + rot_name + '/' + enz_name_1s2r, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('SRT'):
                        scores = line[:-1].split(' ')
                        ts_score_1s2r = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                ts_sheet.write(row, col_1s2r, ts_score_1s2r)
        except:
            pass
    scores_sheet.write(row, lowest_col_1s2r, round(lowest_sc_1s2r, 2), green_style)
    try:
        if preferred_stereoisomer == 3:
            positive_state_lowest_score = lowest_sc_1s2r
            positive_state_partition_function = partition_function
        elif preferred_stereoisomer == 2:
            enantio_negative_state_lowest_score = lowest_sc_1s2r
            enantio_negative_state_partition_function = partition_function
        else:
            diastereo_negative_state_lowest_scores.append(lowest_sc_1s2r)
            diastereo_negative_state_partition_function += partition_function
    except:
        pass
    
    if lowest_sc_1r2r < lowest_sc_1s2s:
        sum_sheets[0].write(row, 2, round(lowest_sc_1r2r, 2), green_style)
        sum_sheets[0].write(row, 3, round(lowest_sc_1s2s, 2), orange_style)
    else:
        sum_sheets[0].write(row, 2, round(lowest_sc_1r2r, 2), orange_style)
        sum_sheets[0].write(row, 3, round(lowest_sc_1s2s, 2), green_style)
    if lowest_sc_1r2s < lowest_sc_1s2r:
        sum_sheets[0].write(row, 4, round(lowest_sc_1r2s, 2), green_style)
        sum_sheets[0].write(row, 5, round(lowest_sc_1s2r, 2), orange_style)
    else:
        sum_sheets[0].write(row, 4, round(lowest_sc_1r2s, 2), orange_style)
        sum_sheets[0].write(row, 5, round(lowest_sc_1s2r, 2), green_style)
    
    if baseline_dGs is not None:
        ddG_fold = enz_dG_fold - baseline_dGs[0]
        positive_state_ddG = positive_state_lowest_score - baseline_dGs[1]
        positive_state_ddG_bind = positive_state_ddG - ddG_fold
    else:
        baseline_dGs = [enz_dG_fold, positive_state_lowest_score]
        ddG_fold = 0
        positive_state_ddG = 0
        positive_state_ddG_bind = 0

    enantioselectivity = positive_state_lowest_score - enantio_negative_state_lowest_score
    diastereo_negative_state_lowest_score = min(diastereo_negative_state_lowest_scores)
    diastereoselectivity = positive_state_lowest_score - diastereo_negative_state_lowest_score
    stereoselectivity = max(enantioselectivity, diastereoselectivity)
    multistate_design_score = positive_state_ddG + stereoselectivity

    sum_sheets[1].write(row, 1, round(multistate_design_score, 2))
    sum_sheets[1].write(row, 2, round(positive_state_ddG, 2))
    sum_sheets[1].write(row, 3, round(ddG_fold, 2))
    sum_sheets[1].write(row, 4, round(positive_state_ddG_bind, 2))
    sum_sheets[1].write(row, 5, round(stereoselectivity, 2))
    sum_sheets[1].write(row, 6, round(diastereoselectivity, 2))
    sum_sheets[1].write(row, 7, round(enantioselectivity, 2))
    
    # if baseline_partition_function is not None:
    #     positive_state_ddG_bind = -math.log(positive_state_partition_function/baseline_partition_function) * temperature
    #     positive_state_ddG = ddG_fold + positive_state_ddG_bind
    # else:
    #     baseline_partition_function = positive_state_partition_function
    #     positive_state_ddG_bind = 0
    #     positive_state_ddG = 0
    
    # enantio_negative_design_fitness = positive_state_partition_function/enantio_negative_state_partition_function
    # enantioselectivity = -math.log(enantio_negative_design_fitness) * temperature
    # diastereo_negative_design_fitness = positive_state_partition_function/diastereo_negative_state_partition_function
    # diastereoselectivity = -math.log(diastereo_negative_design_fitness) * temperature
    # negative_state_partition_function = enantio_negative_state_partition_function + diastereo_negative_state_partition_function
    # negative_design_fitness = positive_state_partition_function/negative_state_partition_function
    # stereoselectivity = -math.log(negative_design_fitness) * temperature
    # minus_log_multistate_design_fitness = positive_state_ddG + stereoselectivity

    # sum_sheets[2].write(row, 1, round(minus_log_multistate_design_fitness, 2))
    # sum_sheets[2].write(row, 2, round(positive_state_ddG, 2))
    # sum_sheets[2].write(row, 3, round(ddG_fold, 2))
    # sum_sheets[2].write(row, 4, round(positive_state_ddG_bind, 2))
    # sum_sheets[2].write(row, 5, round(stereoselectivity, 2))
    # sum_sheets[2].write(row, 6, round(diastereoselectivity, 2))
    # sum_sheets[2].write(row, 7, round(enantioselectivity, 2))

    return baseline_dGs, baseline_partition_function


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline_variant", type=str)
    parser.add_argument("-s", "--preferred_stereoisomer", type=int, choices=[0, 1, 2, 3])
    parser.add_argument("-t", "--temperature", type=float, default=3)
    parser.add_argument("-rm", "--remove_redundant_decoys", action="store_true")
    args = parser.parse_args()
    baseline_variant = args.baseline_variant
    pdb = baseline_variant.split("_")[0]
    preferred_stereoisomer = args.preferred_stereoisomer
    temperature = args.temperature

    stereoisomers = ["1R2R", "1S2S", "1R2S", "1S2R"]
    workbook = xlwt.Workbook(encoding="ascii")
    sum_sheets = list()
    sum_sheets.append(workbook.add_sheet("summary"))
    sum_sheets.append(workbook.add_sheet("minimum"))
    # sum_sheets.append(workbook.add_sheet("Boltzmann factor"))
    scores_sheets = list()
    for i in range(4):
        scores_sheets.append(workbook.add_sheet("scores_" + stereoisomers[i], cell_overwrite_ok=True))
    ts_sheets = list()
    for i in range(4):
        ts_sheets.append(workbook.add_sheet("ts_" + stereoisomers[i], cell_overwrite_ok=True))
    green_style = xlwt.XFStyle()
    green_font = xlwt.Font()
    green_font.colour_index = 17
    green_style.font = green_font
    orange_style = xlwt.XFStyle()
    orange_font = xlwt.Font()
    orange_font.colour_index = 52
    orange_style.font = orange_font
    enz_vars = os.listdir()
    sum_sheets[0].write(0, 0, "unique sequence variants")
    sum_sheets[0].write(0, 1, "unbound")
    sum_sheets[0].write(0, 2, "1R,2R")
    sum_sheets[0].write(0, 3, "1S,2S")
    sum_sheets[0].write(0, 4, "1R,2S")
    sum_sheets[0].write(0, 5, "1S,2R")
    sum_sheets[1].write(0, 0, "unique sequence variants")
    sum_sheets[1].write(0, 1, "MSD score")
    sum_sheets[1].write(0, 2, "PSD score")
    sum_sheets[1].write(0, 3, "ddG fold")
    sum_sheets[1].write(0, 4, "ddG bind")
    sum_sheets[1].write(0, 5, "stereoselectivity")
    sum_sheets[1].write(0, 6, "diastereoselectivity")
    sum_sheets[1].write(0, 7, "enantioselectivity")
    # sum_sheets[2].write(0, 0, "unique sequence variants")
    # sum_sheets[2].write(0, 1, "MSD score")
    # sum_sheets[2].write(0, 2, "PSD score")
    # sum_sheets[2].write(0, 3, "ddG fold")
    # sum_sheets[2].write(0, 4, "ddG bind")
    # sum_sheets[2].write(0, 5, "stereoselectivity")
    # sum_sheets[2].write(0, 6, "diastereoselectivity")
    # sum_sheets[2].write(0, 7, "enantioselectivity")
    baseline_dGs, baseline_partition_function = read_variant_scores(baseline_variant, 1, \
            sum_sheets, scores_sheets, ts_sheets, green_style, orange_style, \
            preferred_stereoisomer, temperature=temperature, \
            remove_redundant_decoys=args.remove_redundant_decoys)
    row = 2
    for enz_var in filter(lambda x: x != baseline_variant and os.path.isdir(x) and \
            x.startswith(pdb + "_"), sorted(enz_vars)):
        _, _ = read_variant_scores(enz_var, row, sum_sheets, scores_sheets, ts_sheets, \
        green_style, orange_style, preferred_stereoisomer, baseline_dGs=baseline_dGs, \
        baseline_partition_function=baseline_partition_function, temperature=temperature,
        remove_redundant_decoys=args.remove_redundant_decoys)
        row = row + 1
    workbook.save(pdb + "_scores.xls")
