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

def read_variant_scores(enz_path, preferred_stereoisomer, baseline_dGs=None, \
        baseline_partition_function=None, temperature=4, remove_redundant_decoys=False):
    enz_var = enz_path.split("/")[-1]
    xls_row_info = list()
    pdb = enz_var.split("_")[0]
    enz_scores_apo = read_scores_from_fasc(enz_path + "/" + pdb)
    if enz_scores_apo:
        lowest_enz_scores = extract_n_decoys(enz_scores_apo)
        for enz_name, enz_dG_fold in lowest_enz_scores.items():
            enz_no = int(enz_name.split("_")[-1][:-4])
            break
        for best_decoy in os.listdir(enz_path + "/" + pdb):
            if best_decoy.endswith("_" + str(enz_no) + ".pdb"):
                enz_name = best_decoy
                break
        if remove_redundant_decoys:
            for decoy in os.listdir(enz_path + "/" + pdb):
                if not decoy.endswith("_" + str(enz_no) + ".pdb") and \
                        not decoy.endswith('.fasc'):
                    os.remove(enz_path + "/" + pdb + '/' + decoy)
    elif remove_redundant_decoys:
        shutil.rmtree(enz_path + "/" + pdb)
        enz_dG_fold = 0
    xls_row_info.append(round(enz_dG_fold, 2)) # [0]

    diastereo_negative_state_lowest_scores = list()
    diastereo_negative_state_partition_function = 0

    lowest_sc_1r2r = 10000
    partition_function = 0
    for col_1r2r in range(1, 9):
        if col_1r2r % 2 == 1:
            rot_name = str(int((col_1r2r + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1r2r / 2)) + '-'
        enz_scores_1r2r = read_scores_from_fasc(enz_path + "/" + pdb + "_1R2R-rot" + rot_name)
        if enz_scores_1r2r:
            lowest_enz_scores_1r2r = extract_n_decoys(enz_scores_1r2r)
            for enz_name_1r2r, enz_score_1r2r in lowest_enz_scores_1r2r.items():
                enz_no_1r2r = int(enz_name_1r2r.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_path + "/" + pdb + "_1R2R-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1r2r) + ".pdb"):
                    enz_name_1r2r = best_decoy
                    break
            xls_row_info.append(round(enz_score_1r2r, 2)) # [1]
            partition_function += math.exp(enz_dG_fold - enz_score_1r2r/temperature)
            if enz_score_1r2r < lowest_sc_1r2r:
                lowest_col_1r2r = col_1r2r
                lowest_sc_1r2r = enz_score_1r2r
                # lowest_no_1r2r = enz_no_1r2r
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_path + "/" + pdb + "_1R2R-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no_1r2r) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_path + "/" + pdb + "_1R2R-rot" + rot_name + '/' + decoy)
        else:
            xls_row_info.append(None) # [1]
            if remove_redundant_decoys:
                shutil.rmtree(enz_path + "/" + pdb + "_1R2R-rot" + rot_name)
        try:
            with open(enz_path + "/" + pdb + "_1R2R-rot" + rot_name + '/' + enz_name_1r2r, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('RRT'):
                        scores = line[:-1].split(' ')
                        ts_score_1r2r = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                xls_row_info.append(ts_score_1r2r) # [2]
        except:
            xls_row_info.append(None) # [2]
    xls_row_info.append(lowest_col_1r2r) # [17]
    xls_row_info.append(round(lowest_sc_1r2r, 2)) # [18]
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
    
    lowest_sc_1s2s = 10000
    partition_function = 0
    for col_1s2s in range(1, 9):
        if col_1s2s % 2 == 1:
            rot_name = str(int((col_1s2s + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1s2s / 2)) + '-'
        enz_scores_1s2s = read_scores_from_fasc(enz_path + "/" + pdb + "_1S2S-rot" + rot_name)
        if enz_scores_1s2s:
            lowest_enz_scores_1s2s = extract_n_decoys(enz_scores_1s2s)
            for enz_name_1s2s, enz_score_1s2s in lowest_enz_scores_1s2s.items():
                enz_no_1s2s = int(enz_name_1s2s.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_path + "/" + pdb + "_1S2S-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1s2s) + ".pdb"):
                    enz_name_1s2s = best_decoy
                    break
            xls_row_info.append(round(enz_score_1s2s, 2)) # [19]
            partition_function += math.exp(enz_dG_fold - enz_score_1s2s/temperature)
            if enz_score_1s2s < lowest_sc_1s2s:
                lowest_col_1s2s = col_1s2s
                lowest_sc_1s2s = enz_score_1s2s
                # lowest_no_1s2s = enz_no_1s2s
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_path + "/" + pdb + "_1S2S-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no_1s2s) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_path + "/" + pdb + "_1S2S-rot" + rot_name + '/' + decoy)
        else:
            xls_row_info.append(None) # [19]
            if remove_redundant_decoys:
                shutil.rmtree(enz_path + "/" + pdb + "_1S2S-rot" + rot_name)
        try:
            with open(enz_path + "/" + pdb + "_1S2S-rot" + rot_name + '/' + enz_name_1s2s, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('SST'):
                        scores = line[:-1].split(' ')
                        ts_score_1s2s = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                xls_row_info.append(ts_score_1s2s) # [20]
        except:
            xls_row_info.append(None) # [20]
    xls_row_info.append(lowest_col_1s2s) # [35]
    xls_row_info.append(round(lowest_sc_1s2s, 2)) # [36]
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
    
    lowest_sc_1r2s = 10000
    partition_function = 0
    for col_1r2s in range(1, 9):
        if col_1r2s % 2 == 1:
            rot_name = str(int((col_1r2s + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1r2s / 2)) + '-'
        enz_scores_1r2s = read_scores_from_fasc(enz_path + "/" + pdb + "_1R2S-rot" + rot_name)
        if enz_scores_1r2s:
            lowest_enz_scores_1r2s = extract_n_decoys(enz_scores_1r2s)
            for enz_name_1r2s, enz_score_1r2s in lowest_enz_scores_1r2s.items():
                enz_no_1r2s = int(enz_name_1r2s.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_path + "/" + pdb + "_1R2S-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1r2s) + ".pdb"):
                    enz_name_1r2s = best_decoy
                    break
            xls_row_info.append(round(enz_score_1r2s, 2)) # [37]
            partition_function += math.exp(enz_dG_fold - enz_score_1r2s/temperature)
            if enz_score_1r2s < lowest_sc_1r2s:
                lowest_col_1r2s = col_1r2s
                lowest_sc_1r2s = enz_score_1r2s
                # lowest_no_1r2s = enz_no_1r2s
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_path + "/" + pdb + "_1R2S-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no_1r2s) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_path + "/" + pdb + "_1R2S-rot" + rot_name + '/' + decoy)
        else:
            xls_row_info.append(None) # [37]
            if remove_redundant_decoys:
                shutil.rmtree(enz_path + "/" + pdb + "_1R2S-rot" + rot_name)
        try:
            with open(enz_path + "/" + pdb + "_1R2S-rot" + rot_name + '/' + enz_name_1r2s, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('RST'):
                        scores = line[:-1].split(' ')
                        ts_score_1r2s = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                xls_row_info.append(ts_score_1r2s) # [38]
        except:
            xls_row_info.append(None) # [38]
    xls_row_info.append(lowest_col_1r2s) # [53]
    xls_row_info.append(round(lowest_sc_1r2s, 2)) # [54]
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
    
    lowest_sc_1s2r = 10000
    partition_function = 0
    for col_1s2r in range(1, 9):
        if col_1s2r % 2 == 1:
            rot_name = str(int((col_1s2r + 1) / 2)) + '+'
        else:
            rot_name = str(int(col_1s2r / 2)) + '-'
        enz_scores_1s2r = read_scores_from_fasc(enz_path + "/" + pdb + "_1S2R-rot" + rot_name)
        if enz_scores_1s2r:
            lowest_enz_scores_1s2r = extract_n_decoys(enz_scores_1s2r)
            for enz_name_1s2r, enz_score_1s2r in lowest_enz_scores_1s2r.items():
                enz_no_1s2r = int(enz_name_1s2r.split("_")[-1][:-4])
                break
            for best_decoy in os.listdir(enz_path + "/" + pdb + "_1S2R-rot" + rot_name):
                if best_decoy.endswith("_" + str(enz_no_1s2r) + ".pdb"):
                    enz_name_1s2r = best_decoy
                    break
            xls_row_info.append(round(enz_score_1s2r, 2)) # [55]
            partition_function += math.exp(enz_dG_fold - enz_score_1s2r/temperature)
            if enz_score_1s2r < lowest_sc_1s2r:
                lowest_col_1s2r = col_1s2r
                lowest_sc_1s2r = enz_score_1s2r
                # lowest_no_1s2r = enz_no_1s2r
            if remove_redundant_decoys:
                for decoy in os.listdir(enz_path + "/" + pdb + "_1S2R-rot" + rot_name):
                    if not decoy.endswith("_" + str(enz_no_1s2r) + ".pdb") and \
                            not decoy.endswith('.fasc'):
                        os.remove(enz_path + "/" + pdb + "_1S2R-rot" + rot_name + '/' + decoy)
        else:
            xls_row_info.append(None) # [55]
            if remove_redundant_decoys:
                shutil.rmtree(enz_path + "/" + pdb + "_1S2R-rot" + rot_name)
        try:
            with open(enz_path + "/" + pdb + "_1S2R-rot" + rot_name + '/' + enz_name_1s2r, 'r') as pf:
                for line in pf.readlines():
                    if line.startswith('SRT'):
                        scores = line[:-1].split(' ')
                        ts_score_1s2r = round(float(scores[-1]) - float(scores[-2]) - float(scores[-13]), 2)
                xls_row_info.append(ts_score_1s2r) # [56]
        except:
            xls_row_info.append(None) # [56]
    xls_row_info.append(lowest_col_1s2r) # [71]
    xls_row_info.append(round(lowest_sc_1s2r, 2)) # [72]
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

    xls_row_info.append(round(multistate_design_score, 2)) # [73]
    xls_row_info.append(round(positive_state_ddG, 2)) # [74]
    xls_row_info.append(round(ddG_fold, 2)) # [75]
    xls_row_info.append(round(positive_state_ddG_bind, 2)) # [76]
    xls_row_info.append(round(stereoselectivity, 2)) # [77]
    xls_row_info.append(round(diastereoselectivity, 2)) # [78]
    xls_row_info.append(round(enantioselectivity, 2)) # [79]
    
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

    # xls_row_info.append(round(minus_log_multistate_design_fitness, 2)) # [80]
    # xls_row_info.append(round(positive_state_ddG, 2)) # [81]
    # xls_row_info.append(round(ddG_fold, 2)) # [82]
    # xls_row_info.append(round(positive_state_ddG_bind, 2)) # [83]
    # xls_row_info.append(round(stereoselectivity, 2)) # [84]
    # xls_row_info.append(round(diastereoselectivity, 2)) # [85]
    # xls_row_info.append(round(enantioselectivity, 2)) # [86]

    return xls_row_info, baseline_dGs, baseline_partition_function

def write_xls_row(row, enz_var, xls_row_info, sum_sheets, scores_sheets, ts_sheets, \
        green_style, orange_style):
    for sum_sheet in sum_sheets:
        sum_sheet.write(row, 0, enz_var)
    for scores_sheet in scores_sheets:
        scores_sheet.write(row, 0, enz_var)
    for ts_sheet in ts_sheets:
        ts_sheet.write(row, 0, enz_var)
    sum_sheets[0].write(row, 1, xls_row_info[0])
    for col_1r2r in range(1, 9):
        if xls_row_info[col_1r2r * 2 - 1] is not None:
            scores_sheets[0].write(row, col_1r2r, xls_row_info[col_1r2r * 2 - 1])
        if xls_row_info[col_1r2r * 2] is not None:
            ts_sheets[0].write(row, col_1r2r, xls_row_info[col_1r2r * 2])
    scores_sheets[0].write(row, xls_row_info[17], xls_row_info[18], green_style)
    for col_1s2s in range(1, 9):
        if xls_row_info[col_1s2s * 2 + 17] is not None:
            scores_sheets[1].write(row, col_1s2s, xls_row_info[col_1s2s * 2 + 17])
        if xls_row_info[col_1s2s * 2 + 18] is not None:
            ts_sheets[1].write(row, col_1s2s, xls_row_info[col_1s2s * 2 + 18])
    scores_sheets[1].write(row, xls_row_info[35], xls_row_info[36], green_style)
    for col_1r2s in range(1, 9):
        if xls_row_info[col_1r2s * 2 + 35] is not None:
            scores_sheets[2].write(row, col_1r2s, xls_row_info[col_1r2s * 2 + 35])
        if xls_row_info[col_1r2s * 2 + 36] is not None:
            ts_sheets[2].write(row, col_1r2s, xls_row_info[col_1r2s * 2 + 36])
    scores_sheets[2].write(row, xls_row_info[53], xls_row_info[54], green_style)
    for col_1s2r in range(1, 9):
        if xls_row_info[col_1s2r * 2 + 53] is not None:
            scores_sheets[3].write(row, col_1s2r, xls_row_info[col_1s2r * 2 + 53])
        if xls_row_info[col_1s2r * 2 + 54] is not None:
            ts_sheets[3].write(row, col_1s2r, xls_row_info[col_1s2r * 2 + 54])
    scores_sheets[3].write(row, xls_row_info[71], xls_row_info[72], green_style)
    if xls_row_info[18] < xls_row_info[36]:
        sum_sheets[0].write(row, 2, round(xls_row_info[18], 2), green_style)
        sum_sheets[0].write(row, 3, round(xls_row_info[36], 2), orange_style)
    else:
        sum_sheets[0].write(row, 2, round(xls_row_info[18], 2), orange_style)
        sum_sheets[0].write(row, 3, round(xls_row_info[36], 2), green_style)
    if xls_row_info[54] < xls_row_info[72]:
        sum_sheets[0].write(row, 4, round(xls_row_info[54], 2), green_style)
        sum_sheets[0].write(row, 5, round(xls_row_info[72], 2), orange_style)
    else:
        sum_sheets[0].write(row, 4, round(xls_row_info[54], 2), orange_style)
        sum_sheets[0].write(row, 5, round(xls_row_info[72], 2), green_style)

    sum_sheets[1].write(row, 1, xls_row_info[73])
    sum_sheets[1].write(row, 2, xls_row_info[74])
    sum_sheets[1].write(row, 3, xls_row_info[75])
    sum_sheets[1].write(row, 4, xls_row_info[76])
    sum_sheets[1].write(row, 5, xls_row_info[77])
    sum_sheets[1].write(row, 6, xls_row_info[78])
    sum_sheets[1].write(row, 7, xls_row_info[79])

    # sum_sheets[2].write(row, 1, xls_row_info[80])
    # sum_sheets[2].write(row, 2, xls_row_info[81])
    # sum_sheets[2].write(row, 3, xls_row_info[82])
    # sum_sheets[2].write(row, 4, xls_row_info[83])
    # sum_sheets[2].write(row, 5, xls_row_info[84])
    # sum_sheets[2].write(row, 6, xls_row_info[85])
    # sum_sheets[2].write(row, 7, xls_row_info[86])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, default=".")
    parser.add_argument("-baseline", "--baseline_variant_path", type=str)
    parser.add_argument("-stereo", "--preferred_stereoisomer", type=int, choices=[0, 1, 2, 3])
    parser.add_argument("-t", "--temperature", type=float, default=3)
    parser.add_argument("-rank", "--rank_fitness", type=str, choices=["m", "p", "f", "b", "s", "e", "d"])
    parser.add_argument("-rm", "--remove_redundant_decoys", action="store_true")
    args = parser.parse_args()
    baseline_variant = args.baseline_variant_path.rstrip("/").split("/")[-1]
    pdb = baseline_variant.split("_")[0]
    preferred_stereoisomer = args.preferred_stereoisomer
    temperature = args.temperature

    workbook = xlwt.Workbook(encoding="ascii")
    sum_sheets = list()
    sum_sheets.append(workbook.add_sheet("summary"))
    sum_sheets.append(workbook.add_sheet("argmin"))
    # sum_sheets.append(workbook.add_sheet("softmin"))
    scores_sheets = list()
    stereoisomers = ["1R2R", "1S2S", "1R2S", "1S2R"]
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
    xls_row_info, baseline_dGs, baseline_partition_function = read_variant_scores(\
            args.baseline_variant_path.rstrip("/"), preferred_stereoisomer, temperature=temperature, \
            remove_redundant_decoys=args.remove_redundant_decoys)
    write_xls_row(1, baseline_variant, xls_row_info, sum_sheets, scores_sheets, ts_sheets, \
            green_style, orange_style)
    scores_dict = dict()
    for enz_var in filter(lambda x: x != baseline_variant and os.path.isdir(os.path.join(args.directory, x)) \
                and x.startswith(pdb + "_"), sorted(os.listdir(args.directory))):
        xls_row_info, _, _ = read_variant_scores(os.path.join(args.directory, enz_var), preferred_stereoisomer, \
                baseline_dGs=baseline_dGs, baseline_partition_function=baseline_partition_function, \
                temperature=temperature, remove_redundant_decoys=args.remove_redundant_decoys)
        scores_dict[enz_var] = xls_row_info
    suffix = "_" + stereoisomers[preferred_stereoisomer]
    if args.rank_fitness:
        if args.rank_fitness == "m":
            i_fitness = 73
            suffix += "_multi-state_design"
        elif args.rank_fitness == "p":
            i_fitness = 74
            suffix += "_positive_design"
        elif args.rank_fitness == "f":
            i_fitness = 75
            suffix += "_ddG_fold"
        elif args.rank_fitness == "b":
            i_fitness = 76
            suffix += "_ddG_bind"
        elif args.rank_fitness == "s":
            i_fitness = 77
            suffix += "_stereoselectivity"
        elif args.rank_fitness == "e":
            i_fitness = 78
            suffix += "_enantioselectivity"
        elif args.rank_fitness == "d":
            i_fitness = 79
            suffix += "_diastereoselectivity"
        scores_dict = dict(sorted(scores_dict.items(), key=lambda x: x[1][i_fitness]))
    for row, enz_var_xls_row_info in enumerate(scores_dict.items()):
        enz_var, xls_row_info = enz_var_xls_row_info
        write_xls_row(row + 2, enz_var, xls_row_info, sum_sheets, scores_sheets, ts_sheets, \
            green_style, orange_style)
    workbook.save(pdb + suffix + ".xls")
