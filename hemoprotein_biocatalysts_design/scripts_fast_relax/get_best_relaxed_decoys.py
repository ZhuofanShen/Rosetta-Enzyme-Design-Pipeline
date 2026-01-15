import argparse
import json
import os
import shutil


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', type=str)
    parser.add_argument('-s', '--step', type=int, choices=[1, 2], default=2)
    args = parser.parse_args()
    if args.step == 1:
        for pdb in os.listdir(args.dir):
            if os.path.isfile(args.dir + '/' + pdb + '/' + pdb + '_relaxed.pdb'):
                continue
            best_score = 100000000
            for i in range(1, 11):
                if os.path.isdir(args.dir + '/' + pdb + '/relax' + str(i)):
                    scores1 = read_scores_from_fasc(args.dir + '/' + pdb + '/relax' + str(i))
                    for name1, score1 in extract_n_decoys(scores1, 1).items():
                        break
                    if score1 < best_score:
                        best_score = score1
                        best_decoy = name1
            shutil.copy(args.dir + '/' + pdb + '/relax' + str(i) + '/' + best_decoy, args.dir + '/' + pdb + '/' + pdb + '_relaxed.pdb')

    elif args.step == 2:
        pr = args.dir.split("_")[0]
        scores = read_scores_from_fasc(args.dir + '/' + pr)
        if scores:
            for name, score in extract_n_decoys(scores, 1).items():
                break
            for new_name in os.listdir(args.dir + '/' + pr):
                if new_name.endswith('_' + name[:-4].split('_')[-1] + '.pdb'):
                    print(args.dir + '/' + pr + '/' + new_name)
                    break
            shutil.copy(args.dir + '/' + pr + '/' + new_name, args.dir + '/' + pr + '.pdb')
        for configuration in ['1R2R', '1S2S', '1R2S', '1S2R']:
            for carbene in ['rot1', 'rot2', 'rot3', 'rot4']:
                for ester in ['+', '-']:
                    scores = read_scores_from_fasc(args.dir + '/' + pr + '_' + configuration + '-' + carbene + ester)
                    if scores:
                        for name, score in extract_n_decoys(scores, 1).items():
                            break
                        for new_name in os.listdir(args.dir + '/' + pr + '_' + configuration + '-' + carbene + ester):
                            if new_name.startswith(pr + '_' + configuration + '-' + carbene + ester + '_' + name[:-4].split('_')[-1] + '.'):
                                print(args.dir + '/' + pr + '_' + configuration + '-' + carbene + ester + '/' + new_name)
                                break
                        shutil.copy(args.dir + '/' + pr + '_' + configuration + '-' + carbene + ester + '/' + new_name, args.dir + '/' + pr + '_' + configuration + '-' + carbene + ester + '.pdb')
