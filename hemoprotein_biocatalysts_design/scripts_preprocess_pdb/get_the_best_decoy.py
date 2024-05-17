import argparse
import json
import os
import shutil


def read_scores_from_fasc(var):
    file_names = os.listdir(var)
    for file_name in file_names:
        if file_name.endswith(".fasc"):
            decoy_scores_str = "["
            with open(var + "/" + file_name) as fasc:
                for line in fasc:
                    decoy_scores_str += line[:-1]
                    decoy_scores_str += ","
            decoy_scores_str = decoy_scores_str[:-1] + "]"
    return json.loads(decoy_scores_str)

def extract_n_decoys(mb_scores, n=1, is_reversed=False):
    mb_scores = sorted(mb_scores, key=lambda entry: entry["total_score"], reverse=is_reversed)
    stable_mb_scores = dict()
    for mb_score in mb_scores[0:n]:
        stable_mb_scores.update({str(mb_score["decoy"]): mb_score["total_score"]})
    return stable_mb_scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', type=str)
    parser.add_argument('-s', '--step', type=int, choices=[1, 2], required=True)
    args = parser.parse_args()
    if args.step == 1:
        for pdb in os.listdir(args.dir):
            if os.path.isfile(args.dir + '/' + pdb + '/' + pdb + '_relax.pdb'):
                continue
            #1
            if os.path.isdir(args.dir + '/' + pdb + '/relax1'):
                if not os.path.isfile(args.dir + '/' + pdb + '/relax1/' + pdb + '_clean_relaxed.fasc') and \
                        not os.path.isfile(args.dir + '/' + pdb + '/relax1/' + pdb + '_mono_relaxed.fasc'):
                    shutil.rmtree(args.dir + '/' + pdb + '/relax1')
                    continue
                scores1 = read_scores_from_fasc(args.dir + '/' + pdb + '/relax1')
                for name1, score1 in extract_n_decoys(scores1, 1).items():
                    break
            else:
                name1 = "None"
                score1 = 100000000
            #2
            if os.path.isdir(args.dir + '/' + pdb + '/relax2'):
                if not os.path.isfile(args.dir + '/' + pdb + '/relax2/' + pdb + '_clean_relaxed.fasc') and \
                        not os.path.isfile(args.dir + '/' + pdb + '/relax2/' + pdb + '_mono_relaxed.fasc'):
                    shutil.rmtree(args.dir + '/' + pdb + '/relax2')
                    continue
                scores2 = read_scores_from_fasc(args.dir + '/' + pdb + '/relax2')
                for name2, score2 in extract_n_decoys(scores2, 1).items():
                    break
            else:
                name2 = "None"
                score2 = 100000000
            #clash
            if (score1 > 0 and score1 != 100000000) or (score2 > 0 and score2 != 100000000):
                with open(args.dir + '/clash', 'a') as pf:
                    pf.write(pdb + '\n')
                continue
            #compare
            if score1 != 100000000 and score2 != 100000000:
                if score1 <= score2:
                    shutil.copy(args.dir + '/' + pdb + '/relax1/' + name1, args.dir + '/' + pdb + '/' + pdb + '_relax.pdb')
                else:
                    shutil.copy(args.dir + '/' + pdb + '/relax2/' + name2, args.dir + '/' + pdb + '/' + pdb + '_relax.pdb')
    else:
        variant = args.dir.split('/')[-1]
        tmp = variant.split('_')
        protein = tmp[0] + '_' + tmp[1]
        for configuration in ['1R2R', '1S2S', '1R2S', '1S2R']:
            for carbene in ['rot1', 'rot2', 'rot3', 'rot4']:
                for ester in ['+', '-']:
                    mb_scores = read_scores_from_fasc(args.dir + '/' + variant + '_' + configuration + '-' + carbene + ester)
                    for name, score in extract_n_decoys(mb_scores, 1).items():
                        break
                    for new_name in os.listdir(args.dir + '/' + variant + '_' + configuration + '-' + carbene + ester):
                        if new_name.startswith(protein + '_' + configuration + '-' + carbene + ester + '_' + name[:-4].split('_')[-1]):
                            print(args.dir + '/' + variant + '_' + configuration + '-' + carbene + ester + '/' + new_name)
                            break
                    shutil.copy(args.dir + '/' + variant + '_' + configuration + '-' + carbene + ester + '/' + new_name, protein + '_' + configuration + '-' + carbene + ester + '.pdb')
