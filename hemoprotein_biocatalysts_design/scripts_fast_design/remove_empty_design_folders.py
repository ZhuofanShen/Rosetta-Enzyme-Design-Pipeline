import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('dir', type=str)
parser.add_argument('-fold', type=str)
args = parser.parse_args()

for pdb in filter(lambda x: os.path.isdir(args.dir + '/' + x + '/cyclopropanation_EDA_styrene') \
    and os.path.isfile(args.dir + '/' + x + '/main') \
    and (not args.fold or os.path.isfile(args.dir + '/' + x + '/' + args.fold + '.fold')), os.listdir(args.dir)):
    for pos in ['ABOVE', 'BELOW']:
        if os.path.isdir(args.dir + '/' + pdb + '/cyclopropanation_EDA_styrene/' + pos + '/' + pdb + '_' + pos) \
                and not os.path.isfile(args.dir + '/' + pdb + '/cyclopropanation_EDA_styrene/' + pos + '/' + pdb + '_' + pos + '/' + pdb + '_' + pos + '_1R2S-rot1+/' + pdb + '_' + pos + '_1R2S-rot1+.fasc'):
            shutil.rmtree(args.dir + '/' + pdb + '/cyclopropanation_EDA_styrene/' + pos + '/' + pdb + '_' + pos)
