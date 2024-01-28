import os
from pymol import *


heme_dir = 'HEM'
for i in set(str(x) for x in range(1, 10)):
    for j in set(str(x) for x in range(0, 10)).union({'A','B','C','D','E','F','G',\
        'H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'}):
        for k in set(str(x) for x in range(0, 10)).union({'A','B','C','D','E','F','G',\
            'H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'}):
            unique_pdbs = dict() # {str: best_pdb: (float: lowest_resolution, {str: unique_pdb, str: pdb1, ...})}
            for l in set(str(x) for x in range(0, 10)).union({'A','B','C','D','E','F','G',\
                'H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'}):
                pdb = i + j + k + l
                if os.path.isfile(heme_dir + '/' + pdb + '/' + pdb + '_mono.pdb'):
                    if os.path.isfile(heme_dir + '/' + pdb + '/resolution'):
                        with open(heme_dir + '/' + pdb + '/resolution', 'r') as pf:
                            resolution = float(pf.read())
                    else:
                        resolution = 10.0
                    cmd.load(heme_dir + '/' + pdb + '/' + pdb + '_mono.pdb', pdb)
                    same_pdb = False
                    for best_pdb, res_pdbs in unique_pdbs.items():
                        lowest_resolution = res_pdbs[0]
                        redundant_pdb_set = res_pdbs[1]
                        rmsd = cmd.align(pdb, best_pdb)[0]
                        if rmsd < 0.8:
                            redundant_pdb_set.add(pdb)
                            if resolution < lowest_resolution:
                                find_better_pdb = True
                            else:
                                find_better_pdb = False
                            same_pdb = True
                            break
                    if not same_pdb:
                        unique_pdbs[pdb] = (resolution, {pdb})
                    elif find_better_pdb:
                        unique_pdbs[pdb] = (resolution, redundant_pdb_set)
                        del unique_pdbs[best_pdb]
            for best_pdb, res_pdbs in unique_pdbs.items():
                redundant_pdb_set = res_pdbs[1]
                with open(heme_dir + '/' + best_pdb + '/main', 'w') as pf:
                    pf.write(best_pdb + '\n' + str(redundant_pdb_set))
                redundant_pdb_set.remove(best_pdb)
                for redundant_pdb in redundant_pdb_set:
                    with open(heme_dir + '/' + redundant_pdb + '/redundant', 'w') as pf:
                        pf.write(best_pdb)

