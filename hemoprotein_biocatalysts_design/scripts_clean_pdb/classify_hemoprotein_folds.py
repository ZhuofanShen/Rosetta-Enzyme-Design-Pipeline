import os

folds = ['GLOBIN', 'MYOGLOBIN', 'HEMOGLOBIN', 'DEHALOPEROXIDASE', 'ERYTHROCRUORIN', 'LEGHEMOGLOBIN', \
    'CYTOGLOBIN', 'NEUROGLOBIN', 'THERMOGLOBIN', 'PROTOGLOBIN', 'FLAVOHEMOGLOBIN', 'TRUNCATED HEMOGLOBIN', 'MINI-HEMOGLOBIN', \
    'CYTOCHROME', 'P450', 'CYP51', 'CYP121', 'CYTOCHROME C', 'CYTOCHROME C PRIME', 'CYTOCHROME B5', 'CYTOCHROME BC1', 'CYTOCHROME F', \
    'PEROXIDASE', 'LACTOPEROXIDASE', 'PEROXINECTIN', 'NITRITE REDUCTASE', \
    'HEME OXYGENASE', 'CYCLOOXYGENASE', 'DIOXYGENASE', 'PEROXYGENASE', 'PRNB', 'NITRIC OXIDE SYNTHASE', \
    'CYSTATHIONINE BETA-SYNTHASE', 'NITROPHORIN', 'DISMUTASE', 'DEHYDROGENASE', \
    'H-NOX', 'HEME-DEGRADING', 'HEMOPHORE', 'HEME TRANSFER PROTEIN', 'SENSOR DOMAIN', 'HMUT', 'HMOB', 'PROTEUS MIRABILIS CATALASE']
for output_dir in ['HEM_monomer']:
    for pdb in os.listdir(output_dir):
        with open(output_dir + '/' + pdb + '/' + pdb + '.pdb', 'r') as pf:
            for line in pf:
                if line.startswith('TITLE'):
                    title = line[5:-1].strip(' ')
                    for fold in folds:
                        if fold in title:
                            with open(output_dir + '/' + pdb + '/' + fold + '.fold', 'w') as _:
                                # break
                                pass
                        # break
                elif line.startswith('COMPND   2 MOLECULE:'):
                    info = line[20:-1].strip(' ').strip(';')
                    for fold in folds:
                        if fold in info:
                            with open(output_dir + '/' + pdb + '/' + fold + '.fold', 'w') as _:
                                # break
                                pass
                        # break
                elif line.startswith('KEYWDS'):
                    kwds = line[6:-1].strip(' ')
                    for fold in folds:
                        if fold in kwds:
                            with open(output_dir + '/' + pdb + '/' + fold + '.fold', 'w') as _:
                                # break
                                pass
                    # break
                elif line.startswith('REMARK'):
                    break
