import argparse
import json
import os
from pyrosetta import *
import xlwt

parser = argparse.ArgumentParser()
parser.add_argument('dir', type=str)
parser.add_argument('-fold', type=str)
args = parser.parse_args()
workbook = xlwt.Workbook(encoding="ascii")
design_sheet = workbook.add_sheet(args.fold if args.fold else "pdb")
design_sheet.write(0, 0, "pdb")
design_sheet.write(0, 1, "molecule")
design_sheet.write(0, 2, "synonym")
design_sheet.write(0, 3, "engineered")
design_sheet.write(0, 4, "organism")
design_sheet.write(0, 5, "expression")
design_sheet.write(0, 6, "keywords")
row = 1
for pdb in filter(lambda x: os.path.isfile(args.dir + '/' + x + '/main') \
    and (not args.fold or os.path.isfile(args.dir + '/' + x + '/' + args.fold + '.fold')) \
    and not os.path.isfile(args.dir + '/' + x + '/MYOGLOBIN.fold') \
    and not os.path.isfile(args.dir + '/' + x + '/HEMOGLOBIN.fold'), sorted(os.listdir(args.dir))):
    molecule = str()
    synonym = str()
    engineered = str()
    organism = str()
    expression = str()
    kwd1 = str()
    kwd2 = str()
    kwd3 = str()
    kwd4 = str()
    kwd5 = str()
    with open(args.dir + '/' + pdb + '/' + pdb + '.pdb', 'r') as pf:
        for line in pf:
            if line.startswith('COMPND   2 MOLECULE:'):
                molecule = line[20:-1].strip(' ').strip(';')
            elif line.startswith('COMPND   4 SYNONYM:'):
                synonym = line[19:-1].strip(' ').strip(';')
            elif line.startswith('COMPND   6 ENGINEERED:'):
                engineered = line[22:-1].strip(' ').strip(';')
            elif line.startswith('SOURCE   2 ORGANISM_SCIENTIFIC:'):
                organism = line[31:-1].strip(' ').strip(';')
            elif line.startswith('SOURCE   4 EXPRESSION_SYSTEM:'):
                expression = line[29:-1].strip(' ').strip(';')
            elif line.startswith('KEYWDS'):
                kwd1 = line[6:-1].strip(' ')
            elif line.startswith('KEYWDS   2'):
                kwd2 = line[10:-1].strip(' ')
            elif line.startswith('KEYWDS   3'):
                kwd3 = line[10:-1].strip(' ')
            elif line.startswith('KEYWDS   4'):
                kwd4 = line[10:-1].strip(' ')
            elif line.startswith('KEYWDS   5'):
                kwd5 = line[10:-1].strip(' ')
    kwds = ', '.join([kwd1, kwd2, kwd3, kwd4, kwd5])
    design_sheet.write(row, 0, pdb)
    design_sheet.write(row, 1, molecule)
    design_sheet.write(row, 2, synonym)
    design_sheet.write(row, 3, engineered)
    design_sheet.write(row, 4, organism)
    design_sheet.write(row, 5, expression)
    design_sheet.write(row, 6, kwds)
    row += 1
workbook.save(args.dir + '.xls')
