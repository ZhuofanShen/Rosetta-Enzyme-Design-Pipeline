#!/usr/bin/python3
import argparse
import os
import paramiko
import xlrd

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='/home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product')
parser.add_argument('-u', '--username', type=str, required=True, help='NetID')
parser.add_argument('-p', '--password', type=str, required=True)
parser.add_argument('-s', '--step', type=int, choices=[0, 1, 2], required=True)
args = parser.parse_args()

scaf_subs = args.directory[args.directory.rfind('/') + 1:]

workbook = xlrd.open_workbook(scaf_subs + '.xls')
directory_sheet = workbook.sheet_by_index(0)
design_sheet = workbook.sheet_by_index(args.step)

if args.step == 1:
    directory = 'design'
elif args.step == 2:
    directory = 'manual_design'

ssh_client = paramiko.SSHClient()
ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh_client.connect(hostname='amarel.hpc.rutgers.edu', port=22, username=args.username, password=args.password)
sftp_client = ssh_client.open_sftp()

for row in range(0, directory_sheet.nrows):
    variant = directory_sheet.cell_value(row, 0)
    print(variant)
    if not os.path.isdir(variant):
        os.mkdir(variant)
    if args.step == 0:
        sftp_client.get(args.directory + '/' + variant + '/' + variant + '.pdb', variant + '/' + variant + '.pdb')
        stdin, stdout, stderr = ssh_client.exec_command('ls ' + args.directory + '/' + variant + '/*.rotlib.pdb')
        rotlib_full_path = stdout.readlines()[0][:-1]
        rotlib = rotlib_full_path.split('/')[-1]
        sftp_client.get(rotlib_full_path, variant + '/' + rotlib)
    else:
        try:
            design = design_sheet.cell_value(row, 0)
            print(design)
        except IndexError:
            print(variant + ' has not been designed yet.')
            continue
        num = design_sheet.cell_value(row, 3)
        if num == '':
            print(variant + ' has not been designed yet.')
            continue
        print(str(row) + ': ' + str(int(num)))
        stdin, stdout, stderr = ssh_client.exec_command('ls ' + args.directory + '/' + variant + '/' + directory + '/*_' + str(int(num)) + '.pdb')
        out = stdout.readlines()
        if len(out) > 0:
            best_decoy_full_path = out[0][:-1]
            sftp_client.get(best_decoy_full_path, variant + '/' + variant + '_' + directory + '.pdb')
        else:
            sftp_client.get(args.directory + '/' + variant + '/' + variant + '_' + directory + '.pdb', variant + '/' + variant + '_' + directory + '.pdb')
        sftp_client.get(args.directory + '/' + variant + '/' + variant + '.pse', variant + '/' + variant + '.pse')
        if args.step == 1:
            sftp_client.get(args.directory + '/' + variant + '/pymol.txt', variant + '/pymol.txt')

sftp_client.close()
ssh_client.close()
