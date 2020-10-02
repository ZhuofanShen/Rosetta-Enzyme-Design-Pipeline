#!/usr/bin/python3
import argparse
import os
import paramiko
import xlrd

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=str, help='/home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product')
parser.add_argument('-u', '--username', type=str, required=True, help='NetID')
parser.add_argument('-p', '--password', type=str, required=True)
args = parser.parse_args()

directory = args.directory[args.directory.rfind('/') + 1:]
os.mkdir(directory)

workbook = xlrd.open_workbook(directory + '.xls')
score_sheet = workbook.sheet_by_index(0)
num_sheet = workbook.sheet_by_index(1)

ssh_client = paramiko.SSHClient()
ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh_client.connect(hostname='amarel.hpc.rutgers.edu', port=22, username=args.username, password=args.password)
sftp_client = ssh_client.open_sftp()

for row in range(0, num_sheet.nrows):
    design = score_sheet.cell_value(row, 0)
    os.mkdir(directory + '/' + design)
    match = num_sheet.cell_value(row, 0)
    # design
    num = num_sheet.cell_value(row, 1)
    stdin, stdout, stderr = ssh_client.exec_command('ls ' + args.directory + '/' + match + '/design/*_' + str(int(num)) + '.pdb')
    best_decoy_full_name = stdout.readlines()[0][:-1]
    sftp_client.get(best_decoy_full_name, directory + '/' + design + '/' + design + '.pdb')
    col = 2
    # revert
    while True:
        try:
            reverted_design = num_sheet.cell_value(row, col)
            reverted_design_pdb_index = score_sheet.cell_value(row, col)
        except:
            break
        col += 1
        num = num_sheet.cell_value(row, col)
        col += 1
        stdin, stdout, stderr = ssh_client.exec_command('ls ' + args.directory + '/' + match + '/' + reverted_design + '/*_' + str(int(num)) + '.pdb')
        best_decoy_full_name = stdout.readlines()[0][:-1]
        sftp_client.get(best_decoy_full_name, directory + '/' + design + '/revert_' + reverted_design_pdb_index + '.pdb')
sftp_client.close()
ssh_client.close()
os.rename(directory + '.xls', directory + '/' + directory + '.xls')
