import subprocess
import os
import string
import numpy as np
import json


def run_hhblits(input_fasta, output_dir, database, e_vals, cpus=4, maxmem=32):
    # Setup the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each e_val and run the command
    file_in = input_fasta.split("/")[-1].split(".")[0]
    # Ensure sci notation is converted to string
    e_val_str = f"{e_vals[0]:.0e}".replace("+", "")
    oa3m_0 = os.path.join(output_dir, f"{file_in}_{e_val_str}.a3m")
    # Base command
    base_command = f"hhblits -o /dev/null -mact 0.35 -maxfilt 20000 -neffmax 20 -all -cpu {cpus} -realign_max 20000 -maxmem {maxmem} -n 4 -d {database}"
    # Perform initial search
    command = f"{base_command} -i {input_fasta} -oa3m {oa3m_0} -e {e_val_str}"
    subprocess.run(command, shell=True)
    # Filter the a3m file
    oa3m_0 = run_hhfilter(oa3m_0, output_dir)
    # Perform iterative search - read in the previous a3m file and perform filtering after each iteration
    oa3m_output_file = oa3m_0
    for e_val in e_vals[1:]:
        # Ensure sci notation is converted to string
        e_val_str = e_val_str + "_" + f"{e_val:.0e}".replace("+", "")
        oa3m_input_file = oa3m_output_file
        oa3m_output_file = os.path.join(output_dir, f"{file_in}_{e_val_str}.a3m")
        command = f"{base_command} -i {oa3m_input_file} -oa3m {oa3m_output_file} -e {e_val_str} -v 0"
        subprocess.run(command, shell=True)
        # Filter the a3m file
        oa3m_output_file = run_hhfilter(oa3m_output_file, output_dir)
    return oa3m_output_file

def run_hhfilter(input_a3m, output_dir, id=90, cov=50, qid=30):
    # Setup the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each e_val and run the command
    output_file = input_a3m.split("/")[-1].split(".")[0] + "_id90cov50qid30.a3m"
    output_file = os.path.join(output_dir, output_file)
    command = f"hhfilter -i {input_a3m} -o {output_file} -id {id} -cov {cov} -qid {qid} -maxseq 10000"
    # Run the command
    subprocess.run(command, shell=True)
    return output_file

def parse_a3m(filename):

    # read A3M file line by line
    lab,seq = [],[] # labels and sequences
    for line in open(filename, "r"):
        if line[0] == '>':
            lab.append(line.split()[0][1:])
            seq.append("")
        else:
            seq[-1] += line.rstrip()

    # parse sequences
    msa,ins = [],[]
    table = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    nrow,ncol = len(seq),len(seq[0])

    for seqi in seq:

        # remove lowercase letters and append to MSA
        msa.append(seqi.translate(table))

        # 0 - match or gap; 1 - insertion
        a = np.array([0 if c.isupper() or c=='-' else 1 for c in seqi])
        i = np.zeros((ncol))

        if np.sum(a) > 0:
            # positions of insertions
            pos = np.where(a==1)[0]

            # shift by occurrence
            a = pos - np.arange(pos.shape[0])

            # position of insertions in the cleaned sequence
            # and their length
            pos,num = np.unique(a, return_counts=True)
            i[pos[pos<ncol]] = num[pos<ncol]

        # append to the matrix of insetions
        ins.append(i)

    # convert letters into numbers
    alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV-"), dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in msa], dtype='|U1')
    msa_num = msa.copy().astype('|S1').view(np.uint8)

    for i in range(alphabet.shape[0]):
        msa_num[msa_num == alphabet[i]] = i

    # treat all unknown characters as gaps
    msa[msa_num > 20] = '-'
    msa_num[msa_num > 20] = 20

    ins = np.array(ins, dtype=np.uint8)

    return {"msa":msa, "msa_num":msa_num, "labels":lab, "insertions":ins}

def get_conserved_positions(input_a3m, structure_name, output_dir, interface_residues=None):
    # Read
    aln = parse_a3m(input_a3m)

    msa = aln['msa']
    msa_num = aln['msa_num']
    L = msa.shape[1]

    counts = np.stack([np.bincount(column,minlength=21) for column in msa_num.T]).T
    max_count = np.max(counts,axis=0)

    freq = counts/msa_num.shape[0]
    freq_norm = freq[:20]/freq[:20].sum(axis=0)
    max_freq_norm = np.max(freq_norm,axis=0) # frequency of most frequent AA, excluding gaps
    max_freq_norm[max_count<10] = 0 # don't choose positions that have too low counts

    frac_conserved_0 = 0.3 # 30% conservation
    frac_conserved_1 = 0.5 # 50% conservation
    frac_conserved_2 = 0.7 # 70% conservation

    conserved_0 = np.sort(np.argsort(max_freq_norm)[::-1][:int(L*frac_conserved_0)]+1) # make 1-indexed to input to MPNN
    conserved_1 = np.sort(np.argsort(max_freq_norm)[::-1][:int(L*frac_conserved_1)]+1) # make 1-indexed to input to MPNN
    conserved_2 = np.sort(np.argsort(max_freq_norm)[::-1][:int(L*frac_conserved_2)]+1) # make 1-indexed to input to MPNN

    # Save the conserved positions to csv without the interface residues
    for percent, conserved in zip([30, 50, 70], [conserved_0, conserved_1, conserved_2]):
        # Write as a simple csv
        output_file_name = input_a3m.split("/")[-1].split(".")[0].split("_")[0] + f'_ionly_cpos_{percent}.csv'
        output_file = os.path.join(output_dir, output_file_name)
        with open(output_file, "w") as f:
            f.write("Conserved Positions\n")
            for pos in conserved:
                f.write(f"{pos}\n")

    # Merge in the interface residues with set then sort
    if interface_residues is not None:
        conserved_0 = sorted(list(set(conserved_0) | set(interface_residues)))
        conserved_1 = sorted(list(set(conserved_1) | set(interface_residues)))
        conserved_2 = sorted(list(set(conserved_2) | set(interface_residues)))
    
    # Write the conserved positions to a CSV file
    for percent, conserved in zip([30, 50, 70], [conserved_0, conserved_1, conserved_2]):
        output_file_name = input_a3m.split("/")[-1].split(".")[0].split("_")[0] + f'_cpos_{percent}.jsonl'
        output_file = os.path.join(output_dir, output_file_name)
        # String of the conserved positions i.e. [1, 2, 3, 4, 5]
        conserved_str = ", ".join([str(pos) for pos in conserved])
        conserved_str = '[' + conserved_str + ']'
        with open(output_file, "w") as f:
            # Structure: "A" : [1, 2, 3, 4, 5]
            json_string = json.dumps({structure_name.split(".")[0]: {"A": "-"}})
            json_string = json_string[:-5] + conserved_str + "}}"
            f.write(json_string + "\n")


e_vals = [1e-50, 1e-30, 1e-10, 1e-4]
database="/projects/f_sdk94_1/Tools/RoseTTAFold-All-Atom/uniclust/UniRef30_2020_06"
input_file = "1ysb.fasta"
a3m_file = run_hhblits(input_file, "output", database, e_vals)
get_conserved_positions(a3m_file, input_file, "output")
