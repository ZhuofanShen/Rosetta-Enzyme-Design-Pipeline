import argparse
import json
import os
import torch
from transformers import AutoTokenizer, EsmForMaskedLM



def esm_scoring(sequence, tokenizer, model_output_logits):
    total_logp = 0
    for ctr, aas in enumerate(sequence):
        position_logits = model_output_logits[0, ctr+1]
        probabilities = torch.nn.functional.softmax(position_logits, dim=-1)
        log_probabilities = torch.log(probabilities)
        wt_token_id = tokenizer.convert_tokens_to_ids(aas)
        wt_log_prob = log_probabilities[wt_token_id].item()
        total_logp += wt_log_prob
    return total_logp

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str)
    parser.add_argument('-lig', "--ligand", type=str, default="O2beY-CYS")
    parser.add_argument('--process_variants', type=str, nargs="*", default=list()) # e.g., 1ysb-A_S6_E126
    parser.add_argument('-od', '--output_directory', default='mpnn')
    parser.add_argument('-esm_n', '--best_n_esm_scored_sequences', type=int, default=10)
    args = parser.parse_args()

    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    model = EsmForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
    for pdb in os.listdir(args.directory):
        pdb_path = os.path.join(args.directory, pdb)
        for pdb_chain_ligand in filter(lambda x: x.startswith(pdb) and os.path.isdir(\
                os.path.join(pdb_path, x)), os.listdir(pdb_path)):
            pdb_chain, uaa_nucleophile = pdb_chain_ligand.split("_")
            pdb_ligand_path = os.path.join(pdb_path, pdb_chain_ligand)
            for variant in filter(lambda x: os.path.isdir(os.path.join(pdb_ligand_path, x)) and x.startswith(pdb_chain), os.listdir(pdb_ligand_path)):
                mutations = variant.split("_")[1:]
                uaa_pose_index = int(mutations[0][1:])
                nucleophile_pose_index = int(mutations[1][1:])
                if len(args.process_variants) > 0 and variant not in args.process_variants:
                    continue
                mpnn_path = os.path.join(pdb_ligand_path, variant, args.output_directory)
                mpnn_seqs_path = os.path.join(mpnn_path, "seqs")
                for mpnn_fasta in os.listdir(mpnn_seqs_path):
                    if mpnn_fasta.endswith(".fa"):
                        mpnn_fasta = os.path.join(mpnn_seqs_path, mpnn_fasta)
                        break
                seqs_list = list()
                total_logp_list = list()
                with open(mpnn_fasta, "r") as pf:
                    i_seq = 0
                    for line in pf:
                        if not line.startswith(">"):
                            seq_tmp = line.split(":")[0]
                            if uaa_pose_index < nucleophile_pose_index:
                                seq_canonical = seq_tmp[:uaa_pose_index-1] + "Y" + seq_tmp[uaa_pose_index-1:]
                                seq_canonical = seq_canonical[:nucleophile_pose_index-1] + "C" + seq_canonical[nucleophile_pose_index-1:]
                                seq = seq_tmp[:uaa_pose_index-1] + "Z" + seq_tmp[uaa_pose_index-1:]
                                seq = seq[:nucleophile_pose_index-1] + "X" + seq[nucleophile_pose_index-1:]
                            else:
                                seq_canonical = seq_tmp[:nucleophile_pose_index-1] + "C" + seq_tmp[nucleophile_pose_index-1:]
                                seq_canonical = seq_canonical[:uaa_pose_index-1] + "Y" + seq_canonical[uaa_pose_index-1:]
                                seq = seq_tmp[:nucleophile_pose_index-1] + "X" + seq_tmp[nucleophile_pose_index-1:]
                                seq = seq[:uaa_pose_index-1] + "Z" + seq[uaa_pose_index-1:]
                            if i_seq == 0:
                                inputs = tokenizer(seq_canonical, return_tensors="pt")
                                with torch.no_grad():
                                    outputs = model(**inputs)
                            seqs_list.append(seq)
                            total_logp_list.append(esm_scoring(seq_canonical, tokenizer, outputs.logits))
                            i_seq += 1
                with open(os.path.join(mpnn_path, "sequences.json"), "w") as pf:
                    pf.write(json.dumps(seqs_list))
                with open(os.path.join(mpnn_path, "esm_scores.json"), "w") as pf:
                    pf.write(json.dumps(total_logp_list))
                with open(os.path.join(pdb_ligand_path, args.output_directory + "_sequences.fasta"), "a") as pf:
                    unique_seq_set = set()
                    for mpnn_variant in sorted(zip(total_logp_list[1:], seqs_list[1:], list(range(1, len(seqs_list) + 1))), key=lambda x: x[0], reverse=True):
                        if mpnn_variant[1] not in unique_seq_set:
                            pf.write("> " + variant + " " + str(mpnn_variant[2]) + " " + str(round(mpnn_variant[0], 2)) + "\n")
                            pf.write(mpnn_variant[1] + "\n")
                            unique_seq_set.add(mpnn_variant[1])
                            if len(unique_seq_set) == args.best_n_esm_scored_sequences:
                                break
                    pf.write("\n\n")
