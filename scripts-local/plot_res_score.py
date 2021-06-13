import argparse
import matplotlib.pyplot as plt


def read_res_sc_from_pdb(pdb_dir):
    with open(pdb_dir, "r") as pdb:
        flag = False
        scores = []
        for line in pdb:
            if not flag:
                if line.startswith("pose"):
                    flag = True
            else:
                if line.startswith("#END_POSE_ENERGIES_TABLE"):
                    break
                else:
                    res_scores = line.split(" ")
                    res_score = float(res_scores[-1][:-1]) - float(res_scores[-2]) - float(res_scores[-13])
                    scores.append(res_score)
    return scores


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", type=str)
    parser.add_argument("-pdbs", type=str)
    parser.add_argument("-diff", type=str)
    args = parser.parse_args()

    if args.dir is not None and args.pdbs is None:
        scores = read_res_sc_from_pdb(args.dir)
        indexes = range(1, len(scores) + 1)
        plt.plot(indexes, scores, color='#0000FF', \
	        marker = "D", markerfacecolor='#00008B', markersize=2)
        #plt.plot(indexes, scores, color='#0000FF')
        #plt.scatter(indexes, scores, marker = "D", c='#0000FF', s=2)
        plt.xticks(range(1, len(scores) + 1, 5))
        #plt.xticks(indexes)
        plt.yticks(range(-24, 15, 2))
        plt.xlabel("residue index")
        plt.ylabel("Rosetta energy unit")
        #plt.title(args.dir)
        plt.show()
    elif args.pdbs is not None and args.dir is None:
        pdbs = args.pdbs.split(",")
        scores1 = read_res_sc_from_pdb(pdbs[0])
        scores2 = read_res_sc_from_pdb(pdbs[1])
        differences = []
        for idx, res_scores1 in enumerate(scores1):
            differences.append(res_scores1 - scores2[idx])
        res_indexes = range(1, len(differences) + 1)
        if args.diff != "diff":
            plt.plot(res_indexes, scores1, color='#FF0000', \
                marker = "D", markerfacecolor='#8B0000', markersize=2)
            plt.plot(res_indexes, scores2, color='#800080', \
                marker = "D", markerfacecolor='#9400D3', markersize=2)
            #plt.plot(res_indexes, scores1, color='#FF0000')
            #plt.plot(res_indexes, scores2, color='#800080')
        if args.diff != "no":
            plt.plot(res_indexes, differences, color='#0000FF', \
                marker = "D", markerfacecolor='#00008B', markersize=2)
            #plt.plot(res_indexes, differences, color='#0000FF')
        plt.xticks(range(1, len(differences) + 1, 5))
        #plt.xticks(res_indexes)
        plt.yticks(range(-24, 15, 2))
        plt.xlabel("residue index")
        plt.ylabel("Rosetta energy unit")
        plt.title(args.pdbs)
        plt.show()
