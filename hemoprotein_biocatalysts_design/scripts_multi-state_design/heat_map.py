import argparse
import plotly.express as px
import numpy as np


argparser = argparse.ArgumentParser()
argparser.add_argument("pssm", type=str, help="pssm file")
argparser.add_argument("-pos", "--positions", type=int, nargs="*", help="residue positions included in the plot")
args = argparser.parse_args()

if args.pssm.endswith(".csv"):
    all_probs_concat = np.loadtxt(args.pssm, delimiter=",", dtype=float)
elif args.pssm.endswith(".npz"):
    mpnn_output = np.load(args.pssm)
    all_probs_concat = mpnn_output["log_p"][0]

if args.positions:
    all_probs_concat = np.transpose(all_probs_concat)
    all_probs_concat = np.transpose(list(filter(lambda x: x is not None, 
            map(lambda x, y: x if y in args.positions else None, 
            all_probs_concat, range(len(all_probs_concat))))))
    fig = px.imshow(all_probs_concat,
                    labels=dict(x="positions", y="amino acids", color="probability"),
                    y=list('ACDEFGHIKLMNPQRSTVWYX'),
                    x=list(str(p) for p in args.positions),
                    template="simple_white"
                )
else:
    fig = px.imshow(all_probs_concat,
                    labels=dict(x="positions", y="amino acids", color="probability"),
                    y=list('ACDEFGHIKLMNPQRSTVWYX'),
                    template="simple_white"
               )
# fig.update_xaxes(side="top")
# fig.show()
fig.write_image(args.pssm.strip(".npz").strip(".csv") + ".jpeg")
