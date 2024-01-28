import argparse
import logomaker
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


title = None
yticks = None

# site_aas_list = {"29": {"I", "V", "N", "D"}, "64": {"V", "T", "A"}, "68": {"L", "F", "H"}}
# ylabel = "multi–state design fitness ($10^{1}$)"
# yticks = [0, 10, 20, 30]
# ytick_labels = ["0", "1", "2", "3"]
# title = "round 1 design"
# ylabel = "positive design fitness"
# ylabel = "stability"

# site_aas_list = {"29": {"I", "V", "N", "D"}, "64": {"V", "T", "A"}, "68": {"L", "F", "H"}}
# ylabel = "multi–state design fitness"
# yticks = [0, 50, 100, 150]
# ytick_labels = ["0", "5", "10", "15"]
# title = "round 2 design"
# ylabel = "positive design fitness ($10^{1}$)"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]
# ylabel = "stability"

# site_aas_list = {"29": {"I", "V", "N", "D"}, "64": {"V", "T", "A"}, "68": {"L", "F", "H"}}
# ylabel = "multi–state design fitness ($10^{1}$)"
# yticks = [0, 50, 100]
# ytick_labels = ["0", "5", "10"]
# title = "round 3 design"
# ylabel = "positive design fitness"
# ylabel = "stability"

# # all sites no F33 L40 F46
# site_aas_list = {"29": {"I", "V", "N", "D"}, "32": {"L"}, "43": {"F"}, "61": {"L"}, \
#                  "64": {"V", "T", "A"}, "65": {"G"}, "68": {"L", "F", "H"}, "69": {"L"}, "107": {"I"}}
# ylabel = "multi–state design fitness ($10^{1}$)"
# yticks = [0, 50, 100, 150]
# ytick_labels = ["0", "5", "10", "15"]
# ylabel = "positive design fitness"

# # F33 L40 F46
# site_aas_list = {"33": {"F", "L", "V", "A"}, "40": {"L"}, "46": {"F"}}
# ylabel = "multi–state design fitness ($10^{3}$)"
# yticks = [0, 1000, 2000]
# ytick_labels = ["0", "1", "2"]
# # ylabel = "multi–state design fitness ($10^{5}$)"
# # yticks = [0, 100000, 200000, 300000]
# # ytick_labels = ["0", "1", "2", "3"]
# ylabel = "positive design fitness ($10^{1}$)"
# yticks = [0, 10, 20]
# ytick_labels = ["0", "1", "2"]

# all sites
# site_aas_list = {"29": {"I", "V", "N", "D"}, "32": {"L"}, "33": {"F", "L", "V", "A"}, "40": {"L"}, "43": {"F"}, "46": {"F"}, "61": {"L"}, \
#                  "64": {"V", "T", "A"}, "65": {"G"}, "68": {"L", "F", "H"}, "69": {"L"}, "107": {"I"}}
# ylabel = "positive design fitness ($10^{1}$)"
# yticks = [0, 10, 20]
# ytick_labels = ["0", "1", "2"]
# ylabel = "stability"

# round 4 design
# site_aas_list = {"29": {"I", "V", "N", "D"}, "32": {"L"}, "33": {"F", "L", "V", "A"}, "40": {"L"}, "43": {"F"}, "46": {"F"}, "61": {"L"}, \
#                  "64": {"V", "T", "A"}, "65": {"G"}, "68": {"L", "F", "H"}, "69": {"L"}, "107": {"I"}}
# ylabel = "multi–state design fitness ($10^{2}$)"
# yticks = [0, 200, 400, 600, 800]
# ytick_labels = ["0", "2", "4", "6", "8"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# round 4 design pCF3
# site_aas_list = {"29": {"I", "V", "N", "D"}, "32": {"L"}, "33": {"F", "L", "V", "A"}, "40": {"L"}, "43": {"F"}, "46": {"F"}, "61": {"L"}, \
#                  "64": {"V", "T", "A"}, "65": {"G"}, "68": {"L", "F", "H"}, "69": {"L"}, "107": {"I"}}
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 50000, 100000, 150000]
# ytick_labels = ["0", "5", "10", "15"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # P450cam
# site_aas_list = {"87": {"F", "L", "I", "M"}, "96": {"F", "Y"}, "101": {"L", "V", "M", "T"}, \
#         "185": {"T"}, "244": {"L"}, "247": {"L"}, "252": {"T"}, "295": {"V"}, "297": {"T"}, \
#         "395": {"I", "L"}, "396": {"I"}}
# # round 1
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000]
# ytick_labels = ["0", "1", "2"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # round 2
# ylabel = "multi–state design fitness ($10^{3}$)"
# yticks = [0, 5000, 10000, 15000]
# ytick_labels = ["0", "5", "10", "15"]
# ylabel = "positive design fitness"
# yticks = [0, 2, 4, 6, 8]
# ytick_labels = ["0", "2", "4", "6", "8"]

# # LFI
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000]
# ytick_labels = ["0", "1", "2"]
# ylabel = "positive design fitness"
# yticks = [0, 2, 4, 6, 8]
# ytick_labels = ["0", "2", "4", "6", "8"]

# # FFI
# ylabel = "multi–state design fitness ($10^{3}$)"
# yticks = [0, 5000, 10000, 15000]
# ytick_labels = ["0", "5", "10", "15"]
# ylabel = "positive design fitness"

# # LFL
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000, 30000]
# ytick_labels = ["0", "1", "2", "3"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # FFL
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000]
# ytick_labels = ["0", "1", "2"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # FFL pCF3
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000, 30000]
# ytick_labels = ["0", "1", "2", "3"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]


# # IDO1
# site_aas_list = {"126": {"Y"}, "129": {"C", "A"}, "163": {"F", "G"}, "164": {"F", "Y", "W", "L", "I", "P"}, \
#         "167": {"V", "I", "T", "C", "A", "S"}, "231": {"L", "R"}, "234": {"L", "F", "M"}, "262": {"G"}, \
#         "263": {"S", "A"}, "264": {"A", "G"}, "265": {"G", "A"}}
# # round1
# ylabel = "multi–state design fitness ($10^{5}$)"
# yticks = [0, 100000, 200000, 300000, 400000]
# ytick_labels = ["0", "1", "2", "3", "4"]
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000, 30000, 40000, 50000]
# ytick_labels = ["0", "1", "2", "3", "4", "5"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10, 15]
# ytick_labels = ["0", "5", "10", "15"]

# # round2
# ylabel = "multi–state design fitness ($10^{5}$)"
# yticks = [0, 100000, 200000, 300000, 400000]
# ytick_labels = ["0", "1", "2", "3", "4"]
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 10000, 20000]
# ytick_labels = ["0", "1", "2"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # round2 CF3
# ylabel = "multi–state design fitness ($10^{6}$)"
# yticks = [0, 1000000, 2000000]
# ytick_labels = ["0", "1", "2"]
# ylabel = "multi–state design fitness ($10^{4}$)"
# yticks = [0, 50000, 100000]
# ytick_labels = ["0", "5", "10"]
# ylabel = "positive design fitness"
# yticks = [0, 5, 10]
# ytick_labels = ["0", "5", "10"]

# # round3
# site_aas_list = {"129": {"C", "A"}, "167": {"V", "I", "T", "C", "A", "S"}, "231": {"L", "R"}}
# ylabel = "multi–state design fitness"
# ylabel = "positive design fitness"


font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

color_scheme = {
    'A': 'black',
    'C': 'black',
    'D': 'black',
    'E': 'black',
    'F': 'black',
    'G': 'black',
    'H': 'black',
    'I': 'black',
    'K': 'black',
    'L': 'black',
    'M': 'black',
    'N': 'black',
    'P': 'black',
    'Q': 'black', 
    'R': 'black', 
    'S': 'black', 
    'T': 'black', 
    'V': 'black', 
    'W': 'black', 
    'Y': 'black'
}

parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str)
parser.add_argument("--add_header", action="store_true")
parser.add_argument("-delim", "--delimiter", type=str, choices=["comma", "space"], default="comma")
parser.add_argument("-pos", "--masked_positions", type=str, nargs="*", help="e.g., -pos 52 168-174")
parser.add_argument("-xcd", "--excluded_positions", type=str, nargs="*", help="e.g., -xcd 52 168-174")
parser.add_argument("-xticks", "--xtick_labels", type=str, nargs="*")
parser.add_argument("-height", "--max_height", type=float, default=5)
args = parser.parse_args()

if args.add_header:
    if args.delimiter == "comma":
        delimiter=","
    elif args.delimiter == "space":
        delimiter=" "
    pssm = np.transpose(np.loadtxt(args.filename, delimiter=delimiter, dtype=float)[:20, :])
    aa_type = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    positions = np.arange(1, pssm.shape[0] + 1)
    df = pd.DataFrame(pssm, columns = aa_type, index = positions)
else:
    if args.delimiter == "comma":
        delim_whitespace=False
    elif args.delimiter == "space":
        delim_whitespace=True
    df = pd.read_csv(args.filename, delim_whitespace=delim_whitespace, index_col=0)

if args.masked_positions:
    positions = np.array([])
    for masked_position in args.masked_positions:
        if "-" in masked_position:
            lower_bound, upper_bound = masked_position.split("-")
            positions = np.append(positions, np.arange(int(lower_bound), int(upper_bound) + 1))
        else:
            positions = np.append(positions, masked_position)
    df = df.loc[positions]
    # df = df[df.index.isin(positions)]
if args.excluded_positions:
    excluded_positions = np.array([])
    for excluded_position in args.excluded_positions:
        if "-" in excluded_position:
            lower_bound, upper_bound = excluded_position.split("-")
            excluded_positions = np.append(excluded_positions, np.arange(int(lower_bound), int(upper_bound) + 1))
        else:
            excluded_positions = np.append(excluded_positions, excluded_position)
    df = df[~df.index.isin(excluded_positions)]

positions = df.index.values

npos = len(positions)
df.index = np.arange(1, npos + 1)

# matplotlib.rc('xtick', labelsize=20)
# matplotlib.rc('ytick', labelsize=20)

fig, ax = plt.subplots(1, 1, figsize=[npos, args.max_height])

logo = logomaker.Logo(df, ax=ax, color_scheme=color_scheme, \
        font_name='Arial Rounded MT Bold')

for i, site_aas in enumerate(site_aas_list.items()):
    for aa in color_scheme.keys() - site_aas[1]:
        logo.style_single_glyph(c=aa, p=i+1, color=[0.5, 0.5, 0.5])

ax.spines.top.set_visible(False)
ax.spines.right.set_visible(False)

logo.ax.set_xticks(np.arange(1, npos + 1))
if args.xtick_labels:
    logo.ax.set_xticklabels(args.xtick_labels)
else:
    logo.ax.set_xticklabels(positions)

if yticks:
    logo.ax.set_yticks(yticks, ytick_labels)

ax.set_ylabel(ylabel)

if title:
    ax.set_title(title, pad=30)

# logo.fig.show()

logo.fig.savefig(args.filename[:args.filename.rfind(".")] + ".png", bbox_inches='tight', dpi=300)
