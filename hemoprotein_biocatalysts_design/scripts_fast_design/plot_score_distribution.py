import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import xlrd, xlwt
from xlutils.copy import copy


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('xls', type=str)
    parser.add_argument('-total', action='store_true')
    return parser.parse_args()

def read_scores(sheet, cal_total):
    scores = np.empty(0)
    if cal_total:
        bins = np.arange(-60, 60, 10)
    else:
        bins = np.arange(-3, 5, 1)
    i = 1 if cal_total else 0
    for col in range(2, 18, 2):
        scores_float = list(filter(lambda x: x != '', sheet.col_values(col + i)))
        scores = np.insert(scores, 0, np.array(scores_float).astype(np.float))
    # score_distribution = np.digitize(scores, bins)
    return bins, scores

def plot(bins, score_distribution, protein):
    plt.hist(score_distribution, bins=bins)
    plt.xlabel('scores')
    plt.ylabel('frequency')
    plt.title(protein)
    # plt.show()
    fig = plt.gcf()
    fig.savefig(protein + '.png')

if __name__ == '__main__':
    args = parse_arguments()
    workbook = xlrd.open_workbook(args.xls)
    sheet = workbook.sheet_by_name('design')
    bins, scores = read_scores(sheet, args.total)
    plot(bins, scores, args.xls[:-4])
