#!/usr/bin/python3
import argparse
import glob, os, shutil

parser = argparse.ArgumentParser()
parser.add_argument('-or', '--or_directories', type=str, nargs='*')
parser.add_argument('-and', '--and_directories', type=str, nargs='*', required=True)
args = parser.parse_args()

if args.or_directories:
    match_intersection = set()
    for directory in args.or_directories:
        for variant in filter(lambda x: os.path.isfile(directory + '/' + x + '/' + x + '.pdb'), os.listdir(directory)):
            match_intersection.add(variant)
    flag = True
else:
    flag = False
for directory in args.and_directories:
    match_set = set()
    for variant in filter(lambda x: os.path.isfile(directory + '/' + x + '/' + x + '.pdb'), os.listdir(directory)):
        match_set.add(variant)
    if flag:
        match_intersection = match_intersection.intersection(match_set)
    else:
        match_intersection = match_set
        flag = True

for directory in args.or_directories:
    for variant in filter(lambda x: os.path.isfile(directory + '/' + x + '/' + x + '.pdb'), os.listdir(directory)):
        if variant not in match_intersection:
            shutil.rmtree(directory + '/' + variant)

for directory in args.and_directories:
    for variant in filter(lambda x: os.path.isfile(directory + '/' + x + '/' + x + '.pdb'), os.listdir(directory)):
        if variant not in match_intersection:
            shutil.rmtree(directory + '/' + variant)
