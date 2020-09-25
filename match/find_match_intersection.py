#!/usr/bin/python3
import argparse
import glob, os, shutil

def load_matches(directories):
    state_match_dict = dict()
    for directory in directories:
        match_set = set()
        for variant in filter(lambda x: os.path.isfile(directory + '/' + x + '/' + x + '.pdb'), os.listdir(directory)):
            match_set.add(variant)
        state_match_dict.update({directory: match_set})
    return state_match_dict

def merge_states(state_match_dict):
    states = state_match_dict.keys()
    for state1 in states:
        scaffold_substrate1 = state1.split('_')
        scaffold1 = scaffold_substrate1[0]
        substrate1 = scaffold_substrate1[1]
        if substrate1.find('-') != substrate1.rfind('-'):
            substrate1_general = substrate1[:substrate1.rfind('-')]
            for state2 in filter(lambda x: x != state1, states):
                scaffold_substrate2 = state2.split('_')
                scaffold2 = scaffold_substrate2[0]
                substrate2 = scaffold_substrate2[1]
                if scaffold2 == scaffold1 and substrate2[:substrate2.rfind('-')] == substrate1_general:
                    state_match_dict[state1].update(state_match_dict[state2])

def intersection(state_match_dict):
    first_intersection = False
    for matches in state_match_dict.values():
        if first_intersection:
            match_intersection.intersection_update(matches)
        else:
            match_intersection = matches
            first_intersection = True
    for state in state_match_dict.keys():
        for match in os.listdir(state):
            if match.startswith('X') and match not in match_intersection:
                shutil.rmtree(state + '/' + match)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directories', type=str, nargs='*')
    args = parser.parse_args()
    state_matches = load_matches(args.directories)
    merge_states(state_matches)
    intersection(state_matches)
