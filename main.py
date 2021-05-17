#!/usr/bin/env python3

from bwaligner import *
import sys

import json

debug = False
show_data_structures = False
use_lower_bound_tree_pruning = True #set this to false (in conjunction with debug=True) to see the full search through the suffix trie
#search parameters
indels_allowed = True # turn off for mismatches only, no insertion or deletions allowed
difference_threshold = 1
insertion_penalty = 1
deletion_penalty = 1
mismatch_penalty = 1


if __name__ == "__main__":
    #index the reference and build up all the necessary datastructures (5 of them; BWT, SA, reference alphabet, C and OCC arrays)
    referencePath = sys.argv[1]
    experimentPath = sys.argv[2]
    testPath = sys.argv[3]
    with open(referencePath, "r") as f:

    	genome = Genome.from_file(f)

    with open(testPath) as json_file:
        test_dict = json.load(json_file)

    aligner = BWAligner(genome, debug=False)
    with open(experimentPath, "r") as f:

    	experiment = FastqExperiment.from_file(f)

    res = aligner.align(experiment, difference_threshold=0)
    acc = []
    for key, value in test_dict.items():
        acc.append(res[key] == value)
    print (sum(acc)/len(acc))