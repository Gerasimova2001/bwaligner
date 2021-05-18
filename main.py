#!/usr/bin/env python3
import sys
import json
import argparse

from functools import reduce

from bwaligner import *

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
    parser = argparse.ArgumentParser(description='BWAligner - a framework for assembling genomes and aligning reads to reference genome')
    parser.add_argument('--test', action="store_true",  default = False, help='Test how the program works using simulated data')
    parser.add_argument("--align",  action="store_true", default = False, help = "Use aligner mode")
    parser.add_argument("--assemble",  action="store_true", default = False, help = "Use genome assemply mode")
    parser.add_argument("-g", "--genomePath",   type=str, help = "Path to reference genome (required for aligner mode)")
    parser.add_argument("-f", "--fastq", type=str, help = "Path to fastq file")
    parser.add_argument("-t", "--threshold", type=int, default = 0, help="Threshold that restricts the difference between a read and a substring in the reference genome")
    parser.add_argument("-k", "--kmer", type=int, default = 12, help = "`length of sequences in the De Brujin graph")

    args = parser.parse_args()

    if args.test:
        if args.align:
            print("Starting Aligner test...")
            print("Generating genome...")
            genome = generateGenome()
            with open("genome.txt", "w") as f:
                f.write(genome.data)
            print("Generating Fastq file...")
            test_dict = makeFastqFile(genome, debug=True)
            with open("experiment.fastq", "r") as f:

                    experiment = FastqExperiment.from_file(f)

            aligner = BWAligner(genome, debug=False)
            aligner.align(experiment, difference_threshold=0)
            res  = 0
            for key, value in test_dict.items():
                if test_dict[key] == value:
                    res += 1

            print(f"Accuracy: {res/len(test_dict)}")


        if args.assemble:
            print("Starting assembly test...")
            assembly_test_res = test_assembly((200, 220), n_iter=5, read_length=20, overlap=10, k=5)
            for key, value in assembly_test_res.items():
                print(f"{key}: mean levenstein distance: {value[0]}, mean accuracy: {value[1]}")


    
    else:
        if args.align and args.assemble:
            print("Can only align or only assemble, exiting")
            sys.exit(1)
        if args.align:
            if args.genomePath is None or args.fastq is None:
                print("Please provide both path to genome and path to fastq file")
                sys.exit(1)

            else:
                referencePath = args.genomePath
                experimentPath = args.fastq
                with open(referencePath, "r") as f:

                    genome = Genome.from_file(f)

                with open(experimentPath, "r") as f:

                    experiment = FastqExperiment.from_file(f)

                aligner = BWAligner(genome, debug=False)
                res = aligner.align(experiment, difference_threshold=args.threshold)
                for key, value in res.items():
                    if key != "start" and key != "stop":
                        print(f"{key}\t{','.join(list(map(lambda x: str(x), value)))}")


        elif args.assemble:
            if args.fastq is None:
                print("Please provide path to the fastq file")
                sys.exit(1)
            experimentPath = args.fastq
            with open(experimentPath, "r") as f:

                    experiment = FastqExperiment.from_file(f)

            genome = BuildGenome(experiment, args.k)
            print(genome)



    


