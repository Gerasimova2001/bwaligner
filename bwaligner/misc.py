import numpy as np
import random
import string
import os



def get_random_str(main_str, substr_len):
    idx = np.random.randint(0, len(main_str) - substr_len + 1)    # Randomly select an "idx" such that "idx + substr_len <= len(main_str)".
    return main_str[idx : (idx+substr_len)]



def generateGenome(length=10000, alphabet=['a', 'c', 'g', 't'], probs=None):
    if probs == None:
        probs = [1/len(alphabet)] * len(alphabet)
    else:
        assert len(alphabet) == len(probs), "Provide probabilities for each letter in alphabet"
    return Genome("".join(np.random.choice(alphabet, size=length, p=probs)))
 


def makeFastqFile(genome, n_reads=1000,read_length=100,base_path="", p_gap = 0.01, p_mut:float = 0.01, gap_length=10, debug=False):
    result = ""
    for i in range(n_reads):
        insert_occurance = np.random.choice([0, 1], p=[1-p_gap, p_gap])
        del_occurance = np.random.choice([0, 1], p=[1-p_gap, p_gap])
        if not del_occurance and not insert_occurance:
            read = get_random_str(genome, read_length)
        else:
            if del_occurance: 
                read = get_random_str(genome, read_length + gap_length)
                start = np.random.randint(0, read_length + 1)
                read = read[:start] + read[start+gap_length:]
                assert len(read) == read_length
                
                
            if insert_occurance: 
                read = get_random_str(genome, read_length - gap_length)
                start = np.random.randint(0, read_length + 1)
                read = read[:start] + generateGenome(gap_length) + read[start:]
        assert len(read) == read_length
        score = "".join(np.random.choice(list(string.printable[:94]), size=read_length))
        result += f"@SEQ_{i}\n"\
                + read + "\n"\
                + "+\n"\
                + score + "\n"
    with open(os.path.join(base_path, "experiment.fastq"), "w") as fastq_file:
        fastq_file.write(result)
        