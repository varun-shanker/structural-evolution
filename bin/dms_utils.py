from collections import Counter
import datetime
from dateutil.parser import parse as dparse
import errno
import numpy as np
import os
import random
import sys
import time
import warnings
from Bio import pairwise2
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio import Seq, SeqIO


np.random.seed(1)
random.seed(1)

AAs = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
]

def tprint(string):
    string = str(string)
    sys.stdout.write(str(datetime.datetime.now()) + ' | ')
    sys.stdout.write(string + '\n')
    sys.stdout.flush()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
        
def deep_mutational_scan(sequence, exclude_noop=True):
    for pos, wt in enumerate(sequence):
        for mt in AAs:
            if exclude_noop and wt == mt:
                continue
            yield (pos, wt, mt)


def make_mutations(seq, mutations):
    mut_seq = [ char for char in seq ]
    for mutation in mutations:
        wt, pos, mt = mutation[0], int(mutation[1:-1]) - 1, mutation[-1]
        assert(seq[pos] == wt)
        mut_seq[pos] = mt
    mut_seq = ''.join(mut_seq).replace('-', '')
    return mut_seq

def find_mutations(seq1, seq2):
    alignment = pairwise2.align.globalms(
        seq1, seq2, 5, -4, -6, -.1, one_alignment_only=True,
    )[0]

    mutation_set = []
    pos1, pos2, pos_map = 0, 0, {}
    for wt, mt in zip(alignment[0], alignment[1]):
        if wt != '-':
            pos1 += 1
        if mt != '-':
            pos2 += 1
        if wt != mt and wt != '-' and mt != '-':
            mut_str = f'{wt}{pos1}{mt}'
            mutation_set.append(mut_str)

    return mutation_set