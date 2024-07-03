import argparse
from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
import torch.nn.functional as F
from tqdm import tqdm
import os

import esm

from dms_utils import deep_mutational_scan
import warnings
import pandas as pd
import util
import multichain_util

import score_log_likelihoods

def get_native_seq(pdbfile, chain):
        structure = util.load_structure(pdbfile, chain)
        _ , native_seq = util.extract_coords_from_structure(structure)
        return native_seq
    
def write_dms_lib(args):
    '''Writes a deep mutational scanning library, including the native/wildtype (wt) of the 
        indicated target chain in the structure to an output Fasta file'''

    sequence = get_native_seq(args.pdbfile, args.chain)
    Path(args.seqpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.seqpath, 'w') as f:
        f.write('>wt\n')
        f.write(sequence+'\n')
        for pos, wt, mt in deep_mutational_scan(sequence):
            assert(sequence[pos] == wt)
            mut_seq = sequence[:pos] + mt + sequence[(pos + 1):]
            f.write('>' + str(wt) + str(pos+1+args.offset) + str(mt) + '\n')
            f.write(mut_seq + '\n')

def get_top_n(args):
    recs, rec_inds = [], []
    scores_df = pd.read_csv(args.outpath).sort_values(by = 'log_likelihood', ascending = False)

    for seqid in scores_df['seqid']:
        res_ind = seqid[1:-1]
        if (rec_inds.count(res_ind) < args.maxrep): 
            if args.upperbound == None or (int(res_ind) < int(args.upperbound)):
                recs.append(seqid)
                rec_inds.append(res_ind)
        if len(recs) == args.n:
            break 

    print(f'\n Chain {args.chain}')
    print(*recs, sep='\n')

def get_model_checkpoint_path(filename):
    # Expanding the user's home directory
    return os.path.expanduser(f"~/.cache/torch/hub/checkpoints/{filename}")
        
def main():
    parser = argparse.ArgumentParser(
        description='Score sequences based on a given structure.'
    )
    parser.add_argument(
        'pdbfile', type=str,
        help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
        '--seqpath', type=str,
        help='filepath where fasta of dms library should be saveda',
    )
    parser.add_argument(
        '--outpath', type=str,
        help='output filepath for scores of variant sequences',
    )
    parser.add_argument(
        '--chain', type=str,
        help='chain id for the chain of interest', default='A',
    )
    parser.set_defaults(multichain_backbone=True)
    parser.add_argument(
        '--multichain-backbone', action='store_true',
        help='use the backbones of all chains in the input for conditioning'
    )
    parser.add_argument(
        '--singlechain-backbone', dest='multichain_backbone',
        action='store_false',
        help='use the backbone of only target chain in the input for conditioning'
    )
    parser.add_argument(
            '--order', type=str, default=None,
            help='for multichain, option to specify order of chains'
    )
    parser.add_argument(
        '--n', type=int,
        help='number of desired predictions to be output', 
        default=10,
    )
    parser.add_argument(
        '--maxrep', type=int,
        help='maximum representation of a single site in the top recommendations  \
              (eg: maxrep = 1 is a unique set where no wildtype residue is mutated more than once)', 
        default=1,
    )
    parser.add_argument(
        '--offset', type=int,
        help='integer offset for labeling of residue indices encoded in the structure',
        default=0,
    )
    parser.add_argument(
        '--upperbound', type=int,
        help='only residue positions less than the user-defined upperbound are considered to be recommended for screening \
                (but all positions are still conditioned for scoring)', 
        default=None,
    )
    parser.add_argument(
        "--nogpu", action="store_true", 
        help="Do not use GPU even if available"
    )

    args = parser.parse_args()

    if args.seqpath is None:
        args.seqpath = f'output/{args.pdbfile[:-4]}-chain{args.chain}_dms.fasta'

    if args.outpath is None:
        args.outpath = f'output/{args.pdbfile[:-4]}-chain{args.chain}_scores.csv'

    #write dms library for target chain 
    write_dms_lib(args)

    model_checkpoint_path = get_model_checkpoint_path('esm_if1_20220410.pt')
    with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            model, alphabet = esm.pretrained.load_model_and_alphabet( \
                model_checkpoint_path \
            )
    model = model.eval()
    
    if args.multichain_backbone:
        score_log_likelihoods.score_multichain_backbone(model, alphabet, args)
    else:
        score_log_likelihoods.score_singlechain_backbone(model, alphabet, args)

    get_top_n(args)


if __name__ == '__main__':
        main()