import argparse
from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
import torch.nn.functional as F
from tqdm import tqdm
import warnings
import os
from dms_utils import deep_mutational_scan
import esm
import pandas as pd
from multichain_util import extract_coords_from_complex, _concatenate_coords, _concatenate_seqs, score_sequence_in_complex
from util import get_sequence_loss, load_structure, load_coords, score_sequence, extract_coords_from_structure

def get_native_seq(pdbfile, chain):
        structure = load_structure(pdbfile, chain)
        _ , native_seq = extract_coords_from_structure(structure)
        return native_seq

def score_singlechain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    coords, native_seq = load_coords(args.pdbfile, args.chain)
    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')
    ll, _ = score_sequence(
            model, alphabet, coords, native_seq)
    print('Native sequence')
    print(f'Log likelihood: {ll:.2f}')
    print(f'Perplexity: {np.exp(-ll):.2f}')
    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqpath)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood\n')
        for header, seq in tqdm(seqs.items()):
            ll, _ = score_sequence(
                    model, alphabet, coords, str(seq))
            fout.write(header + ',' + str(ll) + '\n')
    print(f'Results saved to {args.outpath}') 


def score_multichain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    structure = load_structure(args.pdbfile)
    coords, native_seqs = extract_coords_from_complex(structure)
    target_chain_id = args.chain
    native_seq = native_seqs[target_chain_id]
    order = args.order

    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')

    ll_complex, ll_targetchain = score_sequence_in_complex(
        model,
        alphabet,
        coords,
        native_seqs,
        target_chain_id,
        native_seq,
        order=order,
    )
    print('Native sequence')
    print(f'Log likelihood of complex: {ll_complex:.2f}')
    print(f'Log likelihood of target chain: {ll_targetchain:.2f}')
    print(f'Perplexity: {np.exp(ll_complex):.2f}')

    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqpath)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood, log_likelihood_target\n')
        for header, seq in tqdm(seqs.items()):
            ll_complex, ll_targetchain = score_sequence_in_complex(
                model,
                alphabet,
                coords,
                native_seqs,
                target_chain_id,
                str(seq),
                order=order,
            )
            fout.write(header + ',' + str(ll_complex) + ',' + str(ll_targetchain) + '\n')
    print(f'Results saved to {args.outpath}') 

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
            help='input filepath for variant sequences in a .fasta file',
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
            '--order', type=str, default=None,
            help='for multichain, specify order'
    )
    parser.add_argument(
            '--singlechain-backbone', dest='multichain_backbone',
            action='store_false',
            help='use the backbone of only target chain in the input for conditioning'
    )

    parser.add_argument(
        "--nogpu", action="store_true", 
        help="Do not use GPU even if available"
    )
    args = parser.parse_args()

    if args.outpath is None:
        args.outpath = f'output/{args.pdbfile[:-4]}-chain{args.chain}_scores.csv'

    model_checkpoint_path = get_model_checkpoint_path('esm_if1_20220410.pt')
    with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            model, alphabet = esm.pretrained.load_model_and_alphabet( \
                model_checkpoint_path \
            )
    model = model.eval()


    if args.multichain_backbone:
        score_multichain_backbone(model, alphabet, args)
    else:
        score_singlechain_backbone(model, alphabet, args)


if __name__ == '__main__':
    main()
