import argparse
from dms_utils import deep_mutational_scan
from pathlib import Path
import numpy as np

import esm
from util import load_structure, extract_coords_from_structure
import biotite.structure
from collections import defaultdict


def get_native_seq(pdbfile, chain):
    structure = load_structure(pdbfile, chain)
    _ , native_seq = extract_coords_from_structure(structure)
    return native_seq
    
def write_dms_lib(args):
    '''Writes a deep mutational scanning library, including the native/wildtype (wt) of the 
    indicated target chain in the structure to an output Fasta file'''

    sequence = get_native_seq(args.pdbfile, args.chain)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.dmspath, 'w') as f:
        f.write('>wt\n')
        f.write(sequence+'\n')
        for pos, wt, mt in deep_mutational_scan(sequence):
            assert(sequence[pos] == wt)
            mut_seq = sequence[:pos] + mt + sequence[(pos + 1):]

            f.write('>' + str(wt) + str(pos+1) + str(mt) + '\n')
            f.write(mut_seq + '\n')


def main():
    parser = argparse.ArgumentParser(
            description='Create a DMS library based on target chain in the structure.'
    )
    parser.add_argument(
            'pdbfile', type=str,
            help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
            '--dmspath', type=str,
            help='output filepath for dms library',
    )
    parser.add_argument(
            '--chain', type=str,
            help='chain id for the chain of interest', default='A',
    )

    args = parser.parse_args()

    if args.dmspath is None:
        args.dmspath = f'predictions/{args.pdbfile[:-4]}-{args.chain}_dms.fasta'

    write_dms_lib(args)

if __name__ == '__main__':
    
    main()
