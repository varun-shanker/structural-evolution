#!/bin/bash

run_protein_mpnn() {
    python ../ProteinMPNN/protein_mpnn_run.py \
        --path_to_fasta "$1" \
        --score_only 1 \
        --save_score 1 \
        --pdb_path_chains "$2" \
        --out_folder "$3" \
        --pdb_path "$4"
}

# G6 Light Chain
run_protein_mpnn \
    "data/ab_mutagenesis_expts/g6/g6_2fjg_lc_lib.fasta" \
    "L" \
    "output/ab_mutagenesis_expts/g6/proteinMpnnLC" \
    "data/ab_mutagenesis_expts/g6/2fjg_vlh_fvar.pdb"

# G6 Heavy Chain
run_protein_mpnn \
    "data/ab_mutagenesis_expts/g6/g6_2fjg_hc_lib.fasta" \
    "H" \
    "output/ab_mutagenesis_expts/g6/proteinMpnnHC" \
    "data/ab_mutagenesis_expts/g6/2fjg_vlh_fvar.pdb"

# CR6261 Heavy Chain
run_protein_mpnn \
    "data/ab_mutagenesis_expts/cr6261/cr6261_3gbn_hc_lib.fasta" \
    "H" \
    "output/ab_mutagenesis_expts/cr6261/mpnn" \
    "data/ab_mutagenesis_expts/cr6261/3gbn_ablh_fvar.pdb"

# CR9114 Heavy Chain
run_protein_mpnn \
    "data/ab_mutagenesis_expts/cr9114/cr9114_4fqi_hc_lib.fasta" \
    "H" \
    "output/ab_mutagenesis_expts/cr9114/mpnn" \
    "data/ab_mutagenesis_expts/cr9114/4fqi_ablh_fvar.pdb"