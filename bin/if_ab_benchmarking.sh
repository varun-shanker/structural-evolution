abs=("cr6261" "cr9114" "g6" "g6")
ab_fastas=("cr6261_3gbn_hc_lib.fasta" "cr9114_4fqi_hc_lib.fasta" "g6_2fjg_hc_lib.fasta" "g6_2fjg_lc_lib.fasta")
data_path="data/ab_mutagenesis_expts/"
out_prefix="output/ab_mutagenesis_expts/"

for ((i=0; i<${#abs[@]}; i++)); do
    ab="${abs[i]}"
    ab_fasta="${ab_fastas[i]}"
    ab_dir_path="${data_path}${ab}/"
    struc_list=("${ab_dir_path}"*.pdb)
    ab_out_dir="${out_prefix}${ab}/"

    # Set the default chain value
    chain="H"
    # Special handling for 'g6' antibody
    if [[ "$ab" == "g6" ]]; then
        [[ "$ab_fasta" == *"lc"* ]] && chain="L"
    fi

    # gather pdbs and filter the hc/lc only structure from being scored by library for the other chain
    if [[ "$ab" == "g6" && "$chain" == "L" ]]; then
        struc_list=($(echo "${struc_list[@]}" | tr ' ' '\n' | grep -v '_h_' | tr '\n' ' '))
    elif [[ "$ab" == "g6" && "$chain" == "H" ]]; then
        struc_list=($(echo "${struc_list[@]}" | tr ' ' '\n' | grep -v '_l_' | tr '\n' ' '))
    fi
    
    mkdir -p "$ab_out_dir"

    for struc in "${struc_list[@]}"; do
        out_file="${ab_out_dir}${struc##*/}"
        if [[ "$ab" == "g6" ]]; then
            chain_modeled="$([ "$chain" == "H" ] && echo "hc" || echo "lc")"
            out_file="${out_file%_fvar.pdb}_${chain_modeled}_scores.csv"
        else
            out_file="${out_file%_fvar.pdb}_scores.csv"
        fi

        if [[ ! -f "$out_file" ]]; then
            python bin/score_log_likelihoods.py "$struc" --chain "$chain" --seqpath "${ab_dir_path}${ab_fasta}" --outpath "$out_file"
        else
            echo "$out_file already exists. Skipping..."
        fi
    done
done