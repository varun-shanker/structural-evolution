import json
import pandas as pd
import os

AAs = set('ACDEFGHIKLMNPQRSTVWY')

seqs_abs = {
    'cr6261_vh': 'EVQLVESGAEVKKPGSSVKVSCKASGGPFRSYAISWVRQAPGQGPEWMGGIIPIFGTTKYAPKFQGRVTITADDFAGTVYMELSSLRSEDTAMYYCAKHMGYQVRETMDVWGKGTTVTVSS',
    'cr9114_vh': 'QVQLVQSGAEVKKPGSSVKVSCKSSGGTSNNYAISWVRQAPGQGLDWMGGISPIFGSTAYAQKFQGRVTISADIFSNTAYMELNSLTSEDTAVYFCARHGNYYYYSGMDVWGQGTTVTVSS',
    'g6_vh':     'EVQLVESGGGLVQPGGSLRLSCAASGFTISDYWIHWVRQAPGKGLEWVAGITPAGGYTYYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARFVFFLPYAMDYWGQGTLVTV',  
    'g6_vl':     'DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYTTPPTFGQGTKVEIK',
    'beb_vh':   'QITLKESGPTLVKPTQTLTLTCTFSGFSLSISGVGVGWLRQPPGKALEWLALIYWDDDKRYSPSLKSRLTISKDTSKNQVVLKMTNIDPVDTATYYCAHHSISTIFDHWGQGTLVTVSS',
    'beb_vl' :   'QSALTQPASVSGSPGQSITISCTATSSDVGDYNYVSWYQQHPGKAPKLMIFEVSDRPSGISNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTTSSAVFGGGTKLTVL',
    'sa58_vh':  'QVQLAQSGSELRKPGASVKVSCDTSGHSFTSNAIHWVRQAPGQGLEWMGWINTDTGTPTYAQGFTGRFVFSLDTSARTAYLQISSLKADDTAVFYCARERDYSDYFFDYWGQGTLVTVSS',
    'sa58_vl':   'EVVMTQSPASLSVSPGERATLSCRARASLGISTDLAWYQQRPGQAPRLLIYGASTRATGIPARFSGSGSGTEFTLTISSLQSEDSAVYYCQQYSNWPLTFGGGTKVEIK',
}


def parse_abysis(seq, seq_name):
    fname = f'data/abysis/abysis_counts_{seq_name}.json'
    with open(fname) as f:
        json_data = json.load(f)

    data = []
    for idx, freq_table in enumerate(json_data['frequencies']):
        wt = seq[idx]
        total = freq_table['total']
        for entry in freq_table['counts']:
            if entry['aa'] == seq[idx]:
                wt_frac = entry['c'] / total
        aa_to_freq, seen = {}, set()
        for entry in freq_table['counts']:
            mt = entry['aa']
            seen.add(mt)
            counts = entry['c']
            frac = entry['c'] / total
            ratio = frac / wt_frac
            data.append([ idx + 1, wt, mt, counts, frac, ratio ])
        for mt in AAs - seen:
            data.append([ idx + 1, wt, mt, 0, 0., 0. ])

    df = pd.DataFrame(data, columns=[
        'pos',
        'wt',
        'mt',
        'counts',
        'fraction',
        'likelihood_ratio',
    ])

    if any(prefix in seq_name for prefix in ['cr6261', 'cr9114', 'g6']):
        subdirectory = seq_name.split('_')[0]
    else:
        subdirectory = ''

    output_dir = os.path.join('output', 'ab_mutagenesis_expts', subdirectory)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f'abysis_counts_{seq_name}.txt')
    df.to_csv(output_path, sep='\t')

if __name__ == '__main__':
    for seq_name in seqs_abs:
        seq = seqs_abs[seq_name]
        parse_abysis(seq, seq_name)