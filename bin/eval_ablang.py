import ablang
import numpy as np
import pandas as pd
import scipy.special

heavy_ablang = ablang.pretrained("heavy")
heavy_ablang.freeze()

light_ablang = ablang.pretrained("light")
light_ablang.freeze()

ab_dict = dict(
    cr6261= ('EVQLVESGAEVKKPGSSVKVSCKASGGPFRSYAISWVRQAPGQGPEWMGGIIPIFGTTKYAPKFQGRVTITADDFAGTVYMELSSLRSEDTAMYYCAKHMGYQVRETMDVWGKGTTVTVSS', 
                'QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNDYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEANYYCATWDRRPTAYVVFGGGTKLTVL'), 
    cr9114=  ('QVQLVQSGAEVKKPGSSVKVSCKSSGGTSNNYAISWVRQAPGQGLDWMGGISPIFGSTAYAQKFQGRVTISADIFSNTAYMELNSLTSEDTAVYFCARHGNYYYYSGMDVWGQGTTVTVSS',
                'QSALTQPPAVSGTPGQRVTISCSGSDSNIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLKGAVFGGGTQLTVL'),
    g6 =  ('EVQLVESGGGLVQPGGSLRLSCAASGFTISDYWIHWVRQAPGKGLEWVAGITPAGGYTYYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARFVFFLPYAMDYWGQGTLVTV',
                'DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYTTPPTFGQGTKVEIK')

)

def eval_ablang(s, ab, chain):
    if chain == 'hc':
        fname = 'output/ab_mutagenesis_expts/'+ab+'/' + ab +'_hc_ablangScores.csv'
        log_likelihoods = heavy_ablang(s, mode = 'likelihood')[0][1:-1]
        alphabet = heavy_ablang.tokenizer.vocab_to_aa
    elif chain == 'lc':
        fname = 'output/ab_mutagenesis_expts/'+ab+'/' + ab +'_lc_ablangScores.csv'
        log_likelihoods = light_ablang(s, mode = 'likelihood')[0][1:-1]
        alphabet = light_ablang.tokenizer.vocab_to_aa

    assert (log_likelihoods.shape)[0] == len(s)

    filt_alphabet = {key: value for key, value in alphabet.items() if value.isalpha()}
    log_likelihood_ratio = []
    for i,res_log_likelihoods in enumerate(log_likelihoods):
            wt_res = s[i]
            wt_index = list(filt_alphabet.values()).index(wt_res)
            wt_log_likelihood = res_log_likelihoods[wt_index]
            log_likelihood_ratio.extend(res_log_likelihoods-wt_log_likelihood)

    res_order = [alphabet[key] for key in range(1, 21)] #extract order of residues in likelihood
    mt = res_order * len(s)
    wt = [char for char in s for i in range(len(res_order))]
    pos = [i+1 for i in range(len(s)) for j in range(len(res_order))]
    data = {
        'pos': pos,
        'wt': wt,
        'mt': mt,
        'log_likelihood' : log_likelihoods.flatten(),
        'log_likelihood_ratio' : log_likelihood_ratio,
    } 
    df = pd.DataFrame(data)
    df.to_csv(fname, index = False)


def main():
    
    for ab in ab_dict:
        for s,chain in zip(ab_dict[ab], ('hc','lc')):
              eval_ablang(s, ab, chain)
           

if __name__ == '__main__':
        main()