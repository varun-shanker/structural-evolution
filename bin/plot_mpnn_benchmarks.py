import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.palettes 
from bokeh.transform import factor_cmap
import datashader
import holoviews as hv
import holoviews.operation.datashader
import os
from natsort import natsorted

hv.extension("bokeh")

import warnings
# Suppress FutureWarning messages
warnings.simplefilter(action='ignore',)

cr9114_dict = {
                'ab_name'   : 'CR9114',
                'files'     : {
                                'ESM-IF1':      'output/ab_mutagenesis_expts/cr9114/4fqi_ablh_scores.csv',
                                'ProteinMPNN': 'output/ab_mutagenesis_expts/cr9114/mpnn/score_only/',

                                },
                'dms_df'    : pd.read_csv('data/ab_mutagenesis_expts/cr9114/cr9114_exp_data.csv', 
                                dtype = {'genotype': str},
                                ).rename(columns={'h1_mean': 'H1', 'h3_mean' : 'H3'}),
                'ag_columns': ['H1', 'H3'],
                'expt_type':  'Combinatorial Mutagenesis for Affinity',
                'palette':     bokeh.palettes.Spectral6
}

cr6261_dict = {
                'ab_name'   : 'CR6261',
                'files'     : {
                                'ESM-IF1': 'output/ab_mutagenesis_expts/cr6261/3gbn_ablh_scores.csv',
                                'ProteinMPNN': 'output/ab_mutagenesis_expts/cr6261/mpnn/score_only/'
                                },
                'dms_df'    : pd.read_csv('data/ab_mutagenesis_expts/cr6261/cr6261_exp_data.csv',
                            dtype={'genotype': str},
                            ).rename(columns={'h1_mean': 'H1', 'h9_mean' : 'H9'}),
                'ag_columns': ['H1', 'H9'],
                'expt_type':  'Combinatorial Mutagenesis for Affinity',
                'palette':     bokeh.palettes.Dark2_6,
}


g6LC_dict = {
                'ab_name'   : 'g6',
                'files'     : {
                                'ESM-IF1':'output/ab_mutagenesis_expts/g6/2fjg_vlh_lc_scores.csv',
                                'ProteinMPNN':'output/ab_mutagenesis_expts/g6/proteinMpnnLC/score_only/',
                                
                                },
                'dms_df'    : pd.read_csv('data/ab_mutagenesis_expts/g6/g6_lc_exp_data.csv'),
                'ag_columns': ['norm_binding'],
                'expt_type':  'Deep Mutational Scan for Binding',
                'palette':     bokeh.palettes.Pastel1_6,
                'chain':      'VL'
}

g6HC_dict = {
                'ab_name'   : 'g6',
                'files'     : {
                                'ESM-IF1': 'output/ab_mutagenesis_expts/g6/2fjg_vlh_hc_scores.csv',
                                'ProteinMPNN':'output/ab_mutagenesis_expts/g6/proteinMpnnHC/score_only/',
                                },
                'dms_df'    : pd.read_csv('data/ab_mutagenesis_expts/g6/g6_hc_exp_data.csv'),
                'ag_columns': ['norm_binding'],
                'expt_type':  'Deep Mutational Scan for Binding',
                'palette':     bokeh.palettes.Pastel1_4,
                'chain':      'VH'
}
def combine_mpnn_scores(dir):
    global_scores = []
    npz_files = [f for f in os.listdir(dir) if f.endswith('.npz')]
    npz_files = natsorted(npz_files)[:-1]

    for fname in npz_files:
        file_path = os.path.join(dir, fname)
        data = np.load(file_path)
        global_scores.append(-1* data['global_score'][0])
        data.close()


    return global_scores


#get melted correlations matrix of structure input x target antigen, for a given antibody 
def get_corr(ab_name, files, dms_df, ag_columns,):

    #retreive correlations from inverse fold scores
    for key, filepath in files.items():    
        column_label = key
        if key == 'ProteinMPNN':
            scores = combine_mpnn_scores(filepath)
            dms_df[column_label] = scores
        elif key == 'ESM-IF1':
            df = pd.read_csv(files['ESM-IF1'])  # Read the CSV file into a DataFrame
            log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
            dms_df[column_label] = log_likelihood_column


    conditions = list(files.keys())

    correlations = dms_df[ag_columns + conditions].corr(method='spearman')
    correlations = correlations.drop(conditions, axis= 1)
    correlations = correlations.drop(ag_columns, axis= 0)

    # Melt the correlations DataFrame and rename the index column
    melted_correlations = (
        correlations
        .reset_index()
        .melt(id_vars='index', var_name='Target', value_name='Correlation')
        .rename(columns={'index': 'Input'})
    )

    print(melted_correlations)
    return melted_correlations
    
    return all_ag_melted
    
#get melted correlations matrix of structure input x target antigen, for g6 antibody 
def get_g6_corr(g6Hc_dict, g6LC_dict, dropLOQ = False ):

    ab_name = 'g6'
    ag_columns = g6LC_dict['ag_columns']
    vh_and_vl_dms_df = pd.DataFrame({})

    for d in [g6HC_dict, g6LC_dict]:

        files = d['files']
        dms_df = d['dms_df']

        #retreive correlations from abysis, ablang, and esm1v scores
        for key, filepath in files.items():        
            column_label = key
            if key == 'ProteinMPNN':
                scores = combine_mpnn_scores(filepath)
                dms_df[column_label] = scores
            elif key == 'ESM-IF1':
                df = pd.read_csv(files['ESM-IF1'])  # Read the CSV file into a DataFrame
                log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
                dms_df[column_label] = log_likelihood_column


        vh_and_vl_dms_df = pd.concat([vh_and_vl_dms_df, dms_df], ignore_index= True)

    conditions = list(files.keys())

    correlations = vh_and_vl_dms_df[ag_columns + conditions].corr(method='spearman')
    correlations = correlations.drop(conditions, axis= 1)
    correlations = correlations.drop(ag_columns, axis= 0)

    # Melt the correlations DataFrame and rename the index column
    melted_correlations = (
        correlations
        .reset_index()
        .melt(id_vars='index', var_name='Target', value_name='Correlation')
        .rename(columns={'index': 'Input'})
    )
    melted_correlations['Target'] = ['VEGF-A'] * len(melted_correlations)

    print(melted_correlations)
    return melted_correlations

    
#plot correlation bars for cr antibodies: cr6261  and cr9114
def plot_hbar(title, melted_correlations, palette = bokeh.palettes.Spectral6, ):

    melted_correlations['cats'] = melted_correlations.apply(lambda x: (x["Target"], x["Input"]), axis = 1)
    factors = list(melted_correlations.cats)[::-1]
    
    p = bokeh.plotting.figure(
    height=340,
    width=440,
    x_axis_label="Spearman Correlation",
    x_range=[0, 1],
    y_range=bokeh.models.FactorRange(*factors),
    tools="save",
    title = title
    )
        

    p.hbar(
        source=melted_correlations,
        y="cats",
        right="Correlation",
        height=0.6,
        line_color = 'black',
        fill_color=bokeh.palettes.Dark2_6[0]
             
    )

    labels_df = melted_correlations
    labels_df['corr_str'] = labels_df['Correlation'].apply(lambda x: round(x, 2)).astype(str)
    labels_source = bokeh.models.ColumnDataSource(labels_df)

    labels = bokeh.models.LabelSet(x='Correlation', y='cats', text='corr_str',text_font_size = "10px",
              x_offset=12, y_offset=-5, source=labels_source, render_mode='canvas')

    p.ygrid.grid_line_color = None
    p.y_range.range_padding = 0.1
    p.add_layout(labels)
    p.legend.visible = False

    p.output_backend = "svg"
    return p

if __name__ == '__main__':
    datasets = [cr9114_dict, cr6261_dict, (g6HC_dict, g6LC_dict) ]

    all_corr_plots = []

    for d in datasets:
        if type(d) is tuple:
            g6Hc, g6Lc = d
            title = g6Hc['ab_name'] + ', ' + g6Hc['expt_type']
            g6_combined = get_g6_corr(g6Hc, g6Lc)
        
            all_corr_plots.append(plot_hbar(title, g6_combined, g6Lc['palette']))

        else:
            corr_df = get_corr(d['ab_name'], d['files'], d['dms_df'], d['ag_columns'])
            title = d['ab_name'] + ', ' + d['expt_type']
            all_corr_plots.append(plot_hbar(title, corr_df, d['palette']))


    all_corr_fname = f"output/ab_mutagenesis_expts/mpnn_benchmarks.html"
    bokeh.plotting.output_file(all_corr_fname)
    bokeh.io.show(bokeh.layouts.gridplot(all_corr_plots, ncols = len(all_corr_plots)))
