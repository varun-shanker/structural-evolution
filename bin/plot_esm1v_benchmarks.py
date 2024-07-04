import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.palettes 
from bokeh.transform import factor_cmap
import datashader
import holoviews as hv
import holoviews.operation.datashader
hv.extension("bokeh")

import warnings
# Suppress FutureWarning messages
warnings.simplefilter(action='ignore',)

cr9114_dict = {
                'ab_name'   : 'CR9114',
                'files'     : {
                                'ESM-1v Ab-Ag':      'output/ab_mutagenesis_expts/cr9114/cr9114_esm1vAbAg_exp_data_maskMargLabeled.csv',
                                'ESM-1v Ab only':    'output/ab_mutagenesis_expts/cr9114/cr9114_esm1vbothchains_exp_data_maskMargLabeled.csv',
                                'ESM-1v Ab VH only': 'output/ab_mutagenesis_expts/cr9114/cr9114_exp_data_maskMargLabeled.csv',
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
                                'ESM-1v Ab-Ag':      'output/ab_mutagenesis_expts/cr6261/cr6261_esm1vAbAg_exp_data_maskMargLabeled.csv',
                                'ESM-1v Ab only':    'output/ab_mutagenesis_expts/cr6261/cr6261_esm1vbothchains_exp_data_maskMargLabeled.csv',
                                'ESM-1v Ab VH only': 'output/ab_mutagenesis_expts/cr6261/cr6261_exp_data_maskMargLabeled.csv',
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
                                'LM Ab-Ag':      'output/ab_mutagenesis_expts/g6/g6Lc_esm1vAbAg_exp_data_maskMargLabeled.csv',
                                'LM Ab only':    'output/ab_mutagenesis_expts/g6/g6Lc_esm1vbothchains_exp_data_maskMargLabeled.csv',
                                'LM Ab VH/VL only': 'output/ab_mutagenesis_expts/g6/g6Lc_exp_data_maskMargLabeled.csv',
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
                                'LM Ab-Ag':      'output/ab_mutagenesis_expts/g6/g6Hc_esm1vAbAg_exp_data_maskMargLabeled.csv',
                                'LM Ab only':    'output/ab_mutagenesis_expts/g6/g6Hc_esm1vbothchains_exp_data_maskMargLabeled.csv',
                                'LM Ab VH/VL only': 'output/ab_mutagenesis_expts/g6/g6Hc_exp_data_maskMargLabeled.csv',
                                },
                'dms_df'    : pd.read_csv('data/ab_mutagenesis_expts/g6/g6_hc_exp_data.csv'),
                'ag_columns': ['norm_binding'],
                'expt_type':  'Deep Mutational Scan for Binding',
                'palette':     bokeh.palettes.Pastel1_4,
                'chain':      'VH'
}

def apply_mask_and_average(scores, genotype):
    #scores here should already be a in log space

    if '1' in genotype:
        masked_list = []
        for s, mask_char in zip(scores, genotype):
            if mask_char == '1':
                masked_list.append(s)
        multi_avg = np.mean(masked_list)
    else:
        #if genotype is all zeros = wt
        multi_avg = 0
    
    return multi_avg


def transform_single_to_multi( dms_df, singleMuts_df, condition, sort_col, ascending ):
    multi_scores = []
    #use input parameter as method bc key is kwarg in pd.sort_values
    #method should be either esm1v or abysis or AbLang
    
    #for cr sorts the single mutation dataframe in residue order, ie binary '10000' before '00010'. ascending = False
    #for g6 sorts the single mutation dataframe in position order, ie residue index. ascending = True

    singleMuts_sorted = singleMuts_df.sort_values(sort_col, key=lambda x: x.astype(int), ascending=ascending)
    singles_scores = singleMuts_sorted[condition].to_list()

    for g in dms_df['genotype']:
        multi_scores.append(apply_mask_and_average(singles_scores, g ))

    return multi_scores

#get melted correlations matrix of structure input x target antigen, for a given antibody 
def get_corr(ab_name, files, dms_df, ag_columns, dropLOQ = False ):

    #retreive correlations from abysis and esm1v scores
    for key, filepath in files.items():        
        singleMuts_df = pd.read_csv(filepath)
        # Get the columns for all esm1v models
        esm_columns = [col for col in singleMuts_df.columns if col.startswith('esm')]      
        # Calculate the average for each esm column
        esm_avg_values = singleMuts_df[esm_columns].mean(axis=1)
        singleMuts_df[key] = esm_avg_values
        if ab_name == 'g6':
            #g6 data set is only single mutations
            dms_df[key] = esm_avg_values
        else:
            #cr datasets are multiple mutations
            dms_df[key] = transform_single_to_multi(dms_df, singleMuts_df, key, 'genotype', ascending = False) 


    conditions = list(files.keys())

    if not dropLOQ:
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
        melted_correlations['Method'] = ['LM' for i in melted_correlations['Input']]

        return melted_correlations
    #if we wish to drop the samples on the lower limit of quantitation, drop samples from each ag group individually and compute correlation with IF log_like
    else:
        all_ag_melted = pd.DataFrame({})
        for ag in ag_columns:
            filt_df = dms_df[dms_df[ag] > min(dms_df[ag])]

            correlations = filt_df[[ag] + conditions].corr(method='spearman')
            correlations = correlations.drop(conditions, axis= 1)
            correlations = correlations.drop(ag, axis= 0)

            # Melt the correlations DataFrame and rename the index column
            melted_correlations = (
                correlations
                .reset_index()
                .melt(id_vars='index', var_name='Target', value_name='Correlation')
            ).rename(columns={'index': 'Input'})

            all_ag_melted = pd.concat([all_ag_melted, melted_correlations], ignore_index= True)
        all_ag_melted['Method'] =   ['LM' for i in all_ag_melted['Input']]
        
        return all_ag_melted
    
#get melted correlations matrix of structure input x target antigen, for a given antibody 
def get_g6_corr(g6Hc_dict, g6LC_dict, dropLOQ = False ):

    ab_name = 'g6'
    ag_columns = g6LC_dict['ag_columns']
    vh_and_vl_dms_df = pd.DataFrame({})

    for d in [g6HC_dict, g6LC_dict]:

        files = d['files']
        dms_df = d['dms_df']

        #retreive correlations from abysis, ablang, and esm1v scores
        for key, filepath in files.items():        
            singleMuts_df = pd.read_csv(filepath)
            # Get the columns for all esm1v models
            esm_columns = [col for col in singleMuts_df.columns if col.startswith('esm')]      
            # Calculate the average for each esm column
            esm_avg_values = singleMuts_df[esm_columns].mean(axis=1)
            singleMuts_df[key] = esm_avg_values
            if ab_name == 'g6':
                #g6 data set is only single mutations
                dms_df[key] = esm_avg_values


        vh_and_vl_dms_df = pd.concat([vh_and_vl_dms_df, dms_df], ignore_index= True)

        conditions = list(files.keys())

    #
    if not dropLOQ:
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
        melted_correlations['Method'] = ['LM' for i in melted_correlations['Input']]
        melted_correlations['Target'] = ['VEGF-A'] * len(melted_correlations)

        return melted_correlations
    #if we wish to drop the samples on the lower limit of quantitation, drop samples from each ag group individually and compute correlation with IF log_like
    else:
        all_ag_melted = pd.DataFrame({})
        for ag in ag_columns:
            filt_df = vh_and_vl_dms_df[vh_and_vl_dms_df[ag] > min(vh_and_vl_dms_df[ag])]

            correlations = filt_df[[ag] + conditions].corr(method='spearman')
            correlations = correlations.drop(conditions, axis= 1)
            correlations = correlations.drop(ag, axis= 0)

            # Melt the correlations DataFrame and rename the index column
            melted_correlations = (
                correlations
                .reset_index()
                .melt(id_vars='index', var_name='Target', value_name='Correlation')
            ).rename(columns={'index': 'Input'})

            all_ag_melted = pd.concat([all_ag_melted, melted_correlations], ignore_index= True)
        all_ag_melted['Method'] = ['LM' for i in melted_correlations['Input']]
        melted_correlations['Target'] = ['VEGF-A'] * len(melted_correlations)

        return all_ag_melted    
    
#plot correlation bars for cr antibodies: cr6261  and cr9114
def plot_hbar(title, melted_correlations, palette = bokeh.palettes.Spectral6, compare = 'model'):

    if compare == 'model':
        melted_correlations.replace('Ab-Ag','InverseFolding', inplace = True)
        color_by = 'Input'
    elif compare == 'IF':
        color_by = 'Method'
    else:
        color_by = 'Target'


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
        legend_field = color_by,
        # use the palette to colormap based on the the x[1:2] values
        # fill_color=factor_cmap(
        #                         color_by, 
        #                         palette= palette, 
        #                         factors = list(melted_correlations[color_by].unique()), 
        #                         start=1, 
        #                         end=3)
        fill_color=bokeh.palettes.Dark2_6[1]
             
    )

    #labels_df = melted_correlations[melted_correlations['Input'].isin(['Ab-Ag','ESM-1v', 'abYsis', 'InverseFold'])]
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
        
    #compute with and without points on lower limit of quantitation ignored
    for dropLOQ in [False]:

        all_corr_plots = []

        for d in datasets:
            if type(d) is tuple:
                g6Hc, g6Lc = d
                title = g6Hc['ab_name'] + ', ' + g6Hc['expt_type']
                g6_combined = get_g6_corr(g6Hc, g6Lc, dropLOQ= dropLOQ)
            
                all_corr_plots.append(plot_hbar(title, g6_combined, g6Lc['palette'], 'IF'))

            else:
                corr_df = get_corr(d['ab_name'], d['files'], d['dms_df'], d['ag_columns'], dropLOQ= dropLOQ)
                title = d['ab_name'] + ', ' + d['expt_type']
                all_corr_plots.append(plot_hbar(title, corr_df, d['palette'], 'Target'))

    
        all_corr_fname = f"output/ab_mutagenesis_expts/esm1v_benchmarks{'_dropLOQ' if dropLOQ else ''}.html"
        bokeh.plotting.output_file(all_corr_fname)
        bokeh.io.show(bokeh.layouts.gridplot(all_corr_plots, ncols = len(all_corr_plots)))
