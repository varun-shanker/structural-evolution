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

# Your Pandas code here

if_conds = ['ag+ab', 'ab only', 'ab']

cr9114_dict = {
                'ab_name'   : 'CR9114',
                'files'     : {
                                'Ab-Ag':      'output/ab_mutagenesis_expts/cr9114/4fqi_ablh_scores.csv',
                                'Ab only':    'output/ab_mutagenesis_expts/cr9114/4fqi_lh_scores.csv',
                                'Ab VH only': 'output/ab_mutagenesis_expts/cr9114/4fqi_h_scores.csv',
                                'ESM-1v':      'output/ab_mutagenesis_expts/cr9114/cr9114_exp_data_maskMargLabeled.csv',
                                'AbLang':     'output/ab_mutagenesis_expts/cr9114/cr9114_hc_ablangScores.csv',
                                'abYsis':     'output/ab_mutagenesis_expts/cr9114/abysis_counts_cr9114_vh.txt',
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
                                'Ab-Ag':      'output/ab_mutagenesis_expts/cr6261/3gbn_ablh_scores.csv',
                                'Ab only':    'output/ab_mutagenesis_expts/cr6261/3gbn_lh_scores.csv',
                                'Ab VH only': 'output/ab_mutagenesis_expts/cr6261/3gbn_h_scores.csv',
                                'ESM-1v':      'output/ab_mutagenesis_expts/cr6261/cr6261_exp_data_maskMargLabeled.csv',
                                'AbLang':     'output/ab_mutagenesis_expts/cr6261/cr6261_hc_ablangScores.csv',
                                'abYsis':     'output/ab_mutagenesis_expts/cr6261/abysis_counts_cr6261_vh.txt',
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
                                'Ab-Ag':      'output/ab_mutagenesis_expts/g6/2fjg_vlh_lc_scores.csv',
                                'Ab only':    'output/ab_mutagenesis_expts/g6/2fjg_lh_lc_scores.csv',
                                'Ab VH/VL only': 'output/ab_mutagenesis_expts/g6/2fjg_l_lc_scores.csv',
                                'ESM-1v':      'output/ab_mutagenesis_expts/g6/g6Lc_exp_data_maskMargLabeled.csv',
                                'AbLang':     'output/ab_mutagenesis_expts/g6/g6_lc_ablangScores.csv',
                                'abYsis':     'output/ab_mutagenesis_expts/g6/abysis_counts_g6_vl.txt',
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
                                'Ab-Ag':      'output/ab_mutagenesis_expts/g6/2fjg_vlh_hc_scores.csv',
                                'Ab only':      'output/ab_mutagenesis_expts/g6/2fjg_lh_hc_scores.csv',
                                'Ab VH/VL only': 'output/ab_mutagenesis_expts/g6/2fjg_h_hc_scores.csv',
                                'ESM-1v':      'output/ab_mutagenesis_expts/g6/g6Hc_exp_data_maskMargLabeled.csv',
                                'AbLang':     'output/ab_mutagenesis_expts/g6/g6_hc_ablangScores.csv',
                                'abYsis':     'output/ab_mutagenesis_expts/g6/abysis_counts_g6_vh.txt',
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

def get_abysis_single_scores_df(singleMut_ids, abysis_df):
    single_scores = []
    all_pos = []

    for id in singleMut_ids:
        wt, pos, mt = id[0], int(id[1:-1]), id[-1]
        likelihood_ratio = abysis_df.loc[(abysis_df['pos'] == pos) & (abysis_df['wt'] == wt) & (abysis_df['mt'] == mt)]['likelihood_ratio'].to_list()[0]
        log_like_ratio = np.log10(likelihood_ratio)
        single_scores.append(log_like_ratio)
        all_pos.append(pos)

    singleMuts_df = pd.DataFrame({'pos': all_pos, 'abYsis': single_scores})


    return singleMuts_df

def get_ablang_single_scores_df(singleMut_ids, ablang_df):
    single_scores = []
    all_pos = []

    for id in singleMut_ids:
        wt, pos, mt = id[0], int(id[1:-1]), id[-1]
        log_likelihood_ratio = ablang_df.loc[(ablang_df['pos'] == pos) & (ablang_df['wt'] == wt) & (ablang_df['mt'] == mt)]['log_likelihood_ratio'].to_list()[0]
        single_scores.append(log_likelihood_ratio)
        all_pos.append(pos)

    singleMuts_df = pd.DataFrame({'pos': all_pos, 'AbLang': single_scores})


    return singleMuts_df


#get melted correlations matrix of structure input x target antigen, for a given antibody 
def get_corr(ab_name, files, dms_df, ag_columns, dropLOQ = False ):

    inv_fold_files = {key: value for key, value in files.items() if key not in ['ESM-1v', 'AbLang', 'abYsis']}
    other_files = {key: value for key, value in files.items() if key in ['ESM-1v', 'AbLang', 'abYsis']}

    #retreive correlations from inverse fold scores
    for key, filepath in inv_fold_files.items():
        df = pd.read_csv(filepath)  # Read the CSV file into a DataFrame
        column_label = key
        log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
        dms_df[column_label] = log_likelihood_column

    #retreive correlations from abysis and esm1v scores
    for key, filepath in other_files.items():        
        if key == 'ESM-1v':
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
        elif key == 'abYsis':
            abysis_df = pd.read_csv(filepath, sep = '\t', dtype = {'pos': int})
            if ab_name == 'g6':
                dms_df[key] = get_abysis_single_scores_df(dms_df['mutant'], abysis_df)[key]
            else:
                #cr data set are multiple mutations. use same scoring strategy as esm1v paper
                singleMut_ids = pd.read_csv(files['ESM-1v'])['mutant'].to_list()
                singleMuts_df = get_abysis_single_scores_df(singleMut_ids, abysis_df)
                dms_df[key] = transform_single_to_multi(dms_df, singleMuts_df, key, 'pos', ascending = True)
        elif key == 'AbLang':
            ablang_df = pd.read_csv(filepath, dtype = {'pos': int})
            if ab_name == 'g6':
                dms_df[key] = get_ablang_single_scores_df(dms_df['mutant'], ablang_df)[key]
            else:
                #cr data set are multiple mutations. use same scoring strategy as esm1v paper
                singleMut_ids = pd.read_csv(files['ESM-1v'])['mutant'].to_list()
                singleMuts_df = get_ablang_single_scores_df(singleMut_ids, ablang_df)
                dms_df[key] = transform_single_to_multi(dms_df, singleMuts_df, key, 'pos', ascending = True)

    conditions = list(files.keys())

    #
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
        melted_correlations['Method'] = ['LM' if (('ESM' in i) or ('AbLang' in i))
                                         else 'IF' if 'Ab' in i 
                                        else 'MSA'
                                        for i in melted_correlations['Input']
                                    ]

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
        all_ag_melted['Method'] =   ['LM' if (('ESM' in i) or ('AbLang' in i))
                                         else 'IF' if 'Ab' in i 
                                        else 'MSA'
                                        for i in all_ag_melted['Input']
                                    ]
        
        return all_ag_melted
    

#get melted correlations matrix of structure input x target antigen, for a given antibody 
def get_g6_corr(g6Hc_dict, g6LC_dict, dropLOQ = False ):

    ab_name = 'g6'
    ag_columns = g6LC_dict['ag_columns']
    vh_and_vl_dms_df = pd.DataFrame({})

    for d in [g6HC_dict, g6LC_dict]:

        files = d['files']
        dms_df = d['dms_df']

        inv_fold_files =  {key: value for key, value in files.items() if key not in ['ESM-1v', 'AbLang', 'abYsis']}
        other_files =  {key: value for key, value in files.items() if key in ['ESM-1v', 'AbLang', 'abYsis']}

        #retreive correlations from inverse fold scores
        for key, filepath in inv_fold_files.items():
            df = pd.read_csv(filepath)  # Read the CSV file into a DataFrame
            column_label = key
            log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
            dms_df[column_label] = log_likelihood_column

        #retreive correlations from abysis, ablang, and esm1v scores
        for key, filepath in other_files.items():        
            if key == 'ESM-1v':
                singleMuts_df = pd.read_csv(filepath)
                # Get the columns for all esm1v models
                esm_columns = [col for col in singleMuts_df.columns if col.startswith('esm')]      
                # Calculate the average for each esm column
                esm_avg_values = singleMuts_df[esm_columns].mean(axis=1)
                singleMuts_df[key] = esm_avg_values
                if ab_name == 'g6':
                    #g6 data set is only single mutations
                    dms_df[key] = esm_avg_values
            elif key == 'abYsis':
                abysis_df = pd.read_csv(filepath, sep = '\t', dtype = {'pos': int})
                if ab_name == 'g6':
                    dms_df[key] = get_abysis_single_scores_df(dms_df['mutant'], abysis_df)[key]
            elif key == 'AbLang':
                ablang_df = pd.read_csv(filepath, dtype = {'pos': int})
                if ab_name == 'g6':
                    dms_df[key] = get_ablang_single_scores_df(dms_df['mutant'], ablang_df)[key]


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
        melted_correlations['Method'] = ['LM' if (('ESM' in i) or ('AbLang' in i))
                                         else 'IF' if 'Ab' in i 
                                        else 'MSA'
                                        for i in melted_correlations['Input']
                                    ]
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
        all_ag_melted['Method'] = ['LM' if (('ESM' in i) or ('AbLang' in i))
                                         else 'IF' if 'Ab' in i 
                                        else 'MSA'
                                        for i in all_ag_melted['Input']
                                    ]
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

    if 'CR6261' in title:
        x_min = -0.15
    else:
        x_min = 0
    
    p = bokeh.plotting.figure(
    height=340,
    width=440,
    x_axis_label="Spearman Correlation",
    x_range=[x_min, 1],
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
        fill_color=factor_cmap(
                                color_by, 
                                palette= palette, 
                                factors = list(melted_correlations[color_by].unique()), 
                                start=1, 
                                end=3)
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

    if compare == 'IF':
        p.legend.orientation = "vertical"
        p.legend.location = "bottom_right"
        p.legend.label_text_font_size = "8pt"
    else:
        p.legend.visible = False

    p.output_backend = "svg"
    return p

def plot_g6_scatter(g6Hc_df,  g6Lc_df, ag_columns, dropLOQ):
    ag = ag_columns[0]
    if dropLOQ:
        g6Lc_df = g6Lc_df[g6Lc_df[ag] > np.min(g6Lc_df[ag])]
        g6Hc_df = g6Hc_df[g6Hc_df[ag] > np.min(g6Hc_df[ag])]


    hv.extension("bokeh")

    plots = []
    for (plot_df,title) in [(g6Hc_df, 'g6.31 VH'),(g6Lc_df, 'g6.31 VL')]:

        # Generate HoloViews Points Element
        points = hv.Points(
            data=plot_df,
            kdims=[ag, 'Ab-Ag'],
        )

        # Datashade with spreading of points
        p = hv.operation.datashader.dynspread(
            hv.operation.datashader.rasterize(
                points
            ).opts(
                cmap='Magma',
                cnorm='linear',
            )
            ).opts(
                frame_width=350,
                frame_height=300,
                padding=0.05,
                show_grid=False,
                colorbar = True
            )

        p = hv.render(p)
        p.output_backend = "svg"

        plots.append(p)

    return plots

#plot scatter plot colored by n_muts for CR antibodies
def single_cr_scatterMut(plot_df, ag_columns, x, y, title, ):
    x_axis_label = x + '-logKd' if x in ag_columns else 'log likelihood'
    y_axis_label = y + '-logKd' if y in ag_columns else 'log likelihood'


    #set up figure
    p = bokeh.plotting.figure(
        width=500,
        height=400,
        x_axis_label= x_axis_label,
        y_axis_label= y_axis_label,
        title = title
    )

    #Set up color mapper to color by number of mutations from germline
    mapper = bokeh.transform.linear_cmap(field_name='som_mut', palette=bokeh.palettes.inferno(16) ,low=min(plot_df.som_mut) ,high=max(plot_df.som_mut))
    p.circle(
        source=plot_df,
        x= x,
        y= y,
        alpha = 0.2, 
        line_color=mapper,
        color=mapper,
        size = 4
    )

    #
    p.diamond(
        #germline seqid is all 0's
        source = plot_df[plot_df['genotype'].str.contains('1') == False],
        x= x,
        y= y,
        color = 'dodgerblue',
        size = 12,
        legend_label = 'germline'
    )
    p.star(
        #mature seqid is all 0's
        source = plot_df[plot_df['genotype'].str.contains('0') == False],
        x= x,
        y= y,
        color = 'green',
        size = 12,
        legend_label = 'mature'
    )
    p.legend.location = "top_left"
    p.output_backend = "svg"

    return p, mapper

#plot single scatter plot for CR antibodies - with dynamic spreading and datashading
def single_cr_scatter(plot_df, ag_columns, x, y, title, ):
    x_axis_label = x + '-logKd' if x in ag_columns else 'log likelihood'
    y_axis_label = y + '-logKd' if y in ag_columns else 'log likelihood'

    hv.extension("bokeh")

    # Generate HoloViews Points Element
    points = hv.Points(
        data=plot_df,
        kdims=[x, y],
    )

    # Datashade with spreading of points
    p = hv.operation.datashader.dynspread(
        hv.operation.datashader.rasterize(
            points
        ).opts(
            cmap='Magma',
            cnorm='linear',
        )
        ).opts(
            frame_width=350,
            frame_height=300,
            padding=0.05,
            show_grid=False,
            colorbar = True
        )

    p = hv.render(p)
    
    p.grid.visible = False
    #
    p.diamond(
        #germline seqid is all 0's
        source = plot_df[plot_df['genotype'].str.contains('1') == False],
        x= x,
        y= y,
        color = 'dodgerblue',
        size = 12,
        #legend_label = 'germline'
    )
    p.star(
        #mature seqid is all 0's
        source = plot_df[plot_df['genotype'].str.contains('0') == False],
        x= x,
        y= y,
        color = 'green',
        size = 12,
        #legend_label = 'mature'
    )
    #p.legend.location = "top_left"
    p.output_backend = "svg"

    return p

def plot_CR_scatter_colorMuts(title, files, dms_df, ag_columns, dropLOQ = False):

    agAb_files  = {key: value for key, value in files.items() if key in ['Ab-Ag',]}

    for key, filepath in agAb_files.items():
        df = pd.read_csv(filepath)  # Read the CSV file into a DataFrame
        column_label = key
        log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
        dms_df[column_label] = log_likelihood_column

    plots = []
    #predictions vs experiment
    for ag in ag_columns:

        if dropLOQ:
            plot_df = dms_df[dms_df[ag] > min(dms_df[ag])]
            #add back wt which is on LOQ for visualization   
            plot_df = pd.concat([plot_df, dms_df[dms_df['som_mut'] == 0]], ignore_index = True, axis = 0)
        else:
            plot_df = dms_df

        p, _ = single_cr_scatterMut(plot_df, ag_columns, x = ag, y = 'Ab-Ag', title = f'{title} against {ag}')
        plots.append(p)

    #experimental vs experimental 
    plot_df = dms_df
    if dropLOQ: 
        for ag in ag_columns:
            plot_df = plot_df[plot_df[ag] > min(plot_df[ag])]
        #add back wt which is on LOQ for visualization   
        plot_df = pd.concat([plot_df, dms_df[dms_df['som_mut'] == 0]], ignore_index = True, axis = 0)

    p_exp, mapper = single_cr_scatterMut(plot_df, ag_columns, x= ag_columns[0], y = ag_columns[1], title = 'Experimental Cross-Reactive Binding Landscape')
    p_exp.legend.location = "top_left"
    color_bar = bokeh.models.ColorBar(color_mapper=mapper['transform'], width=8, title = 'Amino Acid Mutations')
    #p_exp.add_layout(color_bar, 'right')

    plots.append(p_exp)
    p.output_backend = "svg"

    return plots

def plot_CR_scatterplots(title, files, dms_df, ag_columns, dropLOQ = False):

    agAb_files  = {key: value for key, value in files.items() if key in ['Ab-Ag',]}

    for key, filepath in agAb_files.items():
        df = pd.read_csv(filepath)  # Read the CSV file into a DataFrame
        column_label = key
        log_likelihood_column = df['log_likelihood']  # Extract the 'log_likelihood' column
        dms_df[column_label] = log_likelihood_column

    plots = []
    #predictions vs experiment
    for ag in ag_columns:

        if dropLOQ:
            plot_df = dms_df[dms_df[ag] > min(dms_df[ag])]
            #add back wt which is on LOQ for visualization   
            plot_df = pd.concat([plot_df, dms_df[dms_df['som_mut'] == 0]], ignore_index = True, axis = 0)
        else:
            plot_df = dms_df

        p = single_cr_scatter(plot_df, ag_columns, x = ag, y = 'Ab-Ag', title = f'{title} against {ag}')
        plots.append(p)

    #experimental vs experimental 
    plot_df = dms_df
    if dropLOQ: 
        for ag in ag_columns:
            plot_df = plot_df[plot_df[ag] > min(plot_df[ag])]
        #add back wt which is on LOQ for visualization   
        plot_df = pd.concat([plot_df, dms_df[dms_df['som_mut'] == 0]], ignore_index = True, axis = 0)

    p_exp = single_cr_scatter(plot_df, ag_columns, x= ag_columns[0], y = ag_columns[1], title = 'Experimental Cross-Reactive Binding Landscape')

    plots.append(p_exp)
    p_exp.output_backend = "svg"

    return plots

#plot correlation bars for cr antibodies: cr6261  and cr9114
def plot_LOQ_comparison(ab_name, unfilt_corrs, filt_corrs, palette = bokeh.palettes.Spectral6,):
    unfilt_corrs['LOQ'] = ['Included'] * len(unfilt_corrs)
    filt_corrs['LOQ'] = ['Excluded'] * len(filt_corrs)
    combined_melt = pd.concat([unfilt_corrs,filt_corrs], ignore_index= True)

    combined_melt['cats'] = combined_melt.apply(lambda x: (x["Target"], x["Input"], x["LOQ"]), axis = 1)
    factors = list(combined_melt.cats)[::-1]
    
    p = bokeh.plotting.figure(
    height=850,
    width=400,
    x_axis_label="Spearman Correlation",
    x_range=[-.25, 1],
    y_range=bokeh.models.FactorRange(*factors),
    tools="save",
    title = f'{ab_name}, Impact of Ignoring Lower LOQ '
    )

    p.hbar(
        source=combined_melt,
        y="cats",
        right="Correlation",
        height=0.7,
        line_color = 'black',
        legend_field = 'LOQ',
        # use the palette to colormap based on the the x[1:2] values
        fill_color=factor_cmap(
                                'LOQ', 
                                palette= palette, 
                                factors = list(combined_melt['LOQ'].unique()), 
                                start=1, 
                                end=2)
    )

    p.ygrid.grid_line_color = None
    p.y_range.range_padding = 0.1
    p.legend.location = "bottom_right"
    p.legend.orientation = "vertical"
    p.yaxis.axis_label = "Model Input"
    p.output_backend = "svg"

    return p



if __name__ == '__main__':

    datasets = [cr9114_dict, cr6261_dict, (g6HC_dict, g6LC_dict) ]
    results = []
    
    #compute with and without points on lower limit of quantitation ignored
    for dropLOQ in [False, True]:

        g6_corrs = []
        all_corr_plots = []

        for d in datasets:
            if type(d) is tuple:
                g6Hc, g6Lc = d
                title = g6Hc['ab_name'] + ', ' + g6Hc['expt_type']
                g6_combined = get_g6_corr(g6Hc, g6Lc, dropLOQ= dropLOQ)

                all_corr_plots.append(plot_hbar(title, g6_combined, g6Lc['palette'], 'IF'))

                #also plot scatter plots
                scatter_fname = f"output/ab_mutagenesis_expts/g6_scatterplot_colorMuts{'_dropLOQ' if dropLOQ else ''}.html"
                bokeh.plotting.output_file(scatter_fname)
                g6_scatterMuts_plots = plot_g6_scatter(g6Hc['dms_df'],  g6Lc['dms_df'], g6Hc['ag_columns'], dropLOQ= dropLOQ)
                bokeh.io.show(bokeh.layouts.gridplot(g6_scatterMuts_plots, ncols = 2))

                if (dropLOQ == False):
                    results.append(g6_combined)

            else:
                corr_df = get_corr(d['ab_name'], d['files'], d['dms_df'], d['ag_columns'], dropLOQ= dropLOQ)
                title = d['ab_name'] + ', ' + d['expt_type']

                all_corr_plots.append(plot_hbar(title, corr_df, d['palette'], 'IF'))

                #also plot scatter plots for CR9114, CR6261
                scatter_fname = f"output/ab_mutagenesis_expts/{d['ab_name']}_scatter_colorMuts{'_dropLOQ' if dropLOQ else ''}.html"
                bokeh.plotting.output_file(scatter_fname)
                cr_scatterMuts_plots = plot_CR_scatter_colorMuts(title,  d['files'], d['dms_df'], d['ag_columns'], dropLOQ= dropLOQ)
                bokeh.io.show(bokeh.layouts.gridplot(cr_scatterMuts_plots, ncols = 3))

                #also scatter plots with datashading/dynamic spreading for CR9114, CR6261
                scatter_fname = f"output/ab_mutagenesis_expts/{d['ab_name']}_scatter{'_dropLOQ' if dropLOQ else ''}.html"
                bokeh.plotting.output_file(scatter_fname)
                cr_scattter_plots = plot_CR_scatterplots(title,  d['files'], d['dms_df'], d['ag_columns'], dropLOQ= dropLOQ)
                bokeh.io.show(bokeh.layouts.gridplot(cr_scattter_plots, ncols = 3))

                if (dropLOQ == False):
                    results.append(corr_df)
     
        all_corr_fname = f"output/ab_mutagenesis_expts/benchmarks{'_dropLOQ' if dropLOQ else ''}.html"
        bokeh.plotting.output_file(all_corr_fname)
        bokeh.io.show(bokeh.layouts.gridplot(all_corr_plots, ncols = len(all_corr_plots)))

    compareLOQ_fname = 'output/ab_mutagenesis_expts/LOQ_comparison.html'
    bokeh.plotting.output_file(compareLOQ_fname)
    compareLOQ_plots = []
    for d in datasets:
        if (type(d) != tuple) and (d['ab_name'].startswith('CR')):  
            corr_withLOQ =  get_corr(d['ab_name'], d['files'], d['dms_df'], d['ag_columns'], )
            corr_withoutLOQ =  get_corr(d['ab_name'], d['files'], d['dms_df'], d['ag_columns'], dropLOQ= dropLOQ)

            compareLOQ_plots.append(plot_LOQ_comparison(d['ab_name'], corr_withLOQ, corr_withoutLOQ , d['palette']))
    bokeh.io.show(bokeh.layouts.gridplot(compareLOQ_plots, ncols = len(compareLOQ_plots)))

    pd.concat(results, ignore_index= True).to_csv( 'output/ab_mutagenesis_expts/results.csv' , index = False)
