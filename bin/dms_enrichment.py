import argparse
from pathlib import Path
from biotite.sequence.io.fasta import FastaFile, get_sequences
from Bio import pairwise2
import pandas as pd
import shutil
import os
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import bokeh.palettes 
from bokeh.transform import factor_cmap

from bokeh.io import export_svg
import iqplot 
from colorcet import glasbey_category10

import subprocess

bla_dict = {
    'protein':'bla', 
    'pdbfile':'bla_1m40_a.pdb', 
    'chain':'A', 
    'dmsfile': 'dms_bla.csv', 
    'fitness_col': 'DMS_amp_2500_(b)', 
    'threshold' : 0.01
}

calm1_dict = {
    'protein':'CALM1', 
    'pdbfile':'calm1_5v03_r.pdb', 
    'chain':'R', 
    'dmsfile': 'dms_calm1.csv', 
    'fitness_col':'DMS', 
    'threshold': 1
}

haeiiim_dict = {
    'protein':'haeIIIM', 
    'pdbfile':'haeiiim_3ubt_b.pdb', 
    'chain':'B', 
    'dmsfile': 'dms_haeiiim.csv', 
    'fitness_col':'DMS_G3', 
    'threshold' : 0.01
}

gal4_dict = {
    'protein':'GAL4', 
    'pdbfile':'gal4_3coq_b.pdb', 
    'chain':'B', 
    'dmsfile': 'dms_gal4.csv', 
    'fitness_col':'DMS_nonsel_24', 
    'threshold': 1
}

hras_dict = {
    'protein':'HRAS', 
    'pdbfile':'hras_2ce2_x.pdb', 
    'chain':'X', 
    'dmsfile': 'dms_hras.csv', 
    'fitness_col':'DMS_unregulated', 
    'threshold': 1
}

mapk1_dict = {
    'protein':'MAPK1', 
    'pdbfile':'mapk1_4zzn_a.pdb', 
    'chain':'A', 
    'dmsfile': 'dms_mapk1.csv', 
    'fitness_col':'DMS_VRT', 
    'threshold': 1
}

tpk1_dict = {
    'protein':'TPK1', 
    'pdbfile':'tpk1_3s4y_a.pdb', 
    'chain':'A', 
    'dmsfile': 'dms_tpk1.csv', 
    'fitness_col':'DMS', 
    'threshold': 1
}

tpmt_dict = {
    'protein':'TPMT', 
    'pdbfile':'tpmt_2bzg_a.pdb', 
    'chain':'A', 
    'dmsfile': 'dms_tpmt.csv', 
    'fitness_col':'DMS', 
    'threshold': 1
}

ube2i_dict = {
    'protein':'UBE2I', 
    'pdbfile':'ube2i_5f6e_a.pdb', 
    'chain':'A', 
    'dmsfile': 'dms_tpmt.csv', 
    'fitness_col':'DMS', 
    'threshold': 1
}

ubi4_dict = {
    'protein':'UBI4', 
    'pdbfile':'ubi4_4q5e_b.pdb', 
    'chain':'B', 
    'dmsfile': 'dms_ubi4.csv', 
    'fitness_col':'DMS_limiting_(b)', 
    'threshold': 1
}

studies = [           
        bla_dict,
        calm1_dict,
        gal4_dict, 
        haeiiim_dict, 
        hras_dict,
        tpmt_dict,
        tpk1_dict,
        mapk1_dict,
        ube2i_dict,
        ubi4_dict,
 ]

def score_if(prot_dict):
    prefix = 'data/dms/'
    fasta_outpath = f"dms_libs/{prot_dict['pdbfile'][:-4]}_dmsLib.fasta"
    score_if_cmd = [
        'python', 
        'bin/score_log_likelihoods.py',
        prefix + f"structures/{prot_dict['pdbfile']}", 
        '--seqpath', prefix + fasta_outpath, 
        '--outpath', f"output/dms/if_results/dms_{prot_dict['protein']}_if.csv"   
        ,'--chain', prot_dict['chain']
    ]
    subprocess.run(score_if_cmd, check=True) 

def get_stats(prot_dict, n, top_per):

    top_per = top_per/100
    prefix = 'output/dms/'
    esm1v_outpath = 'dms_esm1vScored/' + f"dms_{prot_dict['protein']}" + '_maskMarginals.csv'
    if_results_df = pd.read_csv(  prefix + f"if_results/dms_{prot_dict['protein']}_if.csv")
    dms_df = pd.read_csv(prefix + esm1v_outpath)
    
    unprofiled_df = dms_df[dms_df[prot_dict['fitness_col']].isna()]
    unprofiled_variants = unprofiled_df['variant'].to_list()
    
    dms_df = dms_df.dropna(subset = [prot_dict['fitness_col']])
    pop_n = len(dms_df)
    if_results_df = if_results_df[~if_results_df['seqid'].isin(unprofiled_variants)]

    assert (len(if_results_df) == len(dms_df))

    threshold = dms_df[prot_dict['fitness_col']].quantile(1-top_per)
    highFit_subset_variants = dms_df[dms_df[prot_dict['fitness_col']] >= threshold]['variant']
    n_highFit = len(highFit_subset_variants)


    dms_df['esm1v'] = dms_df.loc[:, dms_df.columns.str.startswith('esm1v')].mean(axis=1)
    top_esm1v_df = dms_df.sort_values(by = 'esm1v', ascending = False)[:n]
    top_esm1v_variants = top_esm1v_df['variant'].to_list()

    top_if_df = if_results_df.sort_values(by = 'log_likelihood', ascending = False)[:n]
    top_if_variants = top_if_df['seqid'].to_list()
    
    esm1v_hits_df = top_esm1v_df[top_esm1v_df['variant'].isin(highFit_subset_variants)]
    if_hits_df = top_if_df[top_if_df['seqid'].isin(highFit_subset_variants)]
    
    n_esm1v_hits = len(esm1v_hits_df)
    n_if_hits = len(if_hits_df)

    if_hit_enrich = (n_if_hits / n) / top_per
    esm1v_hit_enrich = (n_esm1v_hits / n) / top_per

    if_hit_rate = (n_if_hits / n)
    esm1v_hit_rate = (n_esm1v_hits / n)

    return (prot_dict['protein'], pop_n, n_highFit, n_if_hits, n_esm1v_hits, if_hit_enrich, esm1v_hit_enrich, if_hit_rate, esm1v_hit_rate )

def plot_bars(melted_df):
    p = bokeh.plotting.figure(
    height=350,
    width=1100,
    y_axis_label="High Fitness \n Prediction Precision",
    x_axis_label = "Functional Percentile Threshold \n for High Fitness Classification",
    x_range=bokeh.models.FactorRange(*factors, group_padding = 1.2),
    tools="save",
    title=''
    )
    p.output_backend = "svg"

    p.vbar(
        source=melted_df,
        x="cats",
        top="Hit Rate",
        width = 1,
        line_color='black',
        alpha = 'alpha',
        legend_field='legend_label',
        fill_color=bokeh.transform.factor_cmap(
            'Method',
            palette=['#999999', '#43a2ca'],
            factors=list(melted_df['Method'].unique()),
            start=1,
            end=2
        )
    )

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.x_range.range_padding = 0.03
    p.legend.location = "top_right"
    p.legend.spacing = 20
    p.legend.label_text_font_size = "11pt"
    p.legend.label_text_color = "black"
    p.legend.orientation = "horizontal"
    p.xaxis.major_label_orientation = 1.2
    p.xaxis.separator_line_alpha = 0
    p.xaxis.group_text_font_size = '11pt'
    p.xaxis.subgroup_text_font_size = '0pt'
    p.xaxis.axis_label_text_font_size = '11pt'
    p.xaxis.major_label_text_font_size = '10pt'
    p.xaxis.axis_label_text_font_style = 'normal'
    p.xaxis.major_label_text_color = 'black'
    p.xaxis.axis_label_text_color = 'black'

    p.yaxis.axis_label_text_font_size = '11pt'
    p.yaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_color = 'black'
    p.yaxis.major_label_text_font_size = '10pt'
    p.xaxis.group_text_color = 'dimgrey'

    return p

def plot_ecdf_hits(prot_dict, n, top_per_lst):

    #use least stringent threshold to plot all hits
    top_per = max(top_per_lst)/100
    prefix = 'output/dms/'

    esm1v_outpath = 'dms_esm1vScored/' + f"dms_{prot_dict['protein']}" + '_maskMarginals.csv'
    if_results_df = pd.read_csv(  prefix + f"if_results/dms_{prot_dict['protein']}_if.csv")
    dms_df = pd.read_csv(prefix + esm1v_outpath)

    
    unprofiled_df = dms_df[dms_df[prot_dict['fitness_col']].isna()]
    unprofiled_variants = unprofiled_df['variant'].to_list()
    
    dms_df = dms_df.dropna(subset = [prot_dict['fitness_col']])
    dms_df['percentile'] = dms_df[prot_dict['fitness_col']].rank(pct=True, method='average')
    dms_df['zscore'] = (dms_df[prot_dict['fitness_col']] - dms_df[prot_dict['fitness_col']].mean())/dms_df[prot_dict['fitness_col']].std()
    if_results_df = if_results_df[~if_results_df['seqid'].isin(unprofiled_variants)]
    
    assert (len(if_results_df) == len(dms_df))

    threshold = dms_df[prot_dict['fitness_col']].quantile(1-top_per)
    highFit_subset_variants = dms_df[dms_df[prot_dict['fitness_col']] >= threshold]['variant']
    n_highFit = len(highFit_subset_variants)

    dms_df['esm1v'] = dms_df.loc[:, dms_df.columns.str.startswith('esm1v')].mean(axis=1)
    top_esm1v_df = dms_df.sort_values(by = 'esm1v', ascending = False)[:n]
    top_esm1v_variants = top_esm1v_df['variant'].to_list()

    top_if_df = if_results_df.sort_values(by = 'log_likelihood', ascending = False)[:n]
    top_if_variants = top_if_df['seqid'].to_list()
    
    esm1v_hits_df = top_esm1v_df[top_esm1v_df['variant'].isin(highFit_subset_variants)]
    if_hits_df = top_if_df[top_if_df['seqid'].isin(highFit_subset_variants)]
    if_hits_df = pd.merge(if_hits_df, dms_df[['variant', 'zscore', 'percentile', prot_dict['fitness_col']]], left_on='seqid', right_on='variant', how='left')

    p_1v = iqplot.ecdf(
            data=dms_df,
            q=prot_dict['fitness_col'],
            style="staircase",
            palette = 'lightgrey',
            line_kwargs= {'line_width':6},
            x_axis_label = prot_dict['protein']+' Fitness Measure',
            #title = prot_dict['protein']
    )
    p_if = iqplot.ecdf(
            data=dms_df,
            q=prot_dict['fitness_col'],
            style="staircase",
            palette = 'lightgrey',
            line_kwargs= {'line_width':6},
            x_axis_label = None, #prot_dict['protein']+'Fitness Measure',
            title = prot_dict['protein']
    )
    colors = ['#1f78b4', '#33a02c', '#ff7f00',]
    for i, top_per in enumerate(top_per_lst):
        top_per = top_per/100
        hline = bokeh.models.Span(location=1-top_per, dimension='width', line_color=colors[i], line_dash='dashed', line_width=4, )
        p_1v.add_layout(hline)
        p_if.add_layout(hline)


    p_1v.circle(
            source = esm1v_hits_df,
            x = prot_dict['fitness_col'],
            y = 'percentile',
            color = 'dimgrey',
            size = 11,
            line_color = 'lightgrey',
            alpha = 0.8,
    )
 
    p_if.circle(
            source = if_hits_df,
            x = prot_dict['fitness_col'],
            y = 'percentile',
            color = '#43a2ca',
            size = 11,
            line_color = 'white',
            alpha = 0.8,
    )

    p_if.title.text_font_style = 'normal'
    p_if.title.text_color = 'black'
    p_if.title.align = 'center'
    p_if.title.text_font_size = '16pt'
    p_1v.xaxis.axis_label_text_font_size = '14pt'
    

    for p in [p_if, p_1v]:
        p.xaxis.axis_label_text_font_style = 'normal'
        p.yaxis.axis_label_text_font_style = 'normal'
        p.xaxis.axis_label_text_font_size = '14pt'
        p.yaxis.axis_label_text_font_size = '14pt'
        p.xaxis.axis_label_text_color = 'black'
        p.yaxis.axis_label_text_color = 'black'
        p.xaxis.major_label_text_font_size = '14pt'
        p.yaxis.major_label_text_font_size = '14pt'
        p.xaxis.major_label_text_color = 'black'
        p.yaxis.major_label_text_color = 'black'
        p.output_backend = "svg"

    return [p_if, p_1v]

if __name__ == '__main__':
    for s in studies:
        score_if(s)
    shutil.copytree("data/dms/dms_esm1vScored", "output/dms/dms_esm1vScored", dirs_exist_ok=True)

    name, n_pop, n_popHits, n_if_hits, n_esm1v_hits, if_hit_enrich, esm1v_hit_enrich, if_hit_rate, esm1v_hit_rate, p_threshold = ([] for _ in range(10))
    data_lists = [name, n_pop, n_popHits, n_if_hits, n_esm1v_hits, if_hit_enrich, esm1v_hit_enrich, if_hit_rate, esm1v_hit_rate, p_threshold]
    for s in studies:
        for p in [5, 10, 20 ]:
            outputs = get_stats(s, 10, p)
            for i in range(len(outputs)):
                data_lists[i].append(outputs[i])
            data_lists[len(data_lists)-1].append(p)
    col_names = ['Protein', 'Total Library Size', 'Library Hits (Variants with Fitness >95th)', 'Inverse Folding Hits (in Top 10)', 'ESM1v Hits (in Top 10)', 'Inverse Folding Hit Enrichment', 'ESM1v Hit Enrichment','Inverse Folding Hit Rate', 'ESM1v Hit Rate','Percentile Threshold' ]
    results_df = pd.DataFrame({col_names[i]: data_lists[i] for i in range(len(data_lists))})
    # Melt the dataframe
    melted_df = results_df.melt(id_vars=['Protein', 'Percentile Threshold'], value_vars=[ 'ESM1v Hit Rate', 'Inverse Folding Hit Rate',],
                        var_name='Method', value_name='Hit Rate')
    melted_df['Method'] = melted_df['Method'].replace({
        'Inverse Folding Hit Rate': 'Inverse Folding',
        'ESM1v Hit Rate': 'Language Model'})
    melted_df['alpha'] = (1/melted_df['Percentile Threshold'].to_numpy())*3 + 0.4
    melted_df['Percentile Threshold'] = melted_df['Percentile Threshold'].astype(str)
    melted_df['legend_label'] = ['Structure-Informed Language Model' if x=='Inverse Folding' else 'Language Model' for x in melted_df['Method']]

    melted_df['cats'] = melted_df.apply(lambda x: (x["Protein"],  x["Method"], x['Percentile Threshold'],), axis = 1)
    factors = list(melted_df.cats)

    p = plot_bars(melted_df)
    comparePrecision_fname = 'output/dms/precision_comparison.html'
    bokeh.plotting.output_file(comparePrecision_fname)
    bokeh.io.show(p)

    ecdf_plots = []
    for s in studies:
        ecdf_plots.extend(plot_ecdf_hits(s, 10, [20,10,5]))
    mid = int(len(ecdf_plots)/2)
    ecdf_fname = 'output/dms/fitness_ecdfs.html'
    bokeh.plotting.output_file(ecdf_fname)
    bokeh.io.show(bokeh.layouts.gridplot([ecdf_plots[:mid:2],
                                        ecdf_plots[1:mid:2],
                                        ecdf_plots[mid::2],
                                        ecdf_plots[mid+1::2]])
    )
