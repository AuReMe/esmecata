#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2023-2024 Pauline Hamon-Giraud and Victor Mataigne - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

#import warnings
import os
import json
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
import plotly.express as px
import plotly.graph_objects as go
from plotly.io import write_json
from plotly.subplots import make_subplots

from ontosunburst.ontosunburst import ontosunburst
from esmecata.report.esmecata_compression import RANK2COL

# from statistics import NormalDist
# import matplotlib.pyplot as plt
# import seaborn as sns
# import plotly.tools as tls
# import math
# import scipy.stats as st


__author__ = "Pauline Hamon-Giraud, Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

# warnings.filterwarnings('ignore')

# ================
# Global variables
# ================

PLOT_BGCOLOR='#e6ffe6'
LINECOLOR = 'DarkSlateGrey'
COLORSCALE = px.colors.qualitative.Dark24
EC_CLASSES_COLORS = {"Ligases": "#636EFA",
    "Oxidoreductases": "#EF553B",
    "Transferases": "#00CC96",
    "Hydrolases": "#AB63FA",
    "Lyases": "#FFA15A",
    "Isomerases": "#19D3F3",
    "Translocases": "#FF6692"} # 'Plotly' default color map
DOWNLOAD_FMT = "svg"
FONTPARAMS = dict(
    family="Courier New, monospace",
    size=14,
    color="black"
)

# ========================
# Code factoring functions
# ========================

def _format_trace(fig):
    '''
    Set visual harmonization of a figure trace
    '''
    fig.update_traces(marker=dict(line=dict(width=1, color=LINECOLOR)))

def _format_axes(fig):
    '''
    Set visual harmonization of a figure axes

    Parameters:
        fig : a plotly object
    '''  
    fig.update_xaxes(
        showline=True, 
        linewidth=2, 
        linecolor=LINECOLOR, 
        mirror=True, 
        showgrid=True, 
        gridwidth=2)
    fig.update_yaxes(
        showline=True, 
        linewidth=2, 
        linecolor=LINECOLOR, 
        mirror=True, 
        showgrid=True, 
        gridwidth=2)
    
# def _create_config(format, filename):
#     config = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d'],
#               'toImageButtonOptions': {'format': 'svg', # one of png, svg, jpeg, webp
#                                        'filename': 'custom_image'
#                 }
#             }
#     return config
    
# =============================
# Data gathering and formatting
# =============================

def _format_taxo(proteome_tax_id):
    '''
    Recovers the full taxonomic lineage of all OTUs/ASVs analyzed by esMeCaTa by parsing proteome_tax_id.tsv with ete3

        Parameters:
            proteome_tax_id (pandas): a DataFrame loaded with the output file proteome_tax_id.tsv of EsMeCaTa proteomes

        Returns:
            df (pandas): DataFrame with the full taxonomic lineage (one column = one rank)
    '''

    ncbi_ranks = ["superkingdom","kingdom","subkingdom","superphylum","phylum",
    "subphylum","infraphylum","superclass","class","subclass",
    "infraclass","cohort","subcohort","superorder","order",
    "suborder","infraorder","parvorder","superfamily","family",
    "subfamily","tribe","subtribe","genus","subgenus",
    "section","subsection","series","subseries","species group",
    "species subgroup","species","forma specialis","subspecies","varietas",
    "subvariety","forma","serogroup","serotype","strain","isolate"]

    ncbi = NCBITaxa()

    ranks = {}
    names = {}
    data = {}

    # ete3 to get lineages
    for index, row in proteome_tax_id.iterrows():
        tmp = ncbi.get_lineage(row['tax_id'])
        ranks[index] = ncbi.get_rank(tmp)
        names[index] = ncbi.get_taxid_translator(tmp)

    # rank:name correspondance
    for k in ranks.keys():
        data[k] = {}
        for k2,v in ranks[k].items():
            data[k].update({v:names[k][k2]})

    df = pd.DataFrame.from_dict(data, orient='index')
    if 'no rank' in df.columns:
        df.drop('no rank', axis=1, inplace=True)

    # order columns according to hierarchy
    ncbi_ranks = [rank for rank in ncbi_ranks if rank in df.columns]
    df = df[ncbi_ranks]

    return df

def post_analysis_config(input_table, results_path):
    '''
    Proceeds to format the inputs and outputs data in order to un all the others plotting functions efficiently

        Parameters:
            input_table (str): the path to the input table of taxa
            results_path (str): the path of the esmecata run

        Returns:
            data (dict) : gathers all useful data for the analysis of an esmecata run (inputs, summary tables, taxonomy ...)
    '''

    # Load and concat summaries of results
    stats_nb_proteomes = pd.read_table(
        os.path.join(results_path,'0_proteomes/stat_number_proteome.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')
    
    stats_nb_clustering = pd.read_table(
        os.path.join(results_path,'1_clustering/stat_number_clustering.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')
    
    stats_nb_annotation = pd.read_table(
        os.path.join(results_path,'2_annotation/stat_number_annotation.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')

    df_stats = stats_nb_proteomes.join([stats_nb_clustering, stats_nb_annotation])
    df_stats = df_stats.dropna()

    # Load input
    input_data = pd.read_table(
        input_table,
        sep='\t',
        header=0,
        index_col=None,
        usecols=['observation_name', 'taxonomic_affiliation'])
    
    # Number of input and output taxa (some are discarded along the pipeline if no results)
    n_out = stats_nb_annotation.shape[0]
    n_in = input_data.shape[0]

    # Load tax ids (outputs)
    proteome_tax_id = pd.read_csv(
        os.path.join(results_path, '0_proteomes/proteome_tax_id.tsv'),
        header=0,
        index_col='observation_name',
        sep='\t')

    proteome_tax_id = proteome_tax_id[proteome_tax_id.index.isin(df_stats.index)]
    proteome_tax_id['n_proteomes'] = proteome_tax_id['proteome'].str.split(pat=',').str.len()

    # Input taxonomy must be formatted ; ete3 is used to find the lineage
    output_taxo = _format_taxo(proteome_tax_id)
    df_stats = df_stats.join(output_taxo)

    input_data.set_index('observation_name', inplace=True)
    discarded = input_data[~input_data.index.isin(df_stats.index)]
    n_discarded = n_in - n_out

    # split discarded taxo str
    discarded = discarded['taxonomic_affiliation'].str.split(";", expand = True)

    # Association tax ids (inputs)
    with open(os.path.join(results_path, '0_proteomes/association_taxon_taxID.json')) as f:
        association_taxon_tax_id = json.load(f)

    data = {
        "INPUT_DATA": input_data,
        "DISCARDED": discarded, 
        "N_DISCARDED": n_discarded, 
        "DF_STATS": df_stats, 
        "N_IN": n_in, 
        "N_OUT": n_out, 
        "PROTEOME_TAX_ID": proteome_tax_id, 
        "ASSOCIATION_PROTEOME_TAX_ID": association_taxon_tax_id, 
    }

    return data

# ==============================
# Proteomes summary page figures
# ==============================

def distributions_by_ranks(df_stats, output_folder, n_bins=15, rank='phylum', savefig=True):
    '''
    For the given rank, plots stacked barplots of how many proteomes belongs to each taxa of this rank.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            output_folder (str) : path to output folder to save fig
            rank (str): taxonomic rank for color code. Must match the input. Default 'phylum'
            savefig (bool): boolean to save fig or not

        Returns:
            fig : a plotly html figure
    '''
    df_stats.sort_values(by=rank, ascending=False, inplace=True)

    fig = px.histogram(
        df_stats,
        x="Number_proteomes",
        color=rank,
        height=600,
        width=500,
        nbins=n_bins,
        color_discrete_sequence=COLORSCALE,
        title=f"Distribution of the number of proteomes found by {rank}")
    
    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(
        xaxis_title="Number of proteomes",
        yaxis_title="Number of taxa", 
        plot_bgcolor=PLOT_BGCOLOR,
        font=FONTPARAMS)
   
    if savefig:
        fig.write_image(os.path.join(output_folder, 'proteomes_distribution_by_ranks.pdf'))
        write_json(
            fig, 
            os.path.join(output_folder, 'proteomes_distribution_by_ranks.json'), 
            pretty=True)

    return fig

def n_prot_ec_go_correlations(df_stats, results_path, rank="phylum", savefig=True):
    '''
    Plots the correlation between the number of GO terms and the number of EC numbers. Each dot is a taxa. Dot sizes are linked to the number of proteomes found.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            results_path (str) : the path of the esmecata run
            rank (str): taxonomic rank for color code. Must match the input. Default 'phylum'

        Returns:
            fig : a plotly figure
    '''

    # TODO : add legend for markers size + colors ?
    sizes = df_stats['Number_proteomes']**0.75+5
    fig = px.scatter(
        data_frame=df_stats, 
        x="Number_ecs",
        y="Number_go_terms",
        size=sizes,
        hover_data=["Number_ecs", "Number_go_terms", "Number_proteomes"],
        color=rank,
        #color_discrete_sequence=COLORSCALE,
        width=800,
        height=600,
        labels={"Number_ecs": "Number of ECs", "Number_go_terms": "Number of GO terms", "Number_proteomes": "Number of proteomes"},
        title="Correlation between the number of EC, GO terms, and proteomes")
    
    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(
        plot_bgcolor=PLOT_BGCOLOR,
        font=FONTPARAMS)
    
    if savefig:
        fig.write_image(os.path.join(results_path, "3_analysis/proteomes_figures/correlation_ec_go_proteomes.pdf"))
        write_json(
            fig, 
            os.path.join(results_path, "3_analysis/proteomes_figures/correlation_ec_go_proteomes.json"), 
            pretty=True)

    return fig

def taxo_ranks_contribution(proteome_tax_id, results_path, n_bins=10, savefig=True):
    '''
    For the given rank, plots stacked barplots of how many proteomes belongs to each taxa of this rank.

        Parameters:
            proteome_tax_id (pandas): a pandas df, output of post_analysis_config()
            n_in (int): number of input taxa for esmecata
            n_out (int): number of output taxa by esmecata
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly figure
    '''

    fig = px.histogram(
        proteome_tax_id,
        x="n_proteomes",
        color="tax_rank",
        height=600,
        width=500,
        nbins=n_bins,
        # labels={"n_proteomes": "Number of proteomes"},
        title="Number of proteomes per taxonomic rank") 

    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(
        xaxis_title="Number of proteomes",
        yaxis_title="Number of inputs", 
        plot_bgcolor=PLOT_BGCOLOR,
        font=FONTPARAMS)

    if savefig:
        fig.write_image(os.path.join(results_path, "3_analysis/proteomes_figures/taxonomic_ranks_contribution.pdf"))
        write_json(
            fig, 
            os.path.join(results_path, "3_analysis/proteomes_figures/taxonomic_ranks_contribution.json"), 
            pretty=True)

    return fig

def compare_ranks_in_out(proteome_tax_id, association_taxon_tax_id, output_folder, savefig=True):
    '''
    Plots a heatmap displaying how many input taxonomic ranks were assigned to the same rank by esmecata and how many were assigned to higher tanks.

        Parameters:
            proteome_tax_id (pandas): a pandas df loaded with esmecata 'proteome_tax_id.tsv' file, output of post_analysis_config()
            association_taxon_tax_id (dict) : a dict of taxa loaded from json by post_analysis_config()
            output_folder (str) : path to output folder

        Returns:
            fig : a plotly figure
    '''
    rank_list = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom', 'clade', 'no rank']

    d_in = dict()
    ncbi = NCBITaxa()

    # Loop on each input (=each key of 0_proteomes/association_taxon_taxID.json)
    for taxa in association_taxon_tax_id.keys():
        # Get the list of rank names, ex: ['Gammaproteobacteria', 'Proteobacteria', 'Bacteria', 'cellular organisms']
        reversed_affiliation_taxa = list(reversed(list(association_taxon_tax_id[taxa].keys()))) # That line is f**** up

        i = 0
        found = False

        # Loop until the lowest known taxa rank is found
        # TODO : Discard taxa with only unknown
        while not found and i < len(reversed_affiliation_taxa):
            if "clade" not in reversed_affiliation_taxa[i] and association_taxon_tax_id[taxa][reversed_affiliation_taxa[i]] != ['not_found']:
                tax_id = association_taxon_tax_id[taxa][reversed_affiliation_taxa[i]]                
                tax_rank = ncbi.get_rank(tax_id) # get_rank takes a list for parameter (here, list of one elem, returns a list of one elem)
                # d_in[taxa] = out_asso_rank[list(tax_rank.values())[0]]
                d_in[taxa] = list(tax_rank.values())[0]
                found = True
            else:
                i += 1        
    
    d_out = dict()
    for ASV in proteome_tax_id.index:
        # rank = out_asso_rank[proteome_tax_id.loc[ASV, 'tax_rank']]
        rank = proteome_tax_id.loc[ASV, 'tax_rank']
        d_out[ASV] = rank

    matrix = pd.DataFrame(index=rank_list, columns=rank_list, dtype=int)
    matrix = matrix.fillna(int(0))    

    for ASV in d_out.keys():
        in_rank = d_in[ASV]
        out_rank = d_out[ASV]
        matrix.loc[in_rank, out_rank] += 1
    
    matrix = matrix.loc[(matrix != 0).any(axis=1)]
    matrix = matrix.loc[:,(matrix != 0).any(axis=0)]
    matrix = matrix.iloc[::-1] # Reverse rows for a more suited viz

    fig = px.imshow(
        matrix,
        text_auto=True,
        height=600,
        width=800,
        labels=dict(x="EsMeCaTa model rank", y="Taxa input rank", color="Number of proteomes"),
        title="Difference between input ranks and EsMeCaTa ranks",
        color_continuous_scale='YlGn')
    
    for i in range(0, matrix.shape[1]-1):
        fig.add_vline(x=i+0.5, line_width=2, line_color=LINECOLOR)
    for i in range(0, matrix.shape[0]-1):
        fig.add_hline(y=i+0.5, line_width=2, line_color=LINECOLOR)

    _format_axes(fig)
    fig.update_layout(font=FONTPARAMS)
    
    if savefig:
        fig.write_image(os.path.join(output_folder, 'input_and_output_ranks.pdf'))
        write_json(
            fig, 
            os.path.join(output_folder, 'input_and_output_ranks.json'), 
            pretty=True)
    
    return fig

# ============================================
# EC numbers and GO terms distribution figures
# ============================================

def create_annot_obs_df(dataset_annotation_file, outpath, content):
    ''' Prepare data of 'frequences of EC numbers / Go terms in taxa' and 'fraction of all EC numbers / Go terms in taxa'

    Args:
        dataset_annotation_file (str) : path to an EC*taxa or GO*taxa matrix, computed by create_dataset_annotation_file()
        outpath (str) : the path of the esmecata run
        content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

    Returns:
        data (dict) : stores two dataframes corresponding to each data, as well as summary numbers
    '''

    df = pd.read_csv(
        dataset_annotation_file,
        sep='\t',
        header=0,
        index_col=0)

    # Frequences of EC/GO terms in taxa
    # ---------------------------------

    # EC/GO*taxa matrix rather than taxa*EC/GO
    dfT = df.T

    # Get data
    dfT["frequency"] = dfT.ne(0).sum(axis=1)/dfT.shape[1]    
    dfT["annot_name"] = dfT.index

    if content == 'EC numbers':
            
        ec_class = {"1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases", 
                "4": "Lyases", "5": "Isomerases", "6": "Ligases", "7": "Translocases"}
        
        dfT["EC_class"] = [ec_class[idx[0]] for idx in dfT.index]

        # "class_color" is unused, kept in case a correct way to apply color ranks with update_traces() is found
        dfT["class_color"] = [EC_CLASSES_COLORS[ecclass] for ecclass in dfT["EC_class"]]
        dfT = dfT[["annot_name", "frequency", "EC_class", "class_color"]]
    elif content == 'GO terms':
        dfT = dfT[["annot_name", "frequency"]]

    dfT.sort_values(by=["frequency", "annot_name"], ascending=False, inplace=True)    

    # Count number of EC/GO in each category (left queue, middle, right queue)
    n_9_df = dfT[dfT["frequency"] > 0.9]
    n_1_df = dfT[dfT["frequency"] < 0.1]
    n91_df = dfT[dfT["frequency"].between(0.1,0.9)]
    n_9_an = n_9_df.shape[0]
    n_1_an = n_1_df.shape[0]
    n91_an = n91_df.shape[0]

    # Save lists
    # with open(os.path.join(outpath, "3_analysis", "ecs_last_decile_freqs"), "w") as f:
    #     for idx in n_9_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "ecs_first_decile_freqs"), "w") as f:
    #     for idx in n_1_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "ecs_inbetween_deciles_freqs"), "w") as f:
    #     for idx in n91_df.index:
    #         f.write(f"{idx}\n")

    # Fraction of all EC/GO per taxa
    # ------------------------------

    df["fraction"] = df.ne(0).sum(axis=1)/df.shape[1]
    df["taxa"] = df.index
    df.sort_values(by="fraction", ascending=False, inplace=True)
    df = df[["taxa", "fraction"]]

    # Get number of taxa in each category
    n_9_df = df[df["fraction"] > 0.9]
    n_1_df = df[df["fraction"] < 0.1]
    n91_df = df[df["fraction"].between(0.1,0.9)]
    n_9_ob = n_9_df.shape[0]
    n_1_ob = n_1_df.shape[0]
    n91_ob = n91_df.shape[0]

    # Save lists
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_last_decile_freqs"), "w") as f:
    #     for idx in n_9_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_first_decile_freqs"), "w") as f:
    #     for idx in n_1_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_inbetween_deciles_freqs"), "w") as f:
    #     for idx in n91_df.index:
    #         f.write(f"{idx}\n")

    data = {
        "df_fractionin_obs": df, 
        "df_annot_frequencies": dfT, 
        "n_9_an": n_9_an, 
        "n_1_an": n_1_an, 
        "n91_an": n91_an, 
        "n_9_ob": n_9_ob, 
        "n_1_ob": n_1_ob, 
        "n91_ob": n91_ob
    }

    return data

def annot_frequencies_in_obs(df, output_folder, content, savefig=True):
    '''
    Plots the frequency of each EC in the taxa (i.e. which EC are in many taxa, and which EC are in few taxa)

        Parameters:
            df (pandas): a pandas df with ECs names, frequencies, and classes as columns (key 'df_ec_frequencies' in the output of create_ec_obs_df())
            output_folder (str) : the path of the output folder
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")
    
    labels={"annot_name": content, "frequency": "Frequence in modelled taxa"}
    title=f"{content} ordered by their frequency in modelled taxa"

    if content == "EC numbers":
        fig = px.scatter(
            df,
            x = "annot_name",
            y = "frequency",
            color="EC_class",
            height=500,
            category_orders={"annot_name": df.index},
            labels=labels,
            title=title)
    else:
        fig = px.scatter(
            df,
            x = "annot_name",
            y = "frequency",
            height=500,
            category_orders={"annot_name": df.index},
            labels=labels,
            title=title)

    _format_axes(fig)
    fig.add_hline(y=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_hline(y=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_hrect(y0=0.9, y1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(
        yaxis_title=f"Proportion of modelled taxa containing the {content}", 
        plot_bgcolor=PLOT_BGCOLOR, 
        yaxis_range=[-0.1,1.1],
        font=FONTPARAMS)

    if savefig:
        fig.write_html(os.path.join(output_folder, f"{content.replace(' ', '_')}_frequencies_in_taxa.html"))
        write_json(
            fig, 
            os.path.join(output_folder, f"{content.replace(' ', '_')}_frequencies_in_taxa.json"), 
            pretty=True)

    return fig

def fraction_of_all_annot_in_obs(df, df_taxo, output_folder, content, rank="phylum", savefig=True):
    '''
    Plots the fraction of all EC of the dataset in each taxa (i.e. which taxa have many ECs, and which taxa have few ECs)

        Parameters:
            df (pandas) : a pandas df with taxa names, and its fraction of EC as columns (key 'df_ec_fraction' in the output of create_ec_obs_df())
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of taxa (key 'DF_STATS' in the output of swf.post_analysis_config())
            results_path (str) : path to the output folder
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")

    df = df.join(df_taxo)
    #df.sort_values(by=rank, ascending=False, inplace=True)

    fig = px.scatter(
        df,
        x = "taxa",
        y = "fraction",
        color=rank,
        height=500,
        category_orders={"taxa": df.index},
        labels={"taxa": "Modelled taxa", "fraction": f"Proportion = ({content} in taxa) / (All {content})"},
        title=f"Modelled taxa richness in {content}")

    _format_axes(fig)
    fig.add_hline(y=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_hline(y=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_hrect(y0=0.9, y1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(#yaxis_title=f"Proportion = ({content} in taxa) / (All {content})", 
                      plot_bgcolor=PLOT_BGCOLOR, 
                      yaxis_range=[-0.1,1.1],
                      font=FONTPARAMS)

    if savefig:
        fig.write_html(os.path.join(output_folder, f"{content.replace(' ', '_')}_fraction_per_taxa.html"))
        write_json(
            fig, 
            os.path.join(output_folder, f"{content.replace(' ', '_')}_fraction_per_taxa.json"), 
            pretty=True)

    return fig

def annot_frequencies_in_obs_hist(df, output_folder, content, n_bins=20, savefig=True):
    '''
    For the given rank, plots stacked barplots of how many EC numbers belongs to each taxa of this rank (i.e. which EC are in many taxa, and which EC are in few taxa).

        Parameters:
            df (pandas): a pandas df with EC/GO as rows, taxa as column, filled with the frequencies of EC/GO in taxa
            output_folder (str) : path to the output folder
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")
    
    labels={"annot_name": content},
    title=f"{content} frequencies in modelled taxa"

    if content == "EC numbers":
        fig = px.histogram(
            df,
            x = "frequency",
            color="EC_class",
            height=500,
            width=500,
            nbins=n_bins,
            labels=labels,
            title=title)
    else:
        fig = px.histogram(
            df,
            x = "frequency",
            height=500,
            width=500,
            nbins=20,
            labels=labels,
            title=title)

    _format_axes(fig)
    _format_trace(fig)
    fig.add_vline(x=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_vline(x=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_vrect(x0=0.9, x1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(
        yaxis_title=f"Number of {content}", 
        xaxis_title=f"Proportion of modelled taxa containing the {content}", 
        plot_bgcolor=PLOT_BGCOLOR, 
        xaxis_range=[-0.1,1.1],
        font=FONTPARAMS)
    
    if savefig:
        fig.write_image(os.path.join(output_folder, f"{content.replace(' ', '_')}_frequencies_in_taxa_hist.pdf"))
        write_json(
            fig, 
            os.path.join(output_folder, f"{content.replace(' ', '_')}_frequencies_in_taxa_hist.json"), 
            pretty=True)

    return fig

def fraction_of_all_annot_in_obs_hist(df, df_taxo, output_folder, content, n_bins=10, rank="phylum", savefig=True):
    '''
    Plots the histogram of the fraction of all EC of the dataset in each taxa (i.e. which taxa have many ECs, and which taxa have few ECs)

        Parameters:
            df (pandas): a pandas df with ECs as rows, taxa as column, filled with the frequencies of ECs in the taxa's clusters
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of taxa (key 'DF_STATS' in the output of swf.post_analysis_config())
            output_folder (str) : path to the output folder
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")

    df = df.join(df_taxo)
    df.sort_values(by=rank, ascending=False, inplace=True)
    
    fig = px.histogram(
        df,
        x = "fraction",
        color=rank,
        color_discrete_sequence=COLORSCALE,
        height=500,
        width=500,
        nbins=n_bins,
        title=f"Distribution of modelled taxa richness in {content}")

    _format_axes(fig)
    _format_trace(fig)
    fig.add_vline(x=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_vline(x=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_vrect(x0=0.9, x1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(
        xaxis_title=f"Proportion = ({content} in taxa) / (All {content})",
        yaxis_title="Number of modelled taxa", 
        plot_bgcolor=PLOT_BGCOLOR, 
        xaxis_range=[-0.1,1.1],
        font=FONTPARAMS)
    
    if savefig:
        fig.write_image(os.path.join(output_folder, f"{content.replace(' ', '_')}_fraction_per_taxa_hist.pdf"))
        write_json(
            fig, 
            os.path.join(output_folder, f"{content.replace(' ', '_')}_fraction_per_taxa_hist.json"), 
            pretty=True)

    return fig

# ==========================
# EC sunburst summary figure
# ==========================

def ec_sunburst(ec_classes, output_folder, savefig=True):
    '''
    Calls ontosunburst to create a figure summarizing EC numbers categories, counts and proportions
    for the whole dataset

        Parameters:
            ec_classes (list): a list of EC numbers
            output_folder (str): path to the output folder
            savefig (bool): if the figure should be saved (in html format). True by default

        Returns:
            fig (plotly) : a plotly figure object 
    '''
    fig = ontosunburst(interest_set=ec_classes,
        ontology='ec',
        output=os.path.join(output_folder, "ec_classes_sunburst"),
        root_cut="total")
    
    fig.update_layout(
        title="EC numbers categories, counts and proportions", 
        height=1000,
        paper_bgcolor="#ffffff", 
        font_color='#111111', 
        font_size=20)

    fig.update_traces(marker=dict(colorscale=px.colors.diverging.RdYlGn, line_color=LINECOLOR, reversescale=True))

    if savefig:
        fig.write_html(os.path.join(output_folder, "ec_classes_sunburst.html"))
        write_json(
            fig, 
            os.path.join(output_folder, "ec_classes_sunburst.json"), 
            pretty=True)

    return fig


def ec_sunburst_per_model(ec_classes, output_folder, taxgroup, savefig=True):
    '''
    Calls ontosunburst to create a figure summarizing EC numbers categories, counts and proportions
    for a model in the dataset

    Parameters
        ec_classes (list) : a list of EC numbers
        output_folder (str): path to the output folder
        taxgroup (int) : an ete3 taxonomic group id
        savefig (bool): if the figure should be saved (in svg format). True by default

    Returns:
        fig (plotly) : a plotly figure object
    '''

    title = f"{taxgroup}: EC numbers categories, counts and proportions"
    figpath = os.path.join(output_folder, f"{taxgroup}_ec_classes_sunburst")
    jsonpath = os.path.join(output_folder, f"{taxgroup}_ec_classes_sunburst.json")

    fig = ontosunburst(interest_set=ec_classes,
        ontology='ec',
        output=figpath,
        root_cut="total")
    
    fig.update_layout(
        title=title, 
        height=1000,
        width=1000,
        paper_bgcolor="#ffffff", 
        font_color='#111111', 
        font_size=20)

    fig.update_traces(marker=dict(colorscale=px.colors.diverging.RdYlGn, line_color=LINECOLOR, reversescale=True))

    if savefig:
        # fig.write_html(figpath)
        fig.write_image(figpath + '.svg')
        write_json(
            fig, 
            os.path.join(jsonpath), 
            pretty=True)

    return fig

# =========================
# Clustering summary figure
# =========================

def _hex_to_rgb(hex, opacity):
    '''
    Converts an hexadecimal color code to its corresponding rgba value

        Parameters:
            hex (str): an hexadecimal color code ('#' included)
            opacity (float): the opacity value, between 0 (transparent) to 1 (solid)

        Returns:
            rgba (str): a rgba color code : 'rgba(x,y,z,o)'    
    '''
    nohashtag = hex[1:]

    rgba = [int(nohashtag[i:i+2], 16) for i in (0, 2, 4)] + [opacity]
    rgba = [str(i) for i in rgba]
    rgba = f"rgba({','.join(rgba)})"

    return rgba

def data_proteome_representativeness(proteome_tax_id, computed_threshold_folder):
    ''' 
    Prepare data for plotting summary of the clustering step (i.e. representativeness ratio), from the computed threshold folder and the proteomes_tax_id file.

        Parameters:
            pandas_dict (dict): dict of pandas dataframe, with pandas dataframe of proteomes taxa id (key "PROTEOME_TAX_ID" of post_analysis_config() output)
            computed_threshold_folder (str): pathname to computed threshold folder.

        Returns:
            statsdf (pandas) : a pandas dataframe summarizing the mean, max, and min number of protein clusters for each clustering trhesholdand each input taxa
    '''
    # From proteomes_tax_id get the taxonomic rank for each observation name.
    tax_rank = proteome_tax_id['tax_rank'].to_dict()

    tax_obs_name = {key: value for key, value in proteome_tax_id['tax_id_name'].to_dict().items()}

    data = []
    for obs_name in tax_rank:
        taxon_name = tax_obs_name[obs_name]
        tsv_file_path = os.path.join(computed_threshold_folder, taxon_name+'.tsv')
        obs_tsv_file_path = os.path.join(computed_threshold_folder, obs_name+'.tsv')
        # Add a condition for run of esmecata with version under 0.4.0: in computed_threshold folder, files were named according to observation name.
        if not os.path.exists(tsv_file_path) and os.path.exists(obs_tsv_file_path):
            tsv_file_path = obs_tsv_file_path

        full_name = f"Input name : {obs_name}; EsMeCaTa name :{proteome_tax_id.loc[obs_name,'name']}"# More precision for hover info if needed

        tmp_df = pd.read_csv(tsv_file_path, sep='\t')
        for tmp_threshold in np.arange(0, 1.01, 0.025):
            # Compute the number of protein clusters associated with representativeness ratio of tmp_threshold.
            nb_protein_cluster_ratio = len(tmp_df[tmp_df['cluster_ratio']>= tmp_threshold])
            data.append([obs_name, full_name, tax_rank[obs_name], tmp_threshold, nb_protein_cluster_ratio])

    df = pd.DataFrame(data, columns=['obs_name', 'full_name', 'rank', 'clust', 'count'])

    return df

def create_proteome_representativeness_lineplot_px(df, clust_threshold, output_folder, savefig=True):
    '''     
    This figure shows the representativeness ratio and the number of associated protein clusters according to the taxonomic rank.

        Parameters:
            df (pandas) : a dataframe summarizing proteomes by cluster and by threshold, output of data_proteome_representativeness()
            clust_threshold (float) : the threshold between 0 and 1 fixed by the user for the clustering process of EsMeCaTa
            output_folder (str): pathname to the output figure file.

        Returns:
            fig : a plotly html figure
    '''
    #statsdf = df.groupby(['rank'])['count'].agg(['mean', 'count', 'std'])
    #statsdfs_detail = df.groupby(['obs_name'])

    # Compute confidence interval
    # h = []
    # l = []
    # for rank, row in statsdf.iterrows():
    #     ltmp, htmp = st.t.interval(0.95, row["count"]-1, loc=row["mean"], scale=row["std"])
    #     h.append(htmp)
    #     l.append(ltmp)    
    # statsdf["ciho"] = h
    # statsdf["cilo"] = l

    fig = go.Figure()

    for rank in df["rank"].unique(): # df["rank"].unique() or statsdf.index
        rank_df = df.loc[df["rank"] == rank]
        rank_df = rank_df.groupby(["clust"])["count"].agg(["mean", "min", "max"])

        # tmp_df["lo"] = tmp_df["count"] - statsdf.loc[rank, "cilo"]
        # tmp_df["hi"] = tmp_df["count"] + statsdf.loc[rank, "ciho"]
        fig.add_traces([
            go.Scatter(
                name=f"{rank} mean" ,
                x=rank_df.index,
                y=rank_df["mean"],
                mode="lines+markers",
                line=dict(color=RANK2COL[rank]),
                hovertemplate = "Threshold: %{x}; Clusters: %{y}"),
            go.Scatter(
                name=f"{rank} min",
                x=rank_df.index,
                y=rank_df["min"],
                mode="lines",
                line=dict(width=0, color=RANK2COL[rank]),
                hoverinfo='skip',
                showlegend=False),
            go.Scatter(
                name=f"{rank} max",
                x=rank_df.index,
                y=rank_df["max"],
                mode="lines",
                fill='tonexty',
                fillcolor=_hex_to_rgb(RANK2COL[rank], 0.2),
                line=dict(width=0, color=RANK2COL[rank]),
                hoverinfo='skip',
                showlegend=False)
        ])

    fig.add_vline(
        x=clust_threshold, 
        line_dash="dash", 
        line_width=1,
        line_color="red", 
        annotation_text="EsMeCaTa's clustering threshold")

    _format_axes(fig)
    fig.update_xaxes(range=[0, 1])  
    fig.update_yaxes(type="log")    

    fig.update_layout(
        title="Count of retained clusters of proteins (split by by taxonomic rank) according to the clustering threshold",
        xaxis_title="Clustering threshold (0 : pan-proteome, 1 : core-proteome)", 
        yaxis_title="Number of clusters of proteomes retained", 
        height=800, 
        plot_bgcolor=PLOT_BGCOLOR,
        legend=dict(x=1, y=1, xanchor='right', yanchor='top'),
        font=FONTPARAMS)

    if savefig:
        fig.write_image(os.path.join(output_folder, 'proteome_representativeness.pdf'))
        write_json(
            fig, 
            os.path.join(output_folder, 'proteome_representativeness.json'), 
            pretty=True)

    return fig

def proteomes_representativeness_details(df, clust_threshold, output_folder, savefig=True):
    '''
    This figure shows details of the representativeness ratio and the number of associated protein clusters according to each cluster. There is one sub-pnale by taxonomic rank.

        Parameters:
            df (pandas) : a dataframe summarizing proteomes by cluster and by threshold, output of data_proteome_representativeness().
            clust_threshold (float) : the threshold between 0 and 1 fixed by the user for the clustering process of EsMeCaTa.
            output_folder (str): pathname to the output figure file.
            savefig (bool): boolean to create output files.

        Returns:
            fig : a plotly html figure
    '''

    ncbi_ranks = ["clade", "superkingdom","kingdom","subkingdom","superphylum","phylum",
        "subphylum","infraphylum","superclass","class","subclass",
        "infraclass","cohort","subcohort","superorder","order",
        "suborder","infraorder","parvorder","superfamily","family",
        "subfamily","tribe","subtribe","genus","subgenus",
        "section","subsection","series","subseries","species group",
        "species subgroup","species","forma specialis","subspecies","varietas",
        "subvariety","forma","serogroup","serotype","strain","isolate"]
    
    ranks = df["rank"].unique()
    ncbi_ranks = [rank for rank in ncbi_ranks if rank in ranks]

    titles=[f"Details for {rank}" for rank in ncbi_ranks]

    n = df["rank"].nunique()
    fig = make_subplots(
        rows=n, 
        cols=1,
        shared_xaxes=False,
        vertical_spacing=0.075,
        subplot_titles=titles)  

    i = 1
    # for rank in df["rank"].unique():
    for rank in ncbi_ranks:

        rank_df = df.loc[df["rank"] == rank]        

        for input_taxa in rank_df["obs_name"].unique():
            rank_df_taxa = rank_df.loc[df["obs_name"] == input_taxa]

            fname = rank_df_taxa["full_name"].iloc[0]

            fig.add_trace(
                go.Scatter(
                    name=fname, # input_taxa
                    x=rank_df_taxa["clust"],
                    y=rank_df_taxa["count"],
                    mode="lines+markers",
                    line=dict(color=RANK2COL[rank]),
                    showlegend=False,
                    hovertemplate = "Threshold: %{x}; Clusters: %{y}"
                ),
                row=i,
                col=1
            )
        i += 1

    fig.add_vline(x=clust_threshold, 
        line_dash="dash", 
        line_width=1, 
        line_color="red", 
        annotation_text="EsMeCaTa's clustering threshold")
        
    _format_axes(fig)    
    fig.update_yaxes(type="log")
    fig.update_xaxes(range=[0, 1])
        
    fig.update_layout(#title=f"Details for {rank} clusters",
        xaxis_title="Clustering threshold (0 : pan-proteome, 1 : core-proteome)", 
        yaxis_title="Number of clusters of proteomes retained", 
        height=n*400,
        width=500,
        plot_bgcolor=PLOT_BGCOLOR,
        font=FONTPARAMS)
        
    i = 1
    for rank in df["rank"].unique():
        fig.update_xaxes(title_text="Clustering threshold (0 : pan-proteome, 1 : core-proteome)", row=i, col=1)
        fig.update_yaxes(title_text="Number of clusters of proteomes retained", row=i, col=1)
        i += 1
        
    if savefig:
        fig.write_html(os.path.join(output_folder, 'proteome_representativeness_details.html'))
        write_json(
            fig, 
            os.path.join(output_folder, 'proteome_representativeness_details.json'),
            pretty=True)

    return fig

# ======================
# Metadata Page function
# ======================

def reproducibility_tokens(outdir):
    '''
    Loads json metadata of all the EsmeCata workflow and converts it to a string.

        Parameters:
            outdir (str) : path to the results of an EsMeCaTa workflow

        Returns:
            metadata (str) : all metadata of the workflow
    '''
    reprometa_proteomes = json.load(open(os.path.join(outdir, '0_proteomes', 'esmecata_metadata_proteomes.json'), 'r'))
    reprometa_clustering = json.load(open(os.path.join(outdir, '1_clustering', 'esmecata_metadata_clustering.json'), 'r'))
    reprometa_annotation = json.load(open(os.path.join(outdir, '2_annotation', 'esmecata_metadata_annotation.json'), 'r'))

    metadata = {"proteomes": reprometa_proteomes,
         "clustering": reprometa_clustering,
         "annotation": reprometa_annotation}
    
    # metadata = json.dumps(metadata, indent=4)

    return metadata
