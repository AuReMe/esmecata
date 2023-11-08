#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import warnings
import os
import json
import numpy as np
import pandas as pd
from ete3 import NCBITaxa
import plotly.express as px
import plotly.graph_objects as go
from  ontosunburst.ontosunburst import ec_ontosunburst
from esmecata_compression import RANK2COL

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
EC_CLASSES_COLORS = {"Ligases": "#636EFA",
                    "Oxidoreductases": "#EF553B",
                     "Transferases": "#00CC96",
                     "Hydrolases": "#AB63FA",
                     "Lyases": "#FFA15A",
                     "Isomerases": "#19D3F3",
                     "Translocases": "#FF6692"} # 'Plotly' default color map

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
    fig.update_xaxes(showline=True, 
                     linewidth=2, 
                     linecolor=LINECOLOR, 
                     mirror=True, 
                     showgrid=True, 
                     gridwidth=2)
    fig.update_yaxes(showline=True, 
                     linewidth=2, 
                     linecolor=LINECOLOR, 
                     mirror=True, 
                     showgrid=True, 
                     gridwidth=2)
    
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
    df.drop("no rank", axis=1, inplace=True)

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
    stats_nb_proteomes = pd.read_table(os.path.join(results_path,'0_proteomes/stat_number_proteome.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')
    
    stats_nb_clustering = pd.read_table(os.path.join(results_path,'1_clustering/stat_number_clustering.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')
    
    stats_nb_annotation = pd.read_table(os.path.join(results_path,'2_annotation/stat_number_annotation.tsv'),
        sep='\t',
        header=0,
        index_col='observation_name')
    
    df_stats = stats_nb_proteomes.join([stats_nb_clustering, stats_nb_annotation])
    df_stats = df_stats.dropna()

    # Load input
    input_data = pd.read_table(input_table,
        sep='\t',
        header=0,
        index_col=None,
        usecols=['observation_name', 'taxonomic_affiliation'])
    
    # Number of input and output taxa (some are discarded along the pipeline if no results)
    n_out = stats_nb_annotation.shape[0]
    n_in = input_data.shape[0]

    # Load tax ids (outputs)
    proteome_tax_id = pd.read_csv(os.path.join(results_path, '0_proteomes/proteome_tax_id.tsv'),
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

    data = {"INPUT_DATA": input_data,
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

def distributions_by_ranks(df_stats, results_path, rank='phylum'):
    '''
    For the given rank, plots stacked barplots of how many proteomes / GO terms / EC numbers belongs to each taxa of this rank.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            rank (str): taxonomic rank for color code. Must match the input. Default 'phylum'
            results_path (str) : the path of the emecata run

        Returns:
            fig : a plotly html figure
    '''

    fig = px.histogram(df_stats,
        x="Number_proteomes",
        color=rank,
        height=600,
        width=500,
        labels={"Number_proteomes": "Number of proteomes"},
        title=f"Distribution of the number of proteomes found by {rank}")
    
    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(yaxis_title="Number of taxa", plot_bgcolor=PLOT_BGCOLOR)
   
    fig.write_html(os.path.join(results_path, "3_analysis/proteomes_distribution_by_ranks.html"))

    return fig

def n_prot_ec_go_correlations(df_stats, results_path, rank="phylum"):
    '''
    Plots the correlation between the number of GO terms and the number of EC numbers. Each dot is a taxa. Dot sizes are linked to the number of proteomes found.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            results_path (str) : the path of the esmecata run
            rank (str): taxonomic rank for color code. Must match the input. Default 'phylum'
    '''

    # TODO : add legend for markers size + colors ?
    sizes = df_stats['Number_proteomes']**0.75+5
    fig = px.scatter(data_frame=df_stats, 
                     x="Number_ecs",
                     y="Number_go_terms",
                     size=sizes,
                     hover_data=["Number_ecs", "Number_go_terms", "Number_proteomes"],
                     color=rank,
                     width=800,
                     height=600,
                     labels={"Number_ecs": "Number of ECs", "Number_go_terms": "Number of GO terms", "Number_proteomes": "Number of proteomes"},
                     title="Correlation between the number of EC, GO terms, and proteomes")
    
    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(plot_bgcolor=PLOT_BGCOLOR) 
    
    fig.write_html(os.path.join(results_path, "3_analysis/correlation_ec_go_proteomes.html"))

    return fig

def taxo_ranks_contribution(proteome_tax_id, results_path):
    '''
    For the given rank, plots stacked barplots of how many proteomes / GO terms / EC numbers belongs to each taxa of this rank.

        Parameters:
            proteome_tax_id (pandas): a pandas df, output of post_analysis_config()
            n_in (int): number of input taxa for esmecata
            n_out (int): number of output taxa by esmecata
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    fig = px.histogram(proteome_tax_id,
        x="n_proteomes",
        color="tax_rank",
        height=600,
        width=500,
        labels={"n_proteomes": "Number of proteomes"},
        title="Number of proteomes per taxonomic rank") 

    _format_trace(fig)
    _format_axes(fig)
    fig.update_layout(yaxis_title="Number of inputs", plot_bgcolor=PLOT_BGCOLOR) 

    fig.write_html(os.path.join(results_path, "3_analysis/taxonomic_ranks_contribution.html"))

    return fig

def compare_ranks_in_out(proteome_tax_id, association_taxon_tax_id, results_path):
    '''
    Plots a heatmap displaying how many input taxonomic ranks were assigned to the same rank by esmecata and how many were assigned to higher tanks.

        Parameters:
            proteome_tax_id (pandas): a pandas df loaded with esmecata 'proteome_tax_id.tsv' file, output of post_analysis_config()
            association_taxon_tax_id (dict) : a dict of taxa loaded from json by post_analysis_config()
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    # out_asso_rank = {'superkingdom': 'k',
    #              'kingdom': 'k',
    #              'phylum': 'p',
    #              'class': 'c',
    #              'order': 'o',
    #              'family': 'f',
    #              'clade': 'f',
    #              'genus': 'g',
    #              'species': 's',
    #              'no rank': 'na'}

    # rank_list = ['s', 'g', 'f', 'o', 'c', 'p', 'k', 'na']
    rank_list = ['species', 'genus', 'clade', 'family', 'order', 'class', 'phylum', 'kingdom', 'no rank']

    d_in = dict()
    ncbi = NCBITaxa()

    for taxa in association_taxon_tax_id.keys():
        reversed_affiliation_taxa = list(reversed(list(association_taxon_tax_id[taxa].keys()))) # That line is f**** up

        i = 0
        found = False

        # loop until the lowest known taxa rank is found
        # TODO : Discard taxa with only unknown
        # TODO : gérer les not found des clade (marche pas)
        while not found and i < len(reversed_affiliation_taxa):
            if association_taxon_tax_id[taxa][reversed_affiliation_taxa[i]] != ['not_found']:
                tax_id = association_taxon_tax_id[taxa][reversed_affiliation_taxa[i]]
                tax_rank = ncbi.get_rank(tax_id)
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

    fig = px.imshow(matrix,
        text_auto=True,
        height=600,
        width=800,
        labels=dict(x="EsMeCaTa rank", y="Input rank", color="Number of proteomes"),
        title="Difference between input ranks and EsMeCaTa ranks",
        color_continuous_scale='YlGn')
    
    for i in range(0, matrix.shape[1]-1):
        fig.add_vline(x=i+0.5, line_width=2, line_color=LINECOLOR)
    for i in range(0, matrix.shape[0]-1):
        fig.add_hline(y=i+0.5, line_width=2, line_color=LINECOLOR)

    _format_axes(fig)
    #fig.update_layout(xaxis={"side": "top"}, title=dict(automargin=True, yref='container'))
    
    fig.write_html(os.path.join(results_path, "3_analysis/input_and_output_ranks.html"))
    
    return fig

# ============================================
# EC numbers and GO terms distribution figures
# ============================================

def create_annot_obs_df(dataset_annotation_file, results_path, content):
    """ Prepare data of 'frequences of EC numbers / Go terms in taxa' and 'fraction of all EC numbers / Go terms in taxa'

    Args:
        dataset_annotation_file (str) : path to an EC*taxa or Go*taxa matrix, computed by create_dataset_annotation_file()
        results_path (str) : the path of the esmecata run
        content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

    Returns:
        data (dict) : stores two dataframes corresponding to each data, as well as summary numbers
    """

    df = pd.read_csv(dataset_annotation_file,
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

    data = {"df_fractionin_obs": df, "df_annot_frequencies": dfT, 
        "n_9_an": n_9_an, "n_1_an": n_1_an, "n91_an": n91_an, 
        "n_9_ob": n_9_ob, "n_1_ob": n_1_ob, "n91_ob": n91_ob}

    return data

def annot_frequencies_in_obs(df, results_path, content):
    '''
    Plots the frequency of each EC in the taxa (i.e. which EC are in many taxa, and which EC are in few taxa)

        Parameters:
            df (pandas): a pandas df with ECs names, frequencies, and classes as columns (key 'df_ec_frequencies' in the output of create_ec_obs_df())
            results_path (str) : the path of the esmecata run
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")
    
    labels={"annot_name": content, "frequency": "Frequence in taxa"}
    title=f"{content} ordered by their frequency in taxa"

    if content == "EC numbers":
        fig = px.scatter(df,
            x = "annot_name",
            y = "frequency",
            color="EC_class",
            height=500,
            category_orders={"annot_name": df.index},
            labels=labels,
            title=title)
    else:
        fig = px.scatter(df,
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

    fig.update_layout(yaxis_title=f"Proportion of taxa containing the {content}", plot_bgcolor=PLOT_BGCOLOR, yaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, f"3_analysis/{content.replace(' ', '_')}_frequencies_in_taxa.html"))

    return fig

def fraction_of_all_annot_in_obs(df, df_taxo, results_path, content, rank="phylum"):
    '''
    Plots the fraction of all EC of the dataset in each taxa (i.e. which taxa have many ECs, and which taxa have few ECs)

        Parameters:
            df (pandas) : a pandas df with taxa names, and its fraction of EC as columns (key 'df_ec_fraction' in the output of create_ec_obs_df())
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of taxa (key 'DF_STATS' in the output of swf.post_analysis_config())
            results_path (str) : the path of the esmecata run
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")

    df = df.join(df_taxo)

    fig = px.scatter(df,
        x = "taxa",
        y = "fraction",
        color=rank,
        height=500,
        category_orders={"taxa": df.index},
        labels={"taxa": "Taxa", "fraction": "Frequence in taxa"},
        title=f"Taxa ordered by the proportion (Number of {content} in an taxa) / (Total number of {content} in the dataset)")

    _format_axes(fig)
    fig.add_hline(y=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_hline(y=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_hrect(y0=0.9, y1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title=f"Proportion of present {content}/total {content}", plot_bgcolor=PLOT_BGCOLOR, yaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, f"3_analysis/{content.replace(' ', '_')}_fraction_per_taxa.html"))

    return fig

def annot_frequencies_in_obs_hist(df, results_path, content):
    '''
    For the given rank, plots stacked barplots of how many EC numbers belongs to each taxa of this rank (i.e. which EC are in many taxa, and which EC are in few taxa).

        Parameters:
            df (pandas): a pandas df with ECs as rows, taxa as column, filled with the frequencies of ECs in the taxa's clusters
            results_path (str) : the path of the esmecata run
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")
    
    labels={"annot_name": content},
    title=f"Histogram of the frequencies of {content} in taxa"

    if content == "EC numbers":
        fig = px.histogram(df,
            x = "frequency",
            color="EC_class",
            height=500,
            width=500,
            nbins=20,
            labels=labels,
            title=title)
    else:
        fig = px.histogram(df,
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

    fig.update_layout(yaxis_title=f"Number of {content}", 
                      xaxis_title=f"Proportion of taxa containing the {content}", 
                      plot_bgcolor=PLOT_BGCOLOR, 
                      xaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, f"3_analysis/{content.replace(' ', '_')}_frequencies_in_taxa_hist.html"))

    return fig

def fraction_of_all_annot_in_obs_hist(df, df_taxo, results_path, content, rank="phylum"):
    '''
    Plots the histogram of the fraction of all EC of the dataset in each taxa (i.e. which taxa have many ECs, and which taxa have few ECs)

        Parameters:
            df (pandas): a pandas df with ECs as rows, taxa as column, filled with the frequencies of ECs in the taxa's clusters
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of taxa (key 'DF_STATS' in the output of swf.post_analysis_config())
            results_path (str) : the path of the esmecata run
            content (str) : 'GO terms' or 'EC numbers' to specifies which data is processed

        Returns:
            fig : a plotly html figure
    '''
    if content not in ["GO terms", "EC numbers"]:
        raise ValueError("param `content` must either be 'GO terms' or 'EC numbers'")

    df = df.join(df_taxo)
    
    fig = px.histogram(df,
        x = "fraction",
        color=rank,
        height=500,
        width=500,
        nbins=10,
        labels={"taxa": "Taxa", "fraction": f"Proportion = ({content} in taxa) / (All {content})"},
        title=f"Histogram of the proportion (Number of {content} in a taxa) / (Total number of {content} in the dataset)")

    _format_axes(fig)
    _format_trace(fig)
    fig.add_vline(x=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_vline(x=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_vrect(x0=0.9, x1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title="Number of taxa", 
                      plot_bgcolor=PLOT_BGCOLOR, 
                      xaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, f"3_analysis/{content.replace(' ', '_')}_fraction_per_taxa_hist.html"))

    return fig

# ==========================
# EC sunburst summary figure
# ==========================

def ec_sunburst(ec_classes, results_path):
    fig = ec_ontosunburst(ec_set=ec_classes, 
                          output=os.path.join(results_path, "3_analysis/taxonomic_ranks_contribution.html"))
    fig.update_layout(title="EC numbers categories, counts and proportions", 
                      height=1000,
                      paper_bgcolor="#ffffff", 
                      font_color='#111111', 
                      font_size=20)

    fig.update_traces(marker=dict(colorscale=px.colors.diverging.RdYlGn, line_color=LINECOLOR, reversescale=True))
    fig.write_html(os.path.join(results_path, "3_analysis/ec_classes_sunburst.html"))

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

def create_proteome_representativeness_lineplot_px(proteome_tax_id, computed_threshold_folder, results_path):
    ''' 
    From the computed threshold folder and the proteomes_tax_id file created a lineplot.
    This figures show the representativeness ratio and the number of associated protein clusters according to the taxonomic rank.

        Parameters:
            proteomes_taxon_id (pandas): pandas dataframe of proteomes taxa id (key "PROTEOME_TAX_ID" of post_analysis_config() output)
            computed_threshold_folder (str): pathname to computed threshold folder.
            results_path (str): pathname to the output figure file.

        Returns:
            fig : a plotly html figure
    '''
    # From proteomes_tax_id get the taxonomic rank for each observation name.
    tax_rank = proteome_tax_id['tax_rank'].to_dict()

    data = []
    for tsv_file in os.listdir(computed_threshold_folder):
        obs_name = os.path.splitext(tsv_file)[0]
        tsv_file_path = os.path.join(computed_threshold_folder, tsv_file)
        tmp_df = pd.read_csv(tsv_file_path, sep='\t')
        for tmp_threshold in np.arange(0, 1.01, 0.025):
            # Compute the number of protein clusters associated with representativeness ratio of tmp_threshold.
            nb_protein_cluster_ratio = len(tmp_df[tmp_df['cluster_ratio']>= tmp_threshold])
            data.append([obs_name, tax_rank[obs_name], tmp_threshold, nb_protein_cluster_ratio])

    df = pd.DataFrame(data, columns=['obs_name', 'rank', 'clust', 'count'])

    statsdf = df.groupby(['rank'])['count'].agg(['mean', 'count', 'std'])

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

    for rank in statsdf.index:
        rank_df = df.loc[df["rank"] == rank]
        rank_df = rank_df.groupby(["clust"])["count"].agg(["mean", "min", "max"])

        # tmp_df["lo"] = tmp_df["count"] - statsdf.loc[rank, "cilo"]
        # tmp_df["hi"] = tmp_df["count"] + statsdf.loc[rank, "ciho"]
        fig.add_traces([
            go.Scatter(name=f"{rank} mean" ,
                       x=rank_df.index,
                       y=rank_df["mean"],
                       mode="lines+markers",
                       line=dict(color=RANK2COL[rank])),                       
            go.Scatter(name=f"{rank} min",
                       x=rank_df.index,
                       y=rank_df["min"],
                       mode="lines",
                       line=dict(width=0),
                       showlegend=False),
            go.Scatter(name=f"{rank} max",
                       x=rank_df.index,
                       y=rank_df["max"],
                       mode="lines",
                       fill='tonexty',
                       fillcolor=_hex_to_rgb(RANK2COL[rank], 0.2),
                       line=dict(width=0),
                       showlegend=False)
        ])

    _format_axes(fig)    
    fig.update_yaxes(type="log")    

    fig.update_layout(title="Count of retained clusters of proteins (split by by taxonomic rank) according to the clustering threshold",
                      xaxis_title="Clustering threshold (0 : pan-proteome, 1 : core-proteome)", 
                      yaxis_title="Number of clusters of proteomes retained", 
                      height=800, 
                      plot_bgcolor=PLOT_BGCOLOR) 

    fig.write_html(os.path.join(results_path, "3_analysis/proteome_representativeness.html"))

    return fig

# ======================
# Metadata Page function
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
    
    metadata = json.dumps(metadata, indent=4)

    return metadata

    #fig.for_each_trace(lambda trace: trace.update(opacity = 0.2) if trace.name in ["max"] else (),)

    # def confidence_interval(data, confidence=0.95):
    #     dist = NormalDist.from_samples(data)
    #     z = NormalDist().inv_cdf((1 + confidence) / 2.)
    #     h = dist.stdev * z / ((len(data) - 1) ** .5)
    #     return h
    
    # CI = confidence_interval(df["count"], 0.95)
    # df["CIlb"] = df["count"] - CI
    # df["CIhb"] = df["count"] + CI


    # fig, ax = plt.subplots(figsize=(9,5))
    # g_results = sns.lineplot(data=df, x="clust", y='count', hue="rank", errorbar='ci')
    # g_results.set(yscale='log')
    # # Use max + 0.25 quantile to set ylim.
    # g_results.set(ylim=[1, df['count'].max()+df['count'].quantile(0.25)], xlim=[0,1])
    # # plt.savefig(os.path.join(results_path, "3_analysis/proteome_representativeness.png"))
    # fig = tls.mpl_to_plotly(fig)
    
    # ci95_hi = []
    # ci95_lo = []
    # for i in statsdf.index:
    #     m, c, s = statsdf.loc[i]
    #     # ci95_hi.append(m + 1.96*s/math.sqrt(c))
    #     # ci95_lo.append(m - 1.96*s/math.sqrt(c))
    #     ci95 = st.norm.interval(confidence=0.95, loc=m, scale=s)
    #     ci95_hi.append(ci95[1])
    #     ci95_lo.append(ci95[0])
    # statsdf['ci95_hi'] = ci95_hi
    # statsdf['ci95_lo'] = ci95_lo

    # print(statsdf)

    # ci = st.norm.interval(confidence=0.95, 
    #     loc=np.mean(df['count']),
    #     scale=st.sem(df['count'])) 

    # df["lower_ci"] = df["count"] - ci[1]
    # df["upper_ci"] = df["count"] + ci[1]

    # fig = go.Figure()
    # fig = px.line(df, x="clust", 
    #     y="count", 
    #     color="rank", 
    #     log_y=True) 
    # #     #error_y=df["upper_ci"], 
    # #     #error_y_minus=df["lower_ci"])
    
    # fig.add_traces([go.Scatter(df, 
    #         x="clust",
    #         y='CIlb',
    #         color="rank",
    #         log_y=True),
    #     px.line(df, 
    #         x="clust", 
    #         y='CIhb', 
    #         color="rank",
    #         log_y=True)])
