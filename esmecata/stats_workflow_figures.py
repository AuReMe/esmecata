#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import warnings
import os
import json
import pandas as pd
from ete3 import NCBITaxa
import plotly.express as px

__author__ = "Pauline Hamon-Giraud, Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

# warnings.filterwarnings('ignore')

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

    # Get ECs frequencies
    ecs_matrix = pd.read_table(os.path.join(results_path,'3_analysis/dataset_annotation_normalised.tsv'),
        sep='\t',
        header=0,
        index_col=0)

    data = {"INPUT_DATA": input_data,
            "DISCARDED": discarded, 
            "N_DISCARDED": n_discarded, 
            "DF_STATS": df_stats, 
            "N_IN": n_in, 
            "N_OUT": n_out, 
            "PROTEOME_TAX_ID": proteome_tax_id, 
            "ASSOCIATION_PROTEOME_TAX_ID": association_taxon_tax_id, 
            "ECS_MATRIX": ecs_matrix}
    
    # return input_data, discarded, n_discarded, df_stats, n_in, n_out, proteome_tax_id, association_taxon_tax_id, ecs_matrix
    return data

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
    
    fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_layout(yaxis_title="Number of observations", plot_bgcolor='#e6ffe6')
   
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
                     #size="Number_proteomes",
                     size=sizes,
                     hover_data=["Number_ecs", "Number_go_terms", "Number_proteomes"],
                     color=rank,
                     width=800,
                     height=600,
                     labels={"Number_ecs": "Number of ECs", "Number_go_terms": "Number of GO terms", "Number_proteomes": "Number of proteomes"},
                     title="Correlation between the number of EC, GO terms, and proteomes")
    
    fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_layout(plot_bgcolor='#e6ffe6') 
    
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
        height=500,
        width=500,
        labels={"n_proteomes": "Number of proteomes"},
        title="Number of proteomes per taxonomic rank") 

    fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_layout(yaxis_title="Number of inputs", plot_bgcolor='#e6ffe6') 

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

    for observation_name in association_taxon_tax_id.keys():
        reversed_affiliation_taxa = list(reversed(list(association_taxon_tax_id[observation_name].keys()))) # That line is f**** up

        i = 0
        found = False

        # loop until the lowest known taxa rank is found
        # TODO : Discard taxa with only unknown
        # TODO : gérer les not found des clade (marche pas)
        while not found and i < len(reversed_affiliation_taxa):
            if association_taxon_tax_id[observation_name][reversed_affiliation_taxa[i]] != ['not_found']:
                tax_id = association_taxon_tax_id[observation_name][reversed_affiliation_taxa[i]]
                tax_rank = ncbi.get_rank(tax_id)
                # d_in[observation_name] = out_asso_rank[list(tax_rank.values())[0]]
                d_in[observation_name] = list(tax_rank.values())[0]
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

    fig = px.imshow(matrix,
        text_auto=True,
        height=500,
        width=500,
        labels=dict(x="EsMeCaTa rank", y="Input rank", color="Number of proteomes"),
        title="Difference between input ranks and EsMeCaTa ranks",
        color_continuous_scale='YlGn')
    
    for i in range(0, matrix.shape[1]-1):
        fig.add_vline(x=i+0.5, line_width=2, line_color="DarkSlateGrey")
    for i in range(0, matrix.shape[0]-1):
        fig.add_hline(y=i+0.5, line_width=2, line_color="DarkSlateGrey")

    # fig.add_vline(x=1+0.5, line_width=2, line_color="DarkSlateGrey")
    
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True)
    
    fig.write_html(os.path.join(results_path, "3_analysis/input_and_output_ranks.html"))
    
    return fig

# ===============================
# EC numbers distribution figures
# ===============================

def create_ec_obs_df(dataset_annotation_file, results_path):
    """ Prepare data of 'frequences of EC numbers in observations' and 'fraction of all EC numbers in observations'

    Args:
        dataset_annotation_file (str) : path to an EC*Observations matrix, computed by create_dataset_annotation_file()
        results_path (str) : the path of the esmecata run

    Returns:
        data (dict) : stores two dataframes corresponding to each data, as well as summary numbers
    """

    df = pd.read_csv(dataset_annotation_file,
        sep='\t',
        header=0,
        index_col=0)

    # Frequences of ECs in observations
    # ---------------------------------

    # Used for color code
    ec_class = {"1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases", 
                "4": "Lyases", "5": "Isomerases", "6": "Ligases", "7": "Translocases"}

    # ECs*observations matrix rather than observations*ECs
    dfT = df.T

    # Get data
    dfT["frequency"] = dfT.ne(0).sum(axis=1)/dfT.shape[1]
    dfT["EC_class"] = [ec_class[idx[0]] for idx in dfT.index]
    dfT["EC_name"] = dfT.index
    dfT.sort_values(by=["frequency", "EC_name"], ascending=False, inplace=True)
    dfT = dfT[["EC_name", "frequency", "EC_class"]]

    # Count number of ECs in each category (left queue, middle, right queue)
    n_9_df = dfT[dfT["frequency"] > 0.9]
    n_1_df = dfT[dfT["frequency"] < 0.1]
    n91_df = dfT[dfT["frequency"].between(0.1,0.9)]
    n_9_ec = n_9_df.shape[0]
    n_1_ec = n_1_df.shape[0]
    n91_ec = n91_df.shape[0]

    # Save EC lists
    # with open(os.path.join(outpath, "3_analysis", "ecs_last_decile_freqs"), "w") as f:
    #     for idx in n_9_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "ecs_first_decile_freqs"), "w") as f:
    #     for idx in n_1_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "ecs_inbetween_deciles_freqs"), "w") as f:
    #     for idx in n91_df.index:
    #         f.write(f"{idx}\n")

    # Fraction of all ECs per observation
    # -----------------------------------

    df["fraction"] = df.ne(0).sum(axis=1)/df.shape[1]
    df["observation"] = df.index
    df.sort_values(by="fraction", ascending=False, inplace=True)
    df = df[["observation", "fraction"]]

    # Get number of observations in each category
    n_9_df = df[df["fraction"] > 0.9]
    n_1_df = df[df["fraction"] < 0.1]
    n91_df = df[df["fraction"].between(0.1,0.9)]
    n_9_ob = n_9_df.shape[0]
    n_1_ob = n_1_df.shape[0]
    n91_ob = n91_df.shape[0]

    # Save EC lists
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_last_decile_freqs"), "w") as f:
    #     for idx in n_9_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_first_decile_freqs"), "w") as f:
    #     for idx in n_1_df.index:
    #         f.write(f"{idx}\n")
    # with open(os.path.join(outpath, "3_analysis", "otu_ec_props_inbetween_deciles_freqs"), "w") as f:
    #     for idx in n91_df.index:
    #         f.write(f"{idx}\n")

    data = {"df_fractionin_obs": df, "df_ec_frequencies": dfT, "n_9_ec": n_9_ec, "n_1_ec": n_1_ec, "n91_ec": n91_ec, "n_9_ob": n_9_ob, "n_1_ob": n_1_ob, "n91_ob": n91_ob}

    return data

def ecs_frequencies_in_obs(df, results_path):
    '''
    Plots the frequency of each EC in the observations (i.e. which EC are in many observations, and which EC are in few observations)

        Parameters:
            df (pandas): a pandas df with ECs names, frequencies, and classes as columns (key 'df_ec_frequencies' in the output of create_ec_obs_df())
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    # Plots
    fig = px.scatter(df,
        x = "EC_name",
        y = "frequency",
        color="EC_class",
        height=500,
        category_orders={"EC_name": df.index},
        labels={"EC_name": "ECs", "frequency": "Frequence in observations"},
        title="ECs ordered by their frequency in observations")   

    # fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.add_hline(y=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_hline(y=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_hrect(y0=0.9, y1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title="Proportion of observations containing the EC", plot_bgcolor='#e6ffe6', yaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, "3_analysis/ecs_frequencies_in_observations.html"))

    return fig

def fraction_of_all_ec_in_obs(df, df_taxo, results_path, rank="phylum"):
    '''
    Plots the fraction of all EC of the dataset in each observation (i.e. which observations have many ECs, and which observations have few ECs)

        Parameters:
            df (pandas) : a pandas df with observations names, and its fraction of EC as columns (key 'df_ec_fraction' in the output of create_ec_obs_df())
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of observations (key 'DF_STATS' in the output of swf.post_analysis_config())
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    df = df.join(df_taxo)

    # Plots
    fig = px.scatter(df,
        x = "observation",
        y = "fraction",
        color=rank,
        height=500,
        category_orders={"observation": df.index},
        labels={"observation": "Observation", "fraction": "Frequence in observations"},
        title="Observations ordered by the proportion (Number of EC numbers in an observation) / (Total number of EC numbers in the dataset)")

    # fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.add_hline(y=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_hline(y=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_hrect(y0=0.9, y1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title="Proportion of present Ecs/total ECs", plot_bgcolor='#e6ffe6', yaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, "3_analysis/ec_fraction_per_observation.html"))

    return fig

def ecs_frequencies_in_obs_hist(df, results_path):
    '''
    For the given rank, plots stacked barplots of how many proteomes / GO terms / EC numbers belongs to each taxa of this rank.

        Parameters:
            df (pandas): a pandas df with ECs as rows, observations as column, filled with the frequencies of ECs in the observation's clusters
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    # Plots
    fig = px.histogram(df,
        x = "frequency",
        color="EC_class",
        height=500,
        width=500,
        nbins=20,
        # hover_data=[normalized_dataset_annotation.index],
        labels={"EC_name": "ECs", "frequency": "Proportion of observations containing the EC"},
        title="Histogram of the frequencies of EC numbers in observations")   

    # fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.add_vline(x=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_vline(x=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_vrect(x0=0.9, x1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title="Number of ECs", plot_bgcolor='#e6ffe6', xaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, "3_analysis/ecs_frequencies_in_observations_hist.html"))

    return fig

def fraction_of_all_ec_in_obs_hist(df, df_taxo, results_path, rank="phylum"):
    '''
    plots the frequencies of each OTU in each EC

        Parameters:
            df (pandas): a pandas df with ECs as rows, observations as column, filled with the frequencies of ECs in the observation's clusters
            df_taxo (pandas) : a pandas dataframe containing the taxonomy of observations (key 'DF_STATS' in the output of swf.post_analysis_config())
            results_path (str) : the path of the esmecata run

        Returns:
            fig : a plotly html figure
    '''

    df = df.join(df_taxo)
    
    # Plots
    fig = px.histogram(df,
        x = "fraction",
        color=rank,
        height=500,
        width=500,
        nbins=10,
        # hover_data=[normalized_dataset_annotation.index],
        labels={"observation": "Observation", "fraction": "Proportion = (ECs in observation) / (All ECs)"},
        title="Histogram of the proportion (Number of EC numbers in an observation) / (Total number of EC numbers in the dataset)")

    # fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')))
    fig.update_xaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='DarkSlateGrey', mirror=True, showgrid=True, gridwidth=2)
    fig.add_vline(x=0.9, line_dash="dash", line_width=1, line_color="red")
    fig.add_vline(x=0.1, line_dash="dash", line_width=1, line_color="red")
    fig.add_vrect(x0=0.9, x1=0.1, line_width=0, fillcolor="red", opacity=0.15)

    fig.update_layout(yaxis_title="Number of observations", plot_bgcolor='#e6ffe6', xaxis_range=[-0.1,1.1])
    fig.write_html(os.path.join(results_path, "3_analysis/ec_fraction_per_observation_hist.html"))

    return fig
