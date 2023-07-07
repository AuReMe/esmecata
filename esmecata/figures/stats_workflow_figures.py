#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os

__author__ = "Pauline Hamon-Giraud, Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

warnings.filterwarnings('ignore')
# nb : Python 3.9.13 was used

def post_analysis_config(input_table, results_path):

    # Load and concat summaries of results
    stats_nb_proteomes = pd.read_table(os.path.join(results_path,'0_proteomes/stat_number_proteome.tsv'),
        sep='\t',
        index_col=0)
    
    stats_nb_clustering = pd.read_table(os.path.join(results_path,'1_clustering/stat_number_clustering.tsv'),
        sep='\t',
        index_col=0)
    
    stats_nb_annotation = pd.read_table(os.path.join(results_path,'2_annotation/stat_number_annotation.tsv'),
        sep='\t',
        index_col=0)
    
    stats = stats_nb_proteomes.join([stats_nb_clustering, stats_nb_annotation])
    stats = stats.dropna()

    # Get number of input and output taxa (some are discarded along the pipeline if no results)
    n_in = stats_nb_proteomes.shape[0]
    n_out = stats.shape[0]

    # Load input (index col : 'observation_name')
    input_data = pd.read_table(input_table,
        sep='\t',
        index_col=0)
    input_data = input_data[input_data.index.isin(stats.index)]

    df_stats = stats.join(input_data)

    # Load tax ids
    proteome_tax_id = pd.read_csv(os.path.join(results_path, 'proteome_tax_id.tsv'),
                              index_col=0,
                              sep='\t')
    
    proteome_tax_id = proteome_tax_id[proteome_tax_id.index.isin(stats.index)]
    proteome_tax_id['n_proteomes'] = proteome_tax_id['proteome'].str.split(pat=',').str.len()

    return input_data, df_stats, n_in, n_out, proteome_tax_id

def distributions_by_ranks(df_stats, n_in, n_out, outpath, rank='Phylum'):
    '''
    For the given rank, plots stacked barplots of how many proteomes / GO terms / EC numbers belongs to each taxa of this rank.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            n_in (int): number of input taxa for esmecata
            n_out (int): number of output taxa by esmecata
            rank (str): taxonomic rank for color code. Must match the input. Default Phylum
            outpath (str) : the path of the output figure
    '''
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    sns.histplot(ax=axes[0], 
                 data=df_stats, 
                 x="Number_proteomes", 
                 bins=15, 
                 multiple="stack", 
                 hue=rank).set(xlabel='Number of proteomes', 
                                   ylabel='Number of taxa')

    sns.histplot(ax=axes[1], 
                 data=df_stats, 
                 x="Number_go_terms", 
                 bins=15, 
                 multiple="stack", 
                 hue=rank).set(xlabel='Number of GO terms')

    sns.histplot(ax=axes[2], 
                 data=df_stats, 
                 x="Number_ecs", 
                 bins=15, 
                 multiple="stack", 
                 hue=rank).set(xlabel='Number of EC numbers')

    fig.suptitle(f'Distributions by {rank} {n_out} taxa (out of {n_in})')
    plt.savefig(os.path.join(outpath, 'distribution_by_ranks.png'))

def n_prot_ec_go_correlations(df_stats, n_in, n_out, outpath):
    '''
    Plots the correlation between the number of GO terms and the number of EC numbers. Each dot is a taxa. Dot sizes are linked to the number of proteomes found.

        Parameters:
            df_stats (pandas): a pandas df, output of post_analysis_config()
            n_in (int): number of input taxa for esmecata
            n_out (int): number of output taxa by esmecata
            outpath (str) : the path of the output figure
    '''
    plt.figure()

    sns.scatterplot(data=df_stats, 
                    x="Number_ecs", 
                    y="Number_go_terms", 
                    size="Number_proteomes", 
                    alpha=0.5, 
                    sizes=(20, 400)).set(title=f'{n_out} taxa (out of {n_in})',
                                         xlabel = 'Number of EC numbers',
                                         ylabel = 'Number of GO terms')

    plt.savefig(os.path.join(outpath, 'n_prot_ec_go_correlations.png'))

def taxo_ranks_contribution(proteome_tax_id, n_in, n_out, outpath):
    '''
    For the given rank, plots stacked barplots of how many proteomes / GO terms / EC numbers belongs to each taxa of this rank.

        Parameters:
            proteome_tax_id (pandas): a pandas df, output of post_analysis_config()
            n_in (int): number of input taxa for esmecata
            n_out (int): number of output taxa by esmecata
            outpath (str) : the path of the output figure
    '''

    plt.figure()
    sns.histplot(data=proteome_tax_id, 
                 x="n_proteomes", 
                 bins=15, 
                 multiple="stack", 
                 hue="tax_rank").set(title=f'Part of each taxonomic rank {n_out} taxa (out of {n_in})',
                                     xlabel='Number of proteomes', 
                                     ylabel='Number of taxa')

    plt.savefig(os.path.join(outpath, 'taxo_ranks_contribution.png'))

def compare_ranks_in_out(input_data, proteome_tax_id, outpath):
    '''
    Plots a heatmap displaying how many input taxonomic ranks were assigned to the same rank by esmecata and how many were assigned to higher tanks.

        Parameters:
            input_data (pandas): a pandas df loaded with esmecata input data, output of post_analysis_config()
            proteome_tax_id (pandas): a pandas df loaded with esmecata 'proteome_tax_id.tsv' file, output of post_analysis_config()
            outpath (str) : the path of the output figure
    '''

    out_asso_rank = {'superkingdom': 'k',
                 'kingdom': 'k',
                 'phylum': 'p',
                 'class': 'c',
                 'order': 'o',
                 'family': 'f',
                 'clade': 'f',
                 'genus': 'g',
                 'species': 's',
                 'no rank': 'na'}

    rank_list = ['s', 'g', 'f', 'o', 'c', 'p', 'k', 'na']

    d_in = dict()
    with open(input_data, 'r') as f:
        rows = list(csv.reader(f, delimiter='\t'))
        for row in rows[1:]:
            tax = row[0]
            ranks = list(reversed(row[1].split(';')))
            for i in range(len(rank_list) - 1):
                if ranks[i] != 'unknown':
                    d_in[tax] = rank_list[i]
                    break

    d_out = dict()
    for ASV in proteome_tax_id.index:
        rank = out_asso_rank[proteome_tax_id.loc[ASV, 'tax_rank']]
        d_out[ASV] = rank

    matrix = pd.DataFrame(index=rank_list, columns=rank_list, dtype=int)
    matrix = matrix.fillna(int(0))

    for ASV in d_out.keys():
        in_rank = d_in[ASV]
        out_rank = d_out[ASV]
        matrix.loc[in_rank, out_rank] += 1

    plt.figure()
    sns.heatmap(matrix, annot=True, fmt='.0f', linewidth=0.5, cmap="crest")
    plt.xlabel('Esmecata Rank')
    plt.ylabel('Input Rank')
    plt.savefig(os.path.join(outpath, "compare_ranks_in_out.png"))
