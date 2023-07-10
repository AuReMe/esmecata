#!/usr/bin/env python
# -*- coding: utf-8 -*-

import seaborn as sns
import warnings
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
from ete3 import NCBITaxa

__author__ = "Pauline Hamon-Giraud, Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

warnings.filterwarnings('ignore')
# nb : Python 3.9.13 was used

def post_analysis_config(input_table, results_path):
    '''
    Proceeds to format the inputs and outputs data in order to un all the others plotting functions efficiently

        Parameters:
            input_table (str): the path to the input table of taxa
            results_path (str): the path of the esmecata run

        Returns:
            input_data (pandas) :
            df_stats (pandas) : 
            n_in (int) : 
            n_out (int) : 
            proteome_tax_id (pandas) :
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
    
    stats = stats_nb_proteomes.join([stats_nb_clustering, stats_nb_annotation])
    stats = stats.dropna()

    # Get number of input and output taxa (some are discarded along the pipeline if no results)
    n_in = stats_nb_proteomes.shape[0]
    n_out = stats.shape[0]

    # Load input
    input_data = pd.read_table(input_table,
        sep='\t',
        header=0,
        index_col='observation_name')
    input_data = input_data[input_data.index.isin(stats.index)]

    df_stats = stats.join(input_data)

    # Load tax ids (outputs)
    proteome_tax_id = pd.read_csv(os.path.join(results_path, '0_proteomes/proteome_tax_id.tsv'),
        header=0,
        index_col='observation_name',
        sep='\t')
    
    proteome_tax_id = proteome_tax_id[proteome_tax_id.index.isin(stats.index)]
    proteome_tax_id['n_proteomes'] = proteome_tax_id['proteome'].str.split(pat=',').str.len()

    # Association tax ids (inputs)
    with open(os.path.join(results_path, '0_proteomes/association_taxon_taxID.json')) as f:
        association_taxon_tax_id = json.load(f)

    return input_data, df_stats, n_in, n_out, proteome_tax_id, association_taxon_tax_id

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

def compare_ranks_in_out(input_data, proteome_tax_id, association_taxon_tax_id, outpath):
    '''
    Plots a heatmap displaying how many input taxonomic ranks were assigned to the same rank by esmecata and how many were assigned to higher tanks.

        Parameters:
            input_data (pandas): a pandas df loaded with esmecata input data, output of post_analysis_config()
            proteome_tax_id (pandas): a pandas df loaded with esmecata 'proteome_tax_id.tsv' file, output of post_analysis_config()
            association_taxon_tax_id (dict) : a dict of taxa loaded from json by post_analysis_config()
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
    ncbi = NCBITaxa()

    for observation_name in association_taxon_tax_id.keys():
        reversed_affiliation_taxa = list(reversed(list(association_taxon_tax_id[observation_name].keys()))) # That line is f**** up

        i = 0
        found = False

        # loop until the lowest known taxa rank is found
        # TODO : Discard taxa with only unknown
        while not found and i < len(reversed_affiliation_taxa):
            if association_taxon_tax_id[observation_name][reversed_affiliation_taxa[i]] != ['not_found']:
                tax_id = association_taxon_tax_id[observation_name][reversed_affiliation_taxa[i]]
                tax_rank = ncbi.get_rank(tax_id)
                d_in[observation_name] = out_asso_rank[list(tax_rank.values())[0]]
                found = True
            else:
                i += 1
        
    
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

def esmecata2krona(run_name):

    def _create_krona_input_in(association_taxon_id, name, output):
        '''
        TODO

            Parameters:
                association_taxon_tax_id (dict) : a dict of taxa loaded from json by post_analysis_config()
                name (str) : ...
                output (str) : the path of ...
        '''
        krona_input = os.path.join(output, f'{name}_krona_input.tsv')
        with open(krona_input, 'w') as f:
            f.write('#queryID\t#taxID')

        with open(krona_input, 'a') as f:
            for observation_name in association_taxon_id:
                reversed_affiliation_taxa = list(reversed(association_taxon_id[observation_name]))
                for tax_name in reversed_affiliation_taxa:
                    if tax_name != 'unknown' and association_taxon_id[observation_name][tax_name] != ['not_found']:
                        tax_id = association_taxon_id[observation_name][tax_name][0]
                        f.write(f'\n{observation_name}\t{tax_id}')
                        break

        krona_output = os.path.join(output, f'{name}_krona_input.html')
        os.system(f'ktImportTaxonomy {krona_input} -o {krona_output}')

    def _create_krona_input_esmecata(proteome_tax_id_tsv, name, output):
        krona_input = os.path.join(output, f'{name}_krona_esmecata.tsv')
        with open(krona_input, 'w') as o_f:
            o_f.write('#queryID\t#taxID')

            with open(proteome_tax_id_tsv, 'r') as i_f:
                lines = list(csv.reader(i_f, delimiter='\t'))
                for l in lines[1:]:
                    o_f.write(f'\n{l[0]}\t{l[2]}')

        krona_output = os.path.join(output, f'{name}_krona_esmecata.html')
        os.system(f'ktImportTaxonomy {krona_input} -o {krona_output}')


    def _create_comp_taxonomy_file(association_taxon_id, proteome_tax_id_tsv, name, output):
        d_tax = dict()
        with open(proteome_tax_id_tsv, 'r') as f:
            lines = list(csv.reader(f, delimiter='\t'))
            for l in lines[1:]:
                observation_name = l[0]
                reversed_affiliation_taxa = list(reversed(association_taxon_id[observation_name]))
                for tax_name in reversed_affiliation_taxa:
                    if tax_name != 'unknown' and association_taxon_id[observation_name][tax_name] != ['not_found']:
                        if tax_name not in d_tax.keys():
                            tax_id = association_taxon_id[observation_name][tax_name][0]
                            tax_rank = ncbi.get_rank([tax_id])
                            if tax_id in tax_rank.keys():
                                tax_rank = tax_rank[tax_id]
                            else:
                                tax_rank = ''
                            d_tax[tax_name] = {'Esmecata ID': l[2], 'Esmecata Rank': l[3], 'Esmecata Name': l[1],
                                            'Input ID': tax_id, 'Input Rank': tax_rank}
                        if 'obs' not in d_tax[tax_name].keys():
                            d_tax[tax_name]['obs'] = list()
                        d_tax[tax_name]['obs'].append(observation_name)
                        break

        output = os.path.join(output, f'{name}_taxonomy_diff.tsv')
        with open(output, 'w') as f:
            f.write('\t'.join(['Input Name', 'Esmecata Name', 'Input Rank', 'Esmecata Rank',
                            'Input ID', 'Esmecata ID', 'Observation Names']))
            for name, ref in d_tax.items():
                f.write('\n' + '\t'.join([name, ref['Esmecata Name'], ref['Input Rank'], ref['Esmecata Rank'],
                                        str(ref['Input ID']), str(ref['Esmecata ID']), ';'.join(ref['obs'])]))
                
    INPUT_RANK_FILE = os.path.join(run_name, '0_proteomes', 'association_taxon_taxID.json')
    ESMECATA_RANK_FILE = os.path.join(run_name, '0_proteomes', 'proteome_tax_id.tsv')
    output = os.path.join(run_name, 'Krona')
    
    if not os.path.exists(output):
        os.mkdir(output)

    _create_krona_input_in(INPUT_RANK_FILE, run_name, output)
    _create_krona_input_esmecata(ESMECATA_RANK_FILE, run_name, output)
    _create_comp_taxonomy_file(INPUT_RANK_FILE, ESMECATA_RANK_FILE, run_name, output)

