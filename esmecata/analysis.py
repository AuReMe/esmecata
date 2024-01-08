# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
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

import csv
import os
import pandas as pd
import logging
import seaborn as sns
import matplotlib.pyplot as plt

from collections import Counter
from ete3 import NCBITaxa

logger = logging.getLogger(__name__)


def create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_file_path):
    """ From annotation reference folder, creates a file resuming the number of ECs for each observation names.

    Args:
        annotation_reference_folder (str): path to annotation reference folder
        dataset_annotation_file_path (str): path to output dataset annotation file

    Returns:
        dataset_annotation (dict): annotation dict: observation_name as key and EC number as value
    """
    dataset_annotation = {}
    total_annotations = []
    for annotation_file in os.listdir(annotation_reference_folder):
        annotation_file_name = os.path.splitext(annotation_file)[0]
        annotation_file_path = os.path.join(annotation_reference_folder, annotation_file)

        annotations = []
        with open(annotation_file_path, 'r') as open_annotation_input_file_path:
            csvreader = csv.DictReader(open_annotation_input_file_path, delimiter='\t')
            for line in csvreader:
                #gos = line['GO'].split(',')
                ecs = line['EC'].split(',')
                #annotations.extend(gos)
                annotations.extend(ecs)
        annotations = [annot for annot in annotations if annot != '']
        total_annotations.extend(annotations)
        dataset_annotation[annotation_file_name] = annotations

    total_annotations = list(set(total_annotations))

    with open(dataset_annotation_file_path, 'w') as dataset_annotation_file:
        csvwriter = csv.writer(dataset_annotation_file, delimiter='\t')
        csvwriter.writerow(['observation_name'] + total_annotations)
        for observation_name in dataset_annotation:
            occurrence_annotations = Counter(dataset_annotation[observation_name])
            observation_name_annotations = [occurrence_annotations[annot] for annot in total_annotations]
            csvwriter.writerow([observation_name] + observation_name_annotations)

    return dataset_annotation


def get_taxon_obs_name(proteome_tax_id_file, selected_taxon_rank='family'):
    """ From proteome tax id file, reads it and extract the taxonomic name associated with the observation name.

    Args:
        proteome_tax_id_file (str): path to the proteome fax id file
        selected_taxon_rank (str): taxonomic rank selected

    Returns:
        taxa_name (dict): annotation dict: taxon_name as key and observation_name as value
    """
    ncbi = NCBITaxa()
    taxa_name = {}
    total_obs_name = []
    selected_obs_name = []
    with open(proteome_tax_id_file, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            total_obs_name.append(observation_name)

            taxon_name = line['name']
            tax_id = line['tax_id']
            tax_rank = line['tax_rank']

            found_taxon_name = None

            if tax_rank == selected_taxon_rank:
                found_taxon_name = taxon_name
            else:
                tax_id_lineages = ncbi.get_lineage(tax_id)
                for tax_id in tax_id_lineages:
                    if ncbi.get_rank([tax_id])[tax_id] == selected_taxon_rank:
                        tmp_lineage_taxon_name = ncbi.get_taxid_translator([tax_id])[tax_id]
                        found_taxon_name = tmp_lineage_taxon_name

            if found_taxon_name is not None:
                if found_taxon_name not in taxa_name:
                    taxa_name[found_taxon_name] = [observation_name]
                    selected_obs_name.append(observation_name)
                else:
                    taxa_name[found_taxon_name].append(observation_name)
                    selected_obs_name.append(observation_name)

    selected_obs_name = set(selected_obs_name)

    logger.info('|EsMeCaTa|analysis| {0} of {1} could been associated with taxon rank {2}.'.format(len(selected_obs_name), len(total_obs_name), selected_taxon_rank))

    return taxa_name


def normalise_dataframe(input_dataframe, taxa_name, nb_digit):
    """ Reduce number of ECs (by limiting the number of digit kept) and taxon by summing their occurrences.
    Then normalise the number of occurrence of EC in each taxon by the number of occurrences of the EC in all taxa.

    Args:
        input_dataframe (pandas DataFrame): annotation dataframe of esmecata
        taxa_name (dict): annotation dict: taxon_name as key and observation_name as value
        nb_digit (int): number of digits to keep for visualisation (by default 3)

    Returns:
        input_dataframe (pandas DataFrame): reduced and normalised dataframe on EC number
    """
    # Get all EC at nb_digit-digit (either 1, 2, 3 or 4 digits).
    ec_numbers = sorted(list(set(['.'.join(ec_number.split('.')[:nb_digit]) for ec_number in input_dataframe.columns])))

    # Sum the occurrence of the 4-digit ECs associated with a EC nb_digit-digit.
    for ec_number in ec_numbers:
        input_dataframe[ec_number] = input_dataframe[[col for col in input_dataframe.columns if col.startswith(ec_number)]].sum(axis=1)
    input_dataframe = input_dataframe[ec_numbers]

    input_dataframe = input_dataframe.transpose()

    # Sum the occurrence of an EC for a taxon by the occurrence of this EC in the organisms of the taxon.
    for taxon_name in taxa_name:
        input_dataframe[taxon_name] = input_dataframe[taxa_name[taxon_name]].sum(axis=1)#.div(len(taxa_name[taxon_name]))
    input_dataframe = input_dataframe[list(taxa_name.keys())]

    # Remove all rows equal to 0 (with the selection of taxonomic rank, it is possible to have EC associated with no selected ranks).
    input_dataframe = input_dataframe.loc[~(input_dataframe==0).all(axis=1)]

    # Divide each value by the sum of the columns.
    # Here, the number of each EC for a taxon is divided by the number of EC in all the taxa.
    input_dataframe = input_dataframe.div(input_dataframe.sum(axis=1), axis=0)

    # Keep only EC with variance higher than the median.
    #input_dataframe = input_dataframe[input_dataframe.T.var() > input_dataframe.T.var().median()]

    return input_dataframe 


def create_visualisation_ec(dataset_annotation_file_path, taxa_name, normalised_dataset_annotation_file_path, output_figure, nb_digit=3):
    """ From the dataset annotation file, creates a clustermap.

    Args:
        dataset_annotation_file_path (str): path to the dataset annotation file
        taxa_name (dict): annotation dict: taxon_name as key and observation_name as value
        normalised_dataset_annotation_file_path (str): path to the output normalised dataset annotation file
        output_figure (str): path to the output clustermap file
        nb_digit (int): number of digits to keep for visualisation (by default 3)
    """
    if nb_digit not in [1, 2, 3, 4]:
        logger.info('|EsMeCaTa|analysis| EC digit for visualisation must be 1, 2, 3 or 4, not {nb_digit}.')

    df = pd.read_csv(dataset_annotation_file_path, sep='\t')

    df.set_index('observation_name', inplace=True)

    normalised_df = normalise_dataframe(df, taxa_name, nb_digit)

    normalised_df.to_csv(normalised_dataset_annotation_file_path, sep='\t')

    nb_obs_names = len(list(normalised_df.columns))
    if nb_obs_names > 1:
        if nb_digit in [1, 2, 3]:
            fig_size = (20, 30)
            linewidth = 1
        else:
            nb_ec = len(normalised_df)
            fig_length = nb_ec * 0.05
            fig_size = (20, fig_length)
            linewidth = 0.2
            sns.set(font_scale=0.35)

        sns.clustermap(normalised_df, row_cluster=True, col_cluster=True, yticklabels=True, figsize=fig_size, mask=(normalised_df==0),
                        linewidths=linewidth, linecolor='black', cmap=sns.cm.rocket_r)

        plt.savefig(output_figure)
    else:
        logger.info('|EsMeCaTa|analysis| Not enough observation names for clustermap: {0}.'.format(','.join(list(normalised_df.columns))))


def create_ec_clustermap(annotation_reference_folder, proteome_tax_id_file, output_folder, taxon_rank='family', nb_digit=3):
    """Create a clustermap showing the occurrences of EC for a dataset.
    Reduce the number of taxon and EC.

    Args:
        annotation_reference_folder (str): path to annotation reference folder
        proteome_tax_id_file (str): path to the proteome fax id file
        output_folder (str): path to output folder
        taxon_rank (str): taxonomic rank selected
        nb_digit (int): number of digits to keep for visualisation (by default 3)
    """
    # Retrive taxon name associated to a specific taxon rank to reduce numer of organisms/taxa to show in clustermap.
    taxa_name = get_taxon_obs_name(proteome_tax_id_file, taxon_rank)

    # Extract EC occurrences from all the annotation files.
    dataset_annotation_file_path = os.path.join(output_folder, 'dataset_annotation.tsv')
    create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_file_path)

    # Normalise the number of EC and create a visualisation.
    normalised_dataset_annotation_file_path = os.path.join(output_folder, 'dataset_annotation_normalised.tsv')
    output_figure = os.path.join(output_folder, 'clustermap.pdf')
    create_visualisation_ec(dataset_annotation_file_path, taxa_name, normalised_dataset_annotation_file_path, output_figure, nb_digit)


def perform_analysis(annotation_input_folder, output_folder, taxon_rank='family', nb_digit=3):
    """ Create analysis file (EC clustermap).

    Args:
        input_folder (str): pathname to mmseqs clustering output folder as it is the input of annotation
        output_folder (str): pathname to the output folder
        taxon_rank (str): taxonomic rank selected
        nb_digit (int): number of digits to keep for visualisation (by default 3)
    """
    annotation_reference_folder = os.path.join(annotation_input_folder, 'annotation_reference')
    proteome_tax_id_file = os.path.join(annotation_input_folder, 'proteome_tax_id.tsv')
    create_ec_clustermap(annotation_reference_folder, proteome_tax_id_file, output_folder, taxon_rank, nb_digit)