# Copyright (C) 2021-2025 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import json
import logging
import os
import time

from esmecata.core.proteomes import retrieve_proteomes, compute_stat_proteomes
from esmecata.core.clustering import make_clustering, compute_stat_clustering
from esmecata.core.annotation import annotate_proteins, compute_stat_annotation
from esmecata.core.eggnog import annotate_with_eggnog

logger = logging.getLogger(__name__)


def compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file=None):
    """Compute stat associated with the number of proteome for each taxonomic affiliations.

    Args:
        proteomes_output_folder (str): pathname to the result folder containing subfolder containing proteomes
        clustering_output_folder (str): pathname to the result folder containing mmseqs results
        annotation_output_folder (str): pathname to the annotation reference folder containing annotations for each cluster
        stat_file (str): pathname to the tsv stat file

    Returns:
        workflow_numbers (dict): dict containing observation names (as key) associated with proteomes, protein clusters, GO Terms and EC (as value)
    """
    workflow_numbers = {}

    proteome_numbers = compute_stat_proteomes(proteomes_output_folder)

    clustering_numbers = compute_stat_clustering(clustering_output_folder)

    annotation_reference_folder = os.path.join(annotation_output_folder, 'annotation_reference')
    annotation_numbers = compute_stat_annotation(annotation_reference_folder)

    all_observation_names = {*proteome_numbers.keys(), *clustering_numbers.keys(), *annotation_numbers.keys()}
    for observation_name in all_observation_names:
        if observation_name in proteome_numbers:
            nb_proteomes, tax_name, lowest_tax_rank, esmecata_name, esmecata_rank, reference_status = proteome_numbers[observation_name]
        else:
            nb_proteomes = 'NA'
        if observation_name in clustering_numbers:
            nb_shared_proteins = len(clustering_numbers[observation_name][1])
        else:
            nb_shared_proteins = 'NA'
        if observation_name in annotation_numbers:
            nb_gos = annotation_numbers[observation_name][0]
            nb_ecs = annotation_numbers[observation_name][1]
        else:
            nb_gos = 'NA'
            nb_ecs = 'NA'
        workflow_numbers[observation_name] = [nb_proteomes, tax_name, lowest_tax_rank, esmecata_name, esmecata_rank, reference_status, nb_shared_proteins, nb_gos, nb_ecs]

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Input_taxon_name', 'Input_rank', 'EsMeCaTa_chosen_name', 'EsMeCaTa_chosen_rank','Number_proteomes', 'Only_reference_proteome', 'Proteins_clusters_0_5' , 'Number_go_terms', 'Number_ecs'])
            for observation_name in workflow_numbers:
                nb_proteomes, tax_name, lowest_tax_rank, esmecata_name, esmecata_rank, reference_status = workflow_numbers[observation_name][:6]
                nb_shared_proteins = workflow_numbers[observation_name][6]
                nb_gos = workflow_numbers[observation_name][7]
                nb_ecs = workflow_numbers[observation_name][8]
                csvwriter.writerow([observation_name, tax_name, lowest_tax_rank, esmecata_name, esmecata_rank, nb_proteomes, reference_status, nb_shared_proteins, nb_gos, nb_ecs])

    return workflow_numbers


def perform_workflow(input_file, output_folder, busco_percentage_keep=80, ignore_taxadb_update=None,
                        all_proteomes=None, uniprot_sparql_endpoint=None, remove_tmp=None,
                        limit_maximal_number_proteomes=99, rank_limit=None,
                        nb_core=1, clust_threshold=1, mmseqs_options=None,
                        linclust=None, propagate_annotation=None, uniref_annotation=None,
                        expression_annotation=None, minimal_number_proteomes=1, annotation_files=None,
                        update_affiliations=None, option_bioservices=None):
    """From the proteomes found by esmecata proteomes, create protein cluster for each taxonomic affiliations.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations
        output_folder (str): pathname to the output folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        ignore_taxadb_update (bool): option to ignore ete4 taxa database update
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        remove_tmp (bool): remove the tmp files
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks)
        nb_core (int): number of CPU-cores to be used by mmseqs
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        mmseqs_options (str): use alternative mmseqs option
        linclust (bool): use linclust
        propagate_annotation (float): float between 0 and 1. It is the ratio of proteins in the cluster that should have the annotation to keep this annotation.
        uniref_annotation (bool): option to use uniref annotation to add annotation.
        expression_annotation (bool): option to add expression annotation from UniProt.
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taxon to be kept.
        annotation_files (str): pathnames to UniProt dat files.
        update_affiliations (str): option to update taxonomic affiliations.
        option_bioservices (bool): use bioservices instead of manual queries.
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|workflow| Begin workflow.')
    workflow_metadata = {}

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteomes_output_folder = os.path.join(output_folder, '0_proteomes')
    retrieve_proteomes(input_file, proteomes_output_folder, busco_percentage_keep,
                        ignore_taxadb_update, all_proteomes, uniprot_sparql_endpoint,
                        limit_maximal_number_proteomes, rank_limit, minimal_number_proteomes,
                        update_affiliations, option_bioservices)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    make_clustering(proteomes_output_folder, clustering_output_folder, nb_core, clust_threshold, mmseqs_options, linclust, remove_tmp)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    annotate_proteins(clustering_output_folder, annotation_output_folder, uniprot_sparql_endpoint, propagate_annotation,
                      uniref_annotation, expression_annotation, annotation_files, option_bioservices)

    stat_file = os.path.join(output_folder, 'stat_number_workflow.tsv')
    compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file)

    proteomes_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(proteomes_metadata_file, 'r') as json_idata:
        proteomes_metadata = json.load(json_idata)
    clustering_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'r') as json_idata:
        clustering_metadata = json.load(json_idata)
    annotation_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation.json')
    with open(annotation_metadata_file, 'r') as json_idata:
        annotation_metadata = json.load(json_idata)

    workflow_metadata['proteomes_metadata'] = proteomes_metadata
    workflow_metadata['clustering_metadata'] = clustering_metadata
    workflow_metadata['annotation_metadata'] = annotation_metadata

    endtime = time.time()
    duration = endtime - starttime
    workflow_metadata['esmecata_workflow_duration'] = duration
    workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow.json')
    if os.path.exists(workflow_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_workflow' in metadata_file]
        workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow_{0}.json'.format(len(metadata_files)))
        with open(workflow_metadata_file, 'w') as ouput_file:
            json.dump(workflow_metadata, ouput_file, indent=4)
    else:
        with open(workflow_metadata_file, 'w') as ouput_file:
            json.dump(workflow_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|workflow| Workflow step complete in {0}s.'.format(duration))



def perform_workflow_eggnog(input_file, output_folder, eggnog_database_path, busco_percentage_keep=80,
                            ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                            remove_tmp=None, limit_maximal_number_proteomes=99, rank_limit=None,
                            nb_core=1, clust_threshold=0.5, mmseqs_options=None,
                            linclust=None, minimal_number_proteomes=5, update_affiliations=None,
                            option_bioservices=None, eggnog_tmp_dir=None, no_dbmem=False):
    """From the proteomes found by esmecata proteomes, create protein cluster for each taxonomic affiliations.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations
        output_folder (str): pathname to the output folder
        eggnog_database_path (str): pathname to eggnog database folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        ignore_taxadb_update (bool): option to ignore ete4 taxa database update
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        remove_tmp (bool): remove the tmp files
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks)
        nb_core (int): number of CPU-cores to be used by mmseqs and eggnog-mapper
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        mmseqs_options (str): use alternative mmseqs option
        linclust (bool): use linclust
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taxon to be kept.
        update_affiliations (str): option to update taxonomic affiliations.
        option_bioservices (bool): use bioservices instead of manual queries.
        eggnog_tmp_dir (str): pathname to eggnog-mapper temporary folder.
        no_dbmem (bool): Boolean to choose to not load eggnog database in memory.
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|workflow| Begin workflow.')
    workflow_metadata = {}

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteomes_output_folder = os.path.join(output_folder, '0_proteomes')
    retrieve_proteomes(input_file, proteomes_output_folder, busco_percentage_keep,
                        ignore_taxadb_update, all_proteomes, uniprot_sparql_endpoint,
                        limit_maximal_number_proteomes, rank_limit, minimal_number_proteomes,
                        update_affiliations, option_bioservices)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    make_clustering(proteomes_output_folder, clustering_output_folder, nb_core, clust_threshold, mmseqs_options, linclust, remove_tmp)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    annotate_with_eggnog(clustering_output_folder, annotation_output_folder, eggnog_database_path, nb_core, eggnog_tmp_dir, no_dbmem)

    stat_file = os.path.join(output_folder, 'stat_number_workflow.tsv')
    compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file)

    proteomes_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(proteomes_metadata_file, 'r') as json_idata:
        proteomes_metadata = json.load(json_idata)
    clustering_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'r') as json_idata:
        clustering_metadata = json.load(json_idata)
    annotation_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation.json')
    with open(annotation_metadata_file, 'r') as json_idata:
        annotation_metadata = json.load(json_idata)

    workflow_metadata['proteomes_metadata'] = proteomes_metadata
    workflow_metadata['clustering_metadata'] = clustering_metadata
    workflow_metadata['annotation_metadata'] = annotation_metadata

    endtime = time.time()
    duration = endtime - starttime
    workflow_metadata['esmecata_workflow_duration'] = duration
    workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow.json')
    if os.path.exists(workflow_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_workflow' in metadata_file]
        workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow_{0}.json'.format(len(metadata_files)))
        with open(workflow_metadata_file, 'w') as ouput_file:
            json.dump(workflow_metadata, ouput_file, indent=4)
    else:
        with open(workflow_metadata_file, 'w') as ouput_file:
            json.dump(workflow_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|workflow| Workflow step complete in {0}s.'.format(duration))
