# Copyright (C) 2021-2022 Arnaud Belcour - Inria Dyliss
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
import sys
import time

from esmecata.proteomes import retrieve_proteomes, compute_stat_proteomes
from esmecata.clustering import make_clustering, compute_stat_clustering
from esmecata.annotation import annotate_proteins, compute_stat_annotation
from esmecata.kegg_metabolism import create_draft_networks, compute_stat_kegg

logger = logging.getLogger(__name__)


def compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file=None, kegg_metabolism_output_folder=None):
    """Compute stat associated with the number of proteome for each taxonomic affiliations.

    Args:
        proteomes_output_folder (str): pathname to the result folder containing subfolder containing proteomes
        clustering_output_folder (str): pathname to the result folder containing mmseqs results
        annotation_output_folder (str): pathname to the annotation reference folder containing annotations for each cluster
        stat_file (str): pathname to the tsv stat file
        kegg_metabolism_output_folder (str): pathname to the sbml folder containing draft metabolic network

    Returns:
        workflow_numbers (dict): dict containing observation names (as key) associated with proteomes, protein clusters, GO Terms and EC (as value)
    """
    workflow_numbers = {}

    result_folder = os.path.join(proteomes_output_folder, 'result')
    proteome_numbers = compute_stat_proteomes(result_folder)

    clustering_numbers = compute_stat_clustering(clustering_output_folder)

    annotation_numbers = compute_stat_annotation(annotation_output_folder)

    if kegg_metabolism_output_folder is not None:
        sbml_kegg_folder = os.path.join(kegg_metabolism_output_folder, 'sbml')
        kegg_numbers = compute_stat_kegg(sbml_kegg_folder)

    if kegg_metabolism_output_folder is not None:
        all_observation_names = {*proteome_numbers.keys(), *clustering_numbers.keys(), *annotation_numbers.keys(), *kegg_numbers.keys()}
    else:
        all_observation_names = {*proteome_numbers.keys(), *clustering_numbers.keys(), *annotation_numbers.keys()}

    for observation_name in all_observation_names:
        if observation_name in proteome_numbers:
            nb_proteomes = proteome_numbers[observation_name]
        else:
            nb_proteomes = 'NA'
        if observation_name in clustering_numbers:
            nb_shared_proteins = clustering_numbers[observation_name]
        else:
            nb_shared_proteins = 'NA'
        if observation_name in annotation_numbers:
            nb_gos = annotation_numbers[observation_name][0]
            nb_ecs = annotation_numbers[observation_name][1]
        else:
            nb_gos = 'NA'
            nb_ecs = 'NA'
        if kegg_metabolism_output_folder is not None:
            if observation_name in kegg_numbers:
                nb_reactions = kegg_numbers[observation_name][0]
                nb_metabolites = kegg_numbers[observation_name][1]
            else:
                nb_reactions = 'NA'
                nb_metabolites = 'NA'
        if kegg_metabolism_output_folder is not None:
            workflow_numbers[observation_name] = [nb_proteomes, nb_shared_proteins, nb_gos, nb_ecs, nb_reactions, nb_metabolites]
        else:
            workflow_numbers[observation_name] = [nb_proteomes, nb_shared_proteins, nb_gos, nb_ecs]

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            if kegg_metabolism_output_folder is not None:
                file_headers = ['observation_name', 'Number_proteomes', 'Number_shared_proteins', 'Number_go_terms', 'Number_ecs', 'Number_reactions', 'Number_metabolites']
            else:
                file_headers = ['observation_name', 'Number_proteomes', 'Number_shared_proteins', 'Number_go_terms', 'Number_ecs']
            csvwriter.writerow(file_headers)
            for observation_name in workflow_numbers:
                csvwriter.writerow([observation_name, *workflow_numbers[observation_name]])

    return workflow_numbers


def perform_workflow(input_file, output_folder, busco_percentage_keep=80, ignore_taxadb_update=None,
                        all_proteomes=None, uniprot_sparql_endpoint=None, remove_tmp=None,
                        limit_maximal_number_proteomes=99, rank_limit=None,
                        nb_cpu=1, clust_threshold=1, mmseqs_options=None,
                        linclust=None, propagate_annotation=None, uniref_annotation=None,
                        expression_annotation=None, minimal_number_proteomes=1,
                        kegg_reconstruction=False, mapping_ko=False, recreate_kegg=None,
                        annotate_with_uniprot=None):
    """From the proteomes found by esmecata proteomes, create protein cluster for each taxonomic affiliations.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations
        output_folder (str): pathname to the output folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        ignore_taxadb_update (bool): option to ignore ete3 taxa database update
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        remove_tmp (bool): remove the tmp files
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        rank_limit (str): rank limit to remove from the data
        nb_cpu (int): number of CPUs to be used by mmseqs
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        mmseqs_options (str): use alternative mmseqs option
        linclust (bool): use linclust
        propagate_annotation (float): float between 0 and 1. It is the ratio of proteins in the cluster that should have the annotation to keep this annotation.
        uniref_annotation (bool): option to use uniref annotation to add annotation
        expression_annotation (bool): option to add expression annotation from uniprot
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kept
        kegg_reconstruction (bool): use kegg_metabolism to reconstruct draft metabolic networks
        mapping_ko (bool): option to use KO ID to retrieve reactions
        recreate_kegg (bool): option to recreate KEGG model by guerying KEGG server
        annotate_with_uniprot (bool): query UniProt to retrieve protein annotations
    """
    starttime = time.time()
    workflow_metadata = {}

    if mapping_ko is True and kegg_reconstruction is not True:
        sys.exit('mapping_ko option (--map-ko) must be used with kegg_reconstruciton option (--kegg).')

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteomes_output_folder = os.path.join(output_folder, '0_proteomes')
    retrieve_proteomes(input_file, proteomes_output_folder, busco_percentage_keep,
                        ignore_taxadb_update, all_proteomes, uniprot_sparql_endpoint,
                        remove_tmp, limit_maximal_number_proteomes, rank_limit,
                        minimal_number_proteomes)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    make_clustering(proteomes_output_folder, clustering_output_folder, nb_cpu,
                        clust_threshold, mmseqs_options, linclust, remove_tmp)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    annotate_proteins(clustering_output_folder, annotation_output_folder, uniprot_sparql_endpoint,
                        propagate_annotation, uniref_annotation, expression_annotation,
                        annotate_with_uniprot)

    if kegg_reconstruction:
        kegg_metabolism_output_folder = os.path.join(output_folder, '3_kegg_metabolism')
        create_draft_networks(annotation_output_folder, kegg_metabolism_output_folder, mapping_ko, recreate_kegg)
    else:
        kegg_metabolism_output_folder = None

    clustering_folder = os.path.join(clustering_output_folder, 'reference_proteins')
    annotation_reference_folder = os.path.join(annotation_output_folder, 'annotation_reference')

    for clust_threshold in os.listdir(clustering_folder):
        clust_threshold_clustering_folder = os.path.join(clustering_folder, clust_threshold)
        clust_annotation_reference_folder = os.path.join(annotation_reference_folder, clust_threshold)
        stat_file = os.path.join(output_folder, 'stat_number_workflow_{0}.tsv'.format(clust_threshold))
        compute_stat_workflow(proteomes_output_folder, clust_threshold_clustering_folder, clust_annotation_reference_folder, stat_file, kegg_metabolism_output_folder)

    proteomes_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(proteomes_metadata_file, 'r') as json_idata:
        proteomes_metadata = json.load(json_idata)
    clustering_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'r') as json_idata:
        clustering_metadata = json.load(json_idata)
    annotation_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation.json')
    with open(annotation_metadata_file, 'r') as json_idata:
        annotation_metadata = json.load(json_idata)
    if kegg_metabolism_output_folder is not None:
        kegg_metadata_file = os.path.join(kegg_metabolism_output_folder, 'esmecata_metadata_kegg.json')
        with open(kegg_metadata_file, 'r') as json_idata:
            kegg_metadata = json.load(json_idata)

    workflow_metadata['proteomes_metadata'] = proteomes_metadata
    workflow_metadata['clustering_metadata'] = clustering_metadata
    workflow_metadata['annotation_metadata'] = annotation_metadata
    if kegg_metabolism_output_folder is not None:
        workflow_metadata['kegg_metadata'] = kegg_metadata

    endtime = time.time()
    duration = endtime - starttime
    workflow_metadata['esmecata_workflow_duration'] = duration
    workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow.json')
    with open(workflow_metadata_file, 'w') as ouput_file:
        json.dump(workflow_metadata, ouput_file, indent=4)
