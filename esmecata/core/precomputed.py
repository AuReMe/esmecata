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
import json
import logging
import os
import pandas as pd
import shutil
import sys
import time
import zipfile

from io import TextIOWrapper

from ete3 import __version__ as ete3_version
from ete3 import NCBITaxa, is_taxadb_up_to_date

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir, send_uniprot_sparql_query
from esmecata.core.proteomes import associate_taxon_to_taxon_id, disambiguate_taxon, filter_rank_limit

from esmecata import __version__ as esmecata_version

logger = logging.getLogger(__name__)

def get_taxon_database(archive):
    """ Extract taxon ID contained in precomputed database.

    Args:
        proteome_tax_id_file (str): open file object of esmecata precomputed database

    Returns:
        database_taxon_ids (dict): list of taxon IDs contained in the database
    """
    database_taxon_ids = []
    proteomes_tax_id_names = {}
    with archive.open('proteome_tax_id.tsv', 'r') as open_database_taxon_file_path:
        csvreader = csv.DictReader(TextIOWrapper(open_database_taxon_file_path), delimiter='\t')
        for line in csvreader:
            database_taxon_ids.append(line['tax_id'])
            proteomes_tax_id_names[line['observation_name']] = line['tax_id_name']

    return database_taxon_ids, proteomes_tax_id_names


def find_proteomes_tax_ids(json_taxonomic_affiliations, database_taxon_ids):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        database_taxon_ids (list): list of taxon IDs contained in the database

    Returns:
        association_taxon_database (dict): observation name (key) associated with tax_name and tax_id
    """
    logger.info('|EsMeCaTa|proteomes| Find observation name present in precomputed database.')

    association_taxon_database = {}
    for observation_name in json_taxonomic_affiliations:
        for tax_name in reversed(json_taxonomic_affiliations[observation_name]):
            tax_id = str(json_taxonomic_affiliations[observation_name][tax_name][0])

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in database_taxon_ids:
                association_taxon_database[observation_name] = (tax_name,  tax_id)
                logger.info('|EsMeCaTa|precomputed| "%s" present in database, %s will be associated with the taxon "%s".', tax_name, observation_name, tax_name)
                break

    return association_taxon_database


def precomputed_parse_affiliation(input_file, database_taxon_file_path, output_folder, rank_limit=None, update_affiliations=None):
    """From a tsv file with taxonomic affiliations find the associated proteomes and download them.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        update_affiliations (str): option to update taxonomic affiliations.
    """
    starttime = time.time()

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    if not os.path.exists(clustering_output_folder):
        os.mkdir(clustering_output_folder)
    computed_threshold_folder = os.path.join(clustering_output_folder, 'computed_threshold')
    if not os.path.exists(computed_threshold_folder):
        os.mkdir(computed_threshold_folder)
    reference_proteins_consensus_fasta_folder = os.path.join(clustering_output_folder, 'reference_proteins_consensus_fasta')
    if not os.path.exists(reference_proteins_consensus_fasta_folder):
        os.mkdir(reference_proteins_consensus_fasta_folder)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    if not os.path.exists(annotation_output_folder):
        os.mkdir(annotation_output_folder)

    annotation_reference_output_folder = os.path.join(annotation_output_folder, 'annotation_reference')
    if not os.path.exists(annotation_reference_output_folder):
        os.mkdir(annotation_reference_output_folder)
    eggnog_output_folder = os.path.join(annotation_output_folder, 'eggnog_output')
    if not os.path.exists(eggnog_output_folder):
        os.mkdir(eggnog_output_folder)
    pathologic_folder = os.path.join(annotation_output_folder, 'pathologic')
    if not os.path.exists(pathologic_folder):
        os.mkdir(pathologic_folder)

    known_extensions = ['.xlsx', '.tsv', '.csv']
    file_name, file_extension = os.path.splitext(input_file)

    if file_extension not in known_extensions:
        logger.critical('|EsMeCaTa|proteomes| ERROR: Extension ({0}) of input file ({1}) is not part of the compatible extensions ({2}).'.format(file_extension, input_file, ','.join(known_extensions)))
        sys.exit(1)

    if '.xlsx' in input_file:
        df = pd.read_excel(input_file)
    if '.tsv' in input_file:
        df = pd.read_csv(input_file, sep='\t')
    if '.csv' in input_file:
        df = pd.read_csv(input_file, sep=',')

    # Set index on the column containing the OTU name.
    df.set_index('observation_name', inplace=True)

    # taxonomic_affiliation is the column containing the taxonomic affiliation separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomic_affiliation']

    proteome_tax_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    archive = zipfile.ZipFile(database_taxon_file_path, 'r')
    database_taxon_ids, proteomes_tax_id_names = get_taxon_database(archive)

    ncbi = NCBITaxa()

    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, update_affiliations, ncbi, output_folder)

    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

    if rank_limit:
        json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

    json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    association_taxon_database = find_proteomes_tax_ids(json_taxonomic_affiliations, database_taxon_ids)

    for observation_name in association_taxon_database:
        taxi_id_name = proteomes_tax_id_names[observation_name]
        clustering_consensus_file = os.path.join(taxi_id_name, taxi_id_name+'.faa')
        output_path_consensus_file = os.path.join(reference_proteins_consensus_fasta_folder, taxi_id_name+'.faa')
        with archive.open(clustering_consensus_file) as zf, open(output_path_consensus_file, 'wb') as f:
            shutil.copyfileobj(zf, f)

        annotation_file = os.path.join(taxi_id_name, taxi_id_name+'.tsv')
        df_annotation = pd.read_csv(archive.open(annotation_file), sep='\t')

        output_computed_threshold_file = os.path.join(computed_threshold_folder, taxi_id_name+'.tsv')
        df_annotation_eggnog = df_annotation[['representative_protein', 'seed_ortholog', 'evalue', 'score', 'COG_category',
                                              'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
                                              'CAZy', 'BiGG_Reaction', 'PFAMs']]
        output_eggnog_annotation_file = os.path.join(eggnog_output_folder, taxi_id_name+'.emapper.annotations')
        df_annotation_eggnog.to_csv(output_eggnog_annotation_file, sep='\t')

        df_annotation[['representative_protein', 'cluster_ratio', 'proteomes']].to_csv(output_computed_threshold_file, sep='\t')

        df_annotation_reference = df_annotation[['representative_protein', 'cluster_members', 'Preferred_name', 'GOs', 'EC', 'KEGG_Reaction']]
        df_annotation_reference.columns = ['protein_cluster', 'cluster_members', 'gene_name', 'GO', 'EC', 'KEGG_reaction']
        output_path_annotation_file = os.path.join(annotation_reference_output_folder, observation_name+'.tsv')
        df_annotation_reference.to_csv(output_path_annotation_file, sep='\t')

    archive.close()

    endtime = time.time()
    duration = endtime - starttime
    logger.info('|EsMeCaTa|proteomes| Proteome step complete in {0}s.'.format(duration))
