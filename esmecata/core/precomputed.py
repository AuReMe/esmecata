# Copyright (C) 2024-2025 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import datetime
import json
import logging
import io
import os
import pandas as pd
import shutil
import sys
import time
import zipfile

from io import TextIOWrapper

from ete4 import __version__ as ete_version
from ete4 import NCBITaxa
from Bio import SeqIO
from Bio import __version__ as biopython_version

from esmecata.utils import is_valid_dir, get_domain_or_superkingdom_from_ncbi_tax_database
from esmecata.core.proteomes import associate_taxon_to_taxon_id, disambiguate_taxon, filter_rank_limit, create_comp_taxonomy_file
from esmecata.core.eggnog import compute_stat_annotation, write_pathologic
from esmecata.core.clustering import get_proteomes_tax_id_name, compute_openness_pan_proteomes
from esmecata.core.annotation import create_dataset_annotation_file

from esmecata import __version__ as esmecata_version

logger = logging.getLogger(__name__)


def get_taxon_database(archive):
    """ Extract taxon ID contained in precomputed database.

    Args:
        proteome_tax_id_file (str): open file object of esmecata precomputed database

    Returns:
        database_taxon_ids (dict): list of taxon IDs contained in the database
        proteomes_tax_id_names (dict): dict of observation name as key and tax_id_name as value
        taxon_data (dict): dict of observation name as key and tax_id, tax_name, tax_rank and associated proteomes as value
    """
    database_taxon_ids = {}
    proteomes_tax_id_names = {}
    taxon_data = {}
    with archive.open('proteome_tax_id.tsv', 'r') as open_database_taxon_file_path:
        csvreader = csv.DictReader(TextIOWrapper(open_database_taxon_file_path), delimiter='\t')
        for line in csvreader:
            database_taxon_ids[line['tax_id']] = line['name']
            proteomes_tax_id_names[line['tax_id']] = line['tax_id_name']
            taxon_data[line['tax_id']] = [line['tax_id'], line['name'], line['tax_rank'], line['proteome']]

    return database_taxon_ids, proteomes_tax_id_names, taxon_data


def find_proteomes_tax_ids_in_precomputed_database(json_taxonomic_affiliations, database_taxon_ids):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        database_taxon_ids (dict): list of taxon IDs contained in the database

    Returns:
        association_taxon_database (dict): observation name (key) associated with tax_name and tax_id
        observation_name_not_founds (list): list of observation name which taxonomic affiliation has no match in esmecata database
    """
    logger.info('|EsMeCaTa|precomputed| Search for observation names present in precomputed database.')

    association_taxon_database = {}
    for observation_name in json_taxonomic_affiliations:
        for tax_name in reversed(json_taxonomic_affiliations[observation_name]):
            tax_id = str(json_taxonomic_affiliations[observation_name][tax_name][0])

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in database_taxon_ids:
                if tax_name != database_taxon_ids[tax_id]:
                    logger.info('|EsMeCaTa|precomputed| "%s" has a matching taxon ID in the database (taxon ID: %s), but taxon names of input file and the one from database do not match. Taxon name from database will be used: "%s" instead of "%s".', observation_name, tax_id, database_taxon_ids[tax_id], tax_name)
                    tax_name = database_taxon_ids[tax_id]

                else:
                    logger.info('|EsMeCaTa|precomputed| "%s" has a matching taxon ID in the database (taxon ID: %s), it will be associated with the taxon "%s".', observation_name, tax_id, tax_name)

                association_taxon_database[observation_name] = (tax_name,  tax_id)
                break

    observation_name_not_founds = []
    for observation_name in json_taxonomic_affiliations:
        if observation_name not in association_taxon_database:
            logger.info('|EsMeCaTa|precomputed| No taxon available in the database has been found for %s. Maybe, try a new run of esmecata with it?', observation_name)
            observation_name_not_founds.append(observation_name)

    return association_taxon_database, observation_name_not_founds


def search_taxon_database(json_taxonomic_affiliations, database_taxon_file_path, clust_threshold,
                          computed_threshold_folder, annotation_reference_output_folder, reference_proteins_consensus_fasta_folder,
                          pathologic_folder):
    """ Search in esmecata precomputed database from taxa.

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified).
        database_taxon_file_path (str): path to esmecata precomputed database.
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster (must be equal or superior to the one use in the creation of the precomputed db).
        computed_threshold_folder (str): path to output computed threshold folder.
        annotation_reference_output_folder (str): path to output annotation reference folder.
        reference_proteins_consensus_fasta_folder (str): path to output reference proteins consensus folder.
        pathologic_folder (str): path to output pathologic folder.

    Returns:
        association_taxon_database (dict): observation name (key) associated with tax_name and tax_id.
        proteomes_tax_id_names (dict): dict of observation name as key and tax_id_name as value.
        taxon_data (dict): dict of observation name as key and tax_id, tax_name, tax_rank and associated proteomes as value.
        observation_name_not_founds (list): list of observation name which taxonomic affiliation has no match in esmecata database.
        only_reference_proteome_used (dict): esmecata selected taxon as key and boolean as value indicating if proteomes are reference proteomes.
        proteomes_data_json (dict): esmecata metadta dictionary for proteomes step.
        json_data (dict): esmecata metadata dictionary from database.
        precomputed_db_version (str): version of precomputed database.
    """
    archive = zipfile.ZipFile(database_taxon_file_path, 'r')
    database_taxon_ids, proteomes_tax_id_names, taxon_data = get_taxon_database(archive)

    with archive.open('esmecata_database_metadata.json', 'r') as open_metadata_json:
        json_data = json.load(open_metadata_json)
    for json_key in json_data:
        if json_key.endswith('_proteomes'):
            proteomes_data_json = json_data[json_key]
        if json_key.endswith('_clustering'):
            clustering_data_json = json_data[json_key]

    precomputed_db_version = json_data['database_version']
    logger.critical('|EsMeCaTa|precomputed| EsMeCaTa is using precomputed database "{0}" version {1}.'.format(database_taxon_file_path, precomputed_db_version))

    # Check compatibility between user threshold and database threshold.
    database_clustering_threhsold = clustering_data_json['tool_options']['clust_threshold']
    if clust_threshold < database_clustering_threhsold:
        logger.critical('|EsMeCaTa|precomputed| Selected threshold (-t) too low ({0}) compared to database threhsold ({1}), select one superior or equal to the one of the database.'.format(clust_threshold, database_clustering_threhsold))
        sys.exit(1)

    association_taxon_database, observation_name_not_founds = find_proteomes_tax_ids_in_precomputed_database(json_taxonomic_affiliations, database_taxon_ids)

    # Get only_reference_proteome_used.
    only_reference_proteome_used = {}
    with archive.open('stat_number_proteome.tsv') as open_stat_proteome_file:
        csvreader = csv.DictReader(TextIOWrapper(open_stat_proteome_file), delimiter='\t')
        for line in csvreader:
            only_reference_proteome_used[line['EsMeCaTa_used_taxon']] = line['only_reference_proteome_used']

    tax_id_obs_names = {}
    for observation_name in association_taxon_database:
        tax_id = association_taxon_database[observation_name][1]
        if tax_id != 'not_found':
            if tax_id not in tax_id_obs_names:
                tax_id_obs_names[tax_id] = [observation_name]
            else:
                tax_id_obs_names[tax_id].append(observation_name)

    # For each line of the input files that has a match in the database, recreate an imitation of esmecata output folder.
    for tax_id in tax_id_obs_names:
        tax_id_name = proteomes_tax_id_names[tax_id]

        # Read annotation file
        annotation_file = os.path.join(tax_id_name, tax_id_name+'.tsv')
        with archive.open(annotation_file) as zf:
            df_annotation = pd.read_csv(zf, sep='\t')

        # Create a computed threhsold file.
        output_computed_threshold_file = os.path.join(computed_threshold_folder, tax_id_name+'.tsv')
        df_annotation[['representative_protein', 'cluster_ratio', 'proteomes']].to_csv(output_computed_threshold_file, sep='\t', index=None)

        # Filter according to selected clust_threshold.
        df_annotation = df_annotation[df_annotation['cluster_ratio'] >= clust_threshold]

        kept_protein_ids = set(df_annotation['representative_protein'].tolist())
        for observation_name in tax_id_obs_names[tax_id]:
            # Create an annotaiton_reference file for the observation name.
            output_path_annotation_file = os.path.join(annotation_reference_output_folder, observation_name+'.tsv')
            if 'KEGG_reaction' in df_annotation.columns:
                selected_df_annotation = df_annotation[['representative_protein', 'cluster_members', 'gene_name', 'GO', 'EC', 'KEGG_reaction']]
                selected_df_annotation.columns = ['protein_cluster', 'cluster_members', 'gene_name', 'GO', 'EC', 'KEGG_reaction']
            else:
                selected_df_annotation = df_annotation[['representative_protein', 'cluster_members', 'gene_name', 'GO', 'EC']]
                selected_df_annotation.columns = ['protein_cluster', 'cluster_members', 'gene_name', 'GO', 'EC']
            if not os.path.exists(output_path_annotation_file):
                selected_df_annotation.to_csv(output_path_annotation_file, sep='\t', index=None)

            # Create pathologic files.
            pathologic_organism_folder = os.path.join(pathologic_folder, observation_name)
            # Add _1 to pathologic file as genetic element cannot have the same name as the organism.
            pathologic_file = os.path.join(pathologic_organism_folder, observation_name+'_1.pf')

            if not os.path.exists(pathologic_file):
                # Replace NA by empty string.
                selected_df_annotation = selected_df_annotation.fillna('')
                reference_proteins_string = selected_df_annotation.set_index('protein_cluster')['cluster_members'].to_dict()
                reference_proteins = {protein: reference_proteins_string[protein].split(',') for protein in reference_proteins_string}
                sub_annotated_proteins = selected_df_annotation.set_index('protein_cluster').to_dict('index')
                # Rename dictionary keys according to the ones expected in write_pathologic.
                sub_annotated_proteins = [(protein, {'Preferred_name': sub_annotated_proteins[protein]['gene_name'],
                                                    'GOs': sub_annotated_proteins[protein]['GO'],
                                                    'EC': sub_annotated_proteins[protein]['EC']})
                                                    for protein in sub_annotated_proteins]
                if len(sub_annotated_proteins) > 0:
                    is_valid_dir(pathologic_organism_folder)
                    write_pathologic(observation_name, sub_annotated_proteins, pathologic_file, reference_proteins)

        # Create a consensus proteomes file.
        clustering_consensus_file = os.path.join(tax_id_name, tax_id_name+'.faa')
        output_path_consensus_file = os.path.join(reference_proteins_consensus_fasta_folder, tax_id_name+'.faa')
        if not os.path.exists(output_path_consensus_file):
            records = []
            with archive.open(clustering_consensus_file) as zf:
                open_zf_text  = io.TextIOWrapper(zf)
                for record in SeqIO.parse(open_zf_text, 'fasta'):
                    if '|' in record.id:
                        protein_id = record.id.split('|')[1]
                    else:
                        protein_id = record.id
                    if protein_id in kept_protein_ids:
                        records.append(record)

            if len(records) > 0:
                logger.critical('|EsMeCaTa|precomputed| {0} protein clusters kept for taxon {1} using threshold {2}.'.format(len(records), tax_id_name, clust_threshold))
                SeqIO.write(records, output_path_consensus_file, 'fasta')
            else:
                logger.critical('|EsMeCaTa|precomputed| 0 protein clusters kept for taxon {0} using threshold {1}, it will not have predictions.'.format(tax_id_name, clust_threshold))

    archive.close()

    return association_taxon_database, proteomes_tax_id_names, taxon_data, observation_name_not_founds, only_reference_proteome_used, proteomes_data_json, json_data, precomputed_db_version


def precomputed_parse_affiliation(input_file, database_taxon_file_path, output_folder, rank_limit=None, update_affiliations=None,
                                  clust_threshold=0.5):
    """From a tsv file with taxonomic affiliations find the associated proteomes and download them.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        update_affiliations (str): option to update taxonomic affiliations.
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster (must be equal or superior to the one use in the creation of the precomputed db).
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|precomputed| Reading input file.')

    # Check how domain/superkingdom are called in NCBI Taxonomy database.
    domain_superkingdom_tax_rank_name = get_domain_or_superkingdom_from_ncbi_tax_database()
    if rank_limit is not None:
        if rank_limit == 'superkingdom' and rank_limit != domain_superkingdom_tax_rank_name and domain_superkingdom_tax_rank_name == 'domain':
            logger.info('|EsMeCaTa|proteomes| Rank limit sets as superkingdom but it is called domain inside NCBI Taxonomy database, use domain instead.')
            rank_limit = 'domain'
        elif rank_limit == 'domain' and rank_limit != domain_superkingdom_tax_rank_name and domain_superkingdom_tax_rank_name == 'superkingdom':
            logger.info('|EsMeCaTa|proteomes| Rank limit sets as domain but it is called superkingdom inside NCBI Taxonomy database, use superkingdom instead.')
            rank_limit = 'superkingdom'

    is_valid_dir(output_folder)

    proteomes_output_folder = os.path.join(output_folder, '0_proteomes')
    is_valid_dir(proteomes_output_folder)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    is_valid_dir(clustering_output_folder)

    computed_threshold_folder = os.path.join(clustering_output_folder, 'computed_threshold')
    is_valid_dir(computed_threshold_folder)

    reference_proteins_consensus_fasta_folder = os.path.join(clustering_output_folder, 'reference_proteins_consensus_fasta')
    is_valid_dir(reference_proteins_consensus_fasta_folder)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    is_valid_dir(annotation_output_folder)

    annotation_reference_output_folder = os.path.join(annotation_output_folder, 'annotation_reference')
    is_valid_dir(annotation_reference_output_folder)

    pathologic_folder = os.path.join(annotation_output_folder, 'pathologic')
    is_valid_dir(pathologic_folder)

    # Metadata of the script.
    options = {'input_file': input_file, 'output_folder': output_folder, 'database_taxon_file_path': database_taxon_file_path,
               'rank_limit': rank_limit, 'update_affiliations': update_affiliations, 'clust_threshold': clust_threshold}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['ete4'] = ete_version
    options['tool_dependencies']['python_package']['pandas'] = pd.__version__
    options['tool_dependencies']['python_package']['biopython'] = biopython_version

    esmecata_metadata = {}
    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')
    esmecata_metadata['access_time'] = date
    esmecata_metadata['tool_options'] = options
    esmecata_metadata['precomputed_database'] = {}

    known_extensions = ['.xlsx', '.tsv', '.csv']
    file_name, file_extension = os.path.splitext(input_file)

    if file_extension not in known_extensions:
        logger.critical('|EsMeCaTa|precomputed| ERROR: Extension ({0}) of input file ({1}) is not part of the compatible extensions ({2}).'.format(file_extension, input_file, ','.join(known_extensions)))
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

    ncbi = NCBITaxa()

    # Parse taxonomic affiliations with ete4 to find matching taxon name.
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, update_affiliations, ncbi, output_folder)

    # Disambiguate taxon name using other taxon from the taxonomic affiliations.
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

    # If a rank limit option has been used, limit the taxonomic ranks used.
    if rank_limit:
        json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

    json_log = os.path.join(proteomes_output_folder, 'association_taxon_taxID.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    predictions_from_which_database = {}
    # If multiple database files are given, iter on each of them to associate taxa with observation name. 
    if ' ' in database_taxon_file_path:
        association_taxon_database = {}
        proteomes_tax_id_names = {}
        taxon_data = {}
        observation_name_not_founds =  []
        only_reference_proteome_used = {}

        json_taxonomic_affiliations_to_research = json_taxonomic_affiliations.copy()
        for single_database_taxon_file_path in database_taxon_file_path.split(' '):
            tmp_database_name = os.path.splitext(os.path.basename(single_database_taxon_file_path))[0]
            tmp_association_taxon_database, tmp_proteomes_tax_id_names, tmp_taxon_data, \
                tmp_observation_name_not_founds, tmp_only_reference_proteome_used, tmp_proteomes_data_json, tmp_json_data, \
                tmp_precomputed_db_version = search_taxon_database(json_taxonomic_affiliations_to_research, single_database_taxon_file_path, clust_threshold,
                                                            computed_threshold_folder, annotation_reference_output_folder, reference_proteins_consensus_fasta_folder,
                                                            pathologic_folder)

            # Extract metadata.
            esmecata_metadata['precomputed_database'][tmp_database_name] = {}
            esmecata_metadata['precomputed_database'][tmp_database_name]['esmecata_query_system'] = tmp_proteomes_data_json['esmecata_query_system']
            esmecata_metadata['precomputed_database'][tmp_database_name]['uniprot_release'] = tmp_proteomes_data_json['uniprot_release']
            esmecata_metadata['precomputed_database'][tmp_database_name]['access_time'] = tmp_proteomes_data_json['access_time']
            esmecata_metadata['precomputed_database'][tmp_database_name]['swissprot_release_number'] = tmp_proteomes_data_json['swissprot_release_number']
            esmecata_metadata['precomputed_database'][tmp_database_name]['swissprot_release_date'] = tmp_proteomes_data_json['swissprot_release_date']
            esmecata_metadata['precomputed_database'][tmp_database_name]['trembl_release_number'] = tmp_proteomes_data_json['trembl_release_number']
            esmecata_metadata['precomputed_database'][tmp_database_name]['trembl_release_date'] = tmp_proteomes_data_json['trembl_release_date']
            esmecata_metadata['precomputed_database'][tmp_database_name]['esmecata_precomputed_db_version'] = tmp_precomputed_db_version

            # Merge all results from the different databases into variables.
            for observation_name in tmp_association_taxon_database:
                predictions_from_which_database[observation_name] = tmp_database_name
                # Remove the observation_name from the ones to be searched with the next database.
                del json_taxonomic_affiliations_to_research[observation_name]
                if observation_name not in association_taxon_database:
                    association_taxon_database[observation_name] = tmp_association_taxon_database[observation_name]

            for tax_id in tmp_proteomes_tax_id_names:
                if tax_id not in proteomes_tax_id_names:
                    proteomes_tax_id_names[tax_id] = tmp_proteomes_tax_id_names[tax_id]

            for tax_id in tmp_taxon_data:
                if tax_id not in taxon_data:
                    taxon_data[tax_id] = tmp_taxon_data[tax_id]

            for observation_name in tmp_observation_name_not_founds:
                if observation_name not in observation_name_not_founds:
                    observation_name_not_founds.append(observation_name)

            for esmecata_taxon in tmp_only_reference_proteome_used:
                if esmecata_taxon not in only_reference_proteome_used:
                    only_reference_proteome_used[esmecata_taxon] = tmp_only_reference_proteome_used[esmecata_taxon]

            # Create metadata file for each database.
            for json_key in tmp_json_data:
                if json_key.endswith('_clustering'):
                    clustering_json_data = tmp_json_data[json_key]
            precomputed_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering_'+tmp_database_name+'.json')
            with open(precomputed_metadata_file, 'w') as ouput_file:
                json.dump(clustering_json_data, ouput_file, indent=4)

            for json_key in tmp_json_data:
                if json_key.endswith('_annotation'):
                    annotation_json_data = tmp_json_data[json_key]
            precomputed_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation_'+tmp_database_name+'.json')
            with open(precomputed_metadata_file, 'w') as ouput_file:
                json.dump(annotation_json_data, ouput_file, indent=4)

        # Ensure that observation names were not found with another database.
        observation_name_not_founds = [observation_name for observation_name in observation_name_not_founds if observation_name not in association_taxon_database]

    else:
        # Run esmecata precomputed on one database.
        association_taxon_database, proteomes_tax_id_names, taxon_data, \
            observation_name_not_founds, only_reference_proteome_used, proteomes_data_json, json_data, \
            precomputed_db_version = search_taxon_database(json_taxonomic_affiliations, database_taxon_file_path, clust_threshold,
                                                        computed_threshold_folder, annotation_reference_output_folder, reference_proteins_consensus_fasta_folder,
                                                        pathologic_folder)

        database_name = os.path.splitext(os.path.basename(database_taxon_file_path))[0]
        for observation_name in association_taxon_database:
            tax_id = association_taxon_database[observation_name][1]
            if tax_id != 'not_found':
                predictions_from_which_database[observation_name] = database_name

        esmecata_metadata['precomputed_database']['esmecata_query_system'] = proteomes_data_json['esmecata_query_system']
        esmecata_metadata['precomputed_database']['uniprot_release'] = proteomes_data_json['uniprot_release']
        esmecata_metadata['precomputed_database']['access_time'] = proteomes_data_json['access_time']
        esmecata_metadata['precomputed_database']['swissprot_release_number'] = proteomes_data_json['swissprot_release_number']
        esmecata_metadata['precomputed_database']['swissprot_release_date'] = proteomes_data_json['swissprot_release_date']
        esmecata_metadata['precomputed_database']['trembl_release_number'] = proteomes_data_json['trembl_release_number']
        esmecata_metadata['precomputed_database']['trembl_release_date'] = proteomes_data_json['trembl_release_date']
        esmecata_metadata['precomputed_database']['esmecata_precomputed_db_version'] = precomputed_db_version

        for json_key in json_data:
            if json_key.endswith('_clustering'):
                clustering_json_data = json_data[json_key]
        precomputed_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering.json')
        with open(precomputed_metadata_file, 'w') as ouput_file:
            json.dump(clustering_json_data, ouput_file, indent=4)

        for json_key in json_data:
            if json_key.endswith('_annotation'):
                annotation_json_data = json_data[json_key]
        precomputed_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation.json')
        with open(precomputed_metadata_file, 'w') as ouput_file:
            json.dump(annotation_json_data, ouput_file, indent=4)

    # Create proteome_tax_id file.
    proteome_tax_id_file = os.path.join(proteomes_output_folder, 'proteome_tax_id.tsv')
    tax_id_obs_names = {}
    proteomes_ids = {}
    with open(proteome_tax_id_file, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        csvwriter.writerow(['observation_name', 'name', 'tax_id', 'tax_id_name', 'tax_rank', 'proteome'])
        for observation_name in association_taxon_database:
            tax_id = association_taxon_database[observation_name][1]
            if tax_id != 'not_found':
                tax_name = association_taxon_database[observation_name][0]
                tax_id_name = proteomes_tax_id_names[tax_id]
                tax_rank = taxon_data[tax_id][2]
                proteome = taxon_data[tax_id][3]
                proteomes_ids[observation_name] = (tax_id, proteome.split(','))
                if tax_id not in tax_id_obs_names:
                    tax_id_obs_names[tax_id] = [observation_name]
                else:
                    tax_id_obs_names[tax_id].append(observation_name)
                csvwriter.writerow([observation_name, tax_name, tax_id, tax_id_name, tax_rank, taxon_data[tax_id][3]])

    create_comp_taxonomy_file(json_taxonomic_affiliations, proteomes_ids, tax_id_names, proteomes_output_folder)

    # Create stat_proteome.
    run_proteome_tax_id_file = os.path.join(proteomes_output_folder, 'stat_number_proteome.tsv')
    with open(run_proteome_tax_id_file, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_proteomes', 'Input_taxon_Name', 'Taxon_rank', 'EsMeCaTa_used_taxon', 'EsMeCaTa_used_rank', 'only_reference_proteome_used'])
        for observation_name in association_taxon_database:
            tax_name = association_taxon_database[observation_name][0]
            tax_id = association_taxon_database[observation_name][1]
            if tax_id != 'not_found':
                tax_rank = taxon_data[tax_id][2]
                nb_tax_proteomes = len(taxon_data[tax_id][3].split(','))
                reversed_affiliation_taxa = list(reversed(list(json_taxonomic_affiliations[observation_name].keys())))
                lowest_tax_rank = None
                for input_tax_name in reversed_affiliation_taxa:
                    if json_taxonomic_affiliations[observation_name][input_tax_name] != ['not_found']:
                        if input_tax_name != 'unknown':
                            lowest_tax_id = json_taxonomic_affiliations[observation_name][input_tax_name][0]
                            lowest_tax_rank = ncbi.get_rank([lowest_tax_id])[lowest_tax_id]
                            break
                csvwriter.writerow([observation_name, nb_tax_proteomes, input_tax_name, lowest_tax_rank, tax_name, tax_rank, only_reference_proteome_used[tax_name]])

    clustering_proteome_tax_id_file = os.path.join(clustering_output_folder, 'proteome_tax_id.tsv')
    shutil.copyfile(proteome_tax_id_file, clustering_proteome_tax_id_file)
    annotation_proteome_tax_id_file = os.path.join(annotation_output_folder, 'proteome_tax_id.tsv')
    shutil.copyfile(proteome_tax_id_file, annotation_proteome_tax_id_file)

    # Create mpwt taxon ID file.
    clustering_taxon_id = {}
    with open(proteome_tax_id_file, 'r') as input_taxon_id_file:
        taxon_id_csvreader = csv.DictReader(input_taxon_id_file, delimiter='\t')
        for line in taxon_id_csvreader:
            clustering_taxon_id[line['observation_name']] = line['tax_id']

    pathologic_taxon_id_file = os.path.join(pathologic_folder, 'taxon_id.tsv')
    with open(pathologic_taxon_id_file, 'w') as taxon_id_file:
        taxon_id_csvwriter = csv.writer(taxon_id_file, delimiter='\t')
        taxon_id_csvwriter.writerow(['species', 'taxon_id'])
        for species in clustering_taxon_id:
            taxon_id_csvwriter.writerow([species, clustering_taxon_id[species]])

    # If some organisms have no match in their taxonomic affiliations with esmecata precomputed database, create an input file for esmecata containing them.
    organism_not_found_in_database_file_path = os.path.join(output_folder, 'organism_not_found_in_database.tsv')
    if len(observation_name_not_founds) > 0:
        df_not_found = df.loc[observation_name_not_founds]
        df_not_found.to_csv(organism_not_found_in_database_file_path, sep='\t')
    else:
        with open(organism_not_found_in_database_file_path, 'w'):
            pass

    function_table_file_path = os.path.join(annotation_output_folder, 'function_table.tsv')
    create_dataset_annotation_file(annotation_reference_output_folder, function_table_file_path, 'all')

    tax_name_clustering_numbers = {}
    for clustering_file in os.listdir(computed_threshold_folder):
        clustering_file_path = os.path.join(computed_threshold_folder, clustering_file)
        with open(clustering_file_path, 'r') as open_clustering_file_path:
            csvreader = csv.reader(open_clustering_file_path, delimiter='\t')
            next(csvreader)
            cluster_0 = []
            selected_threshold_cluster = []
            cluster_0_95 = []
            for line in csvreader:
                protein_cluster_threshold = float(line[1])
                if protein_cluster_threshold >= 0.95:
                    cluster_0_95.append(line[0])
                if protein_cluster_threshold >= clust_threshold:
                    selected_threshold_cluster.append(line[0])
                if protein_cluster_threshold >= 0:
                    cluster_0.append(line[0])

        if len(selected_threshold_cluster) > 0:
            tax_name_clustering_numbers[clustering_file.replace('.tsv', '')] = [cluster_0, selected_threshold_cluster, cluster_0_95]
    proteomes_taxa_id_names = get_proteomes_tax_id_name(clustering_proteome_tax_id_file)

    clustering_numbers = {}
    for observation_name in proteomes_taxa_id_names:
        tax_name = proteomes_taxa_id_names[observation_name]
        if tax_name in tax_name_clustering_numbers:
            clustering_numbers[observation_name] = tax_name_clustering_numbers[tax_name]

    clustering_stat_file = os.path.join(clustering_output_folder, 'stat_number_clustering.tsv')
    with open(clustering_stat_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_protein_clusters_panproteome', 'Number_protein_clusters_kept', 'Number_protein_clusters_coreproteome'])
        for observation_name in clustering_numbers:
            cluster_0 = len(clustering_numbers[observation_name][0])
            selected_threshold_cluster = len(clustering_numbers[observation_name][1])
            cluster_0_95 = len(clustering_numbers[observation_name][2])
            csvwriter.writerow([observation_name, cluster_0, selected_threshold_cluster, cluster_0_95])

    # Compute openness of proteomes.
    proteome_openness_file = os.path.join(clustering_output_folder, 'stat_openness_proteomes.tsv')
    compute_openness_pan_proteomes(computed_threshold_folder, proteome_openness_file, clustering_threhsold=0.5, iteration_nb=100)

    annotation_stat_file = os.path.join(annotation_output_folder, 'stat_number_annotation.tsv')
    compute_stat_annotation(annotation_reference_output_folder, annotation_stat_file)

    endtime = time.time()
    duration = endtime - starttime
    esmecata_metadata['esmecata_precomputed_duration'] = duration
    precomputed_metadata_file = os.path.join(output_folder, 'esmecata_metadata_precomputed.json')
    if os.path.exists(precomputed_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_precomputed' in metadata_file]
        precomputed_metadata_file = os.path.join(output_folder, 'esmecata_metadata_precomputed_{0}.json'.format(len(metadata_files)))
        with open(precomputed_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)
    else:
        with open(precomputed_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)

    stat_precomputed_file = os.path.join(output_folder, 'stat_number_precomputed.tsv')
    df_proteomes = pd.read_csv(run_proteome_tax_id_file, sep='\t')
    df_proteomes.set_index('observation_name', inplace=True)
    df_clustering = pd.read_csv(clustering_stat_file, sep='\t')
    df_clustering.set_index('observation_name', inplace=True)
    df_annotation = pd.read_csv(annotation_stat_file, sep='\t')
    df_annotation.set_index('observation_name', inplace=True)
    df_precomputed_stat = df_proteomes.join([df_clustering, df_annotation])
    df_precomputed_stat['precomputed_db'] = [predictions_from_which_database[obs_name] if obs_name in predictions_from_which_database else '' for obs_name in df_precomputed_stat.index]
    df_precomputed_stat.to_csv(stat_precomputed_file, sep='\t')

    precomputed_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(precomputed_metadata_file, 'w') as ouput_file:
        json.dump(esmecata_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|precomputed| Extraction of data from database completed in {0}s.'.format(duration))
