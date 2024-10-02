# Copyright (C) 2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import os
import pandas as pd
import shutil
import sys
import time
import zipfile

from io import TextIOWrapper

from ete3 import __version__ as ete3_version
from ete3 import NCBITaxa

from esmecata.utils import is_valid_dir
from esmecata.core.proteomes import associate_taxon_to_taxon_id, disambiguate_taxon, filter_rank_limit, create_comp_taxonomy_file
from esmecata.core.eggnog import compute_stat_annotation
from esmecata.core.clustering import get_proteomes_tax_id_name
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
        taxon_proteomes (dict): dict of observation name as key and tax_id, tax_name, tax_rank and associated proteomes as value
    """
    database_taxon_ids = []
    proteomes_tax_id_names = {}
    taxon_data = {}
    with archive.open('proteome_tax_id.tsv', 'r') as open_database_taxon_file_path:
        csvreader = csv.DictReader(TextIOWrapper(open_database_taxon_file_path), delimiter='\t')
        for line in csvreader:
            database_taxon_ids.append(line['name'])
            proteomes_tax_id_names[line['tax_id']] = line['tax_id_name']
            taxon_data[line['tax_id']] = [line['tax_id'], line['name'], line['tax_rank'], line['proteome']]

    return database_taxon_ids, proteomes_tax_id_names, taxon_data


def find_proteomes_tax_ids_in_precomputed_database(json_taxonomic_affiliations, database_taxon_ids):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        database_taxon_ids (list): list of taxon IDs contained in the database

    Returns:
        association_taxon_database (dict): observation name (key) associated with tax_name and tax_id
        observation_name_not_founds (list): list of observation name which taxonomic affiliation has no match in esmecata database
    """
    logger.info('|EsMeCaTa|proteomes| Find observation name present in precomputed database.')

    association_taxon_database = {}
    for observation_name in json_taxonomic_affiliations:
        for tax_name in reversed(json_taxonomic_affiliations[observation_name]):
            tax_id = str(json_taxonomic_affiliations[observation_name][tax_name][0])

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_name in database_taxon_ids:
                association_taxon_database[observation_name] = (tax_name,  tax_id)
                logger.info('|EsMeCaTa|precomputed| "%s" present in database, %s will be associated with the taxon "%s".', tax_name, observation_name, tax_name)
                break

    observation_name_not_founds = []
    for observation_name in json_taxonomic_affiliations:
        if observation_name not in association_taxon_database:
            logger.info('|EsMeCaTa|precomputed| No taxon available in the database has been found for %s. Maybe, try a new run of esmecata with it?', observation_name)
            observation_name_not_founds.append(observation_name)

    return association_taxon_database, observation_name_not_founds


def precomputed_parse_affiliation(input_file, database_taxon_file_path, output_folder, rank_limit=None, update_affiliations=None):
    """From a tsv file with taxonomic affiliations find the associated proteomes and download them.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        update_affiliations (str): option to update taxonomic affiliations.
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|precomputed| Reading input file.')

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
               'rank_limit': rank_limit, 'update_affiliations': update_affiliations}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['ete3'] = ete3_version
    options['tool_dependencies']['python_package']['pandas'] = pd.__version__

    esmecata_metadata = {}
    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')
    esmecata_metadata['access_time'] = date
    esmecata_metadata['tool_options'] = options

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

    archive = zipfile.ZipFile(database_taxon_file_path, 'r')
    database_taxon_ids, proteomes_tax_id_names, taxon_data = get_taxon_database(archive)

    esmecata_metadata['precomputed_database'] = {}

    with archive.open('esmecata_database_metadata.json', 'r') as open_metadata_json:
        json_data = json.load(open_metadata_json)
    for json_key in json_data:
        if json_key.endswith('_proteomes'):
            proteomes_data_json = json_data[json_key]
    esmecata_metadata['precomputed_database']['esmecata_query_system'] = proteomes_data_json['esmecata_query_system']
    esmecata_metadata['precomputed_database']['uniprot_release'] = proteomes_data_json['uniprot_release']
    esmecata_metadata['precomputed_database']['access_time'] = proteomes_data_json['access_time']
    esmecata_metadata['precomputed_database']['swissprot_release_number'] = proteomes_data_json['swissprot_release_number']
    esmecata_metadata['precomputed_database']['swissprot_release_date'] = proteomes_data_json['swissprot_release_date']
    esmecata_metadata['precomputed_database']['trembl_release_number'] = proteomes_data_json['trembl_release_number']
    esmecata_metadata['precomputed_database']['trembl_release_date'] = proteomes_data_json['trembl_release_date']

    ncbi = NCBITaxa()

    # Parse taxonomic affiliations with ete3 to find matching taxon name.
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, update_affiliations, ncbi, output_folder)

    # Disambiguate taxon name using other taxon from the taxonomic affiliations.
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

    # If a rank limit option has been used, limit the taxonomic ranks used.
    if rank_limit:
        json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

    json_log = os.path.join(proteomes_output_folder, 'association_taxon_taxID.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    association_taxon_database, observation_name_not_founds = find_proteomes_tax_ids_in_precomputed_database(json_taxonomic_affiliations, database_taxon_ids)

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
                tax_name = association_taxon_database[observation_name][1]
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

    # Get only_reference_proteome_used.
    only_reference_proteome_used = {}
    with archive.open('stat_number_proteome.tsv') as open_stat_proteome_file:
        csvreader = csv.DictReader(TextIOWrapper(open_stat_proteome_file), delimiter='\t')
        for line in csvreader:
            only_reference_proteome_used[line['EsMeCaTa_used_taxon']] = line['only_reference_proteome_used']

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

    # If some organisms have no match in their taxonomic affiliations with esmecata precomputed database, create an input file for esmecata containing them.
    if len(observation_name_not_founds) > 0:
        df_not_found = df.loc[observation_name_not_founds]
        organism_not_found_in_database_file_path = os.path.join(output_folder, 'organism_not_found_in_database.tsv')
        df_not_found.to_csv(organism_not_found_in_database_file_path, sep='\t')

    # For each line of the input files that has a match in the database, recreate an imitation of esmecata output folder.
    for tax_id in tax_id_obs_names:
        taxi_id_name = proteomes_tax_id_names[tax_id]

        # Create a consensus proteoems file.
        clustering_consensus_file = os.path.join(taxi_id_name, taxi_id_name+'.faa')
        output_path_consensus_file = os.path.join(reference_proteins_consensus_fasta_folder, taxi_id_name+'.faa')
        if not os.path.exists(output_path_consensus_file):
            with archive.open(clustering_consensus_file) as zf, open(output_path_consensus_file, 'wb') as f:
                shutil.copyfileobj(zf, f)

        # Read annotation file
        annotation_file = os.path.join(taxi_id_name, taxi_id_name+'.tsv')
        with archive.open(annotation_file) as zf:
            df_annotation = pd.read_csv(zf, sep='\t')

        # Create a computed threhsold file.
        output_computed_threshold_file = os.path.join(computed_threshold_folder, taxi_id_name+'.tsv')
        df_annotation[['representative_protein', 'cluster_ratio', 'proteomes']].to_csv(output_computed_threshold_file, sep='\t', index=None)

        for observation_name in tax_id_obs_names[tax_id]:
            # Create an annotaiton_reference file for the observation name.
            output_path_annotation_file = os.path.join(annotation_reference_output_folder, observation_name+'.tsv')
            df_annotation.to_csv(output_path_annotation_file, sep='\t', index=None)

    archive.close()

    function_table_file_path = os.path.join(annotation_output_folder, 'function_table.tsv')
    create_dataset_annotation_file(annotation_reference_output_folder, function_table_file_path, 'all')

    tax_name_clustering_numbers = {}
    for clustering_file in os.listdir(computed_threshold_folder):
        clustering_file_path = os.path.join(computed_threshold_folder, clustering_file)
        with open(clustering_file_path, 'r') as open_clustering_file_path:
            csvreader = csv.reader(open_clustering_file_path, delimiter='\t')
            next(csvreader)
            cluster_0 = []
            cluster_0_5 = []
            cluster_0_95 = []
            for line in csvreader:
                protein_cluster_threshold = float(line[1])
                if protein_cluster_threshold >= 0.95:
                    cluster_0_95.append(line[0])
                if protein_cluster_threshold >= 0.5:
                    cluster_0_5.append(line[0])
                if protein_cluster_threshold >= 0:
                    cluster_0.append(line[0])

        tax_name_clustering_numbers[clustering_file.replace('.tsv', '')] = [cluster_0, cluster_0_5, cluster_0_95]
    proteomes_taxa_id_names = get_proteomes_tax_id_name(clustering_proteome_tax_id_file)
    clustering_numbers = {}
    for observation_name in proteomes_taxa_id_names:
        tax_name = proteomes_taxa_id_names[observation_name]
        clustering_numbers[observation_name] = tax_name_clustering_numbers[tax_name]

    clustering_stat_file = os.path.join(clustering_output_folder, 'stat_number_clustering.tsv')
    with open(clustering_stat_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_protein_clusters_panproteome', 'Number_protein_clusters_kept', 'Number_protein_clusters_coreproteome'])
        for observation_name in clustering_numbers:
            cluster_0 = len(clustering_numbers[observation_name][0])
            cluster_0_5 = len(clustering_numbers[observation_name][1])
            cluster_0_95 = len(clustering_numbers[observation_name][2])
            csvwriter.writerow([observation_name, cluster_0, cluster_0_5, cluster_0_95])

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
    df_precomputed_stat.to_csv(stat_precomputed_file, sep='\t')

    precomputed_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(precomputed_metadata_file, 'w') as ouput_file:
        json.dump(esmecata_metadata, ouput_file, indent=4)

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
    logger.info('|EsMeCaTa|precomputed| Extraction of data from database completed in {0}s.'.format(duration))
