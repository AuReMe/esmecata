# Copyright (C) 2024 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
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

import os
import shutil
import pandas as pd
import json
import csv
import zipfile

from esmecata.core.eggnog import get_proteomes_tax_id_name

from multiprocessing import Pool
from collections import Counter


def get_proteomes_tax_id_name(proteomes_tax_id_file_path):
    """ Extract tax_name + tax_id associated with observation name.

    Args:
        proteomes_tax_id_file_path (str): pathname to the proteomes_tax_id file.

    Returns:
        proteomes_taxa_id_names (dict): dict containing observation names (as key) associated with tax name + tax_id used for proteomes (as value)
    """
    proteomes_taxa_id_names = {}
    with open(proteomes_tax_id_file_path, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            tax_id_name = line['tax_id_name']
            proteomes_taxa_id_names[observation_name] = tax_id_name

    return proteomes_taxa_id_names


def create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, database_version, prefix=''):
    """ Copy json metadata file.

    Args:
        esmecata_proteomes_folder (str): path to esmecata proteomes folder.
        esmecata_clustering_folder (str): path to esmecata clustering folder.
        esmecata_annotation_folder (str): path to esmecata anntoation folder.
        output_database_folder (str): path to output folder containing zip database of esmecata.
        prefix (str): add prefix to json file.
    """
    output_database_json = os.path.join(output_database_folder, 'esmecata_database_metadata.json')
    if os.path.exists(output_database_json):
        with open(output_database_json, 'r') as open_json:
            database_json_data = json.load(open_json)
    else:
        database_json_data = {}

    for input_file in os.listdir(esmecata_proteomes_folder):
        if 'metadata' in input_file and input_file.endswith('.json'):
            esmecata_metadata_proteomes_input_path = os.path.join(esmecata_proteomes_folder, input_file)
            with open(esmecata_metadata_proteomes_input_path, 'r') as open_json:
                json_data = json.load(open_json)
                if prefix != '':
                    json_name = prefix + '_' + input_file.replace('.json', '')
                else:
                    json_name = input_file.replace('.json', '')
                database_json_data[json_name] = json_data

    for input_file in os.listdir(esmecata_clustering_folder):
        if 'metadata' in input_file and input_file.endswith('.json'):
            esmecata_clustering_folder_input_path = os.path.join(esmecata_clustering_folder, input_file)
            with open(esmecata_clustering_folder_input_path, 'r') as open_json:
                json_data = json.load(open_json)
                if prefix != '':
                    json_name = prefix + '_' + input_file.replace('.json', '')
                else:
                    json_name = input_file.replace('.json', '')
                database_json_data[json_name] = json_data

    for input_file in os.listdir(esmecata_annotation_folder):
        if 'metadata' in input_file and input_file.endswith('.json'):
            esmecata_annotation_folder_input_path = os.path.join(esmecata_annotation_folder, input_file)
            with open(esmecata_annotation_folder_input_path, 'r') as open_json:
                json_data = json.load(open_json)
                if prefix != '':
                    json_name = prefix + '_' + input_file.replace('.json', '')
                else:
                    json_name = input_file.replace('.json', '')
                database_json_data[json_name] = json_data

    database_json_data['database_version'] = database_version
    with open(output_database_json, 'w') as dumpfile:
        json.dump(database_json_data, dumpfile, indent=4)


def copy_file(annotation_file, proteomes_taxa_names, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder, output_database_folder):
    """ Copy file from esmecata run into esmecata database.

    Args:
        annotation_file (str): path to esmecata annotation reference file for a taxon.
        proteomes_taxa_names (str): path to esmecata proteomes taxa file.
        computed_threshold_folder (str): path to computed threshold folder.
        annotation_reference_folder (str): path to esmecata annotation reference folder.
        consensus_sequence_folder (str): path to esmecata consensus sequence folder.
        output_database_folder (str): path to output folder containing zip database of esmecata.
    """
    observation_name = os.path.splitext(annotation_file)[0]
    taxon_id_name = proteomes_taxa_names[observation_name]

    output_database_taxon_folder = os.path.join(output_database_folder, taxon_id_name)
    if not os.path.exists(output_database_taxon_folder):
        os.mkdir(output_database_taxon_folder)

        computed_threshold_file_path = os.path.join(computed_threshold_folder, taxon_id_name + '.tsv')
        df_computed_threshold = pd.read_csv(computed_threshold_file_path, sep='\t')
        df_computed_threshold.set_index('representative_protein', inplace=True)

        annotation_data_file_path = os.path.join(annotation_reference_folder, annotation_file)
        df_annotation = pd.read_csv(annotation_data_file_path, sep='\t')
        df_annotation.set_index('protein_cluster', inplace=True)
        df_join = df_computed_threshold.join(df_annotation)

        taxon_annotation_file = os.path.join(output_database_taxon_folder, taxon_id_name+'.tsv')
        df_join.to_csv(taxon_annotation_file, sep='\t')

        consensus_sequence_file_path = os.path.join(consensus_sequence_folder, taxon_id_name + '.faa')
        taxon_consensus_sequence_file_path = os.path.join(output_database_taxon_folder, taxon_id_name+'.faa')
        shutil.copyfile(consensus_sequence_file_path, taxon_consensus_sequence_file_path)


def concat_tsv_file(input_tsv_file_path, tsv_file_concat_path):
    """ Concat tsv file for database.

    Args:
        input_tsv_file_path (str): path to input tsv file.
        tsv_file_concat_path (str): path to tsv file containing concatenation.
    """
    if not os.path.exists(tsv_file_concat_path):
        shutil.copyfile(input_tsv_file_path, tsv_file_concat_path)
    else:
        df_proteomes_tax_id_file_path = pd.read_csv(input_tsv_file_path, sep='\t')
        df_tsv_file_concat = pd.read_csv(tsv_file_concat_path, sep='\t')
        df_concat = pd.concat([df_proteomes_tax_id_file_path, df_tsv_file_concat])
        df_concat.to_csv(tsv_file_concat_path, sep='\t', index=False)


def create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, nb_core=1):
    """ Create esmecata database from a run of esmecata.

    Args:
        esmecata_proteomes_folder (str): path to esmecata proteomes folder.
        esmecata_clustering_folder (str): path to esmecata clustering folder.
        esmecata_annotation_folder (str): path to esmecata anntoation folder.
        output_database_folder (str): path to output folder containing zip database of esmecata.
        nb_core (int): number of core to use when creating database.
    """
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    # Create/merge proteome_tax_id file for each taxon level.
    proteomes_tax_id_file_path = os.path.join(esmecata_clustering_folder, 'proteome_tax_id.tsv')
    proteomes_tax_id_file_database_path = os.path.join(output_database_folder, 'proteome_tax_id.tsv')
    concat_tsv_file(proteomes_tax_id_file_path, proteomes_tax_id_file_database_path)

    # Create/merge stat_proteome file for each taxon level.
    stat_number_proteome_file_path = os.path.join(esmecata_proteomes_folder, 'stat_number_proteome.tsv')
    stat_number_proteome_database_path = os.path.join(output_database_folder, 'stat_number_proteome.tsv')
    concat_tsv_file(stat_number_proteome_file_path, stat_number_proteome_database_path)

    # Create/merge stat_number_clustering file for each taxon level.
    stat_number_clustering_file_path = os.path.join(esmecata_clustering_folder, 'stat_number_clustering.tsv')
    stat_number_clustering_database_path = os.path.join(output_database_folder, 'stat_number_clustering.tsv')
    concat_tsv_file(stat_number_clustering_file_path, stat_number_clustering_database_path)

    # Create/merge stat_proteome file for each taxon level.
    stat_number_annotation_file_path = os.path.join(esmecata_annotation_folder, 'stat_number_annotation.tsv')
    stat_number_annotation_database_path = os.path.join(output_database_folder, 'stat_number_annotation.tsv')
    concat_tsv_file(stat_number_annotation_file_path, stat_number_annotation_database_path)

    proteomes_taxa_names = get_proteomes_tax_id_name(proteomes_tax_id_file_path)

    consensus_sequence_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins_consensus_fasta')
    computed_threshold_folder = os.path.join(esmecata_clustering_folder, 'computed_threshold')
    reference_proteins_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins')

    annotation_reference_folder = os.path.join(esmecata_annotation_folder, 'annotation_reference')

    # Use multiprocessing to copy fasta and annotation files.
    database_copy_pool = Pool(nb_core)

    similar_tax_id_names = []
    for observation_name in proteomes_taxa_names:
        similar_tax_id_names.append(proteomes_taxa_names[observation_name])
    similar_tax_id_names = Counter(similar_tax_id_names)

    multiprocessing_data = []
    for annotation_file in os.listdir(annotation_reference_folder):
        observation_name = os.path.splitext(annotation_file)[0]
        taxon_id_name = proteomes_taxa_names[observation_name]
        # If more than 1 observation name associated with a similar tax_id_name, do not use multiprocessing as it leads error.
        if similar_tax_id_names[taxon_id_name] > 1:
            copy_file(annotation_file, proteomes_taxa_names, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder, output_database_folder)
        else:
            multiprocessing_data.append([annotation_file, proteomes_taxa_names, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder, output_database_folder])
    database_copy_pool.starmap(copy_file, multiprocessing_data)

    database_copy_pool.close()
    database_copy_pool.join()


def create_database_from_run(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, nb_core=1):
    """ Create esmecata database from a run of esmecata.

    Args:
        esmecata_proteomes_folder (str): path to esmecata proteomes folder.
        esmecata_clustering_folder (str): path to esmecata clustering folder.
        esmecata_annotation_folder (str): path to esmecata anntoation folder.
        output_database_folder (str): path to output folder containing zip database of esmecata.
        nb_core (int): number of core to use when creating database.
    """
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    # Create/merge proteome_tax_id file for each taxon level.
    proteomes_tax_id_file_path = os.path.join(esmecata_clustering_folder, 'proteome_tax_id.tsv')
    df_proteomes_tax_id = pd.read_csv(proteomes_tax_id_file_path, sep='\t')
    df_proteomes_tax_id.set_index('observation_name', inplace=True)

    # Create/merge stat_proteome file for each taxon level.
    stat_number_proteome_file_path = os.path.join(esmecata_proteomes_folder, 'stat_number_proteome.tsv')
    df_stat_number_proteome = pd.read_csv(stat_number_proteome_file_path, sep='\t')
    df_stat_number_proteome.set_index('observation_name', inplace=True)

    # Create/merge stat_number_clustering file for each taxon level.
    stat_number_clustering_file_path = os.path.join(esmecata_clustering_folder, 'stat_number_clustering.tsv')
    df_stat_number_clustering = pd.read_csv(stat_number_clustering_file_path, sep='\t')
    df_stat_number_clustering.set_index('observation_name', inplace=True)

    # Create/merge stat_proteome file for each taxon level.
    stat_number_annotation_file_path = os.path.join(esmecata_annotation_folder, 'stat_number_annotation.tsv')
    df_stat_number_annotation = pd.read_csv(stat_number_annotation_file_path, sep='\t')
    df_stat_number_annotation.set_index('observation_name', inplace=True)

    proteomes_taxa_names = get_proteomes_tax_id_name(proteomes_tax_id_file_path)

    consensus_sequence_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins_consensus_fasta')
    computed_threshold_folder = os.path.join(esmecata_clustering_folder, 'computed_threshold')
    reference_proteins_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins')

    annotation_reference_folder = os.path.join(esmecata_annotation_folder, 'annotation_reference')

    # Use multiprocessing to copy fasta and annotation files.
    database_copy_pool = Pool(nb_core)

    similar_tax_id_names = []
    for observation_name in proteomes_taxa_names:
        similar_tax_id_names.append(proteomes_taxa_names[observation_name])
    similar_tax_id_names = Counter(similar_tax_id_names)

    proteome_tax_id_data = []
    stat_number_proteome_data = []
    stat_number_clustering_data = []
    stat_number_annotation_data = []
    already_processed_tax_id = []
    multiprocessing_data = []
    for annotation_file in os.listdir(annotation_reference_folder):
        observation_name = os.path.splitext(annotation_file)[0]
        taxon_id_name = proteomes_taxa_names[observation_name]
        if taxon_id_name not in already_processed_tax_id:
            proteome_tax_id_data.append([taxon_id_name, *df_proteomes_tax_id.loc[observation_name].to_list()])
            stat_number_proteome_data.append([taxon_id_name, *df_stat_number_proteome.loc[observation_name].to_list()])
            stat_number_clustering_data.append([taxon_id_name, *df_stat_number_clustering.loc[observation_name].to_list()])
            stat_number_annotation_data.append([taxon_id_name, *df_stat_number_annotation.loc[observation_name].to_list()])

            already_processed_tax_id.append(taxon_id_name)
        # If more than 1 observation name associated with a similar tax_id_name, do not use multiprocessing as it leads error.
        if similar_tax_id_names[taxon_id_name] > 1:
            copy_file(annotation_file, proteomes_taxa_names, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder, output_database_folder)
        else:
            multiprocessing_data.append([annotation_file, proteomes_taxa_names, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder, output_database_folder])
    database_copy_pool.starmap(copy_file, multiprocessing_data)

    database_copy_pool.close()
    database_copy_pool.join()

    proteomes_tax_id_file_database_path = os.path.join(output_database_folder, 'proteome_tax_id.tsv')
    pd.DataFrame(proteome_tax_id_data, columns=['observation_name', *df_proteomes_tax_id.columns]).to_csv(proteomes_tax_id_file_database_path, sep='\t', index=None)

    stat_number_proteome_database_path = os.path.join(output_database_folder, 'stat_number_proteome.tsv')
    pd.DataFrame(stat_number_proteome_data, columns=['observation_name', *df_stat_number_proteome.columns]).to_csv(stat_number_proteome_database_path, sep='\t', index=None)

    stat_number_clustering_database_path = os.path.join(output_database_folder, 'stat_number_clustering.tsv')
    pd.DataFrame(stat_number_clustering_data, columns=['observation_name', *df_stat_number_clustering.columns]).to_csv(stat_number_clustering_database_path, sep='\t', index=None)

    stat_number_annotation_database_path = os.path.join(output_database_folder, 'stat_number_annotation.tsv')
    pd.DataFrame(stat_number_annotation_data, columns=['observation_name', *df_stat_number_annotation.columns]).to_csv(stat_number_annotation_database_path, sep='\t', index=None)


def create_database_from_esmecata_run(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_folder, database_version, nb_core=1):
    """ Extract data from esmecata run and create an esmecata database from these.

    Args:
        esmecata_proteomes_folder (str): path to esmecata proteomes folder.
        esmecata_clustering_folder (str): path to esmecata clustering folder.
        esmecata_annotation_folder (str): path to esmecata anntoation folder.
        database_version (str): database version number.
        output_folder (str): path to output folder containing zip database of esmecata.
        nb_core (int): number of core to use when creating database.
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_database_folder = os.path.join(output_folder, 'database_folder')
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, database_version)
    create_database_from_run(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, nb_core)

    compress_database_file = os.path.join(output_folder, 'esmecata_database')
    if not os.path.exists(compress_database_file):
        os.mkdir(compress_database_file)
    shutil.make_archive(compress_database_file, 'zip', output_database_folder)
    shutil.rmtree(output_database_folder)


def create_database_from_esmecata_workflow_run(esmecata_workflow_folder, output_folder, database_version, nb_core=1):
    """ Extract data from esmecata workflow run and create an esmecata database from these.

    Args:
        esmecata_workflow_folder (str): path to esmecata workflow folder.
        output_folder (str): path to output folder containing zip database of esmecata.
        database_version (str): database version number.
        nb_core (int): number of core to use when creating database.
    """
    esmecata_proteomes_folder = os.path.join(esmecata_workflow_folder, '0_proteomes')
    esmecata_clustering_folder = os.path.join(esmecata_workflow_folder, '1_clustering')
    esmecata_annotation_folder = os.path.join(esmecata_workflow_folder, '2_annotation')
    create_database_from_esmecata_run(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_folder, database_version, nb_core)


def create_database_from_multiple_esmecata_runs(esmecata_proteomes_folders, esmecata_clustering_folders, esmecata_annotation_folders, output_folder, database_version, cpu_number=1):
    """ Extract data from esmecata runs and create an esmecata database from these.

    Args:
        esmecata_proteomes_folders (list): list of path to esmecata proteomes folder.
        esmecata_clustering_folders (list): list of path to esmecata clustering folder.
        esmecata_annotation_folders (list): list of path to esmecata anntoation folder.
        output_folder (str): path to output folder containing zip database of esmecata.
        nb_core (int): number of core to use when creating database.
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_database_folder = os.path.join(output_folder, 'database_folder')
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    for index, esmecata_proteomes_folder in enumerate(esmecata_proteomes_folders):
        taxon_level = os.path.basename(esmecata_proteomes_folder).split('_')[0]
        esmecata_clustering_folder = esmecata_clustering_folders[index]
        esmecata_annotation_folder = esmecata_annotation_folders[index]
        create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, database_version, taxon_level)
        create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, cpu_number)

    compress_database_file = os.path.join(output_folder, 'esmecata_database')
    if not os.path.exists(compress_database_file):
        os.mkdir(compress_database_file)
    shutil.make_archive(compress_database_file, 'zip', output_database_folder)
    shutil.rmtree(output_database_folder)


def merge_db_files(list_db_files, output_folder):
    stat_number_proteome_df = pd.DataFrame(columns=['observation_name', 'Number_proteomes', 'Input_taxon_Name', 'Taxon_rank', 'EsMeCaTa_used_taxon', 'EsMeCaTa_used_rank', 'only_reference_proteome_used'])
    stat_number_clustering_df = pd.DataFrame(columns=['observation_name', 'Number_protein_clusters_panproteome', 'Number_protein_clusters_kept', 'Number_protein_clusters_coreproteome'])
    stat_number_annotation_df = pd.DataFrame(columns=['observation_name', 'Number_go_terms', 'Number_ecs'])
    proteome_tax_id_df_data = []
    already_extracted_proteomes = []
    archive_name = os.path.join(output_folder, 'esmecata_database.zip')
    already_processed_zip_files = []
    with zipfile.ZipFile(archive_name, 'w', compression=zipfile.ZIP_DEFLATED) as open_database:
        for database_file in list_db_files:
            archive = zipfile.ZipFile(database_file, 'r')
            for zip_file in archive.namelist():
                if zip_file == 'proteome_tax_id.tsv':
                    proteome_tmp_df = pd.read_csv(archive.open(zip_file), sep='\t')
                    proteomes_lists = proteome_tmp_df.values.tolist()
                    for proteomes_list in proteomes_lists:
                        # If a tax_id_name will be found in several databases, use only the first occurence.
                        if proteomes_list[0] not in already_extracted_proteomes:
                            proteome_tax_id_df_data.append(proteomes_list)
                            already_extracted_proteomes.append(proteomes_list[0])
                if zip_file == 'stat_number_proteome.tsv':
                    tmp_df = pd.read_csv(archive.open(zip_file), sep='\t')
                    # To avoid adding redundant tax_id_name, remove all the ones already present.
                    tmp_df = tmp_df[~tmp_df['observation_name'].isin(stat_number_proteome_df['observation_name'])]
                    stat_number_proteome_df['only_reference_proteome_used'] = stat_number_proteome_df['only_reference_proteome_used'].astype(bool)
                    stat_number_proteome_df = pd.concat([stat_number_proteome_df, tmp_df])
                if zip_file == 'stat_number_clustering.tsv':
                    tmp_df = pd.read_csv(archive.open(zip_file), sep='\t')
                    tmp_df= tmp_df[~tmp_df['observation_name'].isin(stat_number_clustering_df['observation_name'])]
                    stat_number_clustering_df = pd.concat([stat_number_clustering_df, tmp_df])
                if zip_file == 'stat_number_annotation.tsv':
                    tmp_df = pd.read_csv(archive.open(zip_file), sep='\t')
                    tmp_df = tmp_df[~tmp_df['observation_name'].isin(stat_number_annotation_df['observation_name'])]
                    stat_number_annotation_df = pd.concat([stat_number_annotation_df, tmp_df])
                if not zip_file.endswith('/') and zip_file not in ['proteome_tax_id.tsv', 'stat_number_annotation.tsv', 'stat_number_clustering.tsv', 'stat_number_proteome.tsv']:
                    if zip_file not in already_processed_zip_files:
                        open_database.writestr(zip_file, archive.open(zip_file).read())
                        already_processed_zip_files.append(zip_file)

            archive.close()
        proteome_tax_id_df = pd.DataFrame(proteome_tax_id_df_data, columns=proteome_tmp_df.columns)
        proteome_tax_id_df_str = proteome_tax_id_df.to_csv(sep='\t', index=False)
        open_database.writestr('proteome_tax_id.tsv', proteome_tax_id_df_str)
        stat_number_proteome_df_str = stat_number_proteome_df.to_csv(sep='\t', index=False)
        open_database.writestr('stat_number_proteome.tsv', stat_number_proteome_df_str)
        stat_number_clustering_df_str = stat_number_clustering_df.to_csv(sep='\t', index=False)
        open_database.writestr('stat_number_clustering.tsv', stat_number_clustering_df_str)
        stat_number_annotation_df_str = stat_number_annotation_df.to_csv(sep='\t', index=False)
        open_database.writestr('stat_number_annotation.tsv', stat_number_annotation_df_str)
