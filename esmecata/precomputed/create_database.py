import os
import shutil
import pandas as pd
import zipfile
import json
import csv

from esmecata.core.eggnog import get_proteomes_tax_id_name
from esmecata.core.eggnog import read_annotation
from esmecata.core.annotation import extract_protein_cluster

from multiprocessing import Pool

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


def create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, prefix=''):
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

    with open(output_database_json, 'w') as dumpfile:
        json.dump(database_json_data, dumpfile, indent=4)


def copy_file(annotation_file, proteomes_taxa_names, output_database_folder, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder):
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


def create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, cpu_number=1):
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    proteomes_tax_id_file_path = os.path.join(esmecata_clustering_folder, 'proteome_tax_id.tsv')

    proteomes_tax_id_file_database_path = os.path.join(output_database_folder, 'proteome_tax_id.tsv')
    if not os.path.exists(proteomes_tax_id_file_database_path):
        shutil.copyfile(proteomes_tax_id_file_path, proteomes_tax_id_file_database_path)
    else:
        df_proteomes_tax_id_file_path = pd.read_csv(proteomes_tax_id_file_path, sep='\t')
        df_proteomes_tax_id_file_database_path = pd.read_csv(proteomes_tax_id_file_database_path, sep='\t')
        df_concat = pd.concat([df_proteomes_tax_id_file_path, df_proteomes_tax_id_file_database_path])
        df_concat.to_csv(proteomes_tax_id_file_database_path, sep='\t')

    proteomes_taxa_names = get_proteomes_tax_id_name(proteomes_tax_id_file_path)

    consensus_sequence_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins_consensus_fasta')
    computed_threshold_folder = os.path.join(esmecata_clustering_folder, 'computed_threshold')
    reference_proteins_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins')

    annotation_reference_folder = os.path.join(esmecata_annotation_folder, 'annotation_reference')

    database_copy_pool = Pool(cpu_number)

    multiprocessing_data = []
    for annotation_file in os.listdir(annotation_reference_folder):
        multiprocessing_data.append([annotation_file, proteomes_taxa_names, output_database_folder, computed_threshold_folder, annotation_reference_folder, consensus_sequence_folder])
    database_copy_pool.starmap(copy_file, multiprocessing_data)

    database_copy_pool.close()
    database_copy_pool.join()


def main(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_folder, cpu_number=1):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_database_folder = os.path.join(output_folder, 'database_folder')
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder)
    create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, cpu_number)

    compress_database_file = os.path.join(output_folder, 'esmecata_database')
    if not os.path.exists(compress_database_file):
        os.mkdir(compress_database_file)
    shutil.make_archive(compress_database_file, 'zip', output_database_folder)
    shutil.rmtree(output_database_folder)


def main_multiple_folders(esmecata_proteomes_folders, esmecata_clustering_folders, esmecata_annotation_folders, output_folder, cpu_number=1):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_database_folder = os.path.join(output_folder, 'database_folder')
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    for index, esmecata_proteomes_folder in enumerate(esmecata_proteomes_folders):
        taxon_level = os.path.basename(esmecata_proteomes_folder).split('_')[0]
        esmecata_clustering_folder = esmecata_clustering_folders[index]
        esmecata_annotation_folder = esmecata_annotation_folders[index]
        create_json(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, taxon_level)
        create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder, cpu_number)

    compress_database_file = os.path.join(output_folder, 'esmecata_database')
    if not os.path.exists(compress_database_file):
        os.mkdir(compress_database_file)
    shutil.make_archive(compress_database_file, 'zip', output_database_folder)
    shutil.rmtree(output_database_folder)
