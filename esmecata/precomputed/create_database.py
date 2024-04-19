import os
import shutil
import pandas as pd
import zipfile
import csv

from esmecata.core.eggnog import get_proteomes_tax_id_name
from esmecata.core.eggnog import read_annotation
from esmecata.core.annotation import extract_protein_cluster


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
            tax_id_name = line['name']
            proteomes_taxa_id_names[observation_name] = tax_id_name

    return proteomes_taxa_id_names

def create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder):
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)

    proteomes_tax_id_file_path = os.path.join(esmecata_clustering_folder, 'proteome_tax_id.tsv')
    proteomes_tax_id_file_database_path = os.path.join(output_database_folder, 'proteome_tax_id.tsv')
    shutil.copyfile(proteomes_tax_id_file_path, proteomes_tax_id_file_database_path)
    proteomes_taxa_names = get_proteomes_tax_id_name(proteomes_tax_id_file_path)

    # Copy log.
    esmecata_metadata_proteomes_input_path = os.path.join(esmecata_proteomes_folder, 'esmecata_metadata_proteomes.json')
    esmecata_metadata_proteomes_database_path = os.path.join(output_database_folder, 'esmecata_metadata_proteomes.json')
    shutil.copyfile(esmecata_metadata_proteomes_input_path, esmecata_metadata_proteomes_database_path)
    esmecata_metadata_clustering_input_path = os.path.join(esmecata_clustering_folder, 'esmecata_metadata_clustering.json')
    esmecata_metadata_clustering_database_path = os.path.join(output_database_folder, 'esmecata_metadata_clustering.json')
    shutil.copyfile(esmecata_metadata_clustering_input_path, esmecata_metadata_clustering_database_path)
    esmecata_metadata_annotation_input_path = os.path.join(esmecata_annotation_folder, 'esmecata_metadata_annotation.json')
    esmecata_metadata_annotation_database_path = os.path.join(output_database_folder, 'esmecata_metadata_annotation.json')
    shutil.copyfile(esmecata_metadata_annotation_input_path, esmecata_metadata_annotation_database_path)

    consensus_sequence_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins_consensus_fasta')
    computed_threshold_folder = os.path.join(esmecata_clustering_folder, 'computed_threshold')
    reference_proteins_folder = os.path.join(esmecata_clustering_folder, 'reference_proteins')

    eggnog_folder = True
    annotation_data_folder = os.path.join(esmecata_annotation_folder, 'eggnog_output')
    if not os.path.exists(annotation_data_folder):
        eggnog_folder = False
    annotation_reference_folder = os.path.join(esmecata_annotation_folder, 'annotation_reference')

    if eggnog_folder is False:
        for annotation_file in os.listdir(annotation_reference_folder):
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
    
    elif eggnog_folder is True:
        for annotation_file in os.listdir(annotation_data_folder):
            if annotation_file.endswith('.annotations'):
                taxon_id_name = os.path.splitext(annotation_file)[0].replace('.emapper', '')

                output_database_taxon_folder = os.path.join(output_database_folder, taxon_id_name)
                if not os.path.exists(output_database_taxon_folder):
                    os.mkdir(output_database_taxon_folder)

                    computed_threshold_file_path = os.path.join(computed_threshold_folder, taxon_id_name + '.tsv')
                    df_computed_threshold = pd.read_csv(computed_threshold_file_path, sep='\t')
                    df_computed_threshold.set_index('representative_protein', inplace=True)

                    reference_proteins_file_path = os.path.join(reference_proteins_folder, taxon_id_name+'.tsv')
                    reference_proteins, _ = extract_protein_cluster(reference_proteins_file_path)
                    reference_proteins = [[rep_prot, ','.join(reference_proteins[rep_prot])] for rep_prot in reference_proteins]
                    df_reference_proteins = pd.DataFrame(reference_proteins)
                    df_reference_proteins.columns = ['protein_cluster', 'cluster_members']
                    df_reference_proteins.set_index('protein_cluster', inplace=True)

                    annotation_data_file_path = os.path.join(annotation_data_folder, annotation_file)
                    annotations = read_annotation(annotation_data_file_path)
                    df_annotation = pd.DataFrame({i[0]: i[1] for i in annotations}).T
                    df_annotation.index = [dataframe_index.split('|')[1] for dataframe_index in df_annotation.index]

                    df_join = df_computed_threshold.join(df_annotation)
                    df_join = df_join.join(df_reference_proteins)
                    taxon_annotation_file = os.path.join(output_database_taxon_folder, taxon_id_name+'.tsv')
                    df_join.to_csv(taxon_annotation_file, sep='\t')

                    consensus_sequence_file_path = os.path.join(consensus_sequence_folder, taxon_id_name + '.faa')
                    taxon_consensus_sequence_file_path = os.path.join(output_database_taxon_folder, taxon_id_name+'.faa')
                    shutil.copyfile(consensus_sequence_file_path, taxon_consensus_sequence_file_path)


def main(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_database_folder = os.path.join(output_folder, 'database_folder')
    if not os.path.exists(output_database_folder):
        os.mkdir(output_database_folder)
    create_database(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, output_database_folder)

    compress_database_file = os.path.join(output_folder, 'esmecata_database')
    if not os.path.exists(compress_database_file):
        os.mkdir(compress_database_file)
    shutil.make_archive(compress_database_file, 'zip', output_database_folder)
    shutil.rmtree(output_database_folder)