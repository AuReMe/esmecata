import pandas as pd
import json
import csv
import os
from esmecata.core.proteomes import compute_stat_proteomes
from esmecata.core.clustering import compute_stat_clustering
from esmecata.core.workflow import compute_stat_workflow
import shutil

def convert_old_results_to_0_5_0(esmecata_proteomes_result_folder, esmecata_clustering_result_folder, esmecata_annotation_result_folder,
                                 output_folder):
    """ From old results of esmecata (before 0.5.0) convert some part of the result folder into 0.5.0.

    Args:
        esmecata_proteomes_result_folder (str): pathname to esmecata proteomes output folder
        esmecata_clustering_result_folder (str): pathname to esmecata clustering output folder
        esmecata_annotation_result_folder (str): pathname to esmecata annotation output folder
        output_folder (str): pathanme to ouput folder
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if not os.path.exists(os.path.join(output_folder, '0_proteomes')):
        os.mkdir(os.path.join(output_folder, '0_proteomes'))
    if not os.path.exists(os.path.join(output_folder, '1_clustering')):
        os.mkdir(os.path.join(output_folder, '1_clustering'))
    if not os.path.exists(os.path.join(output_folder, '2_annotation')):
        os.mkdir(os.path.join(output_folder, '2_annotation'))

    df = pd.read_csv(os.path.join(esmecata_proteomes_result_folder, 'proteome_tax_id.tsv'), sep='\t')

    observation_name_tax_ids = {}
    tax_id_names = []
    for index, row in df.iterrows():
        observation_name = row['observation_name']
        tax_name = row['name']
        tax_id = row['tax_id']
        # Convert characters in tax_name into unicode code.
        convert_tax_name = ''.join(c if c.isalnum() or c in '-_' else ('_c' + str(ord(c)) + '_') for c in tax_name)
        tax_id_name = convert_tax_name + '__taxid__' + str(tax_id)
        tax_id_names.append(tax_id_name)
        observation_name_tax_ids[observation_name] = tax_id_name

    df['tax_id_name'] = tax_id_names
    df.to_csv(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), sep='\t', index=None)
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv'))
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '2_annotation', 'proteome_tax_id.tsv'))

    proteome_status = {}
    proteome_description_folder = os.path.join(esmecata_proteomes_result_folder, 'proteomes_description')
    for proteome_description_file in os.listdir(proteome_description_folder):
        with open(os.path.join(proteome_description_folder, proteome_description_file), 'r') as open_proteome_description_file:
                csvreader = csv.DictReader(open_proteome_description_file, delimiter='\t')
                for line in csvreader:
                    if line['reference_proteome'] == 'False':
                        reference_proteome = False
                    else:
                        reference_proteome = True
                    proteome_status[line['proteome_id']] = reference_proteome
    compute_stat_proteomes(os.path.join(esmecata_proteomes_result_folder), stat_file=os.path.join(output_folder, '0_proteomes', 'stat_number_proteome.tsv'))

    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'esmecata_metadata_proteomes.json'), os.path.join(output_folder, '0_proteomes', 'esmecata_metadata_proteomes.json'))
    shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'esmecata_metadata_clustering.json'), os.path.join(output_folder, '1_clustering', 'esmecata_metadata_clustering.json'))
    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'esmecata_metadata_annotation.json'), os.path.join(output_folder, '2_annotation', 'esmecata_metadata_annotation.json'))

    df_stat_clustering = pd.read_csv(os.path.join(esmecata_clustering_result_folder, 'stat_number_clustering.tsv'), sep='\t')

    computed_threshold_folder = os.path.join(output_folder, '1_clustering', 'computed_threshold')
    if not os.path.exists(computed_threshold_folder):
        os.mkdir(computed_threshold_folder)
    reference_proteins_consensus_fasta_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins_consensus_fasta')
    if not os.path.exists(reference_proteins_consensus_fasta_folder):
        os.mkdir(reference_proteins_consensus_fasta_folder)

    already_processed_tax_id_names = []
    for index, row in df_stat_clustering.iterrows():
        observation_name = row['observation_name']
        tax_id_name = observation_name_tax_ids[observation_name]
        if tax_id_name not in  already_processed_tax_id_names:
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'computed_threshold', observation_name+'.tsv'), os.path.join(computed_threshold_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins_consensus_fasta', observation_name+'.faa'), os.path.join(reference_proteins_consensus_fasta_folder, tax_id_name+'.faa'))
            already_processed_tax_id_names.append(tax_id_name)
    compute_stat_clustering(os.path.join(output_folder, '1_clustering'), stat_file=os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv'))

    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'stat_number_annotation.tsv'), os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv'))

    shutil.copytree(os.path.join(esmecata_annotation_result_folder, 'annotation_reference'), os.path.join(output_folder, '2_annotation', 'annotation_reference'))

def convert_old_results_from_workflow_to_0_5_0(esmecata_workflow_results_folder, output_folder):
    """ From old results of esmecata workflow (before 0.5.0) convert some part of the result folder into 0.5.0.

    Args:
        esmecata_workflow_results_folder (str): pathname to esmecata workflow output folder
        output_folder (str): pathanme to ouput folder
    """
    esmecata_proteomes_result_folder = os.path.join(esmecata_workflow_results_folder, '0_proteomes')
    esmecata_clustering_result_folder = os.path.join(esmecata_workflow_results_folder, '1_clustering')
    esmecata_annotation_result_folder = os.path.join(esmecata_workflow_results_folder, '2_annotation')
    convert_old_results_to_0_5_0(esmecata_proteomes_result_folder, esmecata_clustering_result_folder, esmecata_annotation_result_folder,
                                 output_folder)


def convert_old_results_to_0_5_0_all_folder(esmecata_proteomes_result_folder, esmecata_clustering_result_folder, esmecata_annotation_result_folder,
                                 output_folder):
    """ From old results of esmecata (before 0.5.0) convert all the result folder into 0.5.0.

    Args:
        esmecata_proteomes_result_folder (str): pathname to esmecata proteomes output folder
        esmecata_clustering_result_folder (str): pathname to esmecata clustering output folder
        esmecata_annotation_result_folder (str): pathname to esmecata annotation output folder
        output_folder (str): pathanme to ouput folder
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if not os.path.exists(os.path.join(output_folder, '0_proteomes')):
        os.mkdir(os.path.join(output_folder, '0_proteomes'))
    if not os.path.exists(os.path.join(output_folder, '1_clustering')):
        os.mkdir(os.path.join(output_folder, '1_clustering'))
    if not os.path.exists(os.path.join(output_folder, '2_annotation')):
        os.mkdir(os.path.join(output_folder, '2_annotation'))

    df = pd.read_csv(os.path.join(esmecata_proteomes_result_folder, 'proteome_tax_id.tsv'), sep='\t')

    observation_name_tax_ids = {}
    tax_id_names = []
    for index, row in df.iterrows():
        observation_name = row['observation_name']
        tax_name = row['name']
        tax_id = row['tax_id']
        # Convert characters in tax_name into unicode code.
        convert_tax_name = ''.join(c if c.isalnum() or c in '-_' else ('_c' + str(ord(c)) + '_') for c in tax_name)
        tax_id_name = convert_tax_name + '__taxid__' + str(tax_id)
        tax_id_names.append(tax_id_name)
        observation_name_tax_ids[observation_name] = tax_id_name

    df['tax_id_name'] = tax_id_names
    df.to_csv(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), sep='\t', index=None)
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv'))
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '2_annotation', 'proteome_tax_id.tsv'))

    proteome_status = {}
    proteome_description_folder = os.path.join(esmecata_proteomes_result_folder, 'proteomes_description')
    for proteome_description_file in os.listdir(proteome_description_folder):
        with open(os.path.join(proteome_description_folder, proteome_description_file), 'r') as open_proteome_description_file:
                csvreader = csv.DictReader(open_proteome_description_file, delimiter='\t')
                for line in csvreader:
                    if line['reference_proteome'] == 'False':
                        reference_proteome = False
                    else:
                        reference_proteome = True
                    proteome_status[line['proteome_id']] = reference_proteome
    compute_stat_proteomes(os.path.join(esmecata_proteomes_result_folder), stat_file=os.path.join(output_folder, '0_proteomes', 'stat_number_proteome.tsv'))

    shutil.copytree(os.path.join(esmecata_proteomes_result_folder, 'proteomes'), os.path.join(output_folder, '0_proteomes', 'proteomes'))
    shutil.copytree(os.path.join(esmecata_proteomes_result_folder, 'proteomes_description'), os.path.join(output_folder, '0_proteomes', 'proteomes_description'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'association_taxon_taxID.json'), os.path.join(output_folder, '0_proteomes', 'association_taxon_taxID.json'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'tax_comparison_rank.tsv'), os.path.join(output_folder, '0_proteomes', 'tax_comparison_rank.tsv'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'taxonomy_diff.tsv'), os.path.join(output_folder, '0_proteomes', 'taxonomy_diff.tsv'))

    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'esmecata_metadata_proteomes.json'), os.path.join(output_folder, '0_proteomes', 'esmecata_metadata_proteomes.json'))
    shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'esmecata_metadata_clustering.json'), os.path.join(output_folder, '1_clustering', 'esmecata_metadata_clustering.json'))
    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'esmecata_metadata_annotation.json'), os.path.join(output_folder, '2_annotation', 'esmecata_metadata_annotation.json'))

    df_stat_clustering = pd.read_csv(os.path.join(esmecata_clustering_result_folder, 'stat_number_clustering.tsv'), sep='\t')

    computed_threshold_folder = os.path.join(output_folder, '1_clustering', 'computed_threshold')
    if not os.path.exists(computed_threshold_folder):
        os.mkdir(computed_threshold_folder)
    cluster_founds_folder = os.path.join(output_folder, '1_clustering', 'cluster_founds')
    if not os.path.exists(cluster_founds_folder):
        os.mkdir(cluster_founds_folder)
    reference_proteins_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins')
    if not os.path.exists(reference_proteins_folder):
        os.mkdir(reference_proteins_folder)
    reference_proteins_consensus_fasta_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins_consensus_fasta')
    if not os.path.exists(reference_proteins_consensus_fasta_folder):
        os.mkdir(reference_proteins_consensus_fasta_folder)
    reference_proteins_representative_fasta_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins_representative_fasta')
    if not os.path.exists(reference_proteins_representative_fasta_folder):
        os.mkdir(reference_proteins_representative_fasta_folder)

    eggnog_output_fasta_folder = os.path.join(output_folder, '2_annotation', 'eggnog_output')
    if not os.path.exists(eggnog_output_fasta_folder):
        os.mkdir(eggnog_output_fasta_folder)

    already_processed_tax_id_names = []
    for index, row in df_stat_clustering.iterrows():
        observation_name = row['observation_name']
        tax_id_name = observation_name_tax_ids[observation_name]
        if tax_id_name not in  already_processed_tax_id_names:
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'computed_threshold', observation_name+'.tsv'), os.path.join(computed_threshold_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'cluster_founds', observation_name+'.tsv'), os.path.join(cluster_founds_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins', observation_name+'.tsv'), os.path.join(reference_proteins_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins_consensus_fasta', observation_name+'.faa'), os.path.join(reference_proteins_consensus_fasta_folder, tax_id_name+'.faa'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins_representative_fasta', observation_name+'.faa'), os.path.join(reference_proteins_representative_fasta_folder, tax_id_name+'.faa'))

            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.seed_orthologs')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.seed_orthologs'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.seed_orthologs'))
            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.hits')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.hits'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.hits'))
            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.annotations')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', observation_name+'.emapper.annotations'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.annotations'))

            already_processed_tax_id_names.append(tax_id_name)

    compute_stat_clustering(os.path.join(output_folder, '1_clustering'), stat_file=os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv'))

    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'stat_number_annotation.tsv'), os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv'))

    shutil.copytree(os.path.join(esmecata_annotation_result_folder, 'annotation_reference'), os.path.join(output_folder, '2_annotation', 'annotation_reference'))
    shutil.copytree(os.path.join(esmecata_annotation_result_folder, 'pathologic'), os.path.join(output_folder, '2_annotation', 'pathologic'))
    compute_stat_workflow(os.path.join(output_folder, '0_proteomes'), os.path.join(output_folder, '1_clustering'), os.path.join(output_folder, '2_annotation'), os.path.join(output_folder, 'stat_number_workflow.tsv'))


def convert_intermediate_results_to_0_5_0_all_folder(esmecata_proteomes_result_folder, esmecata_clustering_result_folder, esmecata_annotation_result_folder,
                                 output_folder):
    """ From intemerdiary results of esmecata (during 0.5.0) convert all the result folder into 0.5.0.
    The principal changes in these fodler is the fact that the taxon name was used as filename in several intermediary files.

    Args:
        esmecata_proteomes_result_folder (str): pathname to esmecata proteomes output folder
        esmecata_clustering_result_folder (str): pathname to esmecata clustering output folder
        esmecata_annotation_result_folder (str): pathname to esmecata annotation output folder
        output_folder (str): pathanme to ouput folder
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if not os.path.exists(os.path.join(output_folder, '0_proteomes')):
        os.mkdir(os.path.join(output_folder, '0_proteomes'))
    if not os.path.exists(os.path.join(output_folder, '1_clustering')):
        os.mkdir(os.path.join(output_folder, '1_clustering'))
    if not os.path.exists(os.path.join(output_folder, '2_annotation')):
        os.mkdir(os.path.join(output_folder, '2_annotation'))

    df = pd.read_csv(os.path.join(esmecata_proteomes_result_folder, 'proteome_tax_id.tsv'), sep='\t')

    observation_name_tax_ids = {}
    observation_name_tax_names = {}

    tax_id_names = []
    for index, row in df.iterrows():
        observation_name = row['observation_name']
        tax_name = row['name']
        tax_id = row['tax_id']
        # Convert characters in tax_name into unicode code.
        convert_tax_name = ''.join(c if c.isalnum() or c in '-_' else ('_c' + str(ord(c)) + '_') for c in tax_name)
        tax_id_name = convert_tax_name + '__taxid__' + str(tax_id)
        tax_id_names.append(tax_id_name)
        observation_name_tax_ids[observation_name] = tax_id_name
        observation_name_tax_names[tax_name] = tax_id_name

    df['tax_id_name'] = tax_id_names
    df.to_csv(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), sep='\t', index=None)
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv'))
    shutil.copyfile(os.path.join(output_folder, '0_proteomes', 'proteome_tax_id.tsv'), os.path.join(output_folder, '2_annotation', 'proteome_tax_id.tsv'))

    proteome_status = {}
    proteome_description_folder = os.path.join(esmecata_proteomes_result_folder, 'proteomes_description')
    for proteome_description_file in os.listdir(proteome_description_folder):
        with open(os.path.join(proteome_description_folder, proteome_description_file), 'r') as open_proteome_description_file:
                csvreader = csv.DictReader(open_proteome_description_file, delimiter='\t')
                for line in csvreader:
                    if line['reference_proteome'] == 'False':
                        reference_proteome = False
                    else:
                        reference_proteome = True
                    proteome_status[line['proteome_id']] = reference_proteome
    compute_stat_proteomes(os.path.join(esmecata_proteomes_result_folder), stat_file=os.path.join(output_folder, '0_proteomes', 'stat_number_proteome.tsv'))

    shutil.copytree(os.path.join(esmecata_proteomes_result_folder, 'proteomes'), os.path.join(output_folder, '0_proteomes', 'proteomes'))
    shutil.copytree(os.path.join(esmecata_proteomes_result_folder, 'proteomes_description'), os.path.join(output_folder, '0_proteomes', 'proteomes_description'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'association_taxon_taxID.json'), os.path.join(output_folder, '0_proteomes', 'association_taxon_taxID.json'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'tax_comparison_rank.tsv'), os.path.join(output_folder, '0_proteomes', 'tax_comparison_rank.tsv'))
    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'taxonomy_diff.tsv'), os.path.join(output_folder, '0_proteomes', 'taxonomy_diff.tsv'))

    shutil.copyfile(os.path.join(esmecata_proteomes_result_folder, 'esmecata_metadata_proteomes.json'), os.path.join(output_folder, '0_proteomes', 'esmecata_metadata_proteomes.json'))
    shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'esmecata_metadata_clustering.json'), os.path.join(output_folder, '1_clustering', 'esmecata_metadata_clustering.json'))
    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'esmecata_metadata_annotation.json'), os.path.join(output_folder, '2_annotation', 'esmecata_metadata_annotation.json'))

    df_stat_clustering = pd.read_csv(os.path.join(esmecata_clustering_result_folder, 'stat_number_clustering.tsv'), sep='\t')

    computed_threshold_folder = os.path.join(output_folder, '1_clustering', 'computed_threshold')
    if not os.path.exists(computed_threshold_folder):
        os.mkdir(computed_threshold_folder)
    cluster_founds_folder = os.path.join(output_folder, '1_clustering', 'cluster_founds')
    if not os.path.exists(cluster_founds_folder):
        os.mkdir(cluster_founds_folder)
    reference_proteins_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins')
    if not os.path.exists(reference_proteins_folder):
        os.mkdir(reference_proteins_folder)
    reference_proteins_consensus_fasta_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins_consensus_fasta')
    if not os.path.exists(reference_proteins_consensus_fasta_folder):
        os.mkdir(reference_proteins_consensus_fasta_folder)
    reference_proteins_representative_fasta_folder = os.path.join(output_folder, '1_clustering', 'reference_proteins_representative_fasta')
    if not os.path.exists(reference_proteins_representative_fasta_folder):
        os.mkdir(reference_proteins_representative_fasta_folder)

    eggnog_output_fasta_folder = os.path.join(output_folder, '2_annotation', 'eggnog_output')
    if not os.path.exists(eggnog_output_fasta_folder):
        os.mkdir(eggnog_output_fasta_folder)

    already_processed_tax_id_names = []
    for index, row in df_stat_clustering.iterrows():
        observation_name = row['observation_name']
        tax_id_name = observation_name_tax_ids[observation_name]
        tax_name = observation_name_tax_names[observation_name]

        if tax_id_name not in  already_processed_tax_id_names:
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'computed_threshold', tax_name+'.tsv'), os.path.join(computed_threshold_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'cluster_founds', tax_name+'.tsv'), os.path.join(cluster_founds_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins', tax_name+'.tsv'), os.path.join(reference_proteins_folder, tax_id_name+'.tsv'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins_consensus_fasta', tax_name+'.faa'), os.path.join(reference_proteins_consensus_fasta_folder, tax_id_name+'.faa'))
            shutil.copyfile(os.path.join(esmecata_clustering_result_folder, 'reference_proteins_representative_fasta', tax_name+'.faa'), os.path.join(reference_proteins_representative_fasta_folder, tax_id_name+'.faa'))

            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.seed_orthologs')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.seed_orthologs'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.seed_orthologs'))
            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.hits')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.hits'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.hits'))
            if os.path.exists(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.annotations')):
                shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'eggnog_output', tax_name+'.emapper.annotations'), os.path.join(eggnog_output_fasta_folder, tax_id_name+'.emapper.annotations'))

            already_processed_tax_id_names.append(tax_id_name)
    compute_stat_clustering(os.path.join(output_folder, '1_clustering'), stat_file=os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv'))

    shutil.copyfile(os.path.join(esmecata_annotation_result_folder, 'stat_number_annotation.tsv'), os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv'))

    shutil.copytree(os.path.join(esmecata_annotation_result_folder, 'annotation_reference'), os.path.join(output_folder, '2_annotation', 'annotation_reference'))
    shutil.copytree(os.path.join(esmecata_annotation_result_folder, 'pathologic'), os.path.join(output_folder, '2_annotation', 'pathologic'))
    compute_stat_workflow(os.path.join(output_folder, '0_proteomes'), os.path.join(output_folder, '1_clustering'), os.path.join(output_folder, '2_annotation'), os.path.join(output_folder, 'stat_number_workflow.tsv'))
