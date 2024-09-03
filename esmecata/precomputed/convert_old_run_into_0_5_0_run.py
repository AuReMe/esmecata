import pandas as pd
import json
import csv
import os
from esmecata.core.proteomes import compute_stat_proteomes
from esmecata.core.clustering import compute_stat_clustering
import shutil

esmecata_workflow_folder = 'biogas20240130'
output_folder = 'test_new_bigoas_database'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)
if not os.path.exists(os.path.join(output_folder, '0_proteomes')):
    os.mkdir(os.path.join(output_folder, '0_proteomes'))
if not os.path.exists(os.path.join(output_folder, '1_clustering')):
    os.mkdir(os.path.join(output_folder, '1_clustering'))
if not os.path.exists(os.path.join(output_folder, '2_annotation')):
    os.mkdir(os.path.join(output_folder, '2_annotation'))

df = pd.read_csv(os.path.join(esmecata_workflow_folder, '0_proteomes', 'proteome_tax_id.tsv'), sep='\t')

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
proteome_description_folder = os.path.join(esmecata_workflow_folder, '0_proteomes', 'proteomes_description')
for proteome_description_file in os.listdir(proteome_description_folder):
    with open(os.path.join(proteome_description_folder, proteome_description_file), 'r') as open_proteome_description_file:
            csvreader = csv.DictReader(open_proteome_description_file, delimiter='\t')
            for line in csvreader:
                if line['reference_proteome'] == 'False':
                    reference_proteome = False
                else:
                    reference_proteome = True
                proteome_status[line['proteome_id']] = reference_proteome
compute_stat_proteomes(os.path.join(esmecata_workflow_folder, '0_proteomes'), stat_file=os.path.join(output_folder, '0_proteomes', 'stat_number_proteome.tsv'))

shutil.copyfile(os.path.join(esmecata_workflow_folder, '0_proteomes', 'esmecata_metadata_proteomes.json'), os.path.join(output_folder, '0_proteomes', 'esmecata_metadata_proteomes.json'))
shutil.copyfile(os.path.join(esmecata_workflow_folder, '1_clustering', 'esmecata_metadata_clustering.json'), os.path.join(output_folder, '1_clustering', 'esmecata_metadata_clustering.json'))
shutil.copyfile(os.path.join(esmecata_workflow_folder, '2_annotation', 'esmecata_metadata_annotation.json'), os.path.join(output_folder, '2_annotation', 'esmecata_metadata_annotation.json'))

df_stat_clustering = pd.read_csv(os.path.join(esmecata_workflow_folder, '1_clustering', 'stat_number_clustering.tsv'), sep='\t')

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
        shutil.copyfile(os.path.join(esmecata_workflow_folder, '1_clustering', 'computed_threshold', observation_name+'.tsv'), os.path.join(computed_threshold_folder, tax_id_name+'.tsv'))
        shutil.copyfile(os.path.join(esmecata_workflow_folder, '1_clustering', 'reference_proteins_consensus_fasta', observation_name+'.faa'), os.path.join(reference_proteins_consensus_fasta_folder, tax_id_name+'.faa'))
        already_processed_tax_id_names.append(tax_id_name)
compute_stat_clustering(os.path.join(output_folder, '1_clustering'), stat_file=os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv'))

shutil.copyfile(os.path.join(esmecata_workflow_folder, '2_annotation', 'stat_number_annotation.tsv'), os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv'))

shutil.copytree(os.path.join(esmecata_workflow_folder, '2_annotation', 'annotation_reference'), os.path.join(output_folder, '2_annotation', 'annotation_reference'))
