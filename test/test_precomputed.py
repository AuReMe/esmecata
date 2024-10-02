import csv
import os
import shutil
import subprocess

from esmecata.core.precomputed import precomputed_parse_affiliation


EXPECTED_RESULTS = {
    'Cluster_1': {'proteomes': 18, 'protein_clusters': 573, 'GOs': 2557, 'ECs': 367}
}

def test_precomputed_parse_affiliation_offline():
    output_folder = 'output_folder'
    precomputed_parse_affiliation('buchnera_workflow.tsv', 'buchnera_database.zip', output_folder)

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    expected_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name] = {}
            expected_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            expected_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == EXPECTED_RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


def test_precomputed_parse_affiliation_cli_offline():
    output_folder = 'output_folder'
    subprocess.call(['esmecata', 'precomputed', '-i', 'buchnera_workflow.tsv', '-d', 'buchnera_database.zip', '-o', output_folder])

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    expected_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name] = {}
            expected_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            expected_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            expected_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == EXPECTED_RESULTS[observation_name][data]

    shutil.rmtree(output_folder)

if __name__ == "__main__":
    test_precomputed_parse_affiliation_offline()
    test_precomputed_parse_affiliation_cli_offline()
