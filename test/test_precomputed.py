import csv
import os
import shutil
import subprocess
import hashlib

from esmecata.core.precomputed import precomputed_parse_affiliation


EXPECTED_RESULTS = {
    'Cluster_1': {'proteomes': 18, 'protein_clusters': 573, 'GOs': 2557, 'ECs': 367}
}

EXPECTED_RESULTS_TWO_DB = {
    'Cluster_1': {'proteomes': 18, 'protein_clusters': 573, 'GOs': 2557, 'ECs': 367},
    'Cluster_2': {'proteomes': 5, 'protein_clusters': 708, 'GOs': 1439, 'ECs': 296}
}

EXPECTED_MD5_BUCHNERA = "501f1308d5954f796155004b1e6bfd44"
EXPECTED_MD5_SECOND = "0a6331f5ea548b1ba6853cfb4c760602"

def test_precomputed_parse_affiliation_offline():
    output_folder = 'output_folder'
    input_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    esmecata_precomputed_db = os.path.join('test_data', 'buchnera_database.zip')
    # Check md5 database.
    md5_database = hashlib.md5(open(esmecata_precomputed_db,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_BUCHNERA

    precomputed_parse_affiliation(input_file, esmecata_precomputed_db, output_folder)

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    predicted_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name] = {}
            predicted_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            predicted_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS:
        for data in EXPECTED_RESULTS[observation_name]:
            assert predicted_results[observation_name][data] == EXPECTED_RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


def test_precomputed_parse_affiliation_cli_offline():
    output_folder = 'output_folder'
    input_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    esmecata_precomputed_db = os.path.join('test_data', 'buchnera_database.zip')
    # Check md5 database.
    md5_database = hashlib.md5(open(esmecata_precomputed_db,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_BUCHNERA

    subprocess.call(['esmecata', 'precomputed', '-i', input_file, '-d', esmecata_precomputed_db, '-o', output_folder])

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    predicted_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name] = {}
            predicted_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            predicted_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS:
        for data in EXPECTED_RESULTS[observation_name]:
            assert predicted_results[observation_name][data] == EXPECTED_RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


def test_precomputed_parse_affiliation_offline_two_databases():
    output_folder = 'output_folder'
    input_file = os.path.join('test_data', 'two_db_precomputed.tsv')
    esmecata_path_buchnera_database = os.path.join('test_data', 'buchnera_database.zip')
    esmecata_path_second_database = os.path.join('test_data', 'second_database.zip')
    esmecata_precomputed_db = esmecata_path_buchnera_database + ' ' + esmecata_path_second_database
    # Check md5 database.
    md5_database = hashlib.md5(open(esmecata_path_buchnera_database,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_BUCHNERA
    md5_database = hashlib.md5(open(esmecata_path_second_database,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_SECOND

    precomputed_parse_affiliation(input_file, esmecata_precomputed_db, output_folder)

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    predicted_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name] = {}
            predicted_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            predicted_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS_TWO_DB:
        for data in EXPECTED_RESULTS_TWO_DB[observation_name]:
            assert predicted_results[observation_name][data] == EXPECTED_RESULTS_TWO_DB[observation_name][data]

    shutil.rmtree(output_folder)


def test_precomputed_parse_affiliation_two_databases_cli_offline():
    output_folder = 'output_folder'
    input_file = os.path.join('test_data', 'two_db_precomputed.tsv')
    esmecata_path_buchnera_database = os.path.join('test_data', 'buchnera_database.zip')
    esmecata_path_second_database = os.path.join('test_data', 'second_database.zip')
    esmecata_precomputed_db = esmecata_path_buchnera_database + ' ' + esmecata_path_second_database
    # Check md5 database.
    md5_database = hashlib.md5(open(esmecata_path_buchnera_database,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_BUCHNERA
    md5_database = hashlib.md5(open(esmecata_path_second_database,'rb').read()).hexdigest()
    assert md5_database == EXPECTED_MD5_SECOND

    subprocess.call(['esmecata', 'precomputed', '-i', input_file, '-d', esmecata_precomputed_db, '-o', output_folder])

    proteome_tax_id_file = os.path.join(output_folder, '1_clustering', 'proteome_tax_id.tsv')
    clustering_stat_file = os.path.join(output_folder, '1_clustering', 'stat_number_clustering.tsv')
    annotation_stat_file = os.path.join(output_folder, '2_annotation', 'stat_number_annotation.tsv')

    predicted_results = {}
    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name] = {}
            predicted_results[observation_name]['proteomes'] = len(line['proteome'].split(','))

    with open(clustering_stat_file, 'r') as open_clustering_stat_file:
        csvreader = csv.DictReader(open_clustering_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['protein_clusters'] = int(line['Number_protein_clusters_kept'])

    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            predicted_results[observation_name]['GOs'] = int(line['Number_go_terms'])
            predicted_results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS_TWO_DB:
        for data in EXPECTED_RESULTS_TWO_DB[observation_name]:
            assert predicted_results[observation_name][data] == EXPECTED_RESULTS_TWO_DB[observation_name][data]

    organism_not_found_in_database_file_path = os.path.join(output_folder, 'organism_not_found_in_database.tsv')
    missed_observation_names = []
    with open(organism_not_found_in_database_file_path, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            missed_observation_names.append(observation_name)

    # All input observation names should be found.
    assert missed_observation_names == []

    shutil.rmtree(output_folder)


if __name__ == "__main__":
    test_precomputed_parse_affiliation_offline()
    test_precomputed_parse_affiliation_cli_offline()
