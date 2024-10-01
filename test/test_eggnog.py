import shutil
import subprocess
import json
import os
import csv
from esmecata.core.eggnog import annotate_with_eggnog

EXPECTED_RESULTS = {
    'Cluster_1': {'proteomes': 18, 'protein_clusters': 573, 'GOs': 2557, 'ECs': 367}
}

def fake_creation_eggnog(input_path, taxon_name, output_dir, temporary_dir, eggnog_database_path, nb_core, no_dbmem=False, taxon_id_scope=None):
    output_folder = 'output_folder'

    eggnog_folder_expected = os.path.join('annotation_expected', 'eggnog_output')
    for expected_file in os.listdir(eggnog_folder_expected):
        shutil.copyfile(os.path.join(eggnog_folder_expected, expected_file), os.path.join(output_folder, 'eggnog_output', expected_file))
 
    eggnog_fasta_expected = os.path.join('annotation_expected', 'merge_fasta', '2.faa')
    ouput_fasta_expected = os.path.join(output_folder, 'merge_fasta', '2.faa')
    shutil.copyfile(eggnog_fasta_expected, ouput_fasta_expected)


def test_annotate_with_eggnog_mocked(mocker):
    eggnog_version = "2.1.9"
    # Mock call to eggnog-mapper by creating expected files.
    mocker.patch("esmecata.core.eggnog.get_eggnog_version", return_value=eggnog_version)
    mocker.patch("esmecata.core.eggnog.call_to_emapper", wraps=fake_creation_eggnog)
    annotate_with_eggnog('annotation_input', 'output_folder', 'eggnog_database', 1)

    results = {}
    annotation_stat_file = os.path.join('output_folder', 'stat_number_annotation.tsv')
    with open(annotation_stat_file, 'r') as open_annotation_stat_file:
        csvreader = csv.DictReader(open_annotation_stat_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            results[observation_name] = {}
            results[observation_name]['GOs'] = int(line['Number_go_terms'])
            results[observation_name]['ECs'] = int(line['Number_ecs'])

    for observation_name in EXPECTED_RESULTS:
        for data in results[observation_name]:
            assert results[observation_name][data] == EXPECTED_RESULTS[observation_name][data]

    eggnog_metadata_file = os.path.join('output_folder', 'esmecata_metadata_annotation.json')
    with open(eggnog_metadata_file, 'r') as open_json_file:
        json_data = json.load(open_json_file)
        assert json_data['tool_options']['tool_dependencies']['eggnog_mapper'] == eggnog_version
    shutil.rmtree('output_folder')
