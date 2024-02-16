import csv
import os
import subprocess
import shutil

RESULTS = {
    'Cluster_1': {'proteomes': 2, 'protein_clusters': 604, 'GOs': 755, 'ECs': 313}
}


def test_workflow():
    subprocess.call(['esmecata', 'workflow_uniprot', '-i', 'buchnera_workflow.tsv', '-o', 'test_output', '-p', '1', '--minimal-nb-proteomes', '2', '--bioservices'])

    output_stat_file = os.path.join('test_output', 'stat_number_workflow.tsv')

    expected_results = {}
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.reader(stat_file_read, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            expected_results[line[0]] = {}
            expected_results[line[0]]['proteomes'] = int(line[5])
            expected_results[line[0]]['protein_clusters'] = int(line[7])
            expected_results[line[0]]['GOs'] = int(line[8])
            expected_results[line[0]]['ECs'] = int(line[9])

    for observation_name in expected_results:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == RESULTS[observation_name][data]

    shutil.rmtree('test_output')

if __name__ == "__main__":
    test_workflow()
