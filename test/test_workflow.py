import csv
import os
import subprocess
import shutil

RESULTS = {
    'Cluster_1': {'proteomes': 2, 'protein_clusters': 460, 'GOs': 627, 'ECs': 245}
}


def test_workflow():
    subprocess.call(['esmecata', 'workflow', '-i', 'buchnera_workflow.tsv', '-o', 'test_output'])

    output_stat_file = os.path.join('test_output', 'stat_number_workflow.tsv')

    expected_results = {}
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.reader(stat_file_read, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            expected_results[line[0]] = {}
            expected_results[line[0]]['proteomes'] = int(line[1])
            expected_results[line[0]]['protein_clusters'] = int(line[2])
            expected_results[line[0]]['GOs'] = int(line[3])
            expected_results[line[0]]['ECs'] = int(line[4])

    for cluster in expected_results:
        for data in expected_results[cluster]:
            assert expected_results[cluster][data] == RESULTS[cluster][data]

    shutil.rmtree('test_output')

if __name__ == "__main__":
    test_workflow()