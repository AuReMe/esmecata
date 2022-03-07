import csv
import os
import shutil

from esmecata.clustering import make_clustering

RESULTS = {
    'Cluster_1': {'Number_shared_proteins': 460}
}



def test_make_clustering():
    output_folder = 'clustering_output'
    make_clustering('clustering_input', output_folder, nb_cpu=1, clust_threshold=0.95, mmseqs_options=None, linclust=None, remove_tmp=None)

    expected_results = {}
    output_stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.reader(stat_file_read, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            expected_results[line[0]] = {}
            expected_results[line[0]]['Number_shared_proteins'] = int(line[1])

    for observation_name in expected_results:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


if __name__ == "__main__":
    test_make_clustering()
