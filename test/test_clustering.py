import csv
import os
import shutil

from esmecata.clustering import make_clustering, filter_protein_cluster
from esmecata.utils import is_valid_path, is_valid_dir

RESULTS = {
    'Cluster_1': {'Number_shared_proteins': 460}
}

def test_filter_protein_cluster():
    output_folder = 'output'
    observation_name = 'Cluster_1'
    protein_clusters = {'Q89AE4': ['Q89AE4', 'P57473'],
                        'Q89AY7': ['Q89AY7']}
    observation_name_proteomes = ['clustering_input/result/Cluster_1/UP000000601.faa.gz', 'clustering_input/result/Cluster_1/UP000001806.faa.gz']
    clust_threshold = 0.95

    expected_protein = {'Q89AE4'}
    protein_cluster_to_keeps = filter_protein_cluster(observation_name, protein_clusters, observation_name_proteomes, output_folder, clust_threshold)

    assert expected_protein == protein_cluster_to_keeps

    clust_threshold = 0
    expected_protein = {'Q89AE4', 'Q89AY7'}
    protein_cluster_to_keeps = filter_protein_cluster(observation_name, protein_clusters, observation_name_proteomes, output_folder, clust_threshold)

    assert expected_protein == protein_cluster_to_keeps


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
