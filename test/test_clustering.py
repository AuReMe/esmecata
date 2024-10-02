import csv
import os
import shutil
import subprocess

from esmecata.core.clustering import make_clustering, filter_protein_cluster, compute_proteome_representativeness_ratio

RESULTS = {
    'Cluster_1': {'Number_protein_clusters_kept': 603}
}

def test_filter_protein_cluster_offline():
    output_folder = 'output'
    remove_output_folder = False
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        remove_output_folder = True
        os.mkdir(os.path.join(output_folder, 'computed_threshold'))
        os.mkdir(os.path.join(output_folder, 'reference_proteins'))
        os.mkdir(os.path.join(output_folder, 'cluster_founds'))

    observation_name = 'Cluster_1'
    protein_clusters = {'Q89AE4': ['Q89AE4', 'P57473'],
                        'Q89AY7': ['Q89AY7']}
    observation_name_proteomes = [os.path.join('clustering_input', 'proteomes', 'UP000000601.faa.gz'),
                                os.path.join('clustering_input', 'proteomes', 'UP000001806.faa.gz')]
    clust_threshold = 0.95

    expected_protein = {'Q89AE4'}

    number_proteomes, rep_prot_organims, computed_threshold_cluster = compute_proteome_representativeness_ratio(protein_clusters, observation_name_proteomes)
    protein_cluster_to_keeps = filter_protein_cluster(protein_clusters, number_proteomes, rep_prot_organims, computed_threshold_cluster,
                            clust_threshold)
    
    assert expected_protein == protein_cluster_to_keeps

    clust_threshold = 0
    expected_protein = {'Q89AE4', 'Q89AY7'}
    number_proteomes, rep_prot_organims, computed_threshold_cluster = compute_proteome_representativeness_ratio(protein_clusters, observation_name_proteomes)
    protein_cluster_to_keeps = filter_protein_cluster(protein_clusters, number_proteomes, rep_prot_organims, computed_threshold_cluster,
                            clust_threshold)
    assert expected_protein == protein_cluster_to_keeps

    if remove_output_folder is True:
        shutil.rmtree(output_folder)


def test_make_clustering_offline():
    output_folder = 'clustering_output'
    make_clustering('clustering_input', output_folder, nb_core=1, clust_threshold=0.5, mmseqs_options=None, linclust=None, remove_tmp=None)

    expected_results = {}
    output_stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.DictReader(stat_file_read, delimiter='\t')
        for line in csvreader:
            expected_results[line['observation_name']] = {}
            expected_results[line['observation_name']]['Number_protein_clusters_kept'] = int(line['Number_protein_clusters_kept'])

    for observation_name in expected_results:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


def test_clustering_cli_offline():
    output_folder = 'clustering_output'
    subprocess.call(['esmecata', 'clustering', '-i', 'clustering_input', '-o', output_folder, '-c', '1', '-t', '0.5'])
    expected_results = {}
    output_stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.DictReader(stat_file_read, delimiter='\t')
        for line in csvreader:
            expected_results[line['observation_name']] = {}
            expected_results[line['observation_name']]['Number_protein_clusters_kept'] = int(line['Number_protein_clusters_kept'])

    for observation_name in expected_results:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == RESULTS[observation_name][data]

    shutil.rmtree(output_folder)


if __name__ == "__main__":
    test_filter_protein_cluster_offline()
    test_make_clustering_offline()
    test_clustering_cli_offline()
