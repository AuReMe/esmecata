import csv
import os
import shutil
import subprocess

from esmecata.core.clustering import make_clustering, filter_protein_cluster, compute_proteome_representativeness_ratio, heap_law_curve_fitting, compute_openess_pan_proteomes

RESULTS = {
    'Cluster_1': {'Number_protein_clusters_kept': 603}
}

def test_filter_protein_cluster_offline():
    output_folder = 'output'
    clustering_input = os.path.join('test_data', 'clustering_input')

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
    observation_name_proteomes = [os.path.join(clustering_input, 'proteomes', 'UP000000601.faa.gz'),
                                os.path.join(clustering_input, 'proteomes', 'UP000001806.faa.gz')]
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
    clustering_input = os.path.join('test_data', 'clustering_input')
    make_clustering(clustering_input, output_folder, nb_core=1, clust_threshold=0.5, mmseqs_options=None, linclust=None, remove_tmp=None)

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
    clustering_input = os.path.join('test_data', 'clustering_input')
    subprocess.call(['esmecata', 'clustering', '-i', clustering_input, '-o', output_folder, '-c', '1', '-t', '0.5'])
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


def test_heap_law_curve_fitting():
    # Test on pangenome from https://doi.org/10.1371/journal.pgen.0030231 by using Table S1
    import pandas as pd
    import random
    iteration_nb = 10000
    input_file = os.path.join('test_data', 'pangenome', 'Prochlorococcus_pangenome.tsv')

    # Extract number of genomes and number newly found genes.
    df = pd.read_csv(input_file, sep='\t')
    clustering_gene_genomes = [df[col].dropna().index.tolist() for col in df.columns]
    list_proteome_to_iter = [random.sample(range(len(clustering_gene_genomes)), k=len(clustering_gene_genomes)) for nb_iter in range(iteration_nb)]
    nb_proteomes = []
    nb_new_protein_discovered= []
    for proteome_list in list_proteome_to_iter:
        protein_discovered = set()
        nb_proteome = 0
        for proteome in proteome_list:
            # Get the protein cluters of the new proteome.
            protein_clusters_associated = set(clustering_gene_genomes[proteome])
            nb_proteome += 1
            if protein_discovered == set():
                # If it is the first proteome, all its protein clusters correspond to newly found protein clusters.
                protein_discovered = protein_discovered.union(protein_clusters_associated)
                new_protein_discovered = protein_discovered
            else:
                # Take previously found protein clusters from the other proteomes and extract how many new protein clusters are added by the new proteome.
                new_protein_discovered = protein_clusters_associated - protein_discovered
                protein_discovered = protein_discovered.union(protein_clusters_associated)
            nb_proteomes.append(nb_proteome)
            nb_new_protein_discovered.append(len(new_protein_discovered))

    # Fit curve.
    k, alpha = heap_law_curve_fitting(nb_proteomes, nb_new_protein_discovered)
    # Results should be inferior to 0.9 and close to 0.8 (0.8 found in https://doi.org/10.1016/j.mib.2008.09.006)
    assert alpha < 0.9
    assert alpha > 0.79

def test_compute_openess_pan_proteomes():
    esmecata_computed_threshold_folder = os.path.join('test_data', 'computed_threshold')
    output_openess_file = 'output_openess.tsv'
    compute_openess_pan_proteomes(esmecata_computed_threshold_folder, output_openess_file)
    os.remove(output_openess_file)

if __name__ == "__main__":
    test_heap_law_curve_fitting()
    test_filter_protein_cluster_offline()
    test_make_clustering_offline()
    test_clustering_cli_offline()
