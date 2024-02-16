# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import csv
import gzip
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import shutil
import subprocess
import sys
import time

from Bio import SeqIO
from Bio import __version__ as biopython_version
from shutil import which

from esmecata import __version__ as esmecata_version
from esmecata.utils import is_valid_path, is_valid_dir

logger = logging.getLogger(__name__)


def compute_stat_clustering(output_folder, stat_file=None):
    """Compute stat associated to the number of proteome for each taxonomic affiliations.

    Args:
        output_folder (str): pathname to the result folder containing mmseqs results
        stat_file (str): pathname to the tsv stat file

    Returns:
        clustering_numbers (dict): dict containing observation names (as key) associated with the number of protein clusters
    """
    clustering_taxon_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')
    result_folder = os.path.join(output_folder, 'computed_threshold')

    proteomes_taxa_names = get_proteomes_tax_name(clustering_taxon_id_file)

    tax_name_clustering_numbers = {}
    for clustering_file in os.listdir(result_folder):
        clustering_file_path = os.path.join(result_folder, clustering_file)
        num_lines = sum(1 for line in open(clustering_file_path))
        tax_name_clustering_numbers[clustering_file.replace('.tsv', '')] = num_lines

    clustering_numbers = {}
    for observation_name in proteomes_taxa_names:
        tax_name = proteomes_taxa_names[observation_name].replace(' ', '_')
        clustering_numbers[observation_name] = tax_name_clustering_numbers[tax_name]

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_protein_clusters'])
            for observation_name in clustering_numbers:
                csvwriter.writerow([observation_name, clustering_numbers[observation_name]])

    return clustering_numbers


def create_proteome_representativeness_lineplot(proteome_tax_id_file, computed_threshold_folder, output_figure_file):
    """ From the computed threshold folder and the proteomes_tax_id file created a lineplot.
    This figures show the representativeness ratio and the number of associated protein clusters according to the taxonomic rank.

    Args:
        proteomes_taxon_id_file (str): pathname to the proteomes_tax_id file.
        computed_threshold_folder (str): pathname to computed threshold folder.
        output_figure_file (str): pathname to the output figure file.
    """
    # From proteomes_tax_id file get the taxonomic rank for each observation name.
    proteome_df = pd.read_csv(proteome_tax_id_file, sep='\t')
    proteome_df.set_index('observation_name', inplace=True)
    tax_rank = proteome_df['tax_rank'].to_dict()

    data = []
    for tsv_file in os.listdir(computed_threshold_folder):
        obs_name = os.path.splitext(tsv_file)[0]
        tsv_file_path = os.path.join(computed_threshold_folder, tsv_file)
        tmp_df = pd.read_csv(tsv_file_path, sep='\t')
        for tmp_threshold in np.arange(0, 1.01, 0.025):
            # Compute the number of protein clusters associated with representativeness ratio of tmp_threshold.
            nb_protein_cluster_ratio = len(tmp_df[tmp_df['cluster_ratio']>= tmp_threshold])
            data.append([obs_name, obs_name, tmp_threshold, nb_protein_cluster_ratio])

    df = pd.DataFrame(data, columns=['obs_name', 'rank', 'clust', 'count'])

    fig, ax = plt.subplots(figsize=(9,5))
    g_results =sns.lineplot(data=df, x="clust", y='count', hue="rank", errorbar='ci')
    g_results.set(yscale='log')
    # Use max + 0.25 quantile to set ylim.
    g_results.set(ylim=[1, df['count'].max()+df['count'].quantile(0.25)], xlim=[0,1])
    plt.savefig(output_figure_file)


def create_proteome_representativeness_lineplot_per_taxon_rank(proteome_tax_id_file, computed_threshold_folder, output_folder):
    """ From the computed threshold folder and the proteomes_tax_id file created a lineplot.
    This figures show the representativeness ratio and the number of associated protein clusters according to the taxonomic rank.

    Args:
        proteomes_taxon_id_file (str): pathname to the proteomes_tax_id file.
        computed_threshold_folder (str): pathname to computed threshold folder.
        output_folder (str): pathname to the output folder.
    """
    proteome_df = pd.read_csv(proteome_tax_id_file, sep='\t')
    proteome_df.set_index('name', inplace=True)
    tax_ranks = proteome_df['tax_rank'].to_dict()
    tax_ids = proteome_df['tax_id'].to_dict()

    data = []
    for tsv_file in os.listdir(computed_threshold_folder):
        tax_name = os.path.splitext(tsv_file)[0].replace('_', ' ')
        tsv_file_path = os.path.join(computed_threshold_folder, tsv_file)
        tmp_df = pd.read_csv(tsv_file_path, sep='\t')
        for tmp_threshold in np.arange(0, 1.01, 0.025):
            nb_protein_cluster_ratio = len(tmp_df[tmp_df['cluster_ratio']>= tmp_threshold])
            data.append([tax_name, tax_ranks[tax_name], tax_ids[tax_name], tmp_threshold, nb_protein_cluster_ratio])

    df = pd.DataFrame(data, columns=['tax_name', 'tax_rank', 'tax_id', 'clust', 'count'])

    for rank in df['tax_rank'].unique():
        tmp_df = df[df['tax_rank']==rank]
        fig, ax = plt.subplots(figsize=(9,5))
        g_results =sns.lineplot(data=tmp_df, x="clust", y='count', hue="tax_name", errorbar='ci')
        g_results.set(yscale='log')
        g_results.set(ylim=[1, tmp_df['count'].max()+df['count'].quantile(0.25)], xlim=[0,1])
        output_figure = os.path.join(output_folder, 'representativeness_ratio_{0}.svg'.format(rank))
        plt.savefig(output_figure)
        plt.clf()


def get_proteomes_tax_name(proteomes_taxon_id_file):
    """ Extract tax_id associated with observation name.

    Args:
        proteomes_taxon_id_file (str): pathname to the proteomes_tax_id file.

    Returns:
        proteomes_taxa_names (dict): dict containing observation names (as key) associated with tax name used for proteomes (as value)
    """
    proteomes_taxa_names = {}
    with open(proteomes_taxon_id_file, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            tax_name = line['name']
            proteomes_taxa_names[observation_name] = tax_name

    return proteomes_taxa_names

def copy_already_clustered_file(output_folder, already_clustered_obs_name, new_observation_name, file_extension='.tsv'):
    """ Copy files from an observation name with the same tax ID than the new ones.
    Args:
        output_folder (str): path to the ouput result folder
        already_clustered_obs_name (str): name of the already clustered observation name
        new_observation_name (str): name of the new observation name
        file_extension (str): file extension (by default '.tsv')
    """
    already_cluster_output_file = os.path.join(output_folder, already_clustered_obs_name + file_extension)
    cluster_output_file = os.path.join(output_folder, new_observation_name + file_extension)
    shutil.copyfile(already_cluster_output_file, cluster_output_file)


def run_mmseqs(observation_name, observation_name_proteomes, mmseqs_tmp_path, nb_cpu, mmseqs_options, linclust):
    """Run MMseqs2 on proteomes for an observation name

    Args:
        observation_name (str): observation name associated to a taxonomic affiliations
        observation_name_proteomes (list): list of pathname to each proteomes associated to the observation_name
        mmseqs_tmp_path (str): pathname to the folder which will contain mmseqs results
        nb_cpu (int): number of CPU for mmseqs
        mmseqs_options (str): options that will be used by mmseqs (default: '--min-seq-id', '0.3', '-c', '0.8')
        linclust (bool): use linclust for faster clustering

    Returns:
        mmseqs_tmp_clustered_tabulated (str): pathname to the mmseqs tabulated file containg protein cluster
        mmseqs_tmp_representative_fasta (str): pathname to the mmseqs fasta file containing representative protein sequences
        mmseqs_consensus_fasta (str): pathname to the mmseqs fasta file containing consensus protein sequences
    """
    mmseqs_tmp_cluster = os.path.join(mmseqs_tmp_path, observation_name)
    is_valid_dir(mmseqs_tmp_cluster)

    mmseqs_tmp_cluster_output = os.path.join(mmseqs_tmp_cluster, 'cluster')
    mmseqs_tmp_clustered_tabulated = mmseqs_tmp_cluster_output+'_cluster.tsv'
    mmseqs_tmp_representative_fasta = mmseqs_tmp_cluster_output + '_rep_seq.fasta'
    mmseqs_consensus_fasta = mmseqs_tmp_cluster_output + '_con_seq.fasta'
    # Run mmmseqs to find protein clusters.

    # Code using mmseqs database and mmseqs modules to cluster instead of easy-cluster.
    mmseqs_tmp_db = os.path.join(mmseqs_tmp_cluster, 'db')
    mmseqs_tmp_db_clustered = os.path.join(mmseqs_tmp_cluster, 'cluster_db')
    mmseqs_tmp_cluster_tmp = os.path.join(mmseqs_tmp_cluster, 'cluster_tmp')
    mmseqs_seq_db =  os.path.join(mmseqs_tmp_cluster, 'cluster_seq')
    mmseqs_profile =  os.path.join(mmseqs_tmp_cluster, 'cluster_profile')
    mmseqs_consensus =  os.path.join(mmseqs_tmp_cluster, 'cluster_consensus')

    if not os.path.exists(mmseqs_tmp_clustered_tabulated):
        # Create database containing the protein sequences from all the proteomes of a taxon.
        subprocess.call(['mmseqs', 'createdb', *observation_name_proteomes, mmseqs_tmp_db, '-v', '2'])

        # Cluster the protein sequences.
        cluster_cmd = ['mmseqs']
        if linclust:
            cluster_cmd += ['linclust']
        else:
            cluster_cmd += ['cluster']

        cluster_cmd += [mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_tmp_cluster_tmp, '--threads', str(nb_cpu), '-v', '2']

        # If no option given by the user, use the default options: '--min-seq-id', '0.3', '-c', '0.8'.
        # Sequence identity of 30% and coverage of 80%.
        if not mmseqs_options:
            cluster_cmd += ['--min-seq-id', '0.3', '-c', '0.8']
        else:
            cluster_cmd += mmseqs_options.split(' ')

        subprocess.call(cluster_cmd)

        # Create sequence database with representative from the clustered proteins.
        subprocess.call(['mmseqs', 'createsubdb', mmseqs_tmp_db_clustered, mmseqs_tmp_db, mmseqs_seq_db, '-v', '2'])
        # Create the profile from the clustering.
        subprocess.call(['mmseqs', 'result2profile', mmseqs_seq_db, mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_profile, '--threads', str(nb_cpu), '-v', '2'])
        # Create the consensus from the profile.
        subprocess.call(['mmseqs', 'profile2consensus', mmseqs_profile, mmseqs_consensus, '--threads', str(nb_cpu), '-v', '2'])
        # Create the consensus fasta file.
        subprocess.call(['mmseqs', 'convert2fasta', mmseqs_consensus, mmseqs_consensus_fasta, '-v', '2'])
        # Create TSV resulting files to be analysed after.
        subprocess.call(['mmseqs', 'createtsv', mmseqs_tmp_db, mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_tmp_clustered_tabulated, '--threads', str(nb_cpu), '-v', '2'])
        # Create fasta file containing representative proteins (representatives are the first protein of the alignement).
        subprocess.call(['mmseqs', 'convert2fasta', mmseqs_seq_db, mmseqs_tmp_representative_fasta, '-v', '2'])

    # Old clustering method used with easy-cluster.
    """
    if not os.path.exists(mmseqs_tmp_clustered_tabulated):
        # Code using easy-cluster to extract representative protein.
        subprocess.call(['mmseqs', 'easy-cluster', *cluster_fasta_files[cluster], mmseqs_tmp_cluster_output, mmseqs_tmp_cluster, '--threads', str(nb_cpu), '-v', '2', '--min-seq-id', '0.3'])
    """

    return mmseqs_tmp_clustered_tabulated, mmseqs_tmp_representative_fasta, mmseqs_consensus_fasta


def extrat_protein_cluster_from_mmseqs(mmseqs_tmp_clustered_tabulated, cluster_proteomes_output_file):
    """Extract protein cluster from mmseqs into a dictionary

    Args:
        mmseqs_tmp_clustered_tabulated (str): pathname to the mmseqs tabulated file containg protein cluster
        cluster_proteomes_output_file (str): pathname to the output file showing the content of each protein cluster

    Returns:
        protein_clusters (dict): protein clusters found by mmseqs (representative protein as key and all the protein in the cluster as value)
    """
    # Extract protein clusters in a dictionary.
    # The representative protein is the key and all the proteins in the cluster are the values.
    protein_clusters = {}
    with open(mmseqs_tmp_clustered_tabulated, 'r') as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        for row in csvreader:
            if row[0] not in protein_clusters:
                protein_clusters[row[0]] = [row[1]]
            else:
                protein_clusters[row[0]].append(row[1])

    with open(cluster_proteomes_output_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        for rep_protein in protein_clusters:
            csvwriter.writerow([rep_protein, *[prot for prot in protein_clusters[rep_protein]]])

    return protein_clusters


def compute_proteome_representativeness_ratio(protein_clusters, observation_name_proteomes, computed_threshold_file=None):
    """Compute for each protein cluster the ratio of representation of each proteomes in the cluster.

    Args:
        protein_clusters (dict): protein clusters found by mmseqs (representative protein as key and all the protein in the cluster as value)
        observation_name_proteomes (list): list of pathname to each proteomes associated to the observation_name
        computed_threshold_file (str): pathname to the output file containing the computed ratio

    Returns:
        number_proteomes (int): number of the proteomes for the observation_name
        rep_prot_organims (dict): names of each proteomes associated with each protein cluster
        computed_threshold_cluster (dict): ratio of proteomes representativeness for each protein cluster
    """
    # Retrieve protein ID and the corresponding proteome.
    organism_prots = {}
    for fasta_file in observation_name_proteomes:
        with gzip.open(fasta_file, 'rt') as fasta_handle:
            for record in SeqIO.parse(fasta_handle, 'fasta'):
                compressed_filebasename = os.path.basename(fasta_file)
                fasta_filebasename = os.path.splitext(compressed_filebasename)[0]
                filebasename = os.path.splitext(fasta_filebasename)[0]
                organism_prots[record.id.split('|')[1]] = filebasename

    number_proteomes = len(observation_name_proteomes)

    # Compute the ratio between the number of proteome represented by a protein in the cluster and the total number of proteome.
    rep_prot_organims = {}
    computed_threshold_cluster = {}
    for rep_protein in protein_clusters:
        rep_prot_organims[rep_protein] = set([organism_prots[prot] for prot in protein_clusters[rep_protein]])
        computed_threshold_cluster[rep_protein] = len(rep_prot_organims[rep_protein]) / number_proteomes

    if computed_threshold_file:
        # Create a tsv file containing the computer threshold (number of organism in the cluster compared to the total number of organism) for each organisms.
        with open(computed_threshold_file, 'w') as output_file:
            csvwriter = csv.writer(output_file, delimiter='\t')
            csvwriter.writerow(['representative_protein', 'cluster_ratio', 'proteomes'])
            for rep_protein in rep_prot_organims:
                csvwriter.writerow([rep_protein, computed_threshold_cluster[rep_protein], ','.join(rep_prot_organims[rep_protein])])

    return number_proteomes, rep_prot_organims, computed_threshold_cluster


def filter_protein_cluster(protein_clusters, number_proteomes, rep_prot_organims, computed_threshold_cluster,
                           clust_threshold, cluster_proteomes_filtered_output_file=None):
    """Filter protein cluster according to the representation of each proteomes in the cluster.

    Args:
        protein_clusters (dict): protein clusters found by mmseqs (representative protein as key and all the protein in the cluster as value)
        number_proteomes (int): number of the proteomes for the observation_name
        rep_prot_organims (dict): names of each proteomes associated with each protein cluster
        computed_threshold_cluster (dict): ratio of proteomes representativeness for each protein cluster
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        cluster_proteomes_filtered_output_file (str): pathname to the output file showing the protein cluster kept after filtering

    Returns:
        protein_cluster_to_keeps (list): list containing representative protein IDs associated with cluster that are kept
    """
    # Keep a protein cluster according to the ratio and create a list of representative proteins (being the protein cluster to keep).
    protein_cluster_to_keeps = []
    reference_threshold = (clust_threshold * number_proteomes) / number_proteomes

    for rep_protein in rep_prot_organims:
        if computed_threshold_cluster[rep_protein] >= reference_threshold:
            protein_cluster_to_keeps.append(rep_protein)

    if cluster_proteomes_filtered_output_file:
        with open(cluster_proteomes_filtered_output_file, 'w') as filtered_output_file:
            filtered_csvwriter = csv.writer(filtered_output_file, delimiter='\t')
            for rep_protein in protein_cluster_to_keeps:
                filtered_csvwriter.writerow([rep_protein, *[prot for prot in protein_clusters[rep_protein]]])

    # Use set for faster search using 'in'.
    protein_cluster_to_keeps = set(protein_cluster_to_keeps)

    return protein_cluster_to_keeps


def make_clustering(proteome_folder, output_folder, nb_cpu, clust_threshold, mmseqs_options, linclust, remove_tmp):
    """From the proteomes found by esmecata proteomes, create protein cluster for each taxonomic affiliations.

    Args:
        proteome_folder (str): pathname to folder from esmecata folder
        output_folder (str): pathname to the output folder
        nb_cpu (int): number of CPUs to be used by mmseqs
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster
        mmseqs_options (str): use alternative mmseqs option
        linclust (bool): use linclust
        remove_tmp (bool): remove the tmp files
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|clustering| Begin clustering.')

    # Check if mmseqs is in path.
    mmseqs_path = which('mmseqs')
    if not mmseqs_path:
        logger.critical('|EsMeCaTa|clustering| mmseqs not available in path, esmecata will not be able to cluster the proteomes.')
        sys.exit(1)

    if not is_valid_dir(proteome_folder):
        logger.critical('|EsMeCaTa|clustering| Input must be a folder %s.', proteome_folder)
        sys.exit(1)

    # Use the proteomes folder created by retrieve_proteome.py.
    proteome_tax_id_pathname = os.path.join(proteome_folder, 'proteome_tax_id.tsv')

    if not is_valid_path(proteome_tax_id_pathname):
        logger.critical(f"|EsMeCaTa|clustering| Missing output from esmecata proteomes in {proteome_tax_id_pathname}.")
        sys.exit(1)

    reference_proteins_path = os.path.join(output_folder, 'reference_proteins')
    is_valid_dir(reference_proteins_path)

    cluster_founds_path = os.path.join(output_folder, 'cluster_founds')
    is_valid_dir(cluster_founds_path)

    computed_threshold_path = os.path.join(output_folder, 'computed_threshold')
    is_valid_dir(computed_threshold_path)

    already_performed_clustering = [reference_protein_file.replace('.tsv', '') for reference_protein_file in os.listdir(reference_proteins_path)]

    # Create a dictionary with observation_name as key and the pathname to the proteomes associated to this observation_name as value.
    observation_name_fasta_files = {}
    with open(proteome_tax_id_pathname, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            proteomes = line['proteome'].split(',')
            tax_name = line['name']
            proteomes_path = [os.path.join(proteome_folder, 'proteomes', proteome+'.faa.gz') for proteome in proteomes]
            if tax_name not in already_performed_clustering:
                observation_name_fasta_files[tax_name] = proteomes_path
            else:
                logger.info('|EsMeCaTa|clustering| Already performed clustering for %s.', tax_name)

    is_valid_dir(output_folder)

    # Create metadata file.
    clustering_metadata = {}
    clustering_metadata['tool_options'] = {'proteome_folder': proteome_folder, 'output_folder': output_folder, 'nb_cpu':nb_cpu,
                                        'clust_threshold':clust_threshold, 'mmseqs_options': mmseqs_options, 'linclust':linclust,
                                        'remove_tmp': remove_tmp}

    clustering_metadata['tool_dependencies'] = {}
    subprocess_output = subprocess.check_output(['mmseqs', 'version'])
    mmseqs_version = subprocess_output.decode('utf-8')
    clustering_metadata['tool_dependencies']['mmseqs_version'] = mmseqs_version
    clustering_metadata['tool_dependencies']['mmseqs_path'] = mmseqs_path
    clustering_metadata['tool_dependencies']['python_package'] = {}
    clustering_metadata['tool_dependencies']['python_package']['Python_version'] = sys.version
    clustering_metadata['tool_dependencies']['python_package']['biopython'] = biopython_version
    clustering_metadata['tool_dependencies']['python_package']['esmecata'] = esmecata_version

    # Create tmp folder for mmseqs analysis.
    mmseqs_tmp_path = os.path.join(output_folder, 'mmseqs_tmp')
    is_valid_dir(mmseqs_tmp_path)

    # Create output folder containing shared representative proteins.
    reference_proteins_representative_fasta_path = os.path.join(output_folder, 'reference_proteins_representative_fasta')
    is_valid_dir(reference_proteins_representative_fasta_path)

    reference_proteins_consensus_fasta_path = os.path.join(output_folder, 'reference_proteins_consensus_fasta')
    is_valid_dir(reference_proteins_consensus_fasta_path)

    proteome_taxon_id_file = os.path.join(proteome_folder, 'proteome_tax_id.tsv')
    clustering_taxon_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    if os.path.exists(clustering_taxon_id_file):
        if not os.path.samefile(proteome_taxon_id_file, clustering_taxon_id_file):
            os.remove(clustering_taxon_id_file)
            shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)
    else:
        shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)

    proteomes_taxa_names = get_proteomes_tax_name(proteome_taxon_id_file)

    all_tax_names = set(list(proteomes_taxa_names.values()))

    # For each OTU run mmseqs easy-cluster on them to found the clusters that have a protein in each proteome of the OTU.
    # We take the representative protein of a cluster if the cluster contains a protein from all the proteomes of the OTU.
    # If this condition is not satisfied the cluster will be ignored.
    # Then a fasta file containing all the representative proteins for each OTU is written in representative_fasta folder.
    for proteomes_tax_name in all_tax_names:
        # Get proteomes associated with taxon name.
        observation_name_proteomes = observation_name_fasta_files[proteomes_tax_name]

        # Change space with '_' to avoid issue.
        proteomes_tax_name = proteomes_tax_name.replace(' ', '_')
        # If the computed threshold file exists, mmseqs has already been run.
        mmseqs_tmp_cluster = os.path.join(mmseqs_tmp_path, proteomes_tax_name)
        # Run mmseqs on organism.
        # Delete previous mmseqs2 run if it exists to avoid overwritting issues.
        if os.path.exists(mmseqs_tmp_cluster):
            shutil.rmtree(mmseqs_tmp_cluster)
        mmseqs_tmp_clustered_tabulated, mmseqs_tmp_representative_fasta, mmseqs_consensus_fasta = run_mmseqs(proteomes_tax_name, observation_name_proteomes, mmseqs_tmp_path, nb_cpu, mmseqs_options, linclust)

        # Extract protein clusters from mmseqs results.
        cluster_proteomes_output_file = os.path.join(cluster_founds_path, proteomes_tax_name+'.tsv')
        protein_clusters = extrat_protein_cluster_from_mmseqs(mmseqs_tmp_clustered_tabulated, cluster_proteomes_output_file)

        # Compute proteome representativeness ratio.
        computed_threshold_file = os.path.join(computed_threshold_path, proteomes_tax_name+'.tsv')
        number_proteomes, rep_prot_organims, computed_threshold_cluster = compute_proteome_representativeness_ratio(protein_clusters,
                                                                                                                    observation_name_proteomes, computed_threshold_file)

        # Filter protein cluster for each protein cluster.
        cluster_proteomes_filtered_output_file = os.path.join(reference_proteins_path, proteomes_tax_name+'.tsv')
        protein_cluster_to_keeps = filter_protein_cluster(protein_clusters, number_proteomes, rep_prot_organims, computed_threshold_cluster,
                                                        clust_threshold, cluster_proteomes_filtered_output_file)

        logger.info('|EsMeCaTa|clustering| %d protein clusters kept for %s.', len(protein_cluster_to_keeps), proteomes_tax_name)

        # Create BioPython records with the representative proteins kept.
        new_records = [record for record in SeqIO.parse(mmseqs_tmp_representative_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Do not create fasta file when 0 sequences were kept.
        if len(new_records) > 0:
            # Create output proteome file for OTU.
            representative_fasta_file = os.path.join(reference_proteins_representative_fasta_path, proteomes_tax_name+'.faa')
            SeqIO.write(new_records, representative_fasta_file, 'fasta')
        else:
            logger.info('|EsMeCaTa|clustering| 0 protein clusters %s, no fasta created.', proteomes_tax_name)
        del new_records

        # Create BioPython records with the consensus proteins kept.
        consensus_new_records = [record for record in SeqIO.parse(mmseqs_consensus_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Do not create fasta file when 0 sequences were kept.
        if len(consensus_new_records) > 0:
            # Create output proteome file for OTU.
            consensus_fasta_file = os.path.join(reference_proteins_consensus_fasta_path, proteomes_tax_name+'.faa')
            SeqIO.write(consensus_new_records, consensus_fasta_file, 'fasta')
        else:
            logger.info('|EsMeCaTa|clustering| 0 protein clusters %s, no fasta created.', proteomes_tax_name)
        del consensus_new_records

        if remove_tmp:
            shutil.rmtree(mmseqs_tmp_cluster)

    # Compute number of protein clusters kept.
    stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    compute_stat_clustering(output_folder, stat_file)
    output_figure_file = os.path.join(output_folder, 'representativeness_clustering_ratio.svg')
    create_proteome_representativeness_lineplot(clustering_taxon_id_file, computed_threshold_path, output_figure_file)

    proteome_ratio_lineplots_path = os.path.join(output_folder, 'proteome_ratio_lineplots')
    is_valid_dir(proteome_ratio_lineplots_path)
    create_proteome_representativeness_lineplot_per_taxon_rank(clustering_taxon_id_file, computed_threshold_path, proteome_ratio_lineplots_path)

    endtime = time.time()
    duration = endtime - starttime
    clustering_metadata['esmecata_clustering_duration'] = duration
    clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering.json')
    if os.path.exists(clustering_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_clustering' in metadata_file]
        clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering_{0}.json'.format(len(metadata_files)))
        with open(clustering_metadata_file, 'w') as ouput_file:
            json.dump(clustering_metadata, ouput_file, indent=4)
    else:
        with open(clustering_metadata_file, 'w') as ouput_file:
            json.dump(clustering_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|clustering| Clustering complete in {0}s.'.format(duration))
