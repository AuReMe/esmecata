# Copyright (C) 2021-2022 Arnaud Belcour - Inria Dyliss
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
import os
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


def compute_stat_clustering(result_folder, stat_file=None):
    """Compute stat associated to the number of proteome for each taxonomic affiliations.

    Args:
        result_folder (str): pathname to the result folder containing mmseqs results
        stat_file (str): pathname to the tsv stat file

    Returns:
        clustering_numbers (dict): dict containing observation names (as key) associated with the number of protein clusters
    """
    clustering_numbers = {}
    for clustering_file in os.listdir(result_folder):
        clustering_file_path = os.path.join(result_folder, clustering_file)
        num_lines = sum(1 for line in open(clustering_file_path))
        clustering_numbers[clustering_file.replace('.tsv', '')] = num_lines

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_protein_clusters'])
            for observation_name in clustering_numbers:
                csvwriter.writerow([observation_name, clustering_numbers[observation_name]])

    return clustering_numbers


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
    # Run mmmseqs to find rapidly protein clusters.

    # Code using mmseqs database and mmseqs modules to cluster instead of easy-cluster.
    mmseqs_tmp_db = os.path.join(mmseqs_tmp_cluster, 'db')
    #mmseqs_tmp_db_h = os.path.join(mmseqs_tmp_cluster, 'db_h')
    mmseqs_tmp_db_clustered = os.path.join(mmseqs_tmp_cluster, 'cluster_db')
    mmseqs_tmp_cluster_tmp = os.path.join(mmseqs_tmp_cluster, 'cluster_tmp')
    mmseqs_seq_db =  os.path.join(mmseqs_tmp_cluster, 'cluster_seq')
    #mmseqs_seq_db_h =  os.path.join(mmseqs_tmp_cluster, 'cluster_seq_h')
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

        subprocess.call(['mmseqs', 'createsubdb', mmseqs_tmp_db_clustered, mmseqs_tmp_db, mmseqs_seq_db, '-v', '2'])
        #subprocess.call(['mmseqs', 'createsubdb', mmseqs_tmp_db_clustered, mmseqs_tmp_db_h, mmseqs_seq_db_h, '-v', '2'])
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


def extrat_protein_cluster_from_mmseqs(mmseqs_tmp_clustered_tabulated):
    """Extract protein cluster from mmseqs into a dictionary

    Args:
        mmseqs_tmp_clustered_tabulated (str): pathname to the mmseqs tabulated file containg protein cluster

    Returns:
        protein_clusters (dict): protein clusters found by mmseqs (representative protein as key and all the protein in the cluster as value)
    """
    # Extract protein clusters in a dictionary.
    # The representative protein is the key and all the proteins in the cluster are the values.
    protein_clusters = {}
    with open(mmseqs_tmp_clustered_tabulated) as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        for row in csvreader:
            if row[0] not in protein_clusters:
                protein_clusters[row[0]] = [row[1]]
            else:
                protein_clusters[row[0]].append(row[1])

    return protein_clusters


def filter_protein_cluster(observation_name, protein_clusters, observation_name_proteomes, output_folder, clust_threshold):
    """Filter protein cluster according to the representation of each proteomes in the cluster

    Args:
        observation_name (str): observation name associated to a taxonomic affiliations
        protein_clusters (dict): protein clusters found by mmseqs (representative protein as key and all the protein in the cluster as value)
        observation_name_proteomes (list): list of pathname to each proteomes associated to the observation_name
        output_folder (str): pathname to the output folder
        clust_threshold (float): threshold to select protein cluster according to the representation of protein proteome in the cluster

    Returns:
        protein_cluster_to_keeps (list): list containing representative protein IDs associated with cluster that are kept
    """
    reference_proteins_path = os.path.join(output_folder, 'reference_proteins')
    is_valid_dir(reference_proteins_path)

    cluster_founds_path = os.path.join(output_folder, 'cluster_founds')
    is_valid_dir(cluster_founds_path)

    computed_threshold_path = os.path.join(output_folder, 'computed_threshold')
    is_valid_dir(computed_threshold_path)
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

    # Create a tsv file contianing the computer threshold (number of organism in the cluster compared to the total number of organism) for each organisms.
    computed_threshold_file = os.path.join(computed_threshold_path, observation_name+'.tsv')
    with open(computed_threshold_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['representative_protein', 'cluster_ratio', 'proteomes'])
        for rep_protein in rep_prot_organims:
            csvwriter.writerow([rep_protein, computed_threshold_cluster[rep_protein], ','.join(rep_prot_organims[rep_protein])])

    # Keep a protein cluster according to the ratio and create a list of representative proteins (being the protein cluster to keep).
    protein_cluster_to_keeps = []
    cluster_proteomes_filtered_output_file = os.path.join(reference_proteins_path, observation_name+'.tsv')
    cluster_proteomes_output_file = os.path.join(cluster_founds_path, observation_name+'.tsv')
    reference_threshold = (clust_threshold * number_proteomes) / number_proteomes
    with open(cluster_proteomes_filtered_output_file, 'w') as filtered_output_file, open(cluster_proteomes_output_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        filtered_csvwriter = csv.writer(filtered_output_file, delimiter='\t')
        for rep_protein in rep_prot_organims:
            csvwriter.writerow([rep_protein, *[prot for prot in protein_clusters[rep_protein]]])
            if computed_threshold_cluster[rep_protein] >= reference_threshold:
                filtered_csvwriter.writerow([rep_protein, *[prot for prot in protein_clusters[rep_protein]]])
                protein_cluster_to_keeps.append(rep_protein)

    # Use set for faster search using 'in'.
    protein_cluster_to_keeps = set(protein_cluster_to_keeps)

    stat_file = os.path.join(output_folder, 'stat_number_clustering.tsv')
    compute_stat_clustering(reference_proteins_path, stat_file)

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

    # Check if mmseqs is in path.
    mmseqs_path = which('mmseqs')
    if not mmseqs_path:
        logger.critical('|EsMeCaTa|clustering| mmseqs not available in path, esmecata will not be able to cluster the proteomes.')
        sys.exit(1)

    if not is_valid_dir(proteome_folder):
        logger.critical(f"|EsMeCaTa|clustering| Input must be a folder {proteome_folder}.")
        sys.exit(1)

    # Use the result folder created by retrieve_proteome.py.
    result_folder = os.path.join(proteome_folder, 'result')
    if not is_valid_path(result_folder):
        logger.critical(f"|EsMeCaTa|clustering| Missing output from esmecata proteomes in {result_folder}.")
        sys.exit(1)

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

    # Create a dictionary with observation_name as key and the pathname to the proteomes associated to this observation_name as value.
    observation_name_fasta_files = {}
    for observation_name in os.listdir(result_folder):
        result_observation_name_folder = os.path.join(result_folder, observation_name)
        if is_valid_dir(result_observation_name_folder):
            for fasta_file in os.listdir(result_observation_name_folder):
                fasta_filename, compressed_file_extension = os.path.splitext(fasta_file)
                _, file_extension = os.path.splitext(fasta_filename)
                complete_extension = file_extension + compressed_file_extension
                if complete_extension == '.faa.gz':
                    fasta_file_path = os.path.join(result_observation_name_folder, fasta_file)
                    if observation_name not in observation_name_fasta_files:
                        observation_name_fasta_files[observation_name] = [fasta_file_path]
                    else:
                        observation_name_fasta_files[observation_name].append(fasta_file_path)

    # Create tmp folder for mmseqs analysis.
    mmseqs_tmp_path = os.path.join(output_folder, 'mmseqs_tmp')
    is_valid_dir(mmseqs_tmp_path)

    # Create output folder containing shared representative proteins.
    representative_fasta_path = os.path.join(output_folder, 'fasta_representative')
    is_valid_dir(representative_fasta_path)

    consensus_fasta_path = os.path.join(output_folder, 'fasta_consensus')
    is_valid_dir(consensus_fasta_path)

    reference_proteins_representative_fasta_path = os.path.join(output_folder, 'reference_proteins_representative_fasta')
    is_valid_dir(reference_proteins_representative_fasta_path)

    reference_proteins_consensus_fasta_path = os.path.join(output_folder, 'reference_proteins_consensus_fasta')
    is_valid_dir(reference_proteins_consensus_fasta_path)

    logger.info('|EsMeCaTa|clustering| Clustering proteins.')
    # For each OTU run mmseqs easy-cluster on them to found the clusters that have a protein in each proteome of the OTU.
    # We take the representative protein of a cluster if the cluster contains a protein from all the proteomes of the OTU.
    # If this condition is not satisfied the cluster will be ignored.
    # Then a fasta file containing all the representative proteins for each OTU is written in representative_fasta folder.
    for observation_name in observation_name_fasta_files:
        mmseqs_tmp_cluster = os.path.join(mmseqs_tmp_path, observation_name)
        observation_name_proteomes = observation_name_fasta_files[observation_name]
        mmseqs_tmp_clustered_tabulated, mmseqs_tmp_representative_fasta, mmseqs_consensus_fasta = run_mmseqs(observation_name, observation_name_proteomes, mmseqs_tmp_path, nb_cpu, mmseqs_options, linclust)
        protein_clusters = extrat_protein_cluster_from_mmseqs(mmseqs_tmp_clustered_tabulated)
        protein_cluster_to_keeps = filter_protein_cluster(observation_name, protein_clusters, observation_name_proteomes, output_folder, clust_threshold)

        logger.info('|EsMeCaTa|clustering| {0} protein clusters kept for {1}.'.format(len(protein_cluster_to_keeps), observation_name))

        # Copy representative and consensus fasta to subfolder in output_folder.
        shutil.copyfile(mmseqs_tmp_representative_fasta, os.path.join(representative_fasta_path, observation_name+'.faa'))
        shutil.copyfile(mmseqs_consensus_fasta, os.path.join(consensus_fasta_path, observation_name+'.faa'))

        # Create BioPtyhon records with the representative proteins kept.
        new_records = [record for record in SeqIO.parse(mmseqs_tmp_representative_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Create output proteome file for OTU.
        representative_fasta_file = os.path.join(reference_proteins_representative_fasta_path, observation_name+'.faa')
        SeqIO.write(new_records, representative_fasta_file, 'fasta')
        del new_records

        # Create BioPtyhon records with the consensus proteins kept.
        consensus_new_records = [record for record in SeqIO.parse(mmseqs_consensus_fasta, 'fasta') if record.id.split('|')[1] in protein_cluster_to_keeps]

        # Create output proteome file for OTU.
        consensus_fasta_file = os.path.join(reference_proteins_consensus_fasta_path, observation_name+'.faa')
        SeqIO.write(consensus_new_records, consensus_fasta_file, 'fasta')
        del consensus_new_records

        if remove_tmp:
            shutil.rmtree(mmseqs_tmp_cluster)

    proteome_taxon_id_file = os.path.join(proteome_folder, 'proteome_tax_id.tsv')
    clustering_taxon_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    if os.path.exists(clustering_taxon_id_file):
        if not os.path.samefile(proteome_taxon_id_file, clustering_taxon_id_file):
            os.remove(clustering_taxon_id_file)
            shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)
    else:
        shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)

    endtime = time.time()
    duration = endtime - starttime
    clustering_metadata['esmecata_clustering_duration'] = duration
    clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'w') as ouput_file:
        json.dump(clustering_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|clustering| Clustering complete.')
