import csv
import gzip
import json
import os
import shutil
import subprocess
import sys

from Bio import SeqIO
from Bio import __version__ as biopython_version
from shutil import which

from esmecata import __version__ as esmecata_version
from esmecata.utils import is_valid_path, is_valid_dir

def make_clustering(proteome_folder, output_folder, nb_cpu, clust_threshold, mmseqs_options, linclust, remove_tmp):
    mmseqs_path = which('mmseqs')
    if not mmseqs_path:
        print('mmseqs not available in path, esmecata will not be able to cluster the proteomes.')
        sys.exit(1)

    if not is_valid_dir(proteome_folder):
        print(f"Input must be a folder {proteome_folder}.")
        sys.exit(1)

    # Use the result folder created by retrieve_proteome.py.
    result_folder = os.path.join(proteome_folder, 'result')
    if not is_valid_path(result_folder):
        print(f"Missing output from esmecata proteomes in {result_folder}.")
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

    clustering_metadata_file = os.path.join(output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'w') as ouput_file:
        json.dump(clustering_metadata, ouput_file, indent=4)

    cluster_fasta_files = {}
    for cluster in os.listdir(result_folder):
        result_cluster_folder = os.path.join(result_folder, cluster)
        if is_valid_dir(result_cluster_folder):
            for cluster_file in os.listdir(result_cluster_folder):
                fasta_filename, compressed_file_extension = os.path.splitext(cluster_file)
                filename, file_extension = os.path.splitext(fasta_filename)
                complete_extension = file_extension + compressed_file_extension
                if complete_extension == '.faa.gz':
                    cluster_file_path = os.path.join(result_cluster_folder, cluster_file)
                    if cluster not in cluster_fasta_files:
                        cluster_fasta_files[cluster] = [cluster_file_path]
                    else:
                        cluster_fasta_files[cluster].append(cluster_file_path)

    # Create tmp folder for mmseqs analysis.
    mmseqs_tmp = os.path.join(output_folder, 'mmseqs_tmp')
    is_valid_dir(mmseqs_tmp)

    # Create output folder containing shared representative proteins.
    representative_fasta = os.path.join(output_folder, 'fasta_representative')
    is_valid_dir(representative_fasta)

    consensus_fasta = os.path.join(output_folder, 'fasta_consensus')
    is_valid_dir(consensus_fasta)

    reference_proteins = os.path.join(output_folder, 'reference_proteins')
    is_valid_dir(reference_proteins)

    reference_proteins_representative_fasta = os.path.join(output_folder, 'reference_proteins_representative_fasta')
    is_valid_dir(reference_proteins_representative_fasta)

    reference_proteins_consensus_fasta = os.path.join(output_folder, 'reference_proteins_consensus_fasta')
    is_valid_dir(reference_proteins_consensus_fasta)

    cluster_founds = os.path.join(output_folder, 'cluster_founds')
    is_valid_dir(cluster_founds)

    computed_threshold = os.path.join(output_folder, 'computed_threshold')
    is_valid_dir(computed_threshold)

    print('Clustering proteins.')
    # For each OTU run mmseqs easy-cluster on them to found the clusters that have a protein in each proteome of the OTU.
    # We take the representative protein of a cluster if the cluster contains a protein from all the proteomes of the OTU.
    # If this condition is not satisfied the cluster will be ignored.
    # Then a fasta file containing all the representative proteins for each OTU is written in representative_fasta folder.
    for cluster in cluster_fasta_files:
        mmseqs_tmp_cluster = os.path.join(mmseqs_tmp, cluster)
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
            subprocess.call(['mmseqs', 'createdb', *cluster_fasta_files[cluster], mmseqs_tmp_db, '-v', '2'])

            # Cluster the protein sequences.
            cluster_cmd = ['mmseqs']
            if linclust:
                cluster_cmd += ['linclust']
            else:
                cluster_cmd += ['cluster']

            cluster_cmd += [mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_tmp_cluster_tmp, '--threads', str(nb_cpu), '-v', '2']

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

        """
        if not os.path.exists(mmseqs_tmp_clustered_tabulated):
            # Code using easy-cluster to extract representative protein.
            subprocess.call(['mmseqs', 'easy-cluster', *cluster_fasta_files[cluster], mmseqs_tmp_cluster_output, mmseqs_tmp_cluster, '--threads', str(nb_cpu), '-v', '2', '--min-seq-id', '0.3'])
        """

        # Extract protein clusters.
        proteins_representatives = {}
        with open(mmseqs_tmp_clustered_tabulated) as input_file:
            csvreader = csv.reader(input_file, delimiter='\t')
            for row in csvreader:
                if row[0] not in proteins_representatives:
                    proteins_representatives[row[0]] = [row[1]]
                else:
                    proteins_representatives[row[0]].append(row[1])

        # Retrieve protein ID and the corresponding proteome.
        organism_prots = {}
        for fasta_file in cluster_fasta_files[cluster]:
            filename, file_extension = os.path.splitext(fasta_file)
            with gzip.open(fasta_file, 'rt') as fasta_handle:
                for record in SeqIO.parse(fasta_handle, 'fasta'):
                    compressed_filebasename = os.path.basename(fasta_file)
                    fasta_filebasename = os.path.splitext(compressed_filebasename)[0]
                    filebasename = os.path.splitext(fasta_filebasename)[0]
                    organism_prots[record.id.split('|')[1]] = filebasename

        number_proteomes = len(cluster_fasta_files[cluster])

        rep_prot_organims = {}
        computed_threshold_cluster = {}
        for rep_protein in proteins_representatives:
            rep_prot_organims[rep_protein] = set([organism_prots[prot] for prot in proteins_representatives[rep_protein]])
            computed_threshold_cluster[rep_protein] = len(rep_prot_organims[rep_protein]) / number_proteomes

        # Create a tsv file contianing the computer threshold (number of organism in the cluster compared to the total number of organism) for each organisms.
        computed_threshold_file = os.path.join(computed_threshold, cluster+'.tsv')
        with open(computed_threshold_file, 'w') as output_file:
            csvwriter = csv.writer(output_file, delimiter='\t')
            csvwriter.writerow(['representative_protein', 'cluster_ratio', 'proteomes'])
            for rep_protein in rep_prot_organims:
                csvwriter.writerow([rep_protein, computed_threshold_cluster[rep_protein], ','.join(rep_prot_organims[rep_protein])])

        # To keep a cluster, we have to find have at least one protein of each proteome of the OTU (except when using threshold option).
        rep_prot_to_keeps = []
        cluster_proteomes_filtered_output_file = os.path.join(reference_proteins, cluster+'.tsv')
        cluster_proteomes_output_file = os.path.join(cluster_founds, cluster+'.tsv')
        reference_threshold = (clust_threshold * number_proteomes) / number_proteomes
        with open(cluster_proteomes_filtered_output_file, 'w') as filtered_output_file, open(cluster_proteomes_output_file, 'w') as output_file:
            csvwriter = csv.writer(output_file, delimiter='\t')
            filtered_csvwriter = csv.writer(filtered_output_file, delimiter='\t')
            for rep_protein in rep_prot_organims:
                csvwriter.writerow([rep_protein, *[prot for prot in proteins_representatives[rep_protein]]])
                if computed_threshold_cluster[rep_protein] >= reference_threshold:
                    filtered_csvwriter.writerow([rep_protein, *[prot for prot in proteins_representatives[rep_protein]]])
                    rep_prot_to_keeps.append(rep_protein)

        # Use set for faster search using 'in'.
        rep_prot_to_keeps = set(rep_prot_to_keeps)

        shutil.copyfile(mmseqs_tmp_representative_fasta, os.path.join(representative_fasta, cluster+'.faa'))
        shutil.copyfile(mmseqs_consensus_fasta, os.path.join(consensus_fasta, cluster+'.faa'))
        # Create BioPtyhon records with the representative proteins kept.
        new_records = [record for record in SeqIO.parse(mmseqs_tmp_representative_fasta, 'fasta') if record.id.split('|')[1] in rep_prot_to_keeps]

        # Create output proteome file for OTU.
        representative_fasta_file = os.path.join(reference_proteins_representative_fasta, cluster+'.faa')
        SeqIO.write(new_records, representative_fasta_file, 'fasta')
        del new_records

        # Create BioPtyhon records with the consensus proteins kept.
        consensus_new_records = [record for record in SeqIO.parse(mmseqs_consensus_fasta, 'fasta') if record.id.split('|')[1] in rep_prot_to_keeps]

        # Create output proteome file for OTU.
        consensus_fasta_file = os.path.join(reference_proteins_consensus_fasta, cluster+'.faa')
        SeqIO.write(consensus_new_records, consensus_fasta_file, 'fasta')
        del consensus_new_records

        if remove_tmp:
            shutil.rmtree(mmseqs_tmp_cluster)

    proteome_taxon_id_file = os.path.join(proteome_folder, 'proteome_cluster_tax_id.tsv')
    clustering_taxon_id_file = os.path.join(output_folder, 'proteome_cluster_tax_id.tsv')

    if os.path.exists(clustering_taxon_id_file):
        if not os.path.samefile(proteome_taxon_id_file, clustering_taxon_id_file):
            os.remove(clustering_taxon_id_file)
            shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)
    else:
        shutil.copyfile(proteome_taxon_id_file, clustering_taxon_id_file)