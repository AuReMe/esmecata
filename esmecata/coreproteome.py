import os
import csv
import subprocess

from Bio import SeqIO

def create_coreproteome(proteome_folder, output_folder, nb_cpu):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Use the result folder created by retrieve_proteome.py.
    result_folder = os.path.join(proteome_folder, 'result')
    cluster_fasta_files = {}
    for cluster in os.listdir(result_folder):
        result_cluster_folder = os.path.join(result_folder, cluster)
        for cluster_file in os.listdir(result_cluster_folder):
            filename, file_extension = os.path.splitext(cluster_file)
            if file_extension == '.faa':
                cluster_file_path = os.path.join(result_cluster_folder, cluster_file)
                if cluster not in cluster_fasta_files:
                    cluster_fasta_files[cluster] = [cluster_file_path]
                else:
                    cluster_fasta_files[cluster].append(cluster_file_path)

    # Create tmp folder for mmseqs analysis.
    coreproteome_tmp = os.path.join(output_folder, 'coreproteome_tmp')
    if not os.path.exists(coreproteome_tmp):
        os.mkdir(coreproteome_tmp)

    # Create output folder containing shared representative proteins.
    coreproteome = os.path.join(output_folder, 'coreproteome')
    if not os.path.exists(coreproteome):
        os.mkdir(coreproteome)

    print('Creating coreproteome')
    # For each OTU run mmseqs easy-cluster on them to found the clusters that have a protein in each proteome of the OTU.
    # We take the representative protein of a cluster if the cluster contains a protein from all the proteomes of the OTU.
    # If this condition is not satisfied the cluster will be ignored.
    # Then a fasta file containing all the representative proteins for each OTU is written in coreproteome folder.
    for cluster in cluster_fasta_files:
        coreproteome_tmp_cluster = os.path.join(coreproteome_tmp, cluster)
        coreproteome_tmp_cluster_output = os.path.join(coreproteome_tmp_cluster, 'cluster')

        # Run mmmseqs to find rapidly protein clusters.
        if not os.path.exists(coreproteome_tmp_cluster):
            os.mkdir(coreproteome_tmp_cluster)
        subprocess.call(['mmseqs', 'easy-cluster', *cluster_fasta_files[cluster], coreproteome_tmp_cluster_output, coreproteome_tmp_cluster, '--threads', str(nb_cpu)])

        # Extract protein clusters.
        proteins_representatives = {}
        with open(coreproteome_tmp_cluster_output+'_cluster.tsv') as input_file:
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
            for record in SeqIO.parse(fasta_file, 'fasta'):
                organism_prots[record.id.split('|')[1]] = fasta_file

        coreproteome_cluster = os.path.join(output_folder, 'reference_proteins')
        if not os.path.exists(coreproteome_cluster):
            os.mkdir(coreproteome_cluster)

        # To keep a cluster, we have to find have at least one protein of each proteome of the OTU.
        rep_prot_to_keeps = []
        rep_prot_organims = {}
        with open(coreproteome_cluster+'/'+cluster+'_coreproteome.tsv', 'w') as output_file:
            csvwriter = csv.writer(output_file, delimiter='\t')
            for rep_protein in proteins_representatives:
                if len(proteins_representatives[rep_protein]) > 1:
                        rep_prot_organims[rep_protein] = set([organism_prots[prot] for prot in proteins_representatives[rep_protein]])
                        if len(rep_prot_organims[rep_protein]) == len(cluster_fasta_files[cluster]):
                            csvwriter.writerow([rep_protein, *[prot for prot in proteins_representatives[rep_protein]]])
                            rep_prot_to_keeps.append(rep_protein)

        # Create BioPtyhon records with the representative proteins kept.
        new_records = []
        for record in SeqIO.parse(coreproteome_tmp_cluster_output + '_rep_seq.fasta', 'fasta'):
            if record.id.split('|')[1] in rep_prot_to_keeps:
                new_records.append(record)

        # Create output proteome file for OTU.
        coreproteome_fasta_fifle = os.path.join(coreproteome, cluster+'.faa')
        SeqIO.write(new_records, coreproteome_fasta_fifle, 'fasta')
