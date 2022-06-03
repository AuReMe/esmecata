import gzip
import os
import shutil
import urllib.request
import subprocess

from Bio.ExPASy import Enzyme
from Bio import SeqIO


def download_input_file():
    """ Download Expsasy enzyme data file and Swissprot sequences.
    """
    expasy_url = 'https://ftp.expasy.org/databases/enzyme/enzyme.dat'
    urllib.request.urlretrieve(expasy_url, 'enzyme.dat')
    uniprot_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
    urllib.request.urlretrieve(uniprot_url, 'uniprot_sprot.fasta.gz')

def create_consensus_ec(nb_cpu=1):
    """ Create consensus sequences for each EC using MMseqs2.
    """
    # Extract protein sequences from file.
    with gzip.open('uniprot_sprot.fasta.gz', 'rt') as handle:
        uniprot_fasta = {record.id.split('|')[1]: record for record in SeqIO.parse(handle, 'fasta')}

    # Parse enzyme dat file.
    enzyme_uniprots = {}
    with open('enzyme.dat') as infile:
        enzymes = Enzyme.parse(infile)
        # For eac henzyme, get UniProt IDs.
        # Then using the sequences of these IDs, make a clustering with MMseqs2.
        # Retrieve the consensus sequences and create a fasta file.
        for enzyme in enzymes:
            enzyme_id = enzyme['ID']
            print(enzyme_id)
            enzyme_uniprots[enzyme_id] = [ids[0] for ids in enzyme['DR']]
            records = [uniprot_fasta[ids[0]] for ids in enzyme['DR']]
            mmseqs_consensus_fasta = os.path.join('output_fasta', enzyme_id+'.fasta')
            if not os.path.exists(mmseqs_consensus_fasta) and len(records) > 0: 
                SeqIO.write(records, 'test.fasta', 'fasta')

                mmseqs_tmp_cluster = 'mmseqs_tmp'
                if not os.path.exists(mmseqs_tmp_cluster):
                    os.mkdir(mmseqs_tmp_cluster)
                mmseqs_tmp_db = os.path.join(mmseqs_tmp_cluster, 'db')
                mmseqs_tmp_db_clustered = os.path.join(mmseqs_tmp_cluster, 'cluster_db')
                mmseqs_tmp_cluster_tmp = os.path.join(mmseqs_tmp_cluster, 'cluster_tmp')
                mmseqs_seq_db =  os.path.join(mmseqs_tmp_cluster, 'cluster_seq')
                mmseqs_profile =  os.path.join(mmseqs_tmp_cluster, 'cluster_profile')
                mmseqs_consensus =  os.path.join(mmseqs_tmp_cluster, 'cluster_consensus')
                subprocess.call(['mmseqs', 'createdb', 'test.fasta', mmseqs_tmp_db, '-v', '2'])
                
                # Cluster the protein sequences.
                cluster_cmd = ['mmseqs']
                cluster_cmd += ['cluster']
                cluster_cmd += [mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_tmp_cluster_tmp, '--threads', str(nb_cpu), '-v', '2', '-s', '7.5']

                subprocess.call(cluster_cmd)

                subprocess.call(['mmseqs', 'createsubdb', mmseqs_tmp_db_clustered, mmseqs_tmp_db, mmseqs_seq_db, '-v', '2'])
                # Create the profile from the clustering.
                subprocess.call(['mmseqs', 'result2profile', mmseqs_seq_db, mmseqs_tmp_db, mmseqs_tmp_db_clustered, mmseqs_profile, '--threads', str(nb_cpu), '-v', '2'])
                # Create the consensus from the profile.
                subprocess.call(['mmseqs', 'profile2consensus', mmseqs_profile, mmseqs_consensus, '--threads', str(nb_cpu), '-v', '2'])
                # Create the consensus fasta file.
                subprocess.call(['mmseqs', 'convert2fasta', mmseqs_consensus, mmseqs_consensus_fasta, '-v', '2'])
                shutil.rmtree(mmseqs_tmp_cluster)

def merge_consensus():
    """ From the output folder of create_consensus_ec, create a single fasta files containing all consensus sequences.
    """
    records = []
    for fasta_file in os.listdir('output_fasta'):
        fasta_path = os.path.join('output_fasta', fasta_file)
        counter = 0
        for record in SeqIO.parse(fasta_path, 'fasta'):
            record.id = 'EC_number_'+fasta_file.replace('.fasta', '') + '_' + str(counter)
            counter += 1
            records.append(record)
    SeqIO.write(records, 'expasy_consensus.fasta', 'fasta')

download_input_file()
create_consensus_ec()
merge_consensus()