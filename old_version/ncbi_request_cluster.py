import argparse
import subprocess
import urllib 
import gzip
import shutil
import os
import time

import pandas as pa
import re
import time

import csv
import zipfile
from Bio import Entrez


def ncbi_request(organism_identification, email, ncbi_api_key=None):
    nb_requests = 0
    Entrez.email = email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key
    try:
        search_handle = Entrez.esearch(db="assembly", retmax=100000, term= str(organism_identification) + '[Organism] AND "complete genome"[FILTER]', idtype='acc')
    except urllib.error.HTTPError:
        time.spleep(10)
        search_handle = Entrez.esearch(db="assembly", retmax=100000, term= str(organism_identification) + '[Organism] AND "complete genome"[FILTER]', idtype='acc')
    search_record = Entrez.read(search_handle)

    nb_requests += 1
    if nb_requests == 2:
        time.sleep(1)
        nb_requests = 0

    accessions = []
    assembly_names = []
    for assembly_id in search_record['IdList']:
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id, retmax=100000)
        except urllib.error.HTTPError:
            time.spleep(10)
            handle = Entrez.esummary(db="assembly", id=assembly_id, retmax=100000)
        
        record = Entrez.read(handle)
        accessions.append(record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession'])
        assembly_names.append(record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName'].replace(' ', '_'))
        nb_requests += 1
        if nb_requests == 3:
            time.sleep(1)
            nb_requests = 0

    return accessions, assembly_names

def download_genome(index, accession, output_folder):
    ftp_address = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/' + accession[4:7] + '/' + accession[7:10] + '/' + accession[10:13]
    genome_folder_result = urllib.request.urlopen(ftp_address).read().decode('utf-8')
    genome_folder_name = genome_folder_result.split(' ')[-1].replace('\r\n', '')
    ftp_address = ftp_address + '/' + genome_folder_name
    genomes_files_result = urllib.request.urlopen(ftp_address).read()
    genomes_files = genomes_files_result.decode('utf-8').split('\r\n')
    time.sleep(1)
    for genomes_file in genomes_files:
        if '.gbff.gz' in genomes_file:
            genomes_name_file = genomes_file.split(' ')[-1].replace('\r\n', '')
    ftp_address = ftp_address + '/' + genomes_name_file

    if not os.path.exists(output_folder + '/' + accession):
        gz_file = accession + '.gbk.gz'
        urllib.request.urlretrieve(ftp_address, gz_file)
        os.mkdir(output_folder + '/' + accession)
        with gzip.open(gz_file, 'rb') as intput_gz:
            with open(output_folder + '/' + accession + '/' + gz_file.replace('.gz', ''), 'wb') as output_gbk:
                shutil.copyfileobj(intput_gz, output_gbk)
        os.remove(gz_file)
        time.sleep(1)

def parsing_frogs():

    def get_otu_name(otu_rows, selected_otus):
        current_otu = otu_rows[-1]
        not_usable_otus = ['unknown species', 'unknown genus', 'unknown family',
                            'possible family ', 'Multi-affiliation', 'W5', 'Z20', 'unknown order', '(',
                            'metagenome']
        if any(not_usable_otu in current_otu for not_usable_otu in not_usable_otus):
            selected_otus = get_otu_name(otu_rows[:-1], selected_otus)
            return selected_otus
        else:
            selected_otus.append(current_otu)
            return selected_otus

    df = pa.read_excel('frogs_data.xlsx')

    df = df[df['% du total'] > 0.1]

    otus = {}
    with open('stat_otu_assembly.tsv', 'w') as stat_otu:
        csvwriter = csv.writer(stat_otu, delimiter='\t')
        for index, row in df.iterrows():
            if isinstance(row['blast_taxonomy'], str):
                otus_rows = row['blast_taxonomy'].split(';')
                selected_otus = get_otu_name(otus_rows, [])

                otus_rdp_rows = row['#rdp_tax_and_bootstrap'].split(';')[:-1]
                selected_rdp_otus = get_otu_name(otus_rdp_rows, [])

                otus[row['observation_name']] = (selected_otus[0], selected_rdp_otus[0], otus_rows)
                csvwriter.writerow([row['observation_name'], selected_otus[0], selected_rdp_otus[0]])

    return otus

def get_genomes_per_otu(output_folder, otus, email, api_key):
    too_big_otus = ['Bacteria', 'Archaea']
    already_computed_cluster = []

    known_otus = {}
    if os.path.exists(output_folder):
        if os.path.exists(output_folder + '/otu_genomes.tsv'):
            open_mode = 'a'
        else:
            open_mode = 'w'
        for cluster_folder in os.listdir(output_folder):
            already_computed_cluster.append(cluster_folder)
    else:
        open_mode = 'w'
        os.mkdir(output_folder)

    with open(output_folder + '/otu_genomes.tsv', open_mode) as otu_genome:
        csvwriter = csv.writer(otu_genome, delimiter='\t')
        if open_mode == 'w':
            csvwriter.writerow(['OTU', 'taxonomy', 'used_taxonomy', 'number_all_genomes', 'all_genomes_accessions', 'number_valid_genomes', 'valid_genomes_accessions'])
        for cluster in otus:
            if cluster not in already_computed_cluster:
                otu = otus[cluster][0]
                accessions, assembly_names = ncbi_request(otu, email, api_key)
                if len(accessions) == 0:
                    for otu in reversed(otus[cluster][2]):
                        if otu not in too_big_otus:
                            accessions, assembly_names = ncbi_request(otu, email, api_key)
                            if len(accessions) > 0:
                                break
                        else:
                            otu = None
                if otu:
                    if otu not in known_otus:
                        if not os.path.exists(output_folder + '/' + cluster):
                            os.mkdir(output_folder + '/' + cluster)
                        valid_accessions = []
                        cmds = ['./datasets', 'download', 'assembly', *accessions, '--include_gbff', '--filename', 'ncbi_temp.zip']
                        print(' '.join(cmds))
                        datasets_subprocess = subprocess.Popen(cmds, universal_newlines="")
                        datasets_subprocess.wait()
                        try:
                            with zipfile.ZipFile('ncbi_temp.zip', 'r') as zip_ref:
                                zip_ref.extractall('ncbi_temp')
                        except:
                            print('Error no data for ' + cluster)
                        if os.path.exists('ncbi_temp/ncbi_dataset/data'):
                            for root, genome_dirs, files in os.walk('ncbi_temp/ncbi_dataset/data', topdown=False):
                                for genome_dir in genome_dirs:
                                    gbff_filepath = root + '/' + genome_dir + '/' + 'genomic.gbff'
                                    if os.path.exists(gbff_filepath):
                                        if not os.path.exists(output_folder + '/' + genome_dir):
                                            os.mkdir(output_folder + '/' + cluster + '/' + genome_dir)
                                        output_gbk = output_folder + '/' + cluster + '/' + genome_dir + '/' + genome_dir +'.gbk'
                                        shutil.move(gbff_filepath, output_gbk)
                                        valid_accessions.append(genome_dir)
                                    else:
                                        print('Error no gbff for ' + cluster)
                            shutil.rmtree('ncbi_temp')
                        else:
                            print('Error no ncbi_temp/dataset/data for ' + cluster)

                        csvwriter.writerow([cluster, otus[cluster][2], otu, len(accessions), ','.join(accessions), len(valid_accessions), ','.join(valid_accessions)])
                        known_otus[otu] = cluster
                        os.remove('ncbi_temp.zip')
                        time.sleep(1)
                    else:
                        known_path = output_folder + '/' + known_otus[otu]
                        output_path = output_folder + '/' + cluster
                        accessions = []
                        for folder in os.listdir(known_path):
                            genome_path = os.path.join(known_path, folder)
                            genome_new_path = os.path.join(output_path, folder)
                            shutil.copytree(genome_path, genome_new_path)
                            accessions.append(folder)
                        csvwriter.writerow([cluster, otus[cluster][2], otu, len(accessions), ','.join(accessions), len(accessions), ','.join(accessions)])
                else:
                    csvwriter.writerow([cluster, otus[cluster][2], 'too big OTU', '', '', '', ''])


from tqdm import tqdm
def count_accession_per_otus(output_file, otus, email, api_key):
    too_big_otus = ['Bacteria', 'Archaea']
    already_computed_cluster = []
    if os.path.exists(output_file):
        open_mode = 'a'
        with open(output_file, 'r') as otu_genome:
            csvreader = csv.reader(otu_genome, delimiter='\t')
            for row in csvreader:
                already_computed_cluster.append(row[0])
    else:
        open_mode = 'w'

    with open(output_file, open_mode) as otu_genome:
        csvwriter = csv.writer(otu_genome, delimiter='\t')
        if open_mode == 'w':
            csvwriter.writerow(['OTU', 'taxonomy', 'used_taxonomy', 'number_genomes', 'genomes_accessions'])
        for cluster in tqdm(otus):
            if cluster not in already_computed_cluster:
                otu = otus[cluster][0]
                accessions, assembly_names = ncbi_request(otu, email, api_key)
                if len(accessions) == 0:
                    for otu in reversed(otus[cluster][2]):
                        if otu not in too_big_otus:
                            accessions, assembly_names = ncbi_request(otu, email, api_key)
                            if len(accessions) > 0:
                                break
                        else:
                            otu = None
                if otu:
                    csvwriter.writerow([cluster, otus[cluster][2], otu, len(accessions), ','.join(accessions)])
                else:
                    csvwriter.writerow([cluster, otus[cluster][2], 'too big OTU', '', ''])
                time.sleep(1)

parser = argparse.ArgumentParser(description='Request NCBI to get genomes.')
parser.add_argument('--email', help='email address for NCBI request')
parser.add_argument('--key', help='NCBI API key.')

args = parser.parse_args()
email = args.email
if args.key:
    api_key = args.key
else:
    api_key = None

otus=parsing_frogs()

get_genomes_per_otu('data', otus, email, api_key)

