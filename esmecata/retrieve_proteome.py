import pandas as pd
import time
import requests
import os
import csv
import shutil
import gzip

from ete3 import NCBITaxa


df = pd.read_excel('Example.xlsx')

# Set index on the column containing the OTU name.
df.set_index('observation_name', inplace=True)

# taxonomy is the column containing the taxonomy separated by ';': phylum;class;order;family;genus;genus + species
taxonomies = df.to_dict()['taxonomy']

ncbi = NCBITaxa()
tax_ids = []

tax_id_names = {}
cluster_ids = {}
clusters_dicts = {}

# For each taxonomy find the taxonomy ID corresponding to each taxon.
for cluster in taxonomies:
    if isinstance(taxonomies[cluster], str):
        taxons = [taxon for taxon in taxonomies[cluster].split(';')]
        names = ncbi.get_name_translator(taxons)
        invert_names = {v: k for k, vs in names.items() for v in vs }
        tax_id_names.update(invert_names)
        for names_tax_ids in list(names.values()):
            tax_ids.extend(names_tax_ids)
        cluster_ids[cluster] = names
        clusters_dicts[cluster] = [names[tax_name] for tax_name in taxonomies[cluster].split(';') if tax_name in names]

# Query the Uniprot proteomes to find all the proteome IDs associated to taxonomy.
# If there is more thant 100 proteomes we do not keep it because there is too many proteome.
print('Find proteome ID associated to taxonomy')
tax_ids = set(tax_ids)
proteomes_ids = {}
for tax_id in tax_ids:
    http_str = 'https://www.uniprot.org/proteomes/?query=reference:yes+taxonomy:{0}&format=tab'.format(tax_id)

    response = requests.get(http_str)
    csvreader = csv.reader(response.text.splitlines(), delimiter='\t')

    proteomes = []
    # Avoid header.
    next(csvreader)
    for line in csvreader:
        # Check that proteome has busco score.
        if line[4] != '':
            proteome = line[0]
            busco_percentage = float(line[4].split(':')[1].split('%')[0])
            completness = line [6]
            if busco_percentage > 90 and completness == 'full':
                proteomes.append(proteome)

    if len(proteomes) > 0 and len(proteomes) < 100:
        proteomes_ids[tax_id] = proteomes
    time.sleep(1)

# Write for each taxon ID the proteomes found.
with open('proteome_tax_id.tsv', 'w') as out_file:
    csvwriter = csv.writer(out_file, delimiter='\t')
    csvwriter.writerow(['tax_id', 'name' ,'proteomes'])
    for tax_id in proteomes_ids:
        csvwriter.writerow([tax_id , tax_id_names[tax_id], ','.join(proteomes_ids[tax_id])])

proteomes_ids = {}
with open('proteome_tax_id.tsv', 'r') as out_file:
    csvwriter = csv.reader(out_file, delimiter='\t')
    for row in csvwriter:
        proteomes_ids[row[0]] = row[2].split(',')

# Retrieve the proteome to download. Meaning the proteome for the basal taxon in the taxonomy associated with Uniprot proteome.
proteome_to_download = []
cluster_proteomes = {}
cluster_tax_ids = {}
for cluster in clusters_dicts:
    cluster_tax_ids[cluster] = [str(tax_id) for tax_id in clusters_dicts[cluster][-1]]
    for tax_id in clusters_dicts[cluster][-1]:
        if str(tax_id) in proteomes_ids:
            tax_id_proteomes = proteomes_ids[str(tax_id)]
            if cluster not in cluster_proteomes:
                cluster_proteomes[cluster] = tax_id_proteomes
            else:
                cluster_proteomes[cluster].extend(tax_id_proteomes)
            proteome_to_download.extend(tax_id_proteomes)
proteome_to_download = set(proteome_to_download)

# Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
with open('proteome_cluster_tax_id.tsv', 'w') as out_file:
    csvwriter = csv.writer(out_file, delimiter='\t')
    csvwriter.writerow(['cluster', 'tax_id', 'name', 'proteome'])
    for cluster in cluster_tax_ids:
        cluster_tax_id_names = [tax_id_names[int(tax_id)] for tax_id in cluster_tax_ids[cluster]]
        if cluster in cluster_proteomes:
            cluster_tax_id_proteomes = [proteome for proteome in cluster_proteomes[cluster]]
        else:
            cluster_tax_id_proteomes = []
        csvwriter.writerow([cluster , ','.join(cluster_tax_ids[cluster]), ','.join(cluster_tax_id_names), ','.join(cluster_tax_id_proteomes)])

# Download all the proteomes in tmp folder.
print('Downlaod proteome')
output_folder = 'tmp_proteome'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
for proteome in proteome_to_download:
    output_proteome_file = os.path.join(output_folder, proteome+'.faa.gz')
    proteome_response = requests.get('https://www.uniprot.org/uniprot/?query=proteome:{0}&format=fasta&compress=yes'.format(proteome))
    with open(output_proteome_file, 'wb') as f:
        f.write(proteome_response.content)
    time.sleep(1)

print('Creating result folder')
# Create a result folder which contains one sub-folder per OTU.
# Each OTU sub-folder will contain the proteome found.
output_folder = 'result'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
tmp_folder = 'tmp_proteome'
for cluster in cluster_proteomes:
    output_cluster = os.path.join(output_folder, cluster)
    if not os.path.exists(output_cluster):
        os.mkdir(output_cluster)
    for proteome in cluster_proteomes[cluster]:
        output_proteome = os.path.join(output_cluster, proteome+'.faa')
        input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
        with gzip.open(input_proteome_file, 'rb') as f_in:
            with open(output_proteome, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
