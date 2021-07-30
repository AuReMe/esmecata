import pandas as pd
import time
import requests
import os
import csv
import shutil
import gzip
import json

from collections import OrderedDict
from ete3 import NCBITaxa

def retrieve_proteome(input_folder, output_folder, busco_percentage_keep=90):
    busco_percentage_keep = float(busco_percentage_keep)
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    if '.xlsx' in input_folder:
        df = pd.read_excel(input_folder)
    if '.tsv' in input_folder:
        df = pd.read_csv(input_folder, sep='\t')
    if '.csv' in input_folder:
        df = pd.read_csv(input_folder, sep=';')

    # Set index on the column containing the OTU name.
    df.set_index('observation_name', inplace=True)

    # taxonomy is the column containing the taxonomy separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomy']

    ncbi = NCBITaxa()
    tax_ids = []

    tax_id_names = {}
    cluster_ids = {}
    clusters_dicts = {}
    json_cluster_taxons = {}

    # For each taxonomy find the taxonomy ID corresponding to each taxon.
    for cluster in taxonomies:
        if isinstance(taxonomies[cluster], str):
            taxons = [taxon for taxon in taxonomies[cluster].split(';')]
            names = ncbi.get_name_translator(taxons)
            for taxon in taxons:
                if taxon not in names:
                    names[taxon] = ['not_found']
            invert_names = {v: k for k, vs in names.items() for v in vs }
            tax_id_names.update(invert_names)
            for names_tax_ids in list(names.values()):
                tax_ids.extend(names_tax_ids)
            cluster_ids[cluster] = names
            clusters_dicts[cluster] = [names[tax_name] for tax_name in taxonomies[cluster].split(';') if tax_name in names]
            json_cluster_taxons[cluster] = names

    # If there is multiple taxon ID for a taxon, use the taxonomy to find the most relevant taxon.
    # The most relevant taxon is the one with the most overlapping lineage with the taxonomy.
    for cluster in clusters_dicts:
        for index, taxons in enumerate(clusters_dicts[cluster]):
            if len(taxons) > 1:
                relevant_taxon = None
                previous_nb_shared_ids = 0
                for taxon in taxons:
                    taxons_ids = list(next(zip(*clusters_dicts[cluster])))
                    lineage = ncbi.get_lineage(taxon)
                    nb_shared_ids = len(set(taxons_ids).intersection(set(lineage)))
                    if nb_shared_ids > previous_nb_shared_ids:
                        relevant_taxon = taxon
                    previous_nb_shared_ids = nb_shared_ids
                clusters_dicts[cluster][index] = [relevant_taxon]

    json_log = os.path.join(output_folder, 'log.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_cluster_taxons, ouput_file, indent=4)

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
                if busco_percentage > busco_percentage_keep and completness == 'full':
                    proteomes.append(proteome)

        if len(proteomes) > 0 and len(proteomes) < 100:
            proteomes_ids[tax_id] = proteomes
        time.sleep(1)

    # Write for each taxon ID the proteomes found.
    proteome_tax_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')
    with open(proteome_tax_id_file, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        csvwriter.writerow(['tax_id', 'name' ,'proteomes'])
        for tax_id in proteomes_ids:
            csvwriter.writerow([tax_id , tax_id_names[tax_id], ','.join(proteomes_ids[tax_id])])

    proteomes_ids = {}
    with open(proteome_tax_id_file, 'r') as out_file:
        csvwriter = csv.reader(out_file, delimiter='\t')
        next(csvwriter)
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

    proteome_cluster_tax_id_file = os.path.join(output_folder, 'proteome_cluster_tax_id.tsv')
    # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
    with open(proteome_cluster_tax_id_file, 'w') as out_file:
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
    print('Download proteome')
    tmp_folder = os.path.join(output_folder, 'tmp_proteome')
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    for proteome in proteome_to_download:
        output_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
        proteome_response = requests.get('https://www.uniprot.org/uniprot/?query=proteome:{0}&format=fasta&compress=yes'.format(proteome))
        with open(output_proteome_file, 'wb') as f:
            f.write(proteome_response.content)
        time.sleep(1)

    tmp_folder = os.path.join(output_folder, 'tmp_proteome')
    print('Creating result folder')
    # Create a result folder which contains one sub-folder per OTU.
    # Each OTU sub-folder will contain the proteome found.
    result_folder = os.path.join(output_folder, 'result')
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)

    for cluster in cluster_proteomes:
        output_cluster = os.path.join(result_folder, cluster)
        if not os.path.exists(output_cluster):
            os.mkdir(output_cluster)
        for proteome in cluster_proteomes[cluster]:
            output_proteome = os.path.join(output_cluster, proteome+'.faa')
            input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
            with gzip.open(input_proteome_file, 'rb') as f_in:
                with open(output_proteome, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
