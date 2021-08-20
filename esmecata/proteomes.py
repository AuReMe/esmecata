import csv
import gzip
import json
import math
import os
import pandas as pd
import random
import requests
import shutil
import time

from collections import OrderedDict
from ete3 import NCBITaxa, is_taxadb_up_to_date
from esmecata.utils import get_uniprot_release

def associate_taxon_to_taxon_id(taxonomies, ncbi):
    tax_id_names = {}
    clusters_dicts = {}
    json_cluster_taxons = {}

    if taxonomies == {}:
        print('Empty taxonomy dictionary.')
        return

    # For each taxonomy find the taxonomy ID corresponding to each taxon.
    for cluster in taxonomies:
        if isinstance(taxonomies[cluster], str):
            taxons = [taxon for taxon in taxonomies[cluster].split(';')]
            names = ncbi.get_name_translator(taxons)
            taxon_ids = OrderedDict()
            for taxon in taxons:
                taxon_translations = ncbi.get_name_translator([taxon])
                if taxon_translations == {}:
                    print('For {0}, no taxon ID has been found associated to the taxon "{1}" in the NCBI taxonomy of ete3.'.format(cluster, taxon))
                    taxon_translations = {taxon: ['not_found']}
                taxon_ids.update(taxon_translations)
            invert_names = {v: k for k, vs in names.items() for v in vs }
            tax_id_names.update(invert_names)
            clusters_dicts[cluster] = [names[tax_name] for tax_name in taxonomies[cluster].split(';') if tax_name in names]
            json_cluster_taxons[cluster] = taxon_ids

    return tax_id_names, json_cluster_taxons

def filter_taxon(json_cluster_taxons, ncbi):
    # If there is multiple taxon ID for a taxon, use the taxonomy to find the most relevant taxon.
    # The most relevant taxon is the one with the most overlapping lineage with the taxonomy.
    taxon_to_modify = {}
    for cluster in json_cluster_taxons:
        cluster_taxons = json_cluster_taxons[cluster].values()
        for taxon in json_cluster_taxons[cluster]:
            taxon_ids = json_cluster_taxons[cluster][taxon]
            if len(taxon_ids) > 1:
                relevant_taxon = None
                previous_nb_shared_ids = 0
                for taxon_id in taxon_ids:
                    taxons_ids = list(next(zip(*cluster_taxons)))
                    lineage = ncbi.get_lineage(taxon_id)
                    nb_shared_ids = len(set(taxons_ids).intersection(set(lineage)))
                    if nb_shared_ids > previous_nb_shared_ids:
                        relevant_taxon = taxon_id
                    previous_nb_shared_ids = nb_shared_ids
                taxon_to_modify[cluster] = {}
                taxon_to_modify[cluster][taxon] = [relevant_taxon]

    for cluster in taxon_to_modify:
        for taxon in taxon_to_modify[cluster]:
            json_cluster_taxons[cluster][taxon] = taxon_to_modify[cluster][taxon]

    return json_cluster_taxons

def find_proteomes_tax_ids(json_cluster_taxons, ncbi, busco_percentage_keep=None):
    # Query the Uniprot proteomes to find all the proteome IDs associated to taxonomy.
    # If there is more thant 100 proteomes we do not keep it because there is too many proteome.
    print('Find proteome ID associated to taxonomy')
    proteomes_ids = {}
    single_proteomes = {}
    tax_id_not_founds = {}
    tax_id_founds = {}
    for taxon in json_cluster_taxons:
        for tax_name in reversed(json_cluster_taxons[taxon]):
            tax_id = json_cluster_taxons[taxon][tax_name][0]

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in tax_id_founds:
                proteomes_ids[taxon] = (tax_id, tax_id_founds[tax_id])
                if len(tax_id_founds[tax_id]) == 1:
                    single_proteomes[taxon] = (tax_id, proteomes)
                break

            # If tax_id has not been found with a request do not try a new request with the same tax_id.
            if tax_id in tax_id_not_founds:
                continue

            # Find proteomes associated with taxon.
            # Take reference proteomes with "reference:yes".
            # Avoid redundant and excluded proteomes with "redundant%3Ano+excluded%3Ano".
            # Use "format=tab" to easily handle the ouput.
            http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+reference:yes+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)

            response = requests.get(http_str)
            # Raise error if we have a bad request.
            response.raise_for_status()
            if response.text == '':
                time.sleep(1)
                print('{0}: No reference proteomes found for {1} try non-reference proteomes.'.format(taxon, tax_id))
                http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)
                response = requests.get(http_str)
            if response.text != '':
                csvreader = csv.reader(response.text.splitlines(), delimiter='\t')

                proteomes = []
                organism_ids = {}
                # Avoid header.
                next(csvreader)
                for line in csvreader:
                    proteome = line[0]
                    completness = line [6]
                    org_tax_id = line[2]

                    # Check that proteome has busco score.
                    if busco_percentage_keep:
                        if line[4] != '':
                            busco_percentage = float(line[4].split(':')[1].split('%')[0])
                            if busco_percentage >= busco_percentage_keep and completness == 'full':
                                proteomes.append(proteome)
                                if org_tax_id not in organism_ids:
                                    organism_ids[org_tax_id] = [proteome]
                                else:
                                    organism_ids[org_tax_id].append(proteome)
                    else:
                        if completness == 'full':
                            proteomes.append(proteome)
                            if org_tax_id not in organism_ids:
                                organism_ids[org_tax_id] = [proteome]
                            else:
                                organism_ids[org_tax_id].append(proteome)

                if len(proteomes) > 0 and len(proteomes) < 100:
                    print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(proteomes)))
                    proteomes_ids[taxon] = (tax_id, proteomes)
                    tax_id_founds[tax_id] = proteomes
                    if len(proteomes) == 1:
                        single_proteomes[taxon] = (tax_id, proteomes)
                    break

                elif len(proteomes) >= 100:
                    print('More than 99 proteomes are associated to the taxa {0} associated to {1}, esmecata will randomly select around 100 proteomes with respect to the taxonomy proportion.'.format(taxon, tax_name))
                    tree = ncbi.get_topology([org_tax_id for org_tax_id in organism_ids])

                    # For each direct descendant taxon of the tree root (our tax_id), we will look for the proteomes inside these subtaxons.
                    childs = {}
                    for parent_node in tree.get_children():
                        parent_tax_id = str(parent_node.taxid)
                        if parent_node.get_descendants() != []:
                            for child_node in parent_node.get_descendants():
                                child_tax_id = str(child_node.taxid)
                                if parent_tax_id not in childs:
                                    if child_tax_id in organism_ids:
                                        childs[parent_tax_id] = organism_ids[child_tax_id]
                                else:
                                    if child_tax_id in organism_ids:
                                        childs[parent_tax_id].extend(organism_ids[child_tax_id])
                        else:
                            if parent_tax_id in organism_ids:
                                childs[parent_tax_id] = organism_ids[parent_tax_id]
                            else:
                                childs[parent_tax_id].extend(organism_ids[parent_tax_id])

                    # For each direct descendant taxon, compute the number of proteomes in their childrens.
                    elements = [v for k,v in childs.items()]
                    elements_counts = [len(element) for element in elements]
                    # Compute the proportion of proteomes for each direct descendant taxon compare to the total number of proteomes in all the descendant taxons.
                    percentages = [(i/sum(elements_counts))*100  for i in elements_counts]
                    # Can be superior to 100 if there is a lot of data with around 0.xxx percentage.
                    percentages_round = [math.ceil(percentage) if percentage < 1 else math.floor(percentage) for percentage in percentages]

                    # Choose randomly a number of proteomes corresponding to the computed percentage.
                    all_proteomes = []
                    for index, element in enumerate(elements):
                        percentage_to_keep = percentages_round[index]
                        proteomes_to_keep = random.sample(element, percentage_to_keep)
                        all_proteomes.extend(proteomes_to_keep)
                    proteomes_ids[taxon] = (tax_id, all_proteomes)
                    tax_id_founds[tax_id] = all_proteomes
                    print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(all_proteomes)))
                    break

            # Answer is empty no corresponding proteomes to the tax_id.
            else:
                if tax_id not in tax_id_not_founds:
                    tax_id_not_founds[tax_id] = [tax_name]
                else:
                    tax_id_not_founds[tax_id].append(tax_name)
                continue
            time.sleep(1)

    return proteomes_ids, single_proteomes, tax_id_not_founds


def retrieve_proteomes(input_folder, output_folder, busco_percentage_keep=None, ignore_taxadb_update=None):
    if is_taxadb_up_to_date() is False:
        print('''WARNING: ncbi taxonomy database is not up to date with the last NCBI Taxonomy. Update it using:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()''')
        if ignore_taxadb_update is None:
            print('If you want to stil use esmecata with the old taxonomy database use the option --ignore-taxadb-update/ignore_taxadb_update.')
            return
        else:
            print('--ignore-taxadb-update/ignore_taxadb_update option detected, esmecata will continue with this version.')

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    if '.xlsx' in input_folder:
        df = pd.read_excel(input_folder)
    if '.tsv' in input_folder:
        df = pd.read_csv(input_folder, sep='\t')
    if '.csv' in input_folder:
        df = pd.read_csv(input_folder, sep=',')

    # Set index on the column containing the OTU name.
    df.set_index('observation_name', inplace=True)

    # taxonomy is the column containing the taxonomy separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomy']

    ncbi = NCBITaxa()

    tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(taxonomies, ncbi)

    json_cluster_taxons = filter_taxon(json_cluster_taxons, ncbi)

    json_log = os.path.join(output_folder, 'log.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_cluster_taxons, ouput_file, indent=4)

    proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_cluster_taxons, ncbi, busco_percentage_keep)

    proteome_to_download = []
    for proteomes_id in proteomes_ids:
        proteome_to_download.extend(proteomes_ids[proteomes_id][1])
    proteome_to_download = set(proteome_to_download)

    proteome_cluster_tax_id_file = os.path.join(output_folder, 'proteome_cluster_tax_id.tsv')
    # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
    with open(proteome_cluster_tax_id_file, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        csvwriter.writerow(['cluster', 'name', 'tax_id', 'proteome'])
        for cluster in proteomes_ids:
            csvwriter.writerow([cluster, tax_id_names[int(proteomes_ids[cluster][0])], proteomes_ids[cluster][0], ','.join(proteomes_ids[cluster][1])])

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

    # Download Uniprot metadata and create a json file containing them.
    uniprot_releases = get_uniprot_release()

    uniprot_metadata_file = os.path.join(output_folder, 'uniprot_release_metadata.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    print('Creating result folder')
    # Create a result folder which contains one sub-folder per OTU.
    # Each OTU sub-folder will contain the proteome found.
    result_folder = os.path.join(output_folder, 'result')
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)

    result_single_folder = os.path.join(output_folder, 'result_single_proteome')
    if not os.path.exists(result_single_folder):
        os.mkdir(result_single_folder)

    for cluster in proteomes_ids:
        if cluster in single_proteomes:
            output_cluster = os.path.join(result_single_folder, cluster)
        else:
            output_cluster = os.path.join(result_folder, cluster)
        if not os.path.exists(output_cluster):
            os.mkdir(output_cluster)
        for proteome in proteomes_ids[cluster][1]:
            input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
            output_proteome = os.path.join(output_cluster, proteome+'.faa')
            with gzip.open(input_proteome_file, 'rb') as f_in:
                with open(output_proteome, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
