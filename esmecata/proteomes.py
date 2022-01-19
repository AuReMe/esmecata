import csv
import gzip
import json
import math
import os
import pandas as pd
import random
import requests
import shutil
import sys
import time

from collections import OrderedDict

from Bio import __version__ as biopython_version
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ete3 import __version__ as ete3_version
from ete3 import NCBITaxa, is_taxadb_up_to_date

from SPARQLWrapper import __version__ as sparqlwrapper_version

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir, send_uniprot_sparql_query
from esmecata import __version__ as esmecata_version

REQUESTS_HEADERS = {'User-Agent': 'EsMeCaTa proteomes v' + esmecata_version + ', request by requests package v' + requests.__version__ }


def associate_taxon_to_taxon_id(taxonomies, ncbi):
    tax_id_names = {}
    clusters_dicts = {}
    json_cluster_taxons = {}

    if taxonomies == {}:
        print('Empty taxonomy dictionary.')
        return

    # For each taxon find the taxon ID corresponding to each taxon.
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
    # If there is multiple taxon ID for a taxon, use the taxon ID lineage to find the most relevant taxon.
    # The most relevant taxon is the one with the most overlapping lineage with the taxonomic annotation.
    taxon_to_modify = {}

    for cluster in json_cluster_taxons:
        cluster_taxons = list(json_cluster_taxons[cluster].values())

        for index, taxon in enumerate(json_cluster_taxons[cluster]):
            taxon_ids = json_cluster_taxons[cluster][taxon]
            # If a taxon name is associated to more than one taxon ID, search for the one matching with the other taxon in the taxonomic annotation.
            if len(taxon_ids) > 1:
                taxon_shared_ids = {}
                for taxon_id in taxon_ids:
                    # Extract the other taxon present in the taxonomic annotation.
                    data_lineage = [tax_id for tax_id_lists in cluster_taxons[0:index] for tax_id in tax_id_lists]
                    # Extract the known lineage corresponding to one of the taxon ID.
                    lineage = ncbi.get_lineage(taxon_id)
                    # Computes the shared taxon ID between the known lineage and the taxonomic annotation associated to the taxon.
                    nb_shared_ids = len(set(data_lineage).intersection(set(lineage)))
                    taxon_shared_ids[taxon_id] = nb_shared_ids

                # If there is no match with the lineage for all the multiple taxon ID, it is not possible to decipher, stop esmecata.
                if all(shared_id==0 for shared_id in taxon_shared_ids.values()):
                    taxids = ','.join([str(taxid) for taxid in list(taxon_shared_ids.keys())])
                    print('It is not possible to find the taxon ID for the taxon named "{0}" (associated to "{1}") as there is multiple taxID possible ({2}), please add a more detailed taxonomic classification to help finding the correct one.'.format(taxon, cluster, taxids))
                    sys.exit()
                # If there is some matches, take the taxon ID with the most match and it will be the "correct one".
                else:
                    previous_nb_shared_ids = 0
                    for taxon_id in taxon_shared_ids:
                        if taxon_shared_ids[taxon_id] > previous_nb_shared_ids:
                            relevant_taxon = taxon_id
                            taxon_to_modify[cluster] = {}
                            taxon_to_modify[cluster][taxon] = [relevant_taxon]
                        previous_nb_shared_ids = nb_shared_ids

    for cluster in taxon_to_modify:
        for taxon in taxon_to_modify[cluster]:
            json_cluster_taxons[cluster][taxon] = taxon_to_modify[cluster][taxon]

    return json_cluster_taxons


def filter_rank_limit(json_cluster_taxons, ncbi, rank_limit):
    """
    Using teh rank_limit specificied, remove the taxon associated ot this rank.
    For example, if rank_limit == 'superkingdom', Bacteria will be removed.
    """
    # Rank level from Supplementary Table S3 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7408187/
    rank_level = {'superkingdom': 1, 'kingdom': 2, 'subkingdom': 3, 'superphylum': 4,
                    'phylum': 5, 'subphylum': 6, 'infraphylum': 7, 'superclass': 8,
                    'class': 9, 'subclass': 10, 'infraclass': 11, 'cohort': 12, 'subcohort': 13,
                    'superorder': 14, 'order': 15, 'suborder': 16, 'infraorder': 17, 'parvorder': 18,
                    'superfamily': 19, 'family': 20, 'subfamily': 21, 'tribe': 22, 'subtribe': 23,
                    'genus': 24, 'subgenus': 25, 'section': 26, 'subsection': 27, 'series': 28,
                    'subseries': 29, 'species group': 30, 'species subgroup': 31, 'species': 32,
                    'forma specialis': 33, 'subspecies': 34, 'varietas': 35, 'subvariety': 36,
                    'forma': 37, 'serogroup': 38, 'serotype': 39, 'strain': 40, 'isolate': 41}
    unclassified_rank = ['unclassified '+rank for rank in rank_level]
    non_hierarchical_ranks = ['clade', 'environmental samples', 'incertae sedis', 'no rank'] + unclassified_rank

    rank_limit_level = rank_level[rank_limit]
    rank_to_keeps = [rank for rank in rank_level if rank_level[rank] > rank_limit_level]

    for cluster in json_cluster_taxons:
        cluster_taxons = json_cluster_taxons[cluster]
        tax_keep = []
        tax_ranks = {}
        tax_names = {}
        tax_rank_position = {}
        tax_ids = []
        for index, tax_name in enumerate(cluster_taxons):
            tax_id = cluster_taxons[tax_name][0]
            tax_ids.append(tax_id)
            if tax_id != 'not_found':
                tax_rank = ncbi.get_rank([tax_id])[tax_id]
                tax_rank_position[tax_rank] = index
                if tax_rank in rank_to_keeps:
                    tax_ranks[tax_rank] = tax_id
                tax_names[tax_id] = tax_name
                # If the rank is below the rank to remove keep it.
                if tax_rank in rank_to_keeps:
                    tax_keep.append(True)
                elif tax_rank in non_hierarchical_ranks:
                    # If a nonhierarchical rank is below a rank that has been checked as being below the rank to remove, keep it.
                    if any(tax_keep):
                        tax_keep.append(True)
                    # If not remove it.
                    else:
                        tax_keep.append(False)
                else:
                    tax_keep.append(False)
            else:
                tax_keep.append(False)

        # If the rank to remove is in the tax_ranks, keep all the rank below this rank.
        if rank_limit in tax_ranks:
            tax_rank_max_index = tax_rank_position[rank_limit]
            keep_tax_ids = [tax_id[0] for tax_name, tax_id in list(cluster_taxons.items())[tax_rank_max_index+1:]]
        # If not find all the rank below this rank (by using rank_level and rank_to_keeps) in the list and keep them.
        # Use the tax_keep to keep non_hierarchical_ranks below the rank to remove.
        # But they need to be below a hierarchical rank that has been checked as being below the rank to remove.
        else:
            keep_tax_ids = [tax_ids[tax_index] for tax_index, bool_choice in enumerate(tax_keep) if bool_choice is True]

        tax_ids_wo_not_found = [tax_id for tax_id in tax_ids if tax_id != 'not_found']
        tax_id_to_deletes = set(tax_ids_wo_not_found) - set(keep_tax_ids)

        # Delete the remvoed rank and all its superior.
        for tax_id in tax_id_to_deletes:
            tax_name = tax_names[tax_id]
            del json_cluster_taxons[cluster][tax_name]

    return json_cluster_taxons


def rest_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, beta=None):
    proteomes = []
    proteomes_data = []
    organism_ids = {}
    # Find proteomes associated with taxon.
    # Take reference proteomes with "reference:yes".
    # Avoid redundant and excluded proteomes with "redundant%3Ano+excluded%3Ano".
    # Use "format=tab" to easily handle the ouput.
    http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+reference:yes+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)
    beta_httpt_str = 'https://rest.uniprot.org/beta/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A1)&format=json'.format(tax_id)

    if all_proteomes:
        http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)
        beta_httpt_str = 'https://rest.uniprot.org/beta/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A2)&format=json'.format(tax_id)

    # If esmecata does not find proteomes with only reference, search for all poroteomes even if they are not reference.
    all_http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)
    all_beta_httpt_str = 'https://rest.uniprot.org/beta/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A2)&format=json'.format(tax_id)

    response_proteome_status = False

    if not beta:
        with requests.get(http_str, REQUESTS_HEADERS) as proteome_response:
            # Raise error if we have a bad request.
            proteome_response.raise_for_status()
            proteome_response_text = proteome_response.text
            if proteome_response_text != '':
                csvreader = csv.reader(proteome_response_text.splitlines(), delimiter='\t')
                response_proteome_status = True
                # Avoid header.
                next(csvreader)
            else:
                csvreader = []
            reference_proteome = True

        if response_proteome_status is False:
            time.sleep(1)
            print('{0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(taxon, tax_id, tax_name))
            with requests.get(all_http_str, headers=REQUESTS_HEADERS) as proteome_response:
                proteome_response_text = proteome_response.text
                if proteome_response_text != '':
                    csvreader = csv.reader(proteome_response_text.splitlines(), delimiter='\t')
                    response_proteome_status = True
                    # Avoid header.
                    next(csvreader)
                else:
                    csvreader = []
            reference_proteome = False

        for line in csvreader:
            proteome = line[0]
            completness = line[6]
            org_tax_id = line[2]

            if line[4] != '':
                busco_percentage = float(line[4].split(':')[1].split('%')[0])
            else:
                busco_percentage = None

            # Check that proteome has busco score.
            if busco_percentage_keep:
                if busco_percentage and busco_percentage >= busco_percentage_keep and completness == 'full':
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

            proteomes_data.append([proteome, busco_percentage, completness, org_tax_id, reference_proteome])
    else:
        proteome_response = requests.get(url=beta_httpt_str, params=REQUESTS_HEADERS)
        proteome_response.raise_for_status()
        check_code = proteome_response.status_code
        data = proteome_response.json()
        reference_proteome = True
        if len(data['results']) == 0:
            print('{0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(taxon, tax_id, tax_name))
            time.sleep(1)
            proteome_response = requests.get(url=all_beta_httpt_str, params=REQUESTS_HEADERS)
            proteome_response.raise_for_status()
            data = proteome_response.json()
            reference_proteome = False

        for proteome_data in data['results']:
            proteome_id = proteome_data['id']
            if 'proteomeCompletenessReport' in proteome_data:
                if 'buscoReport' in proteome_data['proteomeCompletenessReport']:
                    proteome_busco = proteome_data['proteomeCompletenessReport']['buscoReport']
                else:
                    proteome_busco = None
                if proteome_busco:
                    busco_score = (proteome_busco['complete'] / proteome_busco['total']) * 100
                else:
                    busco_score = None
            else:
                busco_score = None

            if 'level' in proteome_data['genomeAssembly']:
                assembly_level = proteome_data['genomeAssembly']['level']
            else:
                assembly_level = ''
            org_tax_id = proteome_data['taxonomy']['taxonId']
            proteome_type = proteome_data['proteomeType']
            if busco_percentage_keep:
                if busco_score and busco_score >= busco_percentage_keep and assembly_level == 'full':
                    proteomes.append(proteome_id)
                    if org_tax_id not in organism_ids:
                        organism_ids[org_tax_id] = [proteome_id]
                    else:
                        organism_ids[org_tax_id].append(proteome_id)
            else:
                if assembly_level == 'full':
                    proteomes.append(proteome_id)
                    if org_tax_id not in organism_ids:
                        organism_ids[org_tax_id] = [proteome_id]
                    else:
                        organism_ids[org_tax_id].append(proteome_id)
            proteomes_data.append([proteome_id, busco_score, assembly_level, org_tax_id, reference_proteome])

    return proteomes, organism_ids, proteomes_data


def sparql_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql'):
    # SPARQL query to retrieve proteome
    # First FILTER NOT EXISTS to avoid redundant proteomes.
    # Second FILTER NOT EXISTS to avoid excluded proteomes.
    # OPTIONAL allows to retrive reference and non reference proteomes and split them.
    # From: uniprot graph location https://sparql.uniprot.org/.well-known/void%3E
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX busco: <http://busco.ezlab.org/schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

    SELECT DISTINCT ?proteome ?score ?fragmented ?missing ?organism ?completion ?type

    FROM <http://sparql.uniprot.org/proteomes>
    FROM <http://sparql.uniprot.org/taxonomy>

    WHERE
    {{
        ?proteome rdf:type up:Proteome .
        ?proteome up:organism ?organism .
        ?organism rdfs:subClassOf* taxon:{0} .
        OPTIONAL {{
            ?proteome busco:has_score ?busco .
            ?busco busco:complete ?score .
            ?busco busco:fragmented ?fragmented .
            ?busco busco:missing ?missing .
        }}
        OPTIONAL {{
            ?proteome rdfs:seeAlso ?metadataID .
            ?metadataID rdfs:label ?completion .
        }}
        FILTER NOT EXISTS {{
        ?proteome up:redundantTo ?proteome2 .
            }}
        FILTER NOT EXISTS {{
        ?proteome up:exclusionReason ?exclusion .
                }}
        OPTIONAL {{
            ?proteome rdf:type ?type .
        }}
    }}""".format(tax_id)

    csvreader = send_uniprot_sparql_query(uniprot_sparql_query, uniprot_sparql_endpoint)

    proteomes = []
    proteomes_data = []
    organism_ids = {}
    reference_proteomes = []
    other_proteomes = []
    for line in csvreader:
        proteome_id = line[0].split('/')[-1]
        if line[1] != '':
            score = float(line[1].split('^^')[0])
        else:
            score = None
        if line[2] != '':
            fragmented = float(line[2].split('^^')[0])
        else:
            fragmented = None
        if line[3] != '':
            missing = float(line[3].split('^^')[0])
        else:
            missing = None
        org_tax_id = line[4].split('/')[-1]
        completness = line[5]
        if 'Representative_Proteome' in line[6]:
            reference = True
        else:
            reference = False

        if all([score, fragmented, missing]) is True:
            busco_percentage = (score / (score+fragmented+missing)) * 100
        else:
            busco_percentage = None
        # Check that proteome has busco score.
        if busco_percentage_keep and busco_percentage:
            if busco_percentage >= busco_percentage_keep and completness == 'full':
                if reference == True:
                    reference_proteomes.append(proteome_id)
                else:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)
        else:
            if completness == 'full':
                if reference == True:
                    reference_proteomes.append(proteome_id)
                else:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)

        proteomes_data.append([proteome_id, busco_percentage, completness, org_tax_id, reference])
    if all_proteomes is not None:
        proteomes = other_proteomes
    else:
        if len(reference_proteomes) == 0:
            print('{0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(taxon, tax_id, tax_name))
            proteomes = other_proteomes
        else:
            proteomes = reference_proteomes

    return proteomes, organism_ids, proteomes_data


def find_proteomes_tax_ids(json_cluster_taxons, ncbi, proteomes_description_folder,
                        busco_percentage_keep=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        limit_maximal_number_proteomes=99, beta=None):
    # Query the Uniprot proteomes to find all the proteome IDs associated to taxonomic annotation.
    # If there is more than limit_maximal_number_proteomes proteomes a method is applied to extract a subset of the data.
    print('Find proteome ID associated to taxonomic annotation')
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
                    single_proteomes[taxon] = (tax_id, tax_id_founds[tax_id])
                break

            # If tax_id has not been found with a request do not try a new request with the same tax_id.
            if tax_id in tax_id_not_founds:
                continue

            # If tax_id not found because the tax_name has no tax_id associated avoid it.
            if tax_id == 'not_found':
                continue

            if uniprot_sparql_endpoint:
                proteomes, organism_ids, data_proteomes = sparql_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint)
            else:
                proteomes, organism_ids, data_proteomes = rest_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, beta)

            proteomes_description_file = os.path.join(proteomes_description_folder, taxon+'.tsv')
            if os.path.exists(proteomes_description_file):
                proteomes_description_file_writer = 'a'
            else:
                proteomes_description_file_writer = 'w'
            with open(proteomes_description_file, proteomes_description_file_writer) as proteome_output:
                csvwriter = csv.writer(proteome_output, delimiter='\t')
                if proteomes_description_file_writer == 'w':
                    csvwriter.writerow(['tax_id', 'tax_name', 'proteome_id', 'busco_percentage', 'completness', 'org_tax_id', 'reference_proteome'])
                for data_proteome in data_proteomes:
                    csvwriter.writerow([tax_id, tax_name, *data_proteome])

            # Answer is empty no corresponding proteomes to the tax_id.
            if len(proteomes) == 0:
                if tax_id not in tax_id_not_founds:
                    tax_id_not_founds[tax_id] = [tax_name]
                else:
                    tax_id_not_founds[tax_id].append(tax_name)
                continue

            elif len(proteomes) > 0 and len(proteomes) < 100:
                print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(proteomes)))
                proteomes_ids[taxon] = (tax_id, proteomes)
                tax_id_founds[tax_id] = proteomes
                if len(proteomes) == 1:
                    single_proteomes[taxon] = (tax_id, proteomes)
                break

            elif len(proteomes) > limit_maximal_number_proteomes:
                print('More than {0} proteomes are associated to the taxa {1} associated to {2}, esmecata will randomly select around {0} proteomes with respect to the taxonomic diversity.'.format(limit_maximal_number_proteomes, taxon, tax_name))
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
                        if parent_tax_id not in childs:
                            if parent_tax_id in organism_ids:
                                childs[parent_tax_id] = organism_ids[parent_tax_id]
                        else:
                            if parent_tax_id in organism_ids:
                                childs[parent_tax_id].extend(organism_ids[parent_tax_id])

                # For each direct descendant taxon, compute the number of proteomes in their childrens.
                elements = [v for k,v in childs.items()]
                elements_counts = [len(element) for element in elements]
                # Compute the proportion of proteomes for each direct descendant taxon compare to the total number of proteomes in all the descendant taxons.
                percentages = [(i/sum(elements_counts))*100  for i in elements_counts]
                # Can be superior to 100 if there is a lot of data with around 0.xxx percentage.
                percentages_round = [math.ceil(percentage) if percentage < 1 else math.floor(percentage) for percentage in percentages]
                selection_percentages_round = [math.ceil((limit_maximal_number_proteomes*percentage)/100) for percentage in percentages_round]

                # Choose randomly a number of proteomes corresponding to the computed percentage.
                selected_proteomes = []
                for index, element in enumerate(elements):
                    percentage_to_keep = selection_percentages_round[index]
                    if percentage_to_keep > len(element):
                        proteomes_to_keep = element
                    else:
                        proteomes_to_keep = random.sample(element, percentage_to_keep)
                    selected_proteomes.extend(proteomes_to_keep)
                proteomes_ids[taxon] = (tax_id, selected_proteomes)
                tax_id_founds[tax_id] = selected_proteomes
                print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(selected_proteomes)))
                break

            time.sleep(1)

    return proteomes_ids, single_proteomes, tax_id_not_founds

def sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint):
    # Implementation of the rdf:type up:Simple_Sequence to check for canonical sequence?
    # But is the rdf:type up:Simple_Sequence really associated to canonical sequence?
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX proteome: <http://purl.uniprot.org/proteomes/>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

    SELECT ?protein ?name ?isoform ?sequenceaa ?review

    FROM <http://sparql.uniprot.org/uniprot>
    FROM <http://sparql.uniprot.org/proteomes>

    WHERE
    {{
        ?protein a up:Protein ;
                up:mnemonic ?name ;
                up:reviewed ?review ;
                up:proteome ?genomicComponent .
        ?proteome skos:narrower ?genomicComponent .
        ?protein up:sequence ?isoform .
        ?isoform rdf:value ?sequenceaa .
        VALUES (?proteome) {{ (proteome:{0}) }}
    }}""".format(proteome)

    csvreader = send_uniprot_sparql_query(uniprot_sparql_query, uniprot_sparql_endpoint)

    records = []
    already_added_proteins = []
    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        protein_name = line[1]
        protein_isoform = line[2].split('/')[-1]
        protein_seq = line[3]
        prot_review = line[4].split('^^')[0]
        # Check protein review.
        # If review == Swissprot else Trembl.
        if prot_review == 'true':
            prefix_record = 'sp'
        elif prot_review == 'false':
            prefix_record = 'tr'
        else:
            print(protein_id, prot_review)

        # Take only the canonical sequences (ending with -1):
        if protein_isoform.endswith('-1'):
            record_id = f'{prefix_record}|{protein_id}|{protein_name}'
            records.append(SeqRecord(Seq(protein_seq), id=record_id, description=''))
            already_added_proteins.append(protein_id)

    intermediary_file = output_proteome_file[:-3]
    SeqIO.write(records, intermediary_file, 'fasta')

    with open(intermediary_file, 'rb') as input_file:
        with gzip.open(output_proteome_file, 'wb') as output_file:
            shutil.copyfileobj(input_file, output_file)
    os.remove(intermediary_file)


def retrieve_proteomes(input_file, output_folder, busco_percentage_keep=80,
                        ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        remove_tmp=None, limit_maximal_number_proteomes=99, beta=None, rank_limit=None):
    if is_valid_file(input_file) is False:
        print('The input {0} is not a valid file pathname.'.format(input_file))
        sys.exit()

    try:
        ete_taxadb_up_to_date = is_taxadb_up_to_date()
    except:
        ncbi = NCBITaxa()
        ete_taxadb_up_to_date = is_taxadb_up_to_date()

    if ete_taxadb_up_to_date is False:
        print('''WARNING: ncbi taxonomy database is not up to date with the last NCBI Taxonomy. Update it using:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()''')
        if ignore_taxadb_update is None:
            print('If you want to stil use esmecata with the old taxonomy database use the option --ignore-taxadb-update/ignore_taxadb_update.')
            return
        else:
            print('--ignore-taxadb-update/ignore_taxadb_update option detected, esmecata will continue with this version.')

    is_valid_dir(output_folder)

    proteomes_description_folder = os.path.join(output_folder, 'proteomes_description')
    is_valid_dir(proteomes_description_folder)

    if '.xlsx' in input_file:
        df = pd.read_excel(input_file)
    if '.tsv' in input_file:
        df = pd.read_csv(input_file, sep='\t')
    if '.csv' in input_file:
        df = pd.read_csv(input_file, sep=',')

    # Set index on the column containing the OTU name.
    df.set_index('observation_name', inplace=True)

    # taxonomic_annotation is the column containing the taxonomic annotation separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomic_annotation']

    proteome_cluster_tax_id_file = os.path.join(output_folder, 'proteome_cluster_tax_id.tsv')

    result_folder = os.path.join(output_folder, 'result')

    if not os.path.exists(proteome_cluster_tax_id_file):
        ncbi = NCBITaxa()

        tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(taxonomies, ncbi)

        json_cluster_taxons = filter_taxon(json_cluster_taxons, ncbi)

        if rank_limit:
            json_cluster_taxons = filter_rank_limit(json_cluster_taxons, ncbi, rank_limit)

        json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
        with open(json_log, 'w') as ouput_file:
            json.dump(json_cluster_taxons, ouput_file, indent=4)

    if not os.path.exists(proteome_cluster_tax_id_file):
        proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_cluster_taxons, ncbi, proteomes_description_folder,
                                                        busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint, limit_maximal_number_proteomes, beta)

        proteome_to_download = []
        for proteomes_id in proteomes_ids:
            proteome_to_download.extend(proteomes_ids[proteomes_id][1])
        proteome_to_download = set(proteome_to_download)

        # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
        with open(proteome_cluster_tax_id_file, 'w') as out_file:
            csvwriter = csv.writer(out_file, delimiter='\t')
            csvwriter.writerow(['cluster', 'name', 'tax_id', 'tax_rank', 'proteome'])
            for cluster in proteomes_ids:
                tax_id = int(proteomes_ids[cluster][0])
                tax_name = tax_id_names[tax_id]
                tax_rank = ncbi.get_rank([tax_id])[tax_id]
                csvwriter.writerow([cluster, tax_name, tax_id, tax_rank, ','.join(proteomes_ids[cluster][1])])
    else:
        proteome_to_download = []
        proteomes_ids = {}
        single_proteomes = {}
        with open(proteome_cluster_tax_id_file, 'r') as proteome_tax_file:
            csvreader = csv.reader(proteome_tax_file, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                observation_name = line[0]
                name = line[1]
                tax_id = line[2]
                proteomes = line[4].split(',')
                taxon_result_folder = os.path.join(result_folder, observation_name)
                if not os.path.exists(taxon_result_folder):
                    proteome_to_download.extend(proteomes)
                    if len(proteomes) == 1:
                        single_proteomes[observation_name] = (tax_id, proteomes)
                    proteomes_ids[observation_name] = (tax_id, proteomes)
                    print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(observation_name, name, len(proteomes)))

        proteome_to_download = set(proteome_to_download)

    # Download all the proteomes in tmp folder.
    print('Downloading {0} proteomes'.format(str(len(proteome_to_download))))
    tmp_folder = os.path.join(output_folder, 'tmp_proteome')
    is_valid_dir(tmp_folder)

    for proteome in proteome_to_download:
        output_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
        if not os.path.exists(output_proteome_file):
            if uniprot_sparql_endpoint is not None:
                sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint)
            else:
                if beta:
                    http_str = 'https://rest.uniprot.org/beta/uniprotkb/stream?query=proteome:{0}&format=fasta&compressed=true'.format(proteome)
                else:
                    http_str = 'https://www.uniprot.org/uniprot/?query=proteome:{0}&format=fasta&compress=yes'.format(proteome)
                with requests.get(http_str, headers=REQUESTS_HEADERS) as proteome_response:
                    with open(output_proteome_file, 'wb') as f:
                        f.write(proteome_response.content)
        time.sleep(1)

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_file': input_file, 'output_folder': output_folder, 'busco_percentage_keep': busco_percentage_keep,
                        'ignore_taxadb_update': ignore_taxadb_update, 'all_proteomes': all_proteomes, 'uniprot_sparql_endpoint': uniprot_sparql_endpoint,
                        'remove_tmp': remove_tmp, 'limit_maximal_number_proteomes': limit_maximal_number_proteomes}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['biopython'] = biopython_version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['ete3'] = ete3_version
    options['tool_dependencies']['python_package']['pandas'] = pd.__version__
    options['tool_dependencies']['python_package']['requests'] = requests.__version__
    options['tool_dependencies']['python_package']['SPARQLWrapper'] = sparqlwrapper_version

    if uniprot_sparql_endpoint:
        uniprot_releases = get_sparql_uniprot_release(uniprot_sparql_endpoint, options)
    else:
        uniprot_releases = get_rest_uniprot_release(options)

    uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_proteomes.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    print('Creating result folder')
    # Create a result folder which contains one sub-folder per OTU.
    # Each OTU sub-folder will contain the proteome found.
    is_valid_dir(result_folder)

    for cluster in proteomes_ids:
        output_cluster = os.path.join(result_folder, cluster)
        is_valid_dir(output_cluster)
        for proteome in proteomes_ids[cluster][1]:
            input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
            output_proteome = os.path.join(output_cluster, proteome+'.faa.gz')
            shutil.copyfile(input_proteome_file, output_proteome)

    if remove_tmp:
        shutil.rmtree(tmp_folder)
