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

logger = logging.getLogger(__name__)


def taxonomic_affiliation_to_taxon_id(observation_name, taxonomic_affiliation, ncbi=None):
    """ From a taxonomic affiliation (such as cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria) find corresponding taxon ID for each taxon.

    Args:
        observation_name (str): observation name associated with taxonomic affiliation
        taxonomic_affiliation (str): str with taxon from highest taxon (such as kingdom) to lowest (such as species)
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
    Returns:
        tax_ids_to_names (dict): mapping between taxon ID and taxon name
        taxon_ids (dict): mapping between taxon name and taxon ID
    """
    if ncbi is None:
        ncbi = NCBITaxa()

    taxons = [taxon for taxon in taxonomic_affiliation.split(';')]
    tax_names_to_ids = ncbi.get_name_translator(taxons)
    taxon_ids = OrderedDict()
    for taxon in taxons:
        taxon_translations = ncbi.get_name_translator([taxon])
        if taxon_translations == {}:
            logger.critical('|EsMeCaTa|proteomes| For {0}, no taxon ID has been found associated with the taxon "{1}" in the NCBI taxonomy of ete3.'.format(observation_name, taxon))
            taxon_translations = {taxon: ['not_found']}
        taxon_ids.update(taxon_translations)
    tax_ids_to_names = {v: k for k, vs in tax_names_to_ids.items() for v in vs }

    return tax_ids_to_names, taxon_ids


def associate_taxon_to_taxon_id(taxonomic_affiliations, ncbi=None):
    """ From a dictionary containing multiple taxonomic affiliations, find the taxon ID for each.

    Args:
        taxonomic_affiliations (dict): dictionary with observation name as key and taxonomic affiliation as value
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
    Returns:
        tax_id_names (dict): mapping between taxon ID and taxon name
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID
    """
    tax_id_names = {}
    json_taxonomic_affiliations = {}

    if taxonomic_affiliations == {}:
        logger.critical('|EsMeCaTa|proteomes| Empty taxonomy dictionary.')
        return

    # For each taxon find the taxon ID corresponding to each taxon.
    for observation_name in taxonomic_affiliations:
        taxonomic_affiliation = taxonomic_affiliations[observation_name]
        if isinstance(taxonomic_affiliation, str):
            tax_ids_to_names, taxon_ids = taxonomic_affiliation_to_taxon_id(observation_name, taxonomic_affiliation, ncbi)
            tax_id_names.update(tax_ids_to_names)
            json_taxonomic_affiliations[observation_name] = taxon_ids

    return tax_id_names, json_taxonomic_affiliations


def disambiguate_taxon(json_taxonomic_affiliations, ncbi):
    """ From json_taxonomic_affiliations (output of associate_taxon_to_taxon_id), disambiguate taxon name having multiples taxon IDs.
    For example, Yersinia is associated with both taxon ID 629 and 444888. By using the other taxon in the taxonomic affiliation, this function selects the correct taxon ID.

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
    Returns:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID (disambiguated)
    """
    # If there is multiple taxon ID for a taxon, use the taxon ID lineage to find the most relevant taxon.
    # The most relevant taxon is the one with the most overlapping lineage with the taxonomic affiliation.
    taxon_to_modify = {}

    for observation_name in json_taxonomic_affiliations:
        observation_name_taxons = list(json_taxonomic_affiliations[observation_name].values())

        for index, taxon in enumerate(json_taxonomic_affiliations[observation_name]):
            taxon_ids = json_taxonomic_affiliations[observation_name][taxon]
            # If a taxon name is associated with more than one taxon ID, search for the one matching with the other taxon in the taxonomic affiliation.
            if len(taxon_ids) > 1:
                taxon_shared_ids = {}
                for taxon_id in taxon_ids:
                    # Extract the other taxon present in the taxonomic affiliation.
                    data_lineage = [tax_id for tax_id_lists in observation_name_taxons[0:index] for tax_id in tax_id_lists]
                    # Extract the known lineage corresponding to one of the taxon ID.
                    lineage = ncbi.get_lineage(taxon_id)
                    # Computes the shared taxon ID between the known lineage and the taxonomic affiliation associated with the taxon.
                    nb_shared_ids = len(set(data_lineage).intersection(set(lineage)))
                    taxon_shared_ids[taxon_id] = nb_shared_ids

                # If there is no match with the lineage for all the multiple taxon ID, it is not possible to decipher, stop esmecata.
                if all(shared_id==0 for shared_id in taxon_shared_ids.values()):
                    taxids = ','.join([str(taxid) for taxid in list(taxon_shared_ids.keys())])
                    logger.critical('|EsMeCaTa|proteomes| It is not possible to find the taxon ID for the taxon named "{0}" (associated with "{1}") as there is multiple taxID possible ({2}), please add a more detailed taxonomic classification to help finding the correct one.'.format(taxon, observation_name, taxids))
                    sys.exit()
                # If there is some matches, take the taxon ID with the most match and it will be the "correct one".
                else:
                    previous_nb_shared_ids = 0
                    for taxon_id in taxon_shared_ids:
                        if taxon_shared_ids[taxon_id] > previous_nb_shared_ids:
                            relevant_taxon = taxon_id
                            taxon_to_modify[observation_name] = {}
                            taxon_to_modify[observation_name][taxon] = [relevant_taxon]
                        previous_nb_shared_ids = nb_shared_ids

    for observation_name in taxon_to_modify:
        for taxon in taxon_to_modify[observation_name]:
            json_taxonomic_affiliations[observation_name][taxon] = taxon_to_modify[observation_name][taxon]

    return json_taxonomic_affiliations


def filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit):
    """Using the rank_limit specificied, remove the taxon associated with this rank.
    For example, if rank_limit == 'superkingdom', Bacteria will be removed.

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
        rank_limit (str): rank limit to remove from the data
    Returns:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID (with remove rank specified)
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

    for observation_name in json_taxonomic_affiliations:
        observation_name_taxons = json_taxonomic_affiliations[observation_name]
        tax_keep = []
        tax_ranks = {}
        tax_names = {}
        tax_rank_position = {}
        tax_ids = []
        for index, tax_name in enumerate(observation_name_taxons):
            tax_id = observation_name_taxons[tax_name][0]
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
            keep_tax_ids = [tax_id[0] for tax_name, tax_id in list(observation_name_taxons.items())[tax_rank_max_index+1:]]
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
            del json_taxonomic_affiliations[observation_name][tax_name]

    return json_taxonomic_affiliations


def requests_query(http_str, nb_retry=5):
    """Use requests to query an http string.

    Args:
        http_str (str): http address to query
        nb_retry (int): number of retry to perform (default = 5)

    Returns:
        response (requests.models.Response): response returns by requests
    """
    passed = False

    if nb_retry == 0:
        sys.exit('5 retry attempts have been performed but were not successful, so esmecata has been stopped. You may try to relaunch it using the same command esmecata should resume.')
    try:
        response = requests.get(http_str, REQUESTS_HEADERS, timeout=60)
        passed = True
    except requests.exceptions.ConnectionError as error:
        logger.critical('|EsMeCaTa|proteomes| Error {0} occurs for query to "{1}", try to relaunch query.'.format(error, http_str))
        time.sleep(10)
        requests_query(http_str, nb_retry-1)
    except requests.exceptions.Timeout as error:
        logger.critical('|EsMeCaTa|proteomes| Error {0} occurs for query to "{1}", try to relaunch query.'.format(error, http_str))
        time.sleep(10)
        requests_query(http_str, nb_retry-1)

    if passed is True:
        return response


def rest_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes):
    """REST query on UniProt to get the proteomes associated with a taxon.

    Args:
        observation_name (str): observation name associated with the taxonomic affiliation
        tax_id (str): taxon ID from a taxon of the taxonomic affiliation
        tax_name (str): taxon name associated with the tax_id
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)

    Returns:
        proteomes (list): list of proteome IDs associated with the taxon ID
        organism_ids (dict): organism ID (key) associated with each proteomes (values)
        proteomes_data (list): list of lists with proteome_id, busco_score, assembly_level, org_tax_id, reference_proteome
    """
    proteomes = []
    proteomes_data = []
    organism_ids = {}
    # Find proteomes associated with taxon.
    # Take reference proteomes with "reference:yes".
    # Avoid redundant and excluded proteomes with "redundant%3Ano+excluded%3Ano".
    # Use "format=tab" to easily handle the ouput.
    httpt_str = 'https://rest.uniprot.org/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A1)&format=json'.format(tax_id)

    if all_proteomes:
        httpt_str = 'https://rest.uniprot.org/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A2)&format=json'.format(tax_id)

    # If esmecata does not find proteomes with only reference, search for all poroteomes even if they are not reference.
    all_http_str = 'https://rest.uniprot.org/proteomes/stream?query=(taxonomy_id%3A{0})AND(proteome_type%3A2)&format=json'.format(tax_id)

    response_proteome_status = False

    proteome_response = requests_query(httpt_str)
    proteome_response.raise_for_status()
    check_code = proteome_response.status_code
    data = proteome_response.json()
    reference_proteome = True
    if len(data['results']) == 0:
        logger.info('|EsMeCaTa|proteomes| {0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(observation_name, tax_id, tax_name))
        time.sleep(1)
        proteome_response = requests_query(all_http_str)
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


def sparql_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql'):
    """SPARQL query on UniProt to get the proteomes associated with a taxon.

    Args:
        observation_name (str): observation name associated with the taxonomic affiliation
        tax_id (str): taxon ID from a taxon of the taxonomic affiliation
        tax_name (str): taxon name associated with the tax_id
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
    Returns:
        proteomes (list): list of proteome IDs associated with the taxon ID
        organism_ids (dict): organism ID (key) associated with each proteomes (values)
        proteomes_data (list): list of lists with proteome_id, busco_score, assembly_level, org_tax_id, reference_proteome
    """
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
        ?organism rdfs:subClassOf taxon:{0} .
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

        if score is not None and fragmented is not None and missing is not None:
            busco_percentage = (score / (score+fragmented+missing)) * 100
        else:
            busco_percentage = None
        # Check that proteome has busco score.
        if busco_percentage_keep:
            if busco_percentage and busco_percentage >= busco_percentage_keep and completness == 'full':
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
            logger.info('|EsMeCaTa|proteomes| {0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(observation_name, tax_id, tax_name))
            proteomes = other_proteomes
        else:
            proteomes = reference_proteomes

    return proteomes, organism_ids, proteomes_data


def subsampling_proteomes(organism_ids, limit_maximal_number_proteomes, ncbi):
    """If the number of proteomes is superior to limit_maximal_number_proteomes, this funciton will perform a subsampling of the proteomes to select a number around the limit_maximal_number_proteomes.
    It will use the organism ID in organism_ids, to find the taxonomic diversity among the proteomes and keep this diversity.
    The idea is that if I have 10 proteomes from 3 organisms: org_1 with 4 proteomes, org_2 with 4 proteomes and org_3 with 2 proteomes, i would like to keep the same reprensetation with the subsampling.
    For example with a limit_maximal_number_proteomes at 5, subsampling will select 2 proteomes for org_1, 2 for org_2 and 1 for org_3

    Args:
        organism_ids (dict): organism ID as key and the proteome IDs associated with it as value
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        ncbi (ete3.NCBITaxa()): ete3 NCBI database

    Returns:
        selected_proteomes (list): subsample proteomes selected by the methods
    """
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
    # Compute the percentage of proteomes for each direct descendant taxon compare to the total number of proteomes in all the descendant taxons.
    percentages = [(i/sum(elements_counts))*100  for i in elements_counts]
    # Can be superior to the number limit_maximal_number_proteomes if there is a lot of data with around 0.xxx percentage.
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

    return selected_proteomes


def find_proteomes_tax_ids(json_taxonomic_affiliations, ncbi, proteomes_description_folder,
                        busco_percentage_keep=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        limit_maximal_number_proteomes=99, minimal_number_proteomes=1):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID (with remove rank specified)
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
        proteomes_description_folder (str): pathname to the proteomes_description_folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kepp

    Returns:
        proteomes_ids (dict): observation name (key) associated with proteome IDs
        single_proteomes (dict): observation name (key) associated with proteome ID (when there is only one proteome)
        tax_id_not_founds (dict): mapping between taxon ID and taxon name for taxon ID without proteomes
    """
    # Query the Uniprot proteomes to find all the proteome IDs associated with taxonomic affiliation.
    # If there is more than limit_maximal_number_proteomes proteomes a method is applied to extract a subset of the data.
    logger.info('|EsMeCaTa|proteomes| Find proteome ID associated with taxonomic affiliation')
    proteomes_ids = {}
    single_proteomes = {}
    tax_id_not_founds = {}
    tax_id_without_minimal_proteomes_number = {}
    tax_id_founds = {}

    for observation_name in json_taxonomic_affiliations:
        proteomes_descriptions = []
        for tax_name in reversed(json_taxonomic_affiliations[observation_name]):
            tax_id = json_taxonomic_affiliations[observation_name][tax_name][0]

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in tax_id_founds:
                proteomes_ids[observation_name] = (tax_id, tax_id_founds[tax_id])
                if len(tax_id_founds[tax_id]) == 1:
                    single_proteomes[observation_name] = (tax_id, tax_id_founds[tax_id])
                break

            # If tax_id has not been found with a request do not try a new request with the same tax_id.
            if tax_id in tax_id_not_founds or tax_id in tax_id_without_minimal_proteomes_number:
                continue

            # If tax_id not found because the tax_name has no tax_id associated avoid it.
            if tax_id == 'not_found':
                continue

            if uniprot_sparql_endpoint:
                proteomes, organism_ids, data_proteomes = sparql_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint)
            else:
                proteomes, organism_ids, data_proteomes = rest_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes)

            for data_proteome in data_proteomes:
                proteomes_descriptions.append([tax_id, tax_name, *data_proteome])

            # Answer is empty no corresponding proteomes to the tax_id.
            if len(proteomes) == 0:
                if tax_id not in tax_id_not_founds:
                    tax_id_not_founds[tax_id] = [tax_name]
                else:
                    tax_id_not_founds[tax_id].append(tax_name)
                continue

            # Check if the number of proteomes is inferior to the minimal_number_proteomes.
            elif len(proteomes) < minimal_number_proteomes:
                logger.info('|EsMeCaTa|proteomes| Less than {0} proteomes ({1}) are associated with the taxa {2} associated with {3}, esmecata will use a higher taxonomic rank to find proteomes.'.format(minimal_number_proteomes, len(proteomes), tax_name, observation_name))
                if tax_id not in tax_id_without_minimal_proteomes_number:
                    tax_id_without_minimal_proteomes_number[tax_id] = [tax_name]
                else:
                    tax_id_without_minimal_proteomes_number[tax_id].append(tax_name)
                continue

            elif len(proteomes) >= minimal_number_proteomes and len(proteomes) <= limit_maximal_number_proteomes:
                logger.info('|EsMeCaTa|proteomes| {0} will be associated with the taxon "{1}" with {2} proteomes.'.format(observation_name, tax_name, len(proteomes)))
                proteomes_ids[observation_name] = (tax_id, proteomes)
                tax_id_founds[tax_id] = proteomes
                if len(proteomes) == 1:
                    single_proteomes[observation_name] = (tax_id, proteomes)
                break

            elif len(proteomes) > limit_maximal_number_proteomes:
                logger.info('|EsMeCaTa|proteomes| More than {0} proteomes ({1}) are associated with the taxa {2} associated with {3}, esmecata will randomly select around {0} proteomes with respect to the taxonomic diversity.'.format(limit_maximal_number_proteomes, len(proteomes), tax_name, observation_name))
                selected_proteomes = subsampling_proteomes(organism_ids, limit_maximal_number_proteomes, ncbi)
                proteomes_ids[observation_name] = (tax_id, selected_proteomes)
                tax_id_founds[tax_id] = selected_proteomes
                logger.info('|EsMeCaTa|proteomes| {0} will be associated with the taxon "{1}" with {2} proteomes.'.format(observation_name, tax_name, len(selected_proteomes)))
                break

            time.sleep(1)

        proteomes_description_file = os.path.join(proteomes_description_folder, observation_name+'.tsv')

        with open(proteomes_description_file, 'w') as proteome_output:
            csvwriter = csv.writer(proteome_output, delimiter='\t')
            csvwriter.writerow(['tax_id', 'tax_name', 'proteome_id', 'busco_percentage', 'completness', 'org_tax_id', 'reference_proteome'])
            for proteomes_description in proteomes_descriptions:
                csvwriter.writerow(proteomes_description)

    return proteomes_ids, single_proteomes, tax_id_not_founds


def sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql'):
    """ SPARQL query to find the protein sequences associated with a proteome.

    Args:
        proteome (str): proteome ID
        output_proteome_file (str): pathname to the proteome fasta file
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
    """
    # Implementation of the rdf:type up:Simple_Sequence to check for canonical sequence?
    # But is the rdf:type up:Simple_Sequence really associated with canonical sequence?
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


def compute_stat_proteomes(result_folder, stat_file=None):
    """Compute stat associated with the number of proteome for each taxonomic affiliations.

    Args:
        result_folder (str): pathname to the result folder containing subfolder containing proteomes
        stat_file (str): pathname to the tsv stat file

    Returns:
        proteome_numbers (dict): dict containing observation names (as key) associated with the number of proteomes
    """
    proteome_numbers = {}
    for folder in os.listdir(result_folder):
        result_folder_path = os.path.join(result_folder, folder)
        number_proteome = len([proteome for proteome in os.listdir(result_folder_path)])
        proteome_numbers[folder] = number_proteome

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_proteomes'])
            for observation_name in proteome_numbers:
                csvwriter.writerow([observation_name, proteome_numbers[observation_name]])

    return proteome_numbers


def retrieve_proteomes(input_file, output_folder, busco_percentage_keep=80,
                        ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        remove_tmp=None, limit_maximal_number_proteomes=99, rank_limit=None,
                        minimal_number_proteomes=1):
    """From a tsv file with taxonomic affiliations find the associated proteomes.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations
        output_folder (str): pathname to the output folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        ignore_taxadb_update (bool): option to ignore ete3 taxa database update
        remove_tmp (bool): remove the tmp files
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        rank_limit (str): rank limit to remove from the data
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kepp
    """
    starttime = time.time()

    if is_valid_file(input_file) is False:
        logger.critical('|EsMeCaTa|proteomes| The input {0} is not a valid file pathname.'.format(input_file))
        sys.exit()

    try:
        ete_taxadb_up_to_date = is_taxadb_up_to_date()
    except:
        ncbi = NCBITaxa()
        ete_taxadb_up_to_date = is_taxadb_up_to_date()

    if ete_taxadb_up_to_date is False:
        logger.warning('''|EsMeCaTa|proteomes| WARNING: ncbi taxonomy database is not up to date with the last NCBI Taxonomy. Update it using:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()''')
        if ignore_taxadb_update is None:
            logger.info('|EsMeCaTa|proteomes| If you want to stil use esmecata with the old taxonomy database use the option --ignore-taxadb-update/ignore_taxadb_update.')
            return
        else:
            logger.info('|EsMeCaTa|proteomes| --ignore-taxadb-update/ignore_taxadb_update option detected, esmecata will continue with this version.')

    if minimal_number_proteomes >= limit_maximal_number_proteomes:
        logger.critical('|EsMeCaTa|proteomes| Error minimal_number_proteomes ({0}) can not be superior or equal to limit_maximal_number_proteomes ({1}).'.format(minimal_number_proteomes, limit_maximal_number_proteomes))
        sys.exit()

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

    # taxonomic_affiliation is the column containing the taxonomic affiliation separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomic_affiliation']

    proteome_tax_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    result_folder = os.path.join(output_folder, 'result')

    if not os.path.exists(proteome_tax_id_file):
        ncbi = NCBITaxa()

        tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, ncbi)

        json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

        if rank_limit:
            json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

        json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
        with open(json_log, 'w') as ouput_file:
            json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    if not os.path.exists(proteome_tax_id_file):
        proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_taxonomic_affiliations, ncbi, proteomes_description_folder,
                                                        busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint, limit_maximal_number_proteomes, minimal_number_proteomes)

        proteome_to_download = []
        for proteomes_id in proteomes_ids:
            proteome_to_download.extend(proteomes_ids[proteomes_id][1])
        proteome_to_download = set(proteome_to_download)

        # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
        with open(proteome_tax_id_file, 'w') as out_file:
            csvwriter = csv.writer(out_file, delimiter='\t')
            csvwriter.writerow(['observation_name', 'name', 'tax_id', 'tax_rank', 'proteome'])
            for observation_name in proteomes_ids:
                tax_id = int(proteomes_ids[observation_name][0])
                tax_name = tax_id_names[tax_id]
                tax_rank = ncbi.get_rank([tax_id])[tax_id]
                csvwriter.writerow([observation_name, tax_name, tax_id, tax_rank, ','.join(proteomes_ids[observation_name][1])])
    else:
        proteome_to_download = []
        proteomes_ids = {}
        single_proteomes = {}
        with open(proteome_tax_id_file, 'r') as proteome_tax_file:
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
                    logger.info('|EsMeCaTa|proteomes| {0} will be associated with the taxon "{1}" with {2} proteomes.'.format(observation_name, name, len(proteomes)))

        proteome_to_download = set(proteome_to_download)

    # Download all the proteomes in tmp folder.
    logger.info('|EsMeCaTa|proteomes| Downloading {0} proteomes'.format(str(len(proteome_to_download))))
    tmp_folder = os.path.join(output_folder, 'tmp_proteome')
    is_valid_dir(tmp_folder)

    for proteome in proteome_to_download:
        output_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
        if not os.path.exists(output_proteome_file):
            if uniprot_sparql_endpoint is not None:
                sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint)
            else:
                http_str = 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:{0}&format=fasta&compressed=true'.format(proteome)
                proteome_response = requests_query(http_str)
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

    # Create a result folder which contains one sub-folder per OTU.
    # Each OTU sub-folder will contain the proteome found.
    is_valid_dir(result_folder)

    for observation_name in proteomes_ids:
        output_observation_name = os.path.join(result_folder, observation_name)
        is_valid_dir(output_observation_name)
        for proteome in proteomes_ids[observation_name][1]:
            input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
            output_proteome = os.path.join(output_observation_name, proteome+'.faa.gz')
            shutil.copyfile(input_proteome_file, output_proteome)

    if remove_tmp:
        shutil.rmtree(tmp_folder)

    stat_file = os.path.join(output_folder, 'stat_number_proteome.tsv')
    compute_stat_proteomes(result_folder, stat_file)

    endtime = time.time()
    duration = endtime - starttime
    uniprot_releases['esmecata_proteomes_duration'] = duration
    json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
    uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_proteomes.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    logger.info('|EsMeCaTa|proteomes| Download complete.')
