# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
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

import copy
import csv
import gzip
import json
import logging
import math
import os
import pandas as pd
import random
import re
import requests
import shutil
import sys
import time
import seaborn as sns
import matplotlib.pyplot as plt

from collections import OrderedDict
from requests.adapters import HTTPAdapter, Retry

from Bio import __version__ as biopython_version
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ete3 import __version__ as ete3_version
from ete3 import NCBITaxa, is_taxadb_up_to_date

from SPARQLWrapper import __version__ as sparqlwrapper_version

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir, send_uniprot_sparql_query
from esmecata import __version__ as esmecata_version

from urllib.parse import unquote

REQUESTS_HEADERS = {'User-Agent': 'EsMeCaTa proteomes v' + esmecata_version + ', request by requests package v' + requests.__version__ }

logger = logging.getLogger(__name__)

# Rank level from Supplementary Table S3 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7408187/
RANK_LEVEL = OrderedDict({'superkingdom': 1, 'kingdom': 2, 'subkingdom': 3, 'superphylum': 4,
                'phylum': 5, 'subphylum': 6, 'infraphylum': 7, 'superclass': 8,
                'class': 9, 'subclass': 10, 'infraclass': 11, 'cohort': 12, 'subcohort': 13,
                'superorder': 14, 'order': 15, 'suborder': 16, 'infraorder': 17, 'parvorder': 18,
                'superfamily': 19, 'family': 20, 'subfamily': 21, 'tribe': 22, 'subtribe': 23,
                'genus': 24, 'subgenus': 25, 'section': 26, 'subsection': 27, 'series': 28,
                'subseries': 29, 'species group': 30, 'species subgroup': 31, 'species': 32,
                'forma specialis': 33, 'subspecies': 34, 'varietas': 35, 'subvariety': 36,
                'forma': 37, 'serogroup': 38, 'serotype': 39, 'strain': 40, 'isolate': 41})


def get_next_link(headers):
    """ From batch queries, get the next link to request.
    Function from: https://www.uniprot.org/help/id_mapping

    Args:
        headers (dict): headers of batch response.

    Returns:
        str: next link to query in the batch process.
    """
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(session, batch_url):
    """ Batch queries.
    Function from: https://www.uniprot.org/help/id_mapping

    Args:
        session (requests Session object): session used to query UniProt.
        batch_url (str): URL of batch queries.

    Returns:
        batch_response (requests Response object): response to batch query.
    """
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        yield response
        batch_url = get_next_link(response.headers)


def ete3_database_update(ignore_taxadb_update):
    """ Check if ete3 taxa database is up to date, if not update it.

    Args:
        ignore_taxadb_update (bool): if True, ignore the need to update database.
        esmecata_step (str): esmecata step name
    """
    try:
        ete_taxadb_up_to_date = is_taxadb_up_to_date()
    except:
        # If database is not set up, set it up.
        ncbi = NCBITaxa()
        ete_taxadb_up_to_date = is_taxadb_up_to_date()

    if ete_taxadb_up_to_date is False:
        logger.warning('|EsMeCaTa|proteomes| WARNING: ncbi taxonomy database is not up to date with the last NCBI Taxonomy. It will be updated')
        if ignore_taxadb_update is None:
            logger.info('|EsMeCaTa|proteomes| Ete3 taxa database will be updated.')
            ncbi = NCBITaxa()
            ncbi.update_taxonomy_database()
        else:
            logger.info('|EsMeCaTa|proteomes| --ignore-taxadb-update/ignore_taxadb_update option detected, esmecata will continue with this version.')


def update_taxonomy(observation_name, taxonomic_affiliation, ncbi=None):
    """ Update the taxonomic affiliation using ete3 and the lowest available taxonomic name.

    Args:
        observation_name (str): observation name associated with taxonomic affiliation
        taxonomic_affiliation (str): str with taxon from highest taxon (such as kingdom) to lowest (such as species)
        ncbi (ete3.NCBITaxa()): ete3 NCBI database

    Returns:
        new_taxonomic_affiliations (str): str with taxon from highest taxon (such as kingdom) to lowest (such as species)
    """
    if ncbi is None:
        ncbi = NCBITaxa()

    # For each taxon in the affiliation, search for the lineage associated in ete3.
    lineage_all_taxa = []
    taxons = [taxon for taxon in taxonomic_affiliation.split(';')]
    for taxon in reversed(taxons):
        # From taxon name to taxon ID.
        taxon_translations = ncbi.get_name_translator([taxon])

        if taxon in taxon_translations:
            if len(taxon_translations[taxon]) > 1:
                # If there are multiple taxon IDs for a taxon name.
                best_lineage_match = 0
                best_lineage_matches = []

                # Need to have other taxon in affiliations to find the correct one, if not skip this affiliation.
                if len(taxons) == 1:
                    logger.critical('|EsMeCaTa|proteomes| Taxonomy of %s has an ambiguous taxon name but not enoug taxa in its affiliations, skip updating it.', observation_name)
                    new_taxonomic_affiliations = None
                    break
                # Choose the ID matching the lineage of the other taxa in the affiliations.
                for tax_id in taxon_translations[taxon]:
                    tax_id_lineages = ncbi.get_lineage(tax_id)
                    tax_id_translator = ncbi.get_taxid_translator(tax_id_lineages)
                    tax_id_name = [tax_id_translator[tax_lineage] for tax_lineage in tax_id_lineages]
                    if len(set(tax_id_name).intersection(set(taxons))) > best_lineage_match:
                        best_lineage_match = len(set(tax_id_name).intersection(set(taxons)))
                        best_lineage_matches = tax_id_name
                lineage_all_taxa.append(best_lineage_matches)

            # If there is only one taxon ID, get its lineage.
            if len(taxon_translations[taxon]) == 1:
                tax_id = taxon_translations[taxon][0]
                tax_id_lineages = ncbi.get_lineage(tax_id)
                tax_id_translator = ncbi.get_taxid_translator(tax_id_lineages)
                tax_id_name = [tax_id_translator[tax_lineage] for tax_lineage in tax_id_lineages]
                lineage_all_taxa.append(tax_id_name)

            # If no taxon found, no possibility to update
            else:
                new_taxonomic_affiliations = None
        else:
            new_taxonomic_affiliations = None

    if lineage_all_taxa != []:
        # If it founds new affiliations, check that all lineages found are similar.
        check_lineage = []
        for index, lineage_taxon in enumerate(lineage_all_taxa):
            tmp_check_lineage = []
            for index_to_compare in range(index, len(lineage_all_taxa)):
                # Look at the intersection between the current lineage and all the ones found for the other taxon in affiliations.
                intersection_lineage = set(lineage_taxon).intersection(set(lineage_all_taxa[index_to_compare]))
                # If this intersection equals to the lineage we use to compare?
                intersection_lineage_bool = intersection_lineage == set(lineage_all_taxa[index_to_compare])
                tmp_check_lineage.append(intersection_lineage_bool)
            check_lineage.append(tmp_check_lineage)

        for index, check_values in enumerate(check_lineage):
            # If all comparison corresponds, keep the lowest affiliation.
            if all(check_values) is True:
                new_taxonomic_affiliations = ';'.join(lineage_all_taxa[index])
                logger.critical('|EsMeCaTa|proteomes| Taxonomy of %s ("%s") updated into "%s" .', observation_name, taxonomic_affiliation, new_taxonomic_affiliations)
                break

    # If no new taxonomic affiliations, keep the old one. 
    if new_taxonomic_affiliations is None:
        new_taxonomic_affiliations = taxonomic_affiliation

    return new_taxonomic_affiliations


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
            logger.critical('|EsMeCaTa|proteomes| For %s, no taxon ID has been found associated with the taxon "%s" in the NCBI taxonomy of ete3.', observation_name, taxon)
            taxon_translations = {taxon: ['not_found']}
        taxon_ids.update(taxon_translations)
    tax_ids_to_names = {v: k for k, vs in tax_names_to_ids.items() for v in vs }

    return tax_ids_to_names, taxon_ids


def associate_taxon_to_taxon_id(taxonomic_affiliations, update_affiliations=None, ncbi=None, output_folder=None):
    """ From a dictionary containing multiple taxonomic affiliations, find the taxon ID for each.

    Args:
        taxonomic_affiliations (dict): dictionary with observation name as key and taxonomic affiliation as value
        update_affiliations (str): option to update taxonomic affiliations.
        ncbi (ete3.NCBITaxa()): ete3 NCBI database.
        output_folder (str): pathname to output folder.

    Returns:
        tax_id_names (dict): mapping between taxon ID and taxon name
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID
    """
    tax_id_names = {}
    json_taxonomic_affiliations = {}

    if taxonomic_affiliations == {}:
        logger.critical('|EsMeCaTa|proteomes| Empty taxonomy dictionary.')
        return

    new_taxonomic_affiliations = []
    # For each taxon find the taxon ID corresponding to each taxon.
    for observation_name in taxonomic_affiliations:
        taxonomic_affiliation = taxonomic_affiliations[observation_name]
        if isinstance(taxonomic_affiliation, str):
            if update_affiliations is not None:
                taxonomic_affiliation = update_taxonomy(observation_name, taxonomic_affiliation, ncbi)
                new_taxonomic_affiliations.append([observation_name, taxonomic_affiliation])
            tax_ids_to_names, taxon_ids = taxonomic_affiliation_to_taxon_id(observation_name, taxonomic_affiliation, ncbi)
            tax_id_names.update(tax_ids_to_names)
            json_taxonomic_affiliations[observation_name] = taxon_ids

    if update_affiliations:
        # Write new taxonomic file if update-affiliations.
        new_taxon_file = os.path.join(output_folder, 'new_taxonomic_affiliation.tsv')
        with open(new_taxon_file, 'w') as output_file:
            csvwriter = csv.writer(output_file, delimiter='\t')
            csvwriter.writerow(['observation_name', 'taxonomic_affiliation'])
            for new_taxonomic_affiliation in new_taxonomic_affiliations:
                csvwriter.writerow(new_taxonomic_affiliation)

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
                    logger.critical('|EsMeCaTa|proteomes| It is not possible to find the taxon ID for the taxon named "%s" (associated with "%s") as there are multiple taxIDs possible (%s), please add a more detailed taxonomic classification to help finding the correct one.', taxon, observation_name, taxids)
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
    """Using the rank_limit specificied, remove the taxon superior to this rank.
    For example, if rank_limit == 'family', taxa associated to rank superior to family will be removed.

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks)
    Returns:
        output_json_taxonomic_affiliations (dict): observation name and dictionary with mapping betwenn taxon name and taxon ID (with remove rank specified)
    """
    unclassified_rank = ['unclassified '+rank for rank in RANK_LEVEL]
    non_hierarchical_ranks = ['clade', 'environmental samples', 'incertae sedis', 'no rank'] + unclassified_rank

    rank_limit_level = RANK_LEVEL[rank_limit]
    rank_to_keeps = [rank for rank in RANK_LEVEL if RANK_LEVEL[rank] >= rank_limit_level]

    # Create a deep copy of the input dictionary, as we will modified this dictionary.
    output_json_taxonomic_affiliations = copy.deepcopy(json_taxonomic_affiliations)
    for observation_name in output_json_taxonomic_affiliations:
        observation_name_taxons = output_json_taxonomic_affiliations[observation_name]
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

        # If the rank limit is in the tax_ranks, keep all the rank below this rank and this rank.
        if rank_limit in tax_ranks:
            tax_rank_max_index = tax_rank_position[rank_limit]
            keep_tax_ids = [tax_id[0] for _, tax_id in list(observation_name_taxons.items())[tax_rank_max_index:]]
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
            del output_json_taxonomic_affiliations[observation_name][tax_name]

    return output_json_taxonomic_affiliations


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
        logger.critical('|EsMeCaTa|proteomes| Error %s occurs for query to "%s", try to relaunch query.', error, http_str)
        time.sleep(10)
        requests_query(http_str, nb_retry-1)
    except requests.exceptions.Timeout as error:
        logger.critical('|EsMeCaTa|proteomes| Error %s occurs for query to "%s", try to relaunch query.', error, http_str)
        time.sleep(10)
        requests_query(http_str, nb_retry-1)

    if passed is True:
        return response


def rest_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep,
                         all_proteomes, session=None, option_bioservices=None,
                         minimal_number_proteomes=1):
    """REST query on UniProt to get the proteomes associated with a taxon.

    Args:
        observation_name (str): observation name associated with the taxonomic affiliation
        tax_id (str): taxon ID from a taxon of the taxonomic affiliation
        tax_name (str): taxon name associated with the tax_id
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        session: request session object
        option_bioservices (bool): use bioservices instead of manual queries.
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kept.

    Returns:
        proteomes (list): list of proteome IDs associated with the taxon ID
        organism_ids (dict): organism ID (key) associated with each proteomes (values)
        proteomes_data (list): list of lists with proteome_id, busco_score, assembly_level, org_tax_id, reference_proteome, component_elements
    """
    if not session:
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))

    if option_bioservices is None:
        # Find proteomes associated with taxon.
        # Search for both representative and non-representative proteomes (with proteome_type%3A2 for non-representative proteome (or Other proteome) and proteome_type%3A1) for representative proteome).
        # Remove redundant and excluded proteomes.
        httpt_str = 'https://rest.uniprot.org/proteomes/stream?query=(taxonomy_id%3A{0})AND((proteome_type%3A2)OR(proteome_type%3A1))&format=json&size=500'.format(tax_id)

        data = {}
        data['results'] = []
        for batch_reponse in get_batch(session, httpt_str):
            batch_json = batch_reponse.json()
            data['results'].extend(batch_json['results'])
    else:
        import bioservices
        uniprot_bioservices = bioservices.UniProt()
        data = uniprot_bioservices.search(f'(taxonomy_id={tax_id})AND((proteome_type=2)OR(proteome_type=1))',
                                            database='proteomes', frmt='json', progress=False)

    organism_ids = {}
    proteomes_data = []
    representative_proteomes = []
    other_proteomes = []

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
        org_tax_id = str(proteome_data['taxonomy']['taxonId'])

        # Check proteome type for representative or non-representative.
        proteome_type = proteome_data['proteomeType']

        if proteome_type in ['Representative proteome', 'Reference and representative proteome', 'Reference proteome']:
            representative_proteome = True
        elif proteome_type == 'Other proteome':
            representative_proteome = False

        if 'components' in proteome_data:
            component_elements = []
            for component in  proteome_data['components']:
                component_name = component['name']
                component_description = component['description']
                component_elements.append([component_name, component_description])
        else:
            component_elements = None

        if busco_percentage_keep:
            if busco_score and busco_score >= busco_percentage_keep and assembly_level == 'full':
                if representative_proteome is True:
                    representative_proteomes.append(proteome_id)
                elif representative_proteome is False:
                    other_proteomes.append(proteome_id)

                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)
        else:
            if assembly_level == 'full':
                if representative_proteome is True:
                    representative_proteomes.append(proteome_id)
                elif representative_proteome is False:
                    other_proteomes.append(proteome_id)

                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)
        proteomes_data.append([proteome_id, busco_score, assembly_level, org_tax_id, representative_proteome, component_elements])

    # In REST queries, reference proteomes are not associated with non-reference proteome tag (compared to SPARQL query).
    # So to have both of them with all_proteomes otpion, we need to add other_proteomes to representative_proteomes.
    if all_proteomes is not None:
        proteomes = other_proteomes + representative_proteomes
    else:
        if len(representative_proteomes) == 0:
            logger.info('|EsMeCaTa|proteomes| %s: No reference proteomes found for %s (%s) try non-reference proteomes.', observation_name, tax_id, tax_name)
            proteomes = other_proteomes
        elif len(representative_proteomes) < minimal_number_proteomes:
            logger.info('|EsMeCaTa|proteomes| %s: No reference proteomes found for %s (%s) try non-reference proteomes.', observation_name, tax_id, tax_name)
            proteomes = other_proteomes + representative_proteomes
        else:
            proteomes = representative_proteomes

    return proteomes, organism_ids, proteomes_data


def sparql_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep,
                           all_proteomes, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql', minimal_number_proteomes=1):
    """SPARQL query on UniProt to get the proteomes associated with a taxon.

    Args:
        observation_name (str): observation name associated with the taxonomic affiliation
        tax_id (str): taxon ID from a taxon of the taxonomic affiliation
        tax_name (str): taxon name associated with the tax_id
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kept.

    Returns:
        proteomes (list): list of proteome IDs associated with the taxon ID
        organism_ids (dict): organism ID (key) associated with each proteomes (values)
        proteomes_data (list): list of lists with proteome_id, busco_score, assembly_level, org_tax_id, reference_proteome, component_elements
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

    SELECT DISTINCT ?proteome ?score ?fragmented ?missing ?organism ?completion ?type ?component ?componentName

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
        OPTIONAL {{
            ?proteome skos:narrower ?component .
            ?component rdfs:comment ?componentName .
        }}
    }}""".format(tax_id)

    csvreader = send_uniprot_sparql_query(uniprot_sparql_query, uniprot_sparql_endpoint)

    proteomes = []
    proteomes_data = []
    organism_ids = {}
    reference_proteomes = []
    other_proteomes = []

    sparql_proteome_data = {}
    # There can be multiple lines for a same proteomes (especially with reference proteomes being also associated with reference and non-reference proteomes).
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
        if 'Representative_Proteome' in line[6] or 'Reference_Proteome' in line[6]:
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
                if reference is True:
                    reference_proteomes.append(proteome_id)
                elif reference is False:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)
        else:
            if completness == 'full':
                if reference == True:
                    reference_proteomes.append(proteome_id)
                elif reference is False:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)

        # Get the components contained in the proteomes (chromosomes, plasmids).
        if line[7] != '':
            component_name = unquote(line[7]).split('#')[1]
        else:
            component_name = None
        if line[8] != '':
            component_description = line[8]
        else:
            component_description = None
        if component_description is not None or component_name is not None:
            component_elements = [component_name, component_description]

        if proteome_id not in sparql_proteome_data:
            sparql_proteome_data[proteome_id] = [proteome_id, busco_percentage, completness, org_tax_id, reference, [component_elements]]
        else:
            if sparql_proteome_data[proteome_id][4] is False and reference is True:
                sparql_proteome_data[proteome_id][4] = True
            sparql_proteome_data[proteome_id][5].append(component_elements)

    for protein_id in sparql_proteome_data:
        proteomes_data.append(sparql_proteome_data[protein_id])

    # In SPARQL reference proteomes are also labelled as non-reference proteome (with 'http://purl.uniprot.org/core/Proteome').
    # To use both reference and non-reference proteomes, use only other_proteomes.
    other_proteomes = set(other_proteomes)
    reference_proteomes = set(reference_proteomes)

    if all_proteomes is not None:
        proteomes = other_proteomes
    else:
        if len(reference_proteomes) == 0:
            logger.info('|EsMeCaTa|proteomes| %s: No reference proteomes found for %s (%s) try non-reference proteomes.', observation_name, tax_id, tax_name)
            proteomes = other_proteomes
        elif len(reference_proteomes) < minimal_number_proteomes:
            logger.info('|EsMeCaTa|proteomes| %s: No reference proteomes found for %s (%s) try non-reference proteomes.', observation_name, tax_id, tax_name)
            proteomes = other_proteomes
        else:
            proteomes = reference_proteomes

    return proteomes, organism_ids, proteomes_data


def compare_input_taxon_with_esmecata(json_taxonomic_affiliations, proteomes_ids, ncbi):
    """Search the taxonomic affiliations found by ete3 to find the lowest tax rank associated with each observation_name.
    Then look at esmecata to find the taxon rank used to find proteome.

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        proteomes_ids (dict):  observation name (key) associated with proteome IDs
        ncbi (ete3.NCBITaxa()): ete3 NCBI database

    Returns:
        taxon_rank_comparison (dict): lowest taxon name in input affiliation (key) associated with dict containing associated taxon ranks found by esmecata and the corresponding occurrence
    """
    # Create taxon comparison.
    taxon_rank_comparison = {}
    esmecata_found_ranks = []
    for observation_name in json_taxonomic_affiliations:
        reversed_affiliation_taxa = list(reversed(json_taxonomic_affiliations[observation_name]))
        lowest_tax_rank = None
        for tax_name in reversed_affiliation_taxa:
            if json_taxonomic_affiliations[observation_name][tax_name] != ['not_found']:
                if tax_name != 'unknown':
                    lowest_tax_id = json_taxonomic_affiliations[observation_name][tax_name][0]
                    lowest_tax_rank = ncbi.get_rank([lowest_tax_id])[lowest_tax_id]

                    # If lowest_taxon_rank is not in RANK_LEVEL, classified it with the upper tax rank.
                    # This is the case for the unclassified taxa.
                    if lowest_tax_rank not in RANK_LEVEL:
                        higher_tax_id = ncbi.get_lineage(lowest_tax_id)[-2]
                        higher_tax_rank = ncbi.get_rank([higher_tax_id])[higher_tax_id]
                        lowest_tax_rank = higher_tax_rank
                    break

        # If no taxon identified by ete3 is found, classified it as not_found.
        if lowest_tax_rank is None:
            lowest_tax_rank = 'not_found'

        if observation_name in proteomes_ids:
            esmecata_tax_id = int(proteomes_ids[observation_name][0])
            esmecata_tax_rank = ncbi.get_rank([esmecata_tax_id])[esmecata_tax_id]
            if esmecata_tax_rank not in RANK_LEVEL:
                esmecata_higher_tax_id = ncbi.get_lineage(esmecata_tax_id)[-2]
                esmecata_higher_tax_rank = ncbi.get_rank([esmecata_higher_tax_id])[esmecata_higher_tax_id]
                esmecata_tax_rank = esmecata_higher_tax_rank
        else:
            esmecata_tax_rank = 'not_found'

        esmecata_found_ranks.append(esmecata_tax_rank)
        if lowest_tax_rank not in taxon_rank_comparison:
            taxon_rank_comparison[lowest_tax_rank] = {}
            taxon_rank_comparison[lowest_tax_rank][esmecata_tax_rank] = 1
        else:
            if esmecata_tax_rank not in taxon_rank_comparison[lowest_tax_rank]:
                taxon_rank_comparison[lowest_tax_rank][esmecata_tax_rank] = 1
            else:
                taxon_rank_comparison[lowest_tax_rank][esmecata_tax_rank] += 1

    for esmecata_taxon_rank in esmecata_found_ranks:
        if esmecata_taxon_rank not in taxon_rank_comparison:
            taxon_rank_comparison[esmecata_taxon_rank] = {}

    return taxon_rank_comparison


def create_taxon_heatmap(taxon_rank_comparison, output_folder):
    """Create a heatmap showing the lowest taxonomic rank in the input compared to the ones used by esmecata.

    Args:
        taxon_rank_comparison (dict): lowest taxon name in input affiliation (key) associated with dict containing associated taxon ranks found by esmecata and the corresponding occurrence
        output_folder (str): path to the output folder.
    """
    taxon_comparison_dataframe = pd.DataFrame.from_dict(taxon_rank_comparison)

    # Reorder taxonomic ranks in dataframe.
    ordered_ranks = []
    # Add rank not found.
    RANK_LEVEL.update({'not_found':'not_found'})
    RANK_LEVEL.move_to_end('not_found', last=False)
    for tax_rank in RANK_LEVEL:
        if tax_rank in taxon_comparison_dataframe.columns:
            ordered_ranks.append(tax_rank)
    taxon_comparison_dataframe = taxon_comparison_dataframe.reindex(ordered_ranks)
    ordered_ranks = reversed(ordered_ranks)
    taxon_comparison_dataframe = taxon_comparison_dataframe[ordered_ranks]

    # Sum all the apparition
    #import numpy as np
    #Flip matrix
    #matrix_sup_to_diagonal = np.triu(np.fliplr(df), 1)
    # Keep only number above diagonal showing differences betwen input taxon and esmecata taxon.
    #am = np.fliplr( matrix_sup_to_diagonal)[:-1]

    taxon_comparison_dataframe.to_csv(os.path.join(output_folder, 'tax_comparison_rank.tsv'), sep='\t', index=None)
    sns.set('poster', rc={'figure.figsize':(30,20), 'lines.linewidth': 10})
    sns.set_style("white")
    # Create output figure.
    tax_comparison_rank_fig = os.path.join(output_folder, 'tax_comparison_rank.svg')
    sns.heatmap(taxon_comparison_dataframe, annot=True, fmt='.0f')
    plt.xlabel('Lowest taxonomic rank')
    plt.ylabel('Taxonomic rank used by EsMeCaTa')
    plt.savefig(tax_comparison_rank_fig)


def create_taxon_heatmap_from_complete_run(esmecata_proteome_folder):
    """ Create tax_comparison_rank picture file from esmecata proteome folder after complete run.

    Args:
        esmecata_proteome_folder (dict): path to esmecata proteome folder with all data
    """
    proteome_tax_id_file = os.path.join(esmecata_proteome_folder, 'proteome_tax_id.tsv')
    proteomes_ids = {}
    with open(proteome_tax_id_file, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            tax_id = line['tax_id']
            proteomes = line['proteome'].split(',')
            proteomes_ids[observation_name] = (tax_id, proteomes)

    json_log = os.path.join(esmecata_proteome_folder, 'association_taxon_taxID.json')
    with open(json_log, 'r') as input_json_file:
        json_taxonomic_affiliations = json.load(input_json_file)

    ncbi = NCBITaxa()
    # Create heatmap comparing input taxon and taxon used by esmecata to find proteomes.
    taxon_rank_comparison = compare_input_taxon_with_esmecata(json_taxonomic_affiliations, proteomes_ids, ncbi)
    create_taxon_heatmap(taxon_rank_comparison, esmecata_proteome_folder)


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
    try:
        tree = ncbi.get_topology([org_tax_id for org_tax_id in organism_ids])
    except KeyError:
        logger.critical('|EsMeCaTa|proteomes| UniProt tax ID not in ete3 NCBI Taxonomy database. Try to update it with the following command: python3 -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(); ncbi.update_taxonomy_database()".')
        raise KeyError()

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
                        limit_maximal_number_proteomes=99, minimal_number_proteomes=1, session=None,
                        option_bioservices=None):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        ncbi (ete3.NCBITaxa()): ete3 NCBI database
        proteomes_description_folder (str): pathname to the proteomes_description_folder
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold)
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes)
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint)
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taoxn to be kept
        session: request session object
        option_bioservices (bool): use bioservices instead of manual queries.

    Returns:
        proteomes_ids (dict): observation name (key) associated with proteome IDs
        single_proteomes (dict): observation name (key) associated with proteome ID (when there is only one proteome)
        tax_id_not_founds (dict): mapping between taxon ID and taxon name for taxon ID without proteomes
    """
    # Query the Uniprot proteomes to find all the proteome IDs associated with taxonomic affiliation.
    # If there is more than limit_maximal_number_proteomes proteomes a method is applied to extract a subset of the data.
    logger.info('|EsMeCaTa|proteomes| Find proteome ID associated with taxonomic affiliation')
    proteomes_ids = {}
    proteome_data = {}
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
                proteomes_descriptions.append(proteome_data[tax_id])
                if len(tax_id_founds[tax_id]) == 1:
                    single_proteomes[observation_name] = (tax_id, tax_id_founds[tax_id])
                logger.info('|EsMeCaTa|proteomes| "%s" already associated with proteomes, %s will be associated with the taxon "%s" with %d proteomes.', tax_name, observation_name, tax_name, len(tax_id_founds[tax_id]))
                break

            # If tax_id has not been found with a request do not try a new request with the same tax_id.
            if tax_id in tax_id_not_founds or tax_id in tax_id_without_minimal_proteomes_number:
                continue

            # If tax_id not found because the tax_name has no tax_id associated avoid it.
            if tax_id == 'not_found':
                continue

            if uniprot_sparql_endpoint:
                proteomes, organism_ids, data_proteomes = sparql_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint, minimal_number_proteomes)
            else:
                proteomes, organism_ids, data_proteomes = rest_query_proteomes(observation_name, tax_id, tax_name, busco_percentage_keep, all_proteomes, session, option_bioservices, minimal_number_proteomes)

            for data_proteome in data_proteomes:
                proteomes_descriptions.append([tax_id, tax_name, *data_proteome])
            proteome_data[tax_id] = proteomes_descriptions

            # Answer is empty no corresponding proteomes to the tax_id.
            if len(proteomes) == 0:
                if tax_id not in tax_id_not_founds:
                    tax_id_not_founds[tax_id] = [tax_name]
                else:
                    tax_id_not_founds[tax_id].append(tax_name)
                continue

            # Check if the number of proteomes is inferior to the minimal_number_proteomes.
            elif len(proteomes) < minimal_number_proteomes:
                logger.info('|EsMeCaTa|proteomes| Less than %d proteomes (%d) are associated with the taxa %s associated with %s, esmecata will use a higher taxonomic rank to find proteomes.', minimal_number_proteomes, len(proteomes), tax_name, observation_name)
                if tax_id not in tax_id_without_minimal_proteomes_number:
                    tax_id_without_minimal_proteomes_number[tax_id] = [tax_name]
                else:
                    tax_id_without_minimal_proteomes_number[tax_id].append(tax_name)
                continue

            elif len(proteomes) >= minimal_number_proteomes and len(proteomes) <= limit_maximal_number_proteomes:
                logger.info('|EsMeCaTa|proteomes| %s will be associated with the taxon "%s" with %d proteomes.', observation_name, tax_name, len(proteomes))
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
                logger.info('|EsMeCaTa|proteomes| %s will be associated with the taxon "%s" with %d proteomes.', observation_name, tax_name, len(selected_proteomes))
                break

            time.sleep(1)

        proteomes_description_file = os.path.join(proteomes_description_folder, observation_name+'.tsv')

        with open(proteomes_description_file, 'w') as proteome_output:
            csvwriter = csv.writer(proteome_output, delimiter='\t')
            csvwriter.writerow(['tax_id', 'tax_name', 'proteome_id', 'busco_percentage', 'completness', 'org_tax_id', 'reference_proteome', 'components'])
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
        if prot_review == 'true' or prot_review == '1':
            prefix_record = 'sp'
        elif prot_review == 'false' or prot_review == '0':
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


def compute_stat_proteomes(proteomes_folder, stat_file=None):
    """Compute stat associated with the number of proteome for each taxonomic affiliations.

    Args:
        proteome_tax_id_file (str): pathname to the proteome_tax_id file indicating the number of proteomes
        stat_file (str): pathname to the tsv stat file

    Returns:
        proteome_numbers (dict): dict containing observation names (as key) associated with the number of proteomes
    """
    ncbi = NCBITaxa()
    proteome_numbers = {}

    json_log = os.path.join(proteomes_folder, 'association_taxon_taxID.json')
    with open(json_log, 'r') as input_json_file:
        json_taxonomic_affiliations = json.load(input_json_file)

    proteome_status = {}
    proteome_description_folder = os.path.join(proteomes_folder, 'proteomes_description')
    for proteome_description_file in os.listdir(proteome_description_folder):
        with open(os.path.join(proteome_description_folder, proteome_description_file), 'r') as open_proteome_description_file:
             csvreader = csv.DictReader(open_proteome_description_file, delimiter='\t')
             for line in csvreader:
                 if line['reference_proteome'] == 'False':
                     reference_proteome = False
                 else:
                     reference_proteome = True
                 proteome_status[line['proteome_id']] = reference_proteome

    proteome_tax_id_file = os.path.join(proteomes_folder, 'proteome_tax_id.tsv')
    with open(proteome_tax_id_file, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            reversed_affiliation_taxa = list(reversed(list(json_taxonomic_affiliations[observation_name].keys())))
            lowest_tax_rank = None
            for tax_name in reversed_affiliation_taxa:
                if json_taxonomic_affiliations[observation_name][tax_name] != ['not_found']:
                    if tax_name != 'unknown':
                        lowest_tax_id = json_taxonomic_affiliations[observation_name][tax_name][0]
                        lowest_tax_rank = ncbi.get_rank([lowest_tax_id])[lowest_tax_id]
                        break
            esmecata_name = line['name']
            esmecata_rank = line['tax_rank']
            proteomes = line['proteome'].split(',')
            proteomes_status = all([proteome_status[proteome] for proteome in proteomes])

            proteome_numbers[observation_name] = [len(proteomes), tax_name, lowest_tax_rank, esmecata_name, esmecata_rank, proteomes_status]

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_proteomes', 'Input_taxon_Name', 'Taxon_rank', 'EsMeCaTa_used_taxon', 'EsMeCaTa_used_rank', 'only_reference_proteome_used'])
            for observation_name in proteome_numbers:
                csvwriter.writerow([observation_name, proteome_numbers[observation_name]])

    return proteome_numbers


def create_comp_taxonomy_file(association_taxon_id_json, proteomes_ids, tax_id_names, output_dir):
    """ Create taxonomy_diff.tsv file in proteome output folder. Compare Input taxa information to esmecata taxa
    information (taxa name, taxa ID, taxa rank) + precise OTUs associated.

    Args:
        association_taxon_id_json (dict): observation name and dictionary with mapping between taxon name and taxon
            ID (with remove rank specified)
        proteomes_ids (dict): observation name (key) associated with proteome IDs
        tax_id_names (dict): associate tax id to tax name
        output_dir (str): pathname to the output folder
    """
    ncbi = NCBITaxa()

    d_tax = dict()
    for observation_name in proteomes_ids:
        reversed_affiliation_taxa = list(reversed(association_taxon_id_json[observation_name]))
        for tax_name in reversed_affiliation_taxa:
            if tax_name != 'unknown' and association_taxon_id_json[observation_name][tax_name] != ['not_found']:
                if tax_name not in d_tax.keys():
                    tax_id = association_taxon_id_json[observation_name][tax_name][0]
                    tax_rank = ncbi.get_rank([tax_id])
                    if tax_id in tax_rank.keys():
                        tax_rank = tax_rank[tax_id]
                    else:
                        tax_rank = ''
                    esmecata_tax_id = int(proteomes_ids[observation_name][0])
                    esmecata_tax_name = tax_id_names[esmecata_tax_id]
                    esmecata_tax_rank = ncbi.get_rank([esmecata_tax_id])[esmecata_tax_id]
                    d_tax[tax_name] = {'Esmecata ID': esmecata_tax_id, 'Esmecata Rank': esmecata_tax_rank,
                                       'Esmecata Name': esmecata_tax_name, 'Input ID': tax_id, 'Input Rank': tax_rank}
                if 'obs' not in d_tax[tax_name].keys():
                    d_tax[tax_name]['obs'] = list()
                d_tax[tax_name]['obs'].append(observation_name)
                break

    output = os.path.join(output_dir, 'taxonomy_diff.tsv')
    with open(output, 'w') as f:
        f.write('\t'.join(['Input Name', 'Esmecata Name', 'Input Rank', 'Esmecata Rank',
                           'Input ID', 'Esmecata ID', 'Observation Names']))
        for name, ref in d_tax.items():
            f.write('\n' + '\t'.join([name, ref['Esmecata Name'], ref['Input Rank'], ref['Esmecata Rank'],
                                      str(ref['Input ID']), str(ref['Esmecata ID']), ';'.join(ref['obs'])]))


def check_proteomes(input_file, output_folder, busco_percentage_keep=80,
                        ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        limit_maximal_number_proteomes=99, rank_limit=None, minimal_number_proteomes=1,
                        update_affiliations=None, option_bioservices=None):
    """From a tsv file with taxonomic affiliations check the associated proteomes for the taxa.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold).
        ignore_taxadb_update (bool): option to ignore ete3 taxa database update.
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes).
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint).
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taxon to be kept.
        update_affiliations (str): option to update taxonomic affiliations.
        option_bioservices (bool): use bioservices instead of manual queries.
    """
    check_starttime = time.time()

    logger.info('|EsMeCaTa|proteomes| Begin proteomes search.')

    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[429, 500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    if is_valid_file(input_file) is False:
        logger.critical('|EsMeCaTa|proteomes| The input %s is not a valid file pathname.', input_file)
        sys.exit()

    ete3_database_update(ignore_taxadb_update)

    if minimal_number_proteomes >= limit_maximal_number_proteomes:
        logger.critical('|EsMeCaTa|proteomes| Error minimal_number_proteomes (%d) can not be superior or equal to limit_maximal_number_proteomes (%d).', minimal_number_proteomes, limit_maximal_number_proteomes)
        sys.exit()

    is_valid_dir(output_folder)

    proteomes_folder = os.path.join(output_folder, 'proteomes')
    is_valid_dir(proteomes_folder)

    proteomes_description_folder = os.path.join(output_folder, 'proteomes_description')
    is_valid_dir(proteomes_description_folder)

    known_extensions = ['.xlsx', '.tsv', '.csv']
    file_name, file_extension = os.path.splitext(input_file)

    if file_extension not in known_extensions:
        logger.critical('|EsMeCaTa|proteomes| ERROR: Extension ({0}) of input file ({1}) is not part of the compatible extensions ({2}).'.format(file_extension, input_file, ','.join(known_extensions)))
        sys.exit(1)

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

    if not os.path.exists(proteome_tax_id_file):
        ncbi = NCBITaxa()

        tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, update_affiliations, ncbi, output_folder)

        json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

        if rank_limit:
            json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

        json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
        with open(json_log, 'w') as ouput_file:
            json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    if not os.path.exists(proteome_tax_id_file):
        proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_taxonomic_affiliations, ncbi, proteomes_description_folder,
                                                        busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint,
                                                        limit_maximal_number_proteomes, minimal_number_proteomes, session, option_bioservices)

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

        create_comp_taxonomy_file(association_taxon_id_json=json_taxonomic_affiliations,
                                  proteomes_ids=proteomes_ids,
                                  tax_id_names=tax_id_names,
                                  output_dir=output_folder)

    else:
        proteome_to_download = []
        proteomes_ids = {}
        single_proteomes = {}
        already_download_proteomes = os.listdir(proteomes_folder)
        with open(proteome_tax_id_file, 'r') as proteome_tax_file:
            csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
            for line in csvreader:
                observation_name = line['observation_name']
                name = line['name']
                tax_id = line['tax_id']
                proteomes = line['proteome'].split(',')
                for proteome in proteomes:
                    if proteome not in already_download_proteomes:
                        proteome_to_download.append(proteome)

                if len(proteomes) == 1:
                    single_proteomes[observation_name] = (tax_id, proteomes)
                proteomes_ids[observation_name] = (tax_id, proteomes)
                logger.info('|EsMeCaTa|proteomes| %s will be associated with the taxon "%s" with %d proteomes.', observation_name, name, len(proteomes))

        proteome_to_download = set(proteome_to_download)
    # Create heatmap comparing input taxon and taxon used by esmecata to find proteomes.
    create_taxon_heatmap_from_complete_run(output_folder)

    check_endtime = time.time()
    check_duration = check_endtime - check_starttime
    logger.info('|EsMeCaTa|proteomes| Check step complete in {0}s.'.format(check_duration))

    return proteome_to_download, session


def retrieve_proteomes(input_file, output_folder, busco_percentage_keep=80,
                        ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None,
                        limit_maximal_number_proteomes=99, rank_limit=None, minimal_number_proteomes=1,
                        update_affiliations=None, option_bioservices=None):
    """From a tsv file with taxonomic affiliations find the associated proteomes and download them.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        busco_percentage_keep (float): BUSCO score to filter proteomes (proteomes selected will have a higher BUSCO score than this threshold).
        ignore_taxadb_update (bool): option to ignore ete3 taxa database update.
        all_proteomes (bool): Option to select all the proteomes (and not only preferentially reference proteomes).
        uniprot_sparql_endpoint (str): uniprot SPARQL endpoint to query (by default query Uniprot SPARQL endpoint).
        limit_maximal_number_proteomes (int): int threshold after which a subsampling will be performed on the data.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        minimal_number_proteomes (int): minimal number of proteomes required to be associated with a taxon for the taxon to be kept.
        update_affiliations (str): option to update taxonomic affiliations.
        option_bioservices (bool): use bioservices instead of manual queries.
    """
    starttime = time.time()

    proteome_to_download, session = check_proteomes(input_file, output_folder, busco_percentage_keep,
                            ignore_taxadb_update, all_proteomes, uniprot_sparql_endpoint,
                            limit_maximal_number_proteomes, rank_limit, minimal_number_proteomes,
                            update_affiliations, option_bioservices)

    logger.info('|EsMeCaTa|proteomes| Start downloading proteomes.')

    proteomes_folder = os.path.join(output_folder, 'proteomes')
    is_valid_dir(proteomes_folder)

    # Download all the proteomes in proteome folder.
    proteomes_already_downloaded = set([proteome_filename.replace('.faa.gz', '') for proteome_filename in os.listdir(proteomes_folder)])
    proteome_to_download = proteome_to_download - proteomes_already_downloaded

    if len(proteome_to_download) == 0:
        logger.info('|EsMeCaTa|proteomes| All proteomes already downloaded.')
    else:
        logger.info('|EsMeCaTa|proteomes| Downloading %d proteomes', len(proteome_to_download))

    for index, proteome in enumerate(proteome_to_download):
        output_proteome_file = os.path.join(proteomes_folder, proteome+'.faa.gz')
        if not os.path.exists(output_proteome_file):
            if uniprot_sparql_endpoint is not None:
                sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint)
            else:
                if option_bioservices is None:
                    http_str = 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:{0}&format=fasta&compressed=true'.format(proteome)
                    proteome_response = session.get(http_str)
                    with open(output_proteome_file, 'wb') as f:
                        f.write(proteome_response.content)
                else:
                    import bioservices
                    uniprot_bioservices = bioservices.UniProt()
                    data_fasta = uniprot_bioservices.search(f'(proteome:{proteome})', database='uniprot',
                                                            frmt='fasta', compress=True, progress=False)
                    with open(output_proteome_file, 'wb') as f:
                        f.write(data_fasta)

            logger.info('|EsMeCaTa|proteomes| Downloaded %d on %d proteomes',index+1, len(proteome_to_download))
        time.sleep(1)

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_file': input_file, 'output_folder': output_folder, 'busco_percentage_keep': busco_percentage_keep,
                    'ignore_taxadb_update': ignore_taxadb_update, 'all_proteomes': all_proteomes, 'uniprot_sparql_endpoint': uniprot_sparql_endpoint,
                    'limit_maximal_number_proteomes': limit_maximal_number_proteomes, 'rank_limit': rank_limit,
                    'minimal_number_proteomes': minimal_number_proteomes}

    # Collect dependencies metadata.
    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['biopython'] = biopython_version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['ete3'] = ete3_version
    options['tool_dependencies']['python_package']['pandas'] = pd.__version__
    options['tool_dependencies']['python_package']['requests'] = requests.__version__
    options['tool_dependencies']['python_package']['SPARQLWrapper'] = sparqlwrapper_version
    if option_bioservices is True:
        import bioservices
        options['tool_dependencies']['python_package']['bioservices'] = bioservices.version


    if uniprot_sparql_endpoint:
        uniprot_releases = get_sparql_uniprot_release(uniprot_sparql_endpoint, options)
    else:
        uniprot_releases = get_rest_uniprot_release(options)

    stat_file = os.path.join(output_folder, 'stat_number_proteome.tsv')
    compute_stat_proteomes(output_folder, stat_file)

    endtime = time.time()
    duration = endtime - starttime
    uniprot_releases['esmecata_proteomes_duration'] = duration

    uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_proteomes.json')
    if os.path.exists(uniprot_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_proteomes' in metadata_file]
        uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_proteomes_{0}.json'.format(len(metadata_files)))
        with open(uniprot_metadata_file, 'w') as ouput_file:
            json.dump(uniprot_releases, ouput_file, indent=4)
    else:
        with open(uniprot_metadata_file, 'w') as ouput_file:
            json.dump(uniprot_releases, ouput_file, indent=4)

    logger.info('|EsMeCaTa|proteomes| Proteome step complete in {0}s.'.format(duration))
