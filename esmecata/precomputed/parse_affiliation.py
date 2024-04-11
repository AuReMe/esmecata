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

import csv
import json
import logging
import os
import pandas as pd
import sys
import time

from ete3 import __version__ as ete3_version
from ete3 import NCBITaxa, is_taxadb_up_to_date

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir, send_uniprot_sparql_query
from esmecata.core.proteomes import associate_taxon_to_taxon_id, disambiguate_taxon, filter_rank_limit
from esmecata import __version__ as esmecata_version

logger = logging.getLogger(__name__)

def get_taxon_database(database_taxon_file_path):
    """ Extract taxon ID contained in precomputed database.

    Args:
        database_taxon_file_path (str): pathanme to esmecata precomputed database

    Returns:
        database_taxon_ids (dict): list of taxon IDs contained in the database
    """
    database_taxon_ids = []
    with open(database_taxon_file_path, 'r') as open_database_taxon_file_path:
        csvreader = csv.DictReader(open_database_taxon_file_path, delimiter='\t')
        for line in csvreader:
            database_taxon_ids.append(line['taxon_id'])

    return database_taxon_ids


def find_proteomes_tax_ids(json_taxonomic_affiliations, database_taxon_ids):
    """Find proteomes associated with taxonomic affiliations

    Args:
        json_taxonomic_affiliations (dict): observation name and dictionary with mapping between taxon name and taxon ID (with remove rank specified)
        database_taxon_ids (list): list of taxon IDs contained in the database

    Returns:
        association_taxon_database (dict): observation name (key) associated with tax_name and tax_id
    """
    logger.info('|EsMeCaTa|proteomes| Find observation name present in precomputed database.')

    association_taxon_database = {}
    for observation_name in json_taxonomic_affiliations:
        for tax_name in reversed(json_taxonomic_affiliations[observation_name]):
            tax_id = json_taxonomic_affiliations[observation_name][tax_name][0]

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in database_taxon_ids:
                association_taxon_database[observation_name] = (tax_name,  tax_id)
                logger.info('|EsMeCaTa|precomputed| "%s" present in database, %s will be associated with the taxon "%s".', tax_name, observation_name, tax_name)
                break

    return association_taxon_database

def precomputed_parse_affiliation(input_file, output_folder, database_taxon_file_path, rank_limit=None, update_affiliations=None):
    """From a tsv file with taxonomic affiliations find the associated proteomes and download them.

    Args:
        input_file (str): pathname to the tsv input file containing taxonomic affiliations.
        output_folder (str): pathname to the output folder.
        rank_limit (str): rank limit to filter the affiliations (keep this rank and all inferior ranks).
        update_affiliations (str): option to update taxonomic affiliations.
    """
    starttime = time.time()

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

    database_taxon_file_path = os.path.join()

    database_taxon_file_path = os.path.join('database_path', 'summary.tsv')
    database_taxon_ids = get_taxon_database(database_taxon_file_path)
    ncbi = NCBITaxa()

    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(taxonomies, update_affiliations, ncbi, output_folder)

    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)

    if rank_limit:
        json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, rank_limit)

    json_log = os.path.join(output_folder, 'association_taxon_taxID.json')
    with open(json_log, 'w') as ouput_file:
        json.dump(json_taxonomic_affiliations, ouput_file, indent=4)

    association_taxon_database = find_proteomes_tax_ids(json_taxonomic_affiliations, database_taxon_ids, ncbi)

    endtime = time.time()
    duration = endtime - starttime
    logger.info('|EsMeCaTa|proteomes| Proteome step complete in {0}s.'.format(duration))
