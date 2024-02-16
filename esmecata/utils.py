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

import argparse
import csv
import datetime
import logging
import os
import urllib.request
import sys
import time

from io import StringIO
from SPARQLWrapper import SPARQLWrapper, TSV
from socket import timeout

from esmecata import __version__ as esmecata_version

MIN_VAL = 0
MAX_VAL = 1
URLLIB_HEADERS = {'User-Agent': 'EsMeCaTa annotation v' + esmecata_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)


def range_limited_float_type(arg):
    """Type function for argparse - a float within some predefined bounds

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < MIN_VAL or f > MAX_VAL:
        raise argparse.ArgumentTypeError("Argument must be < " + str(MAX_VAL) + " and > " + str(MIN_VAL))
    return f


def limited_integer_type(arg):
    """Type function for argparse - an integer

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer number")
    return f


def is_valid_path(filepath):
    """Return True if filepath is valid
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it)
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True


def urllib_query(request, nb_retry=5):
    """Use urllib to query UniProt.

    Args:
        request (urllib.request.Request): request object
        nb_retry (int): number of retry to perform (default = 5)

    Returns:
        response (http.client.HTTPResponse): response returns by requests
    """
    passed = False

    if nb_retry == 0:
        sys.exit('5 retry attempts have been performed but were not successful, so esmecata has been stopped. You may try to relaunch it using the same command esmecata should resume.')
    try:
        response = urllib.request.urlopen(request)
        passed = True
    except timeout:
        logger.critical('Timeout occurs for query, try to relaunch query.')
        time.sleep(10)
        urllib_query(request, nb_retry-1)

    if passed is True:
        return response


def get_rest_uniprot_release(options):
    """Get the release version and date of Uniprot and Trembl and also the date of the query.

    Args:
        options (dict): dictionary with metadata on the esmecata run

    Returns:
        dict: metadata of Uniprot release
    """
    uniprot_releases = {}

    # Get Uniprot release version
    uniprot_urllib_request = urllib.request.Request('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt',
                                                    headers=URLLIB_HEADERS)
    uniprot_response = urllib_query(uniprot_urllib_request)
    uniprot_lines = uniprot_response.readlines()
    uniprot_release_number = uniprot_lines[0].decode('utf-8').split(' ')[3].replace('\n','')
    swissprot_release_number = uniprot_lines[1].decode('utf-8').split(' ')[2].replace('\n','')
    swissprot_release_date = uniprot_lines[1].decode('utf-8').split(' ')[4].replace('\n','')
    trembl_release_number = uniprot_lines[2].decode('utf-8').split(' ')[2].replace('\n','')
    trembl_release_date = uniprot_lines[2].decode('utf-8').split(' ')[4].replace('\n','')

    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')

    uniprot_releases['esmecata_query_system'] = 'REST queries on Uniprot'
    uniprot_releases['uniprot_release'] = uniprot_release_number
    uniprot_releases['access_time'] = date
    uniprot_releases['swissprot_release_number'] = swissprot_release_number
    uniprot_releases['swissprot_release_date'] = swissprot_release_date
    uniprot_releases['trembl_release_number'] = trembl_release_number
    uniprot_releases['trembl_release_date'] = trembl_release_date
    uniprot_releases['tool_options'] = options

    return uniprot_releases


def get_sparql_uniprot_release(uniprot_sparql_endpoint, options):
    """Get the release version and date of Uniprot and Trembl and also the date of the query.

    Args:
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database
        options (dict): dictionary with metadata on the esmecata run

    Returns:
        dict: metadata of Uniprot release
    """
    uniprot_releases = {}

    uniprot_sparql_query = """SELECT ?version
    WHERE
    {{
        [] <http://www.w3.org/2002/07/owl#versionInfo> ?version .
    }}
    """
    csvreader = send_uniprot_sparql_query(uniprot_sparql_query, uniprot_sparql_endpoint)

    uniprot_release_number = [line[0] for line in csvreader][0]
    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')

    uniprot_releases['esmecata_query_system'] = 'SPARQL queries'
    uniprot_releases['esmecata_sparql_endpoint'] = uniprot_sparql_endpoint
    uniprot_releases['uniprot_release'] = uniprot_release_number
    uniprot_releases['access_time'] = date
    uniprot_releases['tool_options'] = options

    return uniprot_releases


def send_uniprot_sparql_query(sparql_query, uniprot_sparql_endpoint):
    """Get the release version and date of Uniprot and Trembl and also the date of the query.

    Args:
        sparql_query (str): string containing SPARQL query
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database

    Returns:
        dict: metadata of Uniprot release
    """
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)
    sparql.agent = 'EsMeCaTa proteomes v' + esmecata_version + ', wrapper=' + sparql.agent

    sparql.setQuery(sparql_query)
    # Parse output.
    sparql.setReturnFormat(TSV)

    results = sparql.query().convert().decode('utf-8')
    if results != '':
        csvreader = csv.reader(StringIO(results), delimiter='\t')
        # Avoid header.
        next(csvreader)
    else:
        csvreader = []

    return csvreader
