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
import json
import logging
import os
import re
import requests
import time
import sys
import urllib.parse
import urllib.request

from SPARQLWrapper import __version__ as sparqlwrapper_version

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_dir, send_uniprot_sparql_query, urllib_query
from esmecata import __version__ as esmecata_version

URLLIB_HEADERS = {'User-Agent': 'EsMeCaTa annotation v' + esmecata_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)


def rest_query_uniprot_to_retrieve_function(protein_queries):
    """REST query to get annotation from proteins.

    Args:
        protein_queries (list): list of proteins

    Returns:
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    output_dict = {}

    url = 'http://rest.uniprot.org/idmapping/run'

    # Column names can be found at: https://www.uniprot.org/help/uniprotkb_column_names
    # from and to are linked to mapping ID: https://www.uniprot.org/help/api%5Fidmapping

    params = {
        'from': 'UniProtKB_AC-ID',
        'to': 'UniProtKB',
        'ids': protein_queries
    }
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)

    data_js = json.load(urllib_query(req))
    job_id = data_js['jobId']

    http_check_status = 'http://rest.uniprot.org/idmapping/status/{0}'.format(job_id)

    check_status_response = requests.head(http_check_status)
    check_status_response.raise_for_status()

    # Not working because the code is 200 instead of 303...
    # Issue with protein_queries formatting?
    # Wait 5 secodns to check again return code.
    # But is there no better way?
    while check_status_response.status_code != 303:
        time.sleep(5)
        check_status_response = requests.head(http_check_status)

    if check_status_response.status_code == 303:
        http_download = 'http://rest.uniprot.org/idmapping/uniprotkb/results/stream/{0}'.format(job_id)
        response = requests.get(http_download)
        response.raise_for_status()
        data = response.json()

        results = {}
        for result in data['results']:
            protein_id = result['from']
            protein_data = result['to']
            reviewed = protein_data['entryType']
            if 'reviewed' in reviewed:
                review = True
            else:
                review = False
            protein_description = protein_data['proteinDescription']
            if 'recommendedName' in protein_description:
                if 'ecNumbers' in protein_description['recommendedName']:
                    protein_ecs = list(set([ecnumber['value'] for ecnumber in protein_description['recommendedName']['ecNumbers']]))
                else:
                    protein_ecs = []
                protein_fullname = protein_description['recommendedName']['fullName']['value']
            else:
                protein_ecs = []
                protein_fullname = ''

            gene_names = [gene['geneName']['value'] for gene in protein_data['genes'] if 'geneName' in gene]
            if len(gene_names) > 0:
                gene_name = gene_names[0]
            else:
                gene_name = ''
            protein_xrefs = protein_data['uniProtKBCrossReferences']

            rhea_ids = []
            if 'comments' in protein_data:
                protein_comments = protein_data['comments']
                protein_reactions = [comment['reaction'] for comment in protein_comments if comment['commentType'] == 'CATALYTIC ACTIVITY']
                for reaction in protein_reactions:
                    if 'reactionCrossReferences' in reaction:
                        rhea_ids.extend([dbxref['id'] for dbxref in reaction['reactionCrossReferences'] if dbxref['database'] == 'Rhea' and 'RHEA:' in dbxref['id']])
                    if 'ecNumber' in reaction:
                        protein_ecs.append(reaction['ecNumber'])
                rhea_ids = list(set(rhea_ids))
            protein_ecs = list(set(protein_ecs))
            gos = list(set([xref['id'] for xref in protein_xrefs if xref['database'] == 'GO']))
            interpros = list(set([xref['id'] for xref in protein_xrefs if xref['database'] == 'InterPro']))
            results[protein_id] = [protein_fullname, review, gos, protein_ecs, interpros, rhea_ids, gene_name]
            output_dict.update(results)

    return output_dict


def sparql_query_uniprot_to_retrieve_function(proteomes, uniprot_sparql_endpoint):
    """SPARQL query to get annotation from proteomes.

    Args:
        proteomes (list): list of proteomes
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database

    Returns:
        dict: annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    proteomes = ' '.join(['( proteome:'+proteome+' )' for proteome in proteomes])

    uniprot_sparql_function_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX proteome: <http://purl.uniprot.org/proteomes/>

    SELECT ?protein
        (GROUP_CONCAT(DISTINCT ?fullName; separator=";") AS ?name)
        (GROUP_CONCAT(DISTINCT ?goTerm; separator=";") AS ?go)
        (GROUP_CONCAT(DISTINCT ?ecNumber; separator=";") AS ?ec)
        (GROUP_CONCAT(DISTINCT ?ecNumber2; separator=";") AS ?ec2)
        (GROUP_CONCAT(DISTINCT ?interpro; separator=";") AS ?ipr)
        (GROUP_CONCAT(DISTINCT ?rheaReaction; separator=";") AS ?rhea)
        (GROUP_CONCAT(DISTINCT ?reviewed; separator=";") AS ?review)
        (GROUP_CONCAT(DISTINCT ?geneLabel; separator=";") AS ?geneName)
        (GROUP_CONCAT(DISTINCT ?subName; separator=";") AS ?submitName)

    FROM <http://sparql.uniprot.org/uniprot>
    FROM <http://sparql.uniprot.org/proteomes>

    WHERE {{
        ?protein a up:Protein ;
            up:proteome ?genomicComponent .
        ?proteome skos:narrower ?genomicComponent .
        OPTIONAL {{
            ?protein up:reviewed ?reviewed  .
        }}
        OPTIONAL {{
            ?protein up:annotation ?annot_catalytic .
            ?annot_catalytic a up:Catalytic_Activity_Annotation ;
                up:catalyticActivity ?catalyticAct .
            ?catalyticAct up:catalyzedReaction ?rheaReaction .
            ?catalyticAct up:enzymeClass ?ecNumber2 .
        }}
        OPTIONAL {{
            ?protein up:classifiedWith ?goTerm .
            FILTER (regex(str(?goTerm), "GO")) .
        }}
        OPTIONAL {{
            ?protein rdfs:seeAlso ?interpro .
            FILTER (regex(str(?interpro), "interpro")) .
        }}
        OPTIONAL {{
            ?protein up:enzyme ?ecNumber .
        }}
        OPTIONAL {{
            ?protein up:recommendedName ?recommendName .
            ?recommendName up:fullName ?fullName .
        }}
        OPTIONAL {{
            ?protein up:submittedName ?submittedName .
            ?submittedName up:fullName ?subName .
        }}
        OPTIONAL {{
            ?protein up:encodedBy ?gene .
            ?gene skos:prefLabel ?geneLabel .
        }}
    VALUES (?proteome) {{ {0} }}
    }}
    GROUP BY ?protein
    """.format(proteomes)

    csvreader = send_uniprot_sparql_query(uniprot_sparql_function_query, uniprot_sparql_endpoint)

    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        protein_name = line[1].split('^^')[0]
        go_terms = [go_uri.split('obo/')[1].replace('_', ':') for go_uri in line[2].split('^^')[0].split(';') if go_uri != '']
        ec_numbers = [ec_uri.split('enzyme/')[1] for ec_uri in line[3].split('^^')[0].split(';') if ec_uri != '']
        ec_numbers.extend([ec_uri.split('enzyme/')[1] for ec_uri in line[4].split('^^')[0].split(';') if ec_uri != ''])
        ec_numbers = list(set(ec_numbers))
        interpros = [ipr_uri.split('interpro/')[1] for ipr_uri in line[5].split('^^')[0].split(';') if ipr_uri != '']
        rhea_ids = [rhea_uri.split('rdf.rhea-db.org/')[1] for rhea_uri in line[6].split('^^')[0].split(';') if rhea_uri != '']
        review = line[7].split('^^')[0]
        if review == '0':
            review = False
        elif review == '1':
            review = True
        gene_name = line[8].split('^^')[0]

        submitted_name = [submitted_name for submitted_name in line[9].split('^^')[0].split(';') if submitted_name != '']
        if submitted_name != []:
            submitted_name = submitted_name[0]
        else:
            submitted_name = ''

        # If there is no recommendedName, use submittedName instead.
        if protein_name == '':
            if submitted_name != '':
                protein_name = submitted_name

        yield protein_id, [protein_name, review, go_terms, ec_numbers, interpros, rhea_ids, gene_name]


def sparql_query_uniprot_annotation_uniref(proteomes, uniref_output_dict, uniprot_sparql_endpoint):
    """SPARQL query to get annotation from proteomes and uniref.

    Args:
        proteomes (list): list of proteomes
        uniref_output_dict (dict): annotation dict: protein as key and annotation as value
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database

    Returns:
        uniref_output_dict (dict): annotation dict: protein as key and annotation as value
    """
    proteomes = ' '.join(['( proteome:'+proteome+' )' for proteome in proteomes])

    sparql_query_uniref = """PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX proteome: <http://purl.uniprot.org/proteomes/>

        SELECT ?protein
            (GROUP_CONCAT(DISTINCT ?goTerm; separator=";") AS ?go)
            (GROUP_CONCAT(DISTINCT ?ecNumber; separator=";") AS ?ec)
            (GROUP_CONCAT(DISTINCT ?ecNumber2; separator=";") AS ?ec2)
            (GROUP_CONCAT(DISTINCT ?cluster; separator=";") AS ?cl)
            (GROUP_CONCAT(DISTINCT ?member; separator=";") AS ?representativeMember)

        FROM <http://sparql.uniprot.org/uniref>
        FROM <http://sparql.uniprot.org/uniprot>
        FROM <http://sparql.uniprot.org/proteomes>

        WHERE
        {{
            ?cluster up:member/up:sequenceFor ?protein ;
                            up:identity ?identity .
            ?protein a up:Protein ;
                up:proteome ?genomicComponent .
            ?proteome skos:narrower ?genomicComponent .

            ?member up:representativeFor ?cluster.

            OPTIONAL {{
                ?member up:annotation ?annot_catalytic .
                ?annot_catalytic a up:Catalytic_Activity_Annotation ;
                    up:catalyticActivity ?catalyticAct .
                ?catalyticAct up:catalyzedReaction ?rheaReaction .
                ?catalyticAct up:enzymeClass ?ecNumber2 .
            }}
            OPTIONAL {{
                ?member up:classifiedWith ?goTerm .
                FILTER (regex(str(?goTerm), "GO")) .
            }}
            OPTIONAL {{
                ?member up:enzyme ?ecNumber .
            }}
            FILTER (?identity > 0.9)
            VALUES (?proteome) {{ {0} }}
        }}
    GROUP BY ?protein
    """.format(proteomes)

    csvreader = send_uniprot_sparql_query(sparql_query_uniref, uniprot_sparql_endpoint)

    results = {}
    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        go_terms = [go_uri.split('obo/')[1].replace('_', ':') for go_uri in line[1].split('^^')[0].split(';') if go_uri != '']
        ec_numbers = [ec_uri.split('enzyme/')[1] for ec_uri in line[2].split('^^')[0].split(';') if ec_uri != '']
        ec_numbers.extend([ec_uri.split('enzyme/')[1] for ec_uri in line[3].split('^^')[0].split(';') if ec_uri != ''])
        ec_numbers = list(set(ec_numbers))

        cluster_id = line[4].split('/')[-1]
        representative_member = line[5].split('/')[-1]

        results[protein_id] = [go_terms, ec_numbers, cluster_id, representative_member]

    uniref_output_dict.update(results)

    return uniref_output_dict


def sparql_query_uniprot_expression(proteomes, expression_output_dict, uniprot_sparql_endpoint):
    """SPARQL query to get expression data from proteomes.

    Args:
        proteomes (list): list of proteomes
        expression_output_dict (dict): annotation dict: protein as key and expression data as value.
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database

    Returns:
        expression_output_dict (dict): annotation dict: protein as key and expression data as value.
    """
    proteomes = ' '.join(['( proteome:'+proteome+' )' for proteome in proteomes])

    sparql_query_expression = """PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX proteome: <http://purl.uniprot.org/proteomes/>

        SELECT ?protein
            (GROUP_CONCAT(DISTINCT ?induct_comment; separator=";") AS ?induction)
            (GROUP_CONCAT(DISTINCT ?tissue_spec_comment; separator=";") AS ?tissue_specificity)
            (GROUP_CONCAT(DISTINCT ?disruption_comment; separator=";") AS ?disruption)

        FROM <http://sparql.uniprot.org/uniprot>
        FROM <http://sparql.uniprot.org/uniref>
        FROM <http://sparql.uniprot.org/proteomes>

        WHERE {{
            ?protein a up:Protein ;
                up:proteome ?genomicComponent .
            ?proteome skos:narrower ?genomicComponent .

            OPTIONAL {{
                ?protein up:annotation ?annot_induction .
                ?annot_induction rdf:type up:Induction_Annotation .
                ?annot_induction rdfs:comment ?induct_comment .
            }}
            OPTIONAL {{
                ?protein up:annotation ?annot_tissue .
                ?annot_tissue rdf:type up:Tissue_Specificity_Annotation .
                ?annot_tissue rdfs:comment ?tissue_spec_comment .
            }}
            OPTIONAL {{
                ?protein up:annotation ?annot_phenotype .
                ?annot_phenotype rdf:type up:Disruption_Phenotype_Annotation .
                ?annot_phenotype rdfs:comment ?disruption_comment .
            }}
        VALUES (?proteome) {{ {0} }}
        }}
        GROUP BY ?protein
    """.format(proteomes)

    csvreader = send_uniprot_sparql_query(sparql_query_expression, uniprot_sparql_endpoint)

    results = {}
    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        induction = line[1].split('^^')[0]
        tissue_specificity = line[2].split('^^')[0]
        disruption = line[3].split('^^')[0]

        results[protein_id] = [induction, tissue_specificity, disruption]

    expression_output_dict.update(results)

    return expression_output_dict


def chunks(elements, n):
    """Yield successive n-sized chunks from list.
    Form: https://stackoverflow.com/a/312464

    Args:
        elements (list): list of elements (proteins or proteomes) to be split in chunks of size n
        n (int): size of the chunks

    Returns:
        list: list of elements of size n
    """
    for i in range(0, len(elements), n):
        yield elements[i:i + n]


def create_pathologic(base_filename, annotated_protein_to_keeps, pathologic_output_file):
    """Create pathologic files.

    Args:
        base_filename (str): observation name associated with the taxonomic affiliations and the annotation
        annotated_protein_to_keeps (dict): annotated protein associated with the annotation of their cluster.
        pathologic_output_file (str): pathname to the pathologic write.
    """
    with open(pathologic_output_file, 'w', encoding='utf-8') as element_file:
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        element_file.write(';; ' + base_filename + '\n')
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        for protein in annotated_protein_to_keeps:
            element_file.write('ID\t' + protein + '\n')
            if annotated_protein_to_keeps[protein][3] != '':
                element_file.write('NAME\t' + annotated_protein_to_keeps[protein][3] + '\n')
            else:
                element_file.write('NAME\t' + protein + '\n')
            if annotated_protein_to_keeps[protein][0] != '':
                element_file.write('FUNCTION\t' + annotated_protein_to_keeps[protein][0] + '\n')
            element_file.write('PRODUCT-TYPE\tP' + '\n')
            element_file.write('PRODUCT-ID\tprot ' + protein + '\n')
            element_file.write('DBLINK\tUNIPROT:' + protein + '\n')
            for go in annotated_protein_to_keeps[protein][1]:
                element_file.write('GO\t' + go + '\n')
            for ec in annotated_protein_to_keeps[protein][2]:
                element_file.write('EC\t' + ec + '\n')
            element_file.write('//\n\n')


def compute_stat_annotation(annotation_reference_folder, stat_file=None):
    """Compute stat associated with the number of proteome for each taxonomic affiliations.

    Args:
        annotation_reference_folder (str): pathname to the annotation reference folder containing annotations for each cluster
        stat_file (str): pathname to the tsv stat file

    Returns:
        annotation_numbers (dict): dict containing observation names (as key) associated with GO Terms and EC (as value)
    """
    annotation_numbers = {}
    for infile in os.listdir(annotation_reference_folder):
        if '.tsv' in infile:
            annotation_input_file_path = os.path.join(annotation_reference_folder, infile)
            infile_gos = []
            infile_ecs = []
            with open(annotation_input_file_path, 'r') as open_annotation_input_file_path:
                csvreader = csv.reader(open_annotation_input_file_path, delimiter='\t')
                next(csvreader)
                for line in csvreader:
                    gos = line[3].split(',')
                    ecs = line[4].split(',')
                    infile_gos.extend(gos)
                    infile_ecs.extend(ecs)
            infile_gos = set([go for go in infile_gos if go != ''])
            infile_ecs = set([ec for ec in infile_ecs if ec != ''])
            annotation_numbers[infile.replace('.tsv','')] = (len(infile_gos), len(infile_ecs))

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_go_terms', 'Number_ecs'])
            for observation_name in annotation_numbers:
                csvwriter.writerow([observation_name, annotation_numbers[observation_name][0], annotation_numbers[observation_name][1]])

    return annotation_numbers


def extract_protein_cluster(reference_protein_pathname):
    """Extract proteins from mmseqs tabulated file

    Args:
        reference_protein_pathname (str): pathname to mmseqs protein cluster tabulated file

    Returns:
        reference_proteins (dict): dict containing representative protein IDs (as key) associated with proteins of the cluster
        set_proteins (set): set all of proteins present in the mmseqs tabulated file
    """
    proteins = []

    reference_proteins = {}
    with open(reference_protein_pathname, 'r') as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        for line in csvreader:
            if line != '':
                proteins.extend(line)
                reference_proteins[line[0]] = line[1:]

    set_proteins = set(proteins)

    return reference_proteins, set_proteins


def search_already_annotated_protein_in_file(set_proteins, already_annotated_proteins_in_file, annotation_folder, output_dict):
    """Use protein annotations already retrieved from UniProt to avoid new queries.
    This option read annotation file to retrived annotation, it is slower thant the second one (search_already_annotated_protein)) but less heavy on memory.

    Args:
        set_proteins (set): set all of proteins present in the mmseqs tabulated file
        already_annotated_proteins_in_file (dict): the already annotated protein ID (as key) and the annotation file containing its annotation (as value)
        annotation_folder (str): pathname to the annotation folder to retrieve annotation file
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])

    Returns:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    # Retrieve already annotated protein using the annotation files from the annotation_folder.
    alreay_annotated_set = set(already_annotated_proteins_in_file.keys())
    reference_files = {}
    for protein in set_proteins:
        if protein in alreay_annotated_set:
            reference_file = already_annotated_proteins_in_file[protein]
            if reference_file not in reference_files:
                reference_files[reference_file] = [protein]
            else:
                reference_files[reference_file].append(protein)

    protein_to_remove = []
    for reference_file in reference_files:
        already_annotated_file = os.path.join(annotation_folder, reference_file+'.tsv')
        with open(already_annotated_file) as already_process_file:
            csvreader = csv.reader(already_process_file, delimiter='\t')
            for line in csvreader:
                if line[0] in set(reference_files[reference_file]):
                    protein_name = line[1]
                    review = line[2]
                    gos = line[3].split(',')
                    ecs = line[4].split(',')
                    interpros = line[5].split(',')
                    rhea_ids = line[6].split(',')
                    gene_name = line[7]
                    output_dict[line[0]] = [protein_name, review, gos, ecs, interpros, rhea_ids, gene_name]
                    protein_to_remove.append(line[0])

    protein_to_search_on_uniprots = set_proteins.difference(set(protein_to_remove))

    return protein_to_search_on_uniprots, output_dict


def search_already_annotated_protein(set_proteins, already_annotated_proteins, output_dict):
    """Use protein annotations already retrieved from UniProt to avoid new queries.
    This option is not used as I fear that it takes a lot of memory to keep each proteins found.

    Args:
        set_proteins (set): set all of proteins present in the mmseqs tabulated file
        already_annotated_proteins (dict): the already annotated protein ID (as key) and the annotation associated (as value)
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])

    Returns:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    alreay_annotated_set = set(already_annotated_proteins.keys())
    protein_to_remove = []
    for protein in set_proteins:
        if protein in alreay_annotated_set:
            output_dict[protein] = already_annotated_proteins[protein]
            protein_to_remove.append(protein)

    protein_to_search_on_uniprots = set_proteins.difference(set(protein_to_remove))

    return protein_to_search_on_uniprots, output_dict


def query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict):
    """Query UniProt with REST to find protein annotations

    Args:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])

    Returns:
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    # The limit of 10 000 proteins per query comes from the help of Uniprot (inferior to 20 000):
    # https://www.uniprot.org/help/uploadlists
    if len(protein_to_search_on_uniprots) < 15000:
        protein_queries = ','.join(protein_to_search_on_uniprots)
        tmp_output_dict = rest_query_uniprot_to_retrieve_function(protein_queries)
        output_dict.update(tmp_output_dict)
        time.sleep(1)
    else:
        protein_chunks = chunks(list(protein_to_search_on_uniprots), 15000)
        for chunk in protein_chunks:
            protein_queries = ','.join(chunk)
            tmp_output_dict = rest_query_uniprot_to_retrieve_function(protein_queries)
            output_dict.update(tmp_output_dict)
            time.sleep(1)

    return output_dict


def query_uniprot_annotation_sparql(proteomes, uniprot_sparql_endpoint, output_dict):
    """Query UniProt with SPARQL to find protein annotations

    Args:
        proteomes (dict): list of proteomes to use as query
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])

    Returns:
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
    """
    # Query Uniprot to get the annotation of each proteins.
    # Add a limit of 100 proteomes per query.
    if len(proteomes) > 100:
        proteomes_chunks = chunks(list(proteomes), 100)
        for proteome_chunk in proteomes_chunks:
            tmp_output_dict = dict(sparql_query_uniprot_to_retrieve_function(proteome_chunk, uniprot_sparql_endpoint))
            output_dict.update(tmp_output_dict)
            time.sleep(1)
    else:
        tmp_output_dict = dict(sparql_query_uniprot_to_retrieve_function(proteomes, uniprot_sparql_endpoint))
        output_dict.update(tmp_output_dict)
        time.sleep(1)

    return output_dict


def write_annotation_file(output_dict, annotation_file):
    """Write annotation file from annotated protein retrieverd from UniProt

    Args:
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
        annotation_file (str): pathname to output tabulated file
    """
    with open(annotation_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        csvwriter.writerow(['protein_id', 'protein_name', 'review', 'gos', 'ecs', 'interpros', 'rhea_ids', 'gene_name'])
        for protein in output_dict:
            protein_name = output_dict[protein][0]
            protein_review_satus = str(output_dict[protein][1])
            go_terms = ','.join(sorted(output_dict[protein][2]))
            ec_numbers = ','.join(sorted(output_dict[protein][3]))
            interpros = ','.join(sorted(output_dict[protein][4]))
            rhea_ids = ','.join(sorted(output_dict[protein][5]))
            gene_name = output_dict[protein][6]

            csvwriter.writerow([protein, protein_name, protein_review_satus, go_terms, ec_numbers, interpros, rhea_ids, gene_name])


def retrieve_annotation_from_uniref(proteomes, uniprot_sparql_endpoint, uniref_annotation_file):
    """Query UniProt uniref with SPARQL to find protein annotations

    Args:
        proteomes (dict): list of proteomes to use as query
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database
        uniref_annotation_file (str): tabulated file containing annotation

    Returns:
        uniref_output_dict (dict): annotation dict: protein as key and annotation as value
    """
    uniref_output_dict = {}
    uniref_output_dict = sparql_query_uniprot_annotation_uniref(proteomes, uniref_output_dict, uniprot_sparql_endpoint)

    with open(uniref_annotation_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        csvwriter.writerow(['protein_id', 'gos', 'ecs', 'uniref_cluster', 'representative_member'])
        for protein in uniref_output_dict:
            go_terms = ','.join(sorted(uniref_output_dict[protein][0]))
            ec_numbers = ','.join(sorted(uniref_output_dict[protein][1]))
            cluster_id = uniref_output_dict[protein][2]
            representative_member = uniref_output_dict[protein][3]
            csvwriter.writerow([protein, go_terms, ec_numbers, cluster_id, representative_member])

    return uniref_output_dict


def retrieve_expresion_from_uniprot(proteomes, uniprot_sparql_endpoint, expression_annotation_file):
    """Query UniProt uniref with SPARQL to find expression data associated with protein

    Args:
        proteomes (dict): list of proteomes to use as query
        uniprot_sparql_endpoint (str): SPARQL endpoint to uniprot database
        expression_annotation_file (str): tabulated file containing expression data

    Returns:
        expression_output_dict (dict): expression dict: protein as key and expression as value
    """
    expression_output_dict = {}
    expression_output_dict = sparql_query_uniprot_expression(proteomes, expression_output_dict, uniprot_sparql_endpoint)
    with open(expression_annotation_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        csvwriter.writerow(['protein_id', 'Induction_Annotation', 'Tissue_Specificity_Annotation', 'Disruption_Phenotype_Annotation'])
        for protein in expression_output_dict:
            induction = expression_output_dict[protein][0]
            tissue_specificity = expression_output_dict[protein][1]
            disruption = expression_output_dict[protein][2]
            csvwriter.writerow([protein, induction, tissue_specificity, disruption])

    return expression_output_dict


def propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict):
    """Propagate the annotation in the cluster. If propagate_annotation is None, will only return annotations from representative proteins.
    If propagate_annotation is 0 meaning all the annotations from all the protein in the cluster will be used (union of annotations).
    If propagate_annotation is 1 an annotation is kept only if it occurs in all protein of the clsuter (intersection of annotation).

    Args:
        output_dict (dict): annotation dict: protein as key and annotation as value ([function_name, review_status, [go_terms], [ec_numbers], [interpros], [rhea_ids], gene_name])
        reference_proteins (dict): dict containing representative protein IDs (as key) associated with proteins of the cluster
        propagate_annotation (float): float between 0 and 1. It is the ratio of proteins in the cluster that should have the annotation to keep this annotation.
        uniref_output_dict (dict): annotation dict: protein as key and annotation as value

    Returns:
        protein_annotations (dict): annotation dict: protein as key and annotation as value ([function_name, [go_terms], [ec_numbers], gene_name])
    """
    # For each reference protein, get the annotation of the other proteins clustered with it and add to its annotation.
    protein_annotations = {}
    for reference_protein in reference_proteins:
        if propagate_annotation is not None:
            gos = []
            ecs = []
            gene_names = []
            protein_names = []
            for protein in reference_proteins[reference_protein]:
                if protein in output_dict:
                    if output_dict[protein][0] != '':
                        protein_names.append(output_dict[protein][0])
                    if output_dict[protein][2] != []:
                        gos.append(list(set(output_dict[protein][2])))
                    if output_dict[protein][3] != []:
                        ecs.append(list(set(output_dict[protein][3])))
                    if output_dict[protein][6] != '':
                        gene_names.append(output_dict[protein][6])

            keep_gos = []
            keep_ecs = []
            keep_gene_names = ''
            keep_protein_names = ''
            # Propagate all the annotations (GO, EC) that have been found to occur in at least X proteins of the cluster. Where X is computed using the ratio given by the user and
            # the number of proteins in the cluster.
            if len(gos) > 0:
                all_gos = set([go for subgos in gos for go in subgos])
                keep_gos = [go for go in all_gos if sum(subgos.count(go) for subgos in gos) >= propagate_annotation * len(reference_proteins[reference_protein])]
            if len(ecs) > 0:
                all_ecs = set([ec for subecs in ecs for ec in subecs])
                keep_ecs = [ec for ec in all_ecs if sum(subecs.count(ec) for subecs in ecs) >= propagate_annotation * len(reference_proteins[reference_protein])]
            # Retrieve the gene and protein names by using the gene/protein name which the maximum occurrence in the proteins of the cluster.
            if len(gene_names) > 0:
                keep_gene_names = max(gene_names,key=gene_names.count)
            if len(protein_names) > 0:
                keep_protein_names = max(protein_names,key=protein_names.count)

            protein_annotations[reference_protein] = [keep_protein_names, keep_gos, keep_ecs, keep_gene_names]
        else:
            # Use only annotation from representative proteins.
            protein_name = output_dict[reference_protein][0]
            gos = output_dict[reference_protein][2]
            ecs = output_dict[reference_protein][3]
            gene_name = output_dict[reference_protein][6]
            protein_annotations[reference_protein] = [protein_name, gos, ecs, gene_name]

        # Propagate anntoations from Uniref
        if uniref_output_dict is not None:
            uniref_gos = uniref_output_dict[reference_protein][0]
            uniref_ecs = uniref_output_dict[reference_protein][1]

            protein_annotations[reference_protein][1].extend(uniref_gos)
            protein_annotations[reference_protein][2].extend(uniref_ecs)

            protein_annotations[reference_protein][1] = set(protein_annotations[reference_protein][1])
            protein_annotations[reference_protein][2] = set(protein_annotations[reference_protein][2])

    return protein_annotations


def write_annotation_reference(protein_annotations, annotation_reference_file, expression_output_dict):
    """Write the annotation associated with a cluster after propagation step.

    Args:
        protein_annotations (dict): annotation dict: protein as key and annotation as value ([function_name, [go_terms], [ec_numbers], gene_name])
        annotation_reference_file (str): pathname to output tabulated file
        expression_output_dict (dict): expression dict: protein as key and expression as value
    """
    with open(annotation_reference_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        if expression_output_dict:
            csvwriter.writerow(['protein', 'protein_name', 'gene_name', 'GO', 'EC', 'Induction', 'Tissue_Specificity', 'Disruption_Phenotype'])
        else:
            csvwriter.writerow(['protein', 'protein_name', 'gene_name', 'GO', 'EC'])
        for protein in protein_annotations:
            protein_name = protein_annotations[protein][0]
            gene_name = protein_annotations[protein][3]
            gos = ','.join(sorted(list(protein_annotations[protein][1])))
            ecs = ','.join(sorted(list(protein_annotations[protein][2])))
            if expression_output_dict:
                induction = expression_output_dict[protein][0]
                tissue_specificity = expression_output_dict[protein][1]
                disruption = expression_output_dict[protein][2]
                csvwriter.writerow([protein, protein_name, gene_name, gos, ecs, induction, tissue_specificity, disruption])
            else:
                csvwriter.writerow([protein, protein_name, gene_name, gos, ecs])


def write_pathologic_file(protein_annotations, pathologic_folder, base_filename, set_proteins):
    """Write the annotation associated with a cluster after propagation step into pathologic file for run on Pathway Tools.

    Args:
        protein_annotations (dict): annotation dict: protein as key and annotation as value ([function_name, [go_terms], [ec_numbers], gene_name])
        pathologic_folder (str): pathname to output pathologic folder
        base_filename (str): observation name
        set_proteins (set): set of proteins in the mmseqs tabulated file
    """
    # Create PathoLogic file and folder for each input.
    annotated_protein_to_keeps = {protein: protein_annotations[protein] for protein in protein_annotations if protein in set_proteins}
    if len(annotated_protein_to_keeps) > 0:
        pathologic_organism_folder = os.path.join(pathologic_folder, base_filename)
        is_valid_dir(pathologic_organism_folder)
        pathologic_file = os.path.join(pathologic_organism_folder, base_filename+'.pf')
        create_pathologic(base_filename, annotated_protein_to_keeps, pathologic_file)
    elif len(annotated_protein_to_keeps) == 0:
        logger.critical('|EsMeCaTa|annotation| No reference proteins for {0}, esmecata will not create a pathologic folder for it.'.format(base_filename))


def annotate_proteins(input_folder, output_folder, uniprot_sparql_endpoint, propagate_annotation, uniref_annotation, expression_annotation):
    """Write the annotation associated with a cluster after propagation step into pathologic file for run on Pathway Tools.

    Args:
        input_folder (str): pathname to mmseqs clustering output folder as it is the input of annotation
        output_folder (str): pathname to the output folder
        propagate_annotation (float): float between 0 and 1. It is the ratio of proteins in the cluster that should have the annotation to keep this annotation.
        uniref_annotation (bool): option to use uniref annotation to add annotation
        expression_annotation (bool): option to add expression annotation from uniprot
    """
    starttime = time.time()

    if uniprot_sparql_endpoint is None and uniref_annotation is not None:
        logger.critical('|EsMeCaTa|annotation| At this moment, --uniref option needs to be used with --sparql option.')
        sys.exit()

    if uniprot_sparql_endpoint is None and expression_annotation is not None:
        logger.critical('|EsMeCaTa|annotation| At this moment, --expression option needs to be used with --sparql option.')
        sys.exit()

    is_valid_dir(output_folder)

    annotation_folder = os.path.join(output_folder, 'annotation')
    is_valid_dir(annotation_folder)

    annotation_reference_folder = os.path.join(output_folder, 'annotation_reference')
    is_valid_dir(annotation_reference_folder)

    pathologic_folder = os.path.join(output_folder, 'pathologic')
    is_valid_dir(pathologic_folder)

    if uniref_annotation:
        uniref_protein_path = os.path.join(output_folder, 'uniref_annotation')
        is_valid_dir(uniref_protein_path)

    if expression_annotation:
        expression_protein_path = os.path.join(output_folder, 'expression_annotation')
        is_valid_dir(expression_protein_path)

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder, 'uniprot_sparql_endpoint': uniprot_sparql_endpoint,
                'propagate_annotation': propagate_annotation, 'uniref_annotation': uniref_annotation, 'expression_annotation': expression_annotation}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['SPARQLWrapper'] = sparqlwrapper_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__

    if uniprot_sparql_endpoint:
        uniprot_releases = get_sparql_uniprot_release(uniprot_sparql_endpoint, options)
    else:
        uniprot_releases = get_rest_uniprot_release(options)

    proteome_tax_id_file = os.path.join(input_folder, 'proteome_tax_id.tsv')

    if uniprot_sparql_endpoint:
        input_proteomes = {}
        with open(proteome_tax_id_file, 'r') as tax_id_file:
            csvreader = csv.reader(tax_id_file, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                input_proteomes[line[0]] = line[4].split(',')

    reference_protein_path = os.path.join(input_folder, 'reference_proteins')
    # There is 2 ways to handle already annotated proteins:
    # one keeping all the proteins in a dicitonary during all the analysis (I fear an issue with the memory).
    already_annotated_proteins = {}
    # a second that use the annotation from the annotation files associated with the proteins.
    # It is slower as it reads the annotation file to retrieve the annotation, but I think it is less heavy on the memory.
    already_annotated_proteins_in_file = {}
    for input_file in os.listdir(reference_protein_path):
        base_file = os.path.basename(input_file)
        base_filename = os.path.splitext(base_file)[0]

        output_dict = {}
        reference_protein_pathname = os.path.join(reference_protein_path, input_file)
        reference_proteins, set_proteins = extract_protein_cluster(reference_protein_pathname)

        # First method to handle protein already annotated
        #protein_to_search_on_uniprots, output_dict = search_already_annotated_protein(set_proteins, already_annotated_proteins, output_dict)
        # Second method to handle protein already annotated
        protein_to_search_on_uniprots, output_dict = search_already_annotated_protein_in_file(set_proteins, already_annotated_proteins_in_file, annotation_folder, output_dict)

        if uniprot_sparql_endpoint:
            proteomes = input_proteomes[base_filename]
            output_dict = query_uniprot_annotation_sparql(proteomes, base_filename, uniprot_sparql_endpoint, output_dict)
        else:
            output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict)

        # First method to handle protein already annotated
        #already_annotated_proteins.update({protein: output_dict[protein] for protein in output_dict if protein not in already_annotated_proteins})
        # Second method to handle protein already annotated
        for protein in output_dict:
            already_annotated_proteins_in_file[protein] = base_filename

        annotation_file = os.path.join(annotation_folder, base_filename+'.tsv')
        write_annotation_file(output_dict, annotation_file)
        if uniref_annotation and uniprot_sparql_endpoint:
            uniref_annotation_file = os.path.join(uniref_protein_path, base_filename+'.tsv')
            uniref_output_dict = retrieve_annotation_from_uniref(proteomes, uniprot_sparql_endpoint, uniref_annotation_file)
        else:
            uniref_output_dict = None
        if expression_annotation and uniprot_sparql_endpoint:
            expression_annotation_file = os.path.join(expression_protein_path, base_filename+'.tsv')
            expression_output_dict = retrieve_expresion_from_uniprot(proteomes, uniprot_sparql_endpoint, expression_annotation_file)
        else:
            expression_output_dict = None
        protein_annotations = propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict)
        annotation_reference_file = os.path.join(annotation_reference_folder, base_filename+'.tsv')
        write_annotation_reference(protein_annotations, annotation_reference_file, expression_output_dict)
        write_pathologic_file(protein_annotations, pathologic_folder, base_filename, set_proteins)

        gos = [go for protein in protein_annotations for go in protein_annotations[protein][1]]
        unique_gos = set(gos)
        ecs = [ec for protein in protein_annotations for ec in protein_annotations[protein][2]]
        unique_ecs = set(ecs)
        logger.info('|EsMeCaTa|annotation| {0} Go Terms (with {1} unique GO Terms) and {2} EC numbers (with {3} unique EC) associated with {4}.'.format(len(gos),
                                                                                                len(unique_gos), len(ecs), len(unique_ecs), base_filename))

    # Create mpwt taxon ID file.
    clustering_taxon_id = {}
    with open(proteome_tax_id_file, 'r') as input_taxon_id_file:
        taxon_id_csvreader = csv.reader(input_taxon_id_file, delimiter='\t')
        next(taxon_id_csvreader)
        for line in taxon_id_csvreader:
            clustering_taxon_id[line[0]] = line[2]

    pathologic_taxon_id_file = os.path.join(pathologic_folder, 'taxon_id.tsv')
    with open(pathologic_taxon_id_file, 'w') as taxon_id_file:
        taxon_id_csvwriter = csv.writer(taxon_id_file, delimiter='\t')
        taxon_id_csvwriter.writerow(['species', 'taxon_id'])
        for species in clustering_taxon_id:
            taxon_id_csvwriter.writerow([species, clustering_taxon_id[species]])

    stat_file = os.path.join(output_folder, 'stat_number_annotation.tsv')
    compute_stat_annotation(annotation_reference_folder, stat_file)

    endtime = time.time()
    duration = endtime - starttime
    uniprot_releases['esmecata_annotation_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_annotation.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    logger.info('|EsMeCaTa|annotation| Annotation complete.')
