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
from bioservices import KEGG, UniProt
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model, read_sbml_model

from esmecata.utils import is_valid_dir, urllib_query, chunks
from esmecata import __version__ as esmecata_version

URLLIB_HEADERS = {'User-Agent': 'EsMeCaTa annotation v' + esmecata_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

# Create KEGG instance of bioservices.KEEG.
kegg = KEGG()
uniprot = UniProt(verbose=False)

def rest_query_uniprot_to_retrieve_kegg_gene(protein_queries, beta=None):
    """REST query to get KEGG gene ID from proteins.

    Args:
        protein_queries (list): list of proteins
        beta (bool): option to use the new API of UniProt (in beta can be unstable)

    Returns:
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value
    """
    output_dict = {}

    # Column names can be found at: https://www.uniprot.org/help/uniprotkb_column_names
    # from and to are linked to mapping ID: https://www.uniprot.org/help/api%5Fidmapping
    if not beta:
        results = uniprot.mapping("ACC", "KEGG_ID", query=protein_queries)
        output_dict.update(results)
    else:
        url = 'http://rest.uniprot.org/beta/idmapping/run'

        # Column names can be found at: https://www.uniprot.org/help/uniprotkb_column_names
        # from and to are linked to mapping ID: https://www.uniprot.org/help/api%5Fidmapping

        params = {
            'from': 'UniProtKB_AC-ID',
            'to': 'KEGG',
            'ids': protein_queries
        }
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)

        data_js = json.load(urllib_query(req))
        job_id = data_js['jobId']

        http_check_status = 'http://rest.uniprot.org/beta/idmapping/status/{0}'.format(job_id)

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
            http_download = 'http://rest.uniprot.org/beta/idmapping/uniprotkb/results/stream/{0}'.format(job_id)
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


def query_uniprot_kegg_rest(protein_to_search_on_uniprots, beta, output_dict):
    """Query UniProt with REST to find KEGG gene ID

    Args:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        beta (bool): option to use the new API of UniProt (in beta can be unstable)
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value

    Returns:
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value
    """
    # The limit of 15 000 proteins per query comes from the help of Uniprot (inferior to 20 000):
    # https://www.uniprot.org/help/uploadlists
    if len(protein_to_search_on_uniprots) < 15000:
        if not beta:
            protein_queries = ' '.join(protein_to_search_on_uniprots)
        else:
            protein_queries = ','.join(protein_to_search_on_uniprots)
        tmp_output_dict = rest_query_uniprot_to_retrieve_kegg_gene(protein_queries, beta)
        output_dict.update(tmp_output_dict)
        time.sleep(1)
    else:
        protein_chunks = chunks(list(protein_to_search_on_uniprots), 15000)
        for chunk in protein_chunks:
            if not beta:
                protein_queries = ' '.join(chunk)
            else:
                protein_queries = ','.join(chunk)
            tmp_output_dict = rest_query_uniprot_to_retrieve_kegg_gene(protein_queries, beta)
            output_dict.update(tmp_output_dict)
            time.sleep(1)

    return output_dict


def retrieve_reactions(reaction_folder):
    """Using bioservices.KEGG retrieve all keg files associated to KEGG reactions

    Args:
        reaction_folder (str): output folder which will contains reaction file
    """
    is_valid_dir(reaction_folder)

    response_text = kegg.list('reaction')

    reaction_ids = [] 
    for line in response_text.splitlines():
        reaction_id = line.split('\t')[0].split(':')[1].strip()
        reaction_ids.append(reaction_id)

    for reaction_id in reaction_ids:
        response_text = kegg.get(reaction_id)
        reaction_file = os.path.join(reaction_folder, reaction_id+'.keg')
        with open(reaction_file, 'w') as output_file:
            output_file.write(response_text)


def extract_reaction(reaction_id, equation_text):
    """Using the EQUATION in the reaction keg file, extract the compounds from the reaction formula.

    Args:
        reaction_id (str): Id of the reaction
        equation_text (str): string containing the reaction formula

    Returns:
        left_compounds (list): compounds at the left of the reaction formula
        right_compounds (list): compounds at the right of the reaction formula
    """
    equation_pattern = r'(?P<stoechiometry>\d|\w|\([\w\d\+]*\))*\ *(?P<compound>[CG]\d{5})|(?P<symbol><*=>*)'
    left_compounds = []
    right_compounds = []
    equation_symbol_passed = False
    for m in re.finditer(equation_pattern, equation_text):
        stoechiometry = m.groupdict()['stoechiometry']
        if stoechiometry is None:
            stoechiometry = 1
        else:
            try:
                stoechiometry = int(stoechiometry)
            except:
                print('Stochiometry is not an int for {0} ignore it.'.format(reaction_id))
                return None, None
        compound = m.groupdict()['compound']
        symbol = m.groupdict()['symbol']
        if symbol is not None:
            equation_symbol_passed = True
        if equation_symbol_passed is False and compound is not None:
            left_compounds.append((compound, stoechiometry))
        elif equation_symbol_passed is True and compound is not None:
            right_compounds.append((compound, stoechiometry))

    return left_compounds, right_compounds


def get_compound_names(compound_file):
    """Using bioservices.KEGG find the name associated to compound and glycan.

    Args:
        compound_file (str): output file which will contains compound ID and name
    """
    response_text = kegg.list('compound')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    compounds = {}
    for line in csvreader:
        cpd_id = line[0].replace('cpd:', '')
        cpd_name = line[1].split('; ')[0]
        compounds[cpd_id] = cpd_name

    response_text = kegg.list('glycan')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []

    for line in csvreader:
        cpd_id = line[0].replace('gl:', '')
        cpd_name = line[1].split('; ')[0]
        compounds[cpd_id] = cpd_name

    with open(compound_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['compound_id', 'compound_name'])
        for cpd_id in compounds:
            csvwriter.writerow([cpd_id, compounds[cpd_id]])


def create_sbml_model_from_kegg_file(reaction_folder, compound_file, output_sbml, output_tsv):
    """Using the reaction keg files (from retrieve_reactions), the compound file (from get_compound_names),
    create a SBML file (containing all reactions of KEGG) and a tsv file (used to map EC and/or KO to reaction ID)

    Args:
        reaction_folder (str): path to the folder containing keg reaction files
        compound_file (str): path to the tsv file containing compound ID and name
        output_sbml (str): path to the sbml output file
        output_tsv (str): path to an output tsv file mapping reaction ID, with KO and EC
    """
    compounds = {}
    with open(compound_file, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compounds[line[0]] = line[1]

    reaction_ecs = {}
    model = Model('KEGG')
    sbml_reactions = []
    reactions = {}
    for reaction_file in os.listdir(reaction_folder):
        reaction_file_path = os.path.join(reaction_folder, reaction_file)
        with open(reaction_file_path) as open_reaction_file_path:
            reaction_data = kegg.parse(open_reaction_file_path.read())

        reaction_id = reaction_data['ENTRY'].split(' ')[0]
        if 'NAME' in reaction_data:
            reaction_name = reaction_data['NAME'][0]
        else:
            reaction_name = None
        left_compounds, right_compounds = extract_reaction(reaction_id, reaction_data['EQUATION'])
        if 'ORTHOLOGY' in reaction_data:
            kegg_orthologs = reaction_data['ORTHOLOGY'].keys()
        else:
            kegg_orthologs = []
        if 'ENZYME' in reaction_data:
            kegg_ecs = reaction_data['ENZYME']
        else:
            kegg_ecs = []
        reaction_ecs[reaction_id] = [kegg_orthologs, kegg_ecs]
        if left_compounds is None:
            continue
        reactions[reaction_id] = {}
        for stochiometry_metabolite in left_compounds:
            metabolite_id = stochiometry_metabolite[0]
            metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=compounds[metabolite_id])
            reactions[reaction_id][metabolite_id_sbml] = - stochiometry_metabolite[1]
        for stochiometry_metabolite in right_compounds:
            metabolite_id = stochiometry_metabolite[0]
            metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=compounds[metabolite_id])
            reactions[reaction_id][metabolite_id_sbml] = stochiometry_metabolite[1]

        reaction = Reaction(reaction_id, upper_bound=1000, lower_bound=-1000)
        if reaction_name:
            reaction.name = reaction_name
        if kegg_orthologs != []:
            reaction.gene_reaction_rule = '( ' + ' or '.join([ko for ko in kegg_orthologs]) + ' )'
        reaction.add_metabolites(reactions[reaction_id])

        sbml_reactions.append(reaction)

    model.add_reactions(sbml_reactions)
    # Create sbml file.
    write_sbml_model(model, output_sbml)

    with open(output_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['reaction_id', 'kegg_orthologs', 'kegg_ec'])
        for reaction_id in reaction_ecs:
            csvwriter.writerow([reaction_id, ','.join(reaction_ecs[reaction_id][0]), ','.join(reaction_ecs[reaction_id][1])])


def map_protein_to_KO_id(proteins_ids):
    """From a list of UniProt protein IDs, retrieved the corresponding KEGG Orthologs.

    Args:
        proteins_ids (list): List of Uniprot protein IDs

    Returns:
        protein_ko_mapping (dict): mapping between UniProt protein ID as key and KO ID as value
    """
    protein_mapping = {}
    beta = None
    protein_mapping = query_uniprot_kegg_rest(proteins_ids, beta, protein_mapping)
    kegg_genes = [gene for prot_id in protein_mapping for gene in protein_mapping[prot_id]]

    chunk_kegg_genes = chunks(kegg_genes, 99)

    mapping_gene_kegg_ko = {}
    for chunk_kegg_gene in chunk_kegg_genes:
        str_chunk_kegg_gene = '+'.join(chunk_kegg_gene)
        response_text = kegg.link('ko', str_chunk_kegg_gene)
        if response_text != '':
            csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
        else:
            csvreader = []
        for line in csvreader:
            kegg_gene = line[0]
            ko = line[1]
            mapping_gene_kegg_ko[kegg_gene] = ko

    protein_ko_mapping = {}
    for protein in protein_mapping:
        kegg_genes = protein_mapping[protein]
        for kegg_gene in kegg_genes:
            if kegg_gene in mapping_gene_kegg_ko:
                ko = mapping_gene_kegg_ko[kegg_gene]
                protein_ko_mapping[protein] = ko

    return protein_ko_mapping

def retrieve_mapping_dictonaries(kegg_rxn_mapping_path):
    """From the KEGG tsv mapping file (creating at create_sbml_model_from_kegg_file)
    retrieve 2 mapping dictionaries, one KO ID to kegg_reaction and another
    EC_number to kegg_reaction

    Args:
        kegg_rxn_mapping_path (str): path to the KEGG tsv file for mapping reaction and EC/KO

    Returns:
        ko_to_reactions (dict): mapping between KO ID as key and kegg reaction as value
        ec_to_reactions (dict): mapping between EC number as key and kegg reaction as value
    """
    ko_to_reactions = {}
    ec_to_reactions = {}
    with open(kegg_rxn_mapping_path, 'r') as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            reaction_id = line[0]
            if line[1] != '':
                ko_ids = line[1].split(',')
            else:
                ko_ids = []
            if line[2] != '':
                ec_ids = line[2].split(',')
            else:
                ec_ids = []
            for ko_id in ko_ids:
                if ko_id not in ko_to_reactions:
                    ko_to_reactions[ko_id] = [reaction_id]
                else:
                    ko_to_reactions[ko_id].append(reaction_id)

            for ec_id in ec_ids:
                if ec_id not in ec_to_reactions:
                    ec_to_reactions[ec_id] = [reaction_id]
                else:
                    ec_to_reactions[ec_id].append(reaction_id)

    return ko_to_reactions, ec_to_reactions


def create_draft_networks(input_folder, output_folder, mapping_ko=False, beta=None):
    """From the output folder of 'esmecata annotation' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        input_folder (str): path to the output folder of esmecata annotation
        output_folder (str): path to the output folder
        mapping_ko (bool): option to use KO ID to retrieve reactions
        beta (bool): option to use the new API of UniProt (in beta can be unstable)
    """
    starttime = time.time()

    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(output_folder, 'kegg_model')
    is_valid_dir(kegg_model_path)

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')

    if not os.path.exists(kegg_sbml_model_path) and not os.path.exists(kegg_rxn_mapping_path):
        if not os.path.exists(kegg_reactions_folder_path):
            logger.info('|EsMeCaTa|kegg| Retrieve reactions from KEGG.')
            retrieve_reactions(kegg_reactions_folder_path)
        if not os.path.exists(compound_file_path):
            logger.info('|EsMeCaTa|kegg| Retrieve compound IDs and names from KEGG.')
            get_compound_names(compound_file_path)
        logger.info('|EsMeCaTa|kegg| Create KEGG reference SBML and mapping tsv file.')
        create_sbml_model_from_kegg_file(kegg_reactions_folder_path, compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path)

    # Create SBML output folder.
    sbml_output_folder_path = os.path.join(output_folder, 'sbml')
    is_valid_dir(sbml_output_folder_path)

    ko_to_reactions, ec_to_reactions = retrieve_mapping_dictonaries(kegg_rxn_mapping_path)

    # Retrieve proteins and annotations from esmecata annotation folder.
    #annotation_folder_path = os.path.join(input_folder, 'annotation')
    annotation_reference_folder_path = os.path.join(input_folder, 'annotation_reference')

    for annot_file in os.listdir(annotation_reference_folder_path):
        annot_file_path = os.path.join(annotation_reference_folder_path, annot_file)

        base_file = os.path.basename(annot_file_path)
        base_filename = os.path.splitext(base_file)[0]

        # Extract protein IDs and EC number from anntotation reference folder.
        protein_ec_numbers = {}
        with open(annot_file_path, 'r') as open_annot_file_path:
            csvreader = csv.reader(open_annot_file_path, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                protein_id = line[0]
                ec_numbers = line[4].split(',')
                protein_ec_numbers[protein_id] = ec_numbers

        taxon_reactions = {}
        # If mapping KO option is used, search for Ko terms associated to proteins.
        if mapping_ko is True:
            reference_protein = protein_ec_numbers.keys()
            protein_ko_mapping = map_protein_to_KO_id(reference_protein) 

            for protein in protein_ko_mapping:
                ko_id = protein_ko_mapping[protein].replace('ko:', '')
                if ko_id in ko_to_reactions:
                    reaction_ids = ko_to_reactions[ko_id]
                    for reaction_id in reaction_ids:
                        if reaction_id not in taxon_reactions:
                            taxon_reactions[reaction_id] = [protein]
                        else:
                            taxon_reactions[reaction_id].append(protein)

        # Use EC found to be associated to reference protein to retrieve KEEG reaction.
        for protein in protein_ec_numbers:
            ec_ids = protein_ec_numbers[protein]
            for ec_id in ec_ids:
                if ec_id in ec_to_reactions:
                    reaction_ids = ec_to_reactions[ec_id]
                    for reaction_id in reaction_ids:
                        if reaction_id not in taxon_reactions:
                            taxon_reactions[reaction_id] = [protein]
                        else:
                            taxon_reactions[reaction_id].append(protein)

        kegg_model = read_sbml_model(kegg_sbml_model_path)

        species_model = Model(base_filename)

        # Map KEGG reaction from KEGG SBML model to taxon SBML.
        sbml_reactions = []
        for reaction in kegg_model.reactions:
            reaction_id = reaction.id.replace('R_','')
            if reaction_id in taxon_reactions:
                reaction.gene_reaction_rule = '( ' + ' or '.join(taxon_reactions[reaction_id]) + ' )'
                sbml_reactions.append(reaction)

        species_model.add_reactions(sbml_reactions)

        # Create SBML file.
        sbml_output_file_path = os.path.join(sbml_output_folder_path, base_filename+'.sbml')
        write_sbml_model(species_model, sbml_output_file_path)

    endtime = time.time()
    duration = endtime - starttime
    logger.info('|EsMeCaTa|kegg| Draft networks creation complete.')
