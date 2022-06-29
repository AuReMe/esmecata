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
import zipfile

from bioservices import __version__ as bioservices_version
from bioservices import KEGG, UniProt
from cobra import __version__ as cobra_version
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model, read_sbml_model

from esmecata.utils import is_valid_dir, urllib_query, chunks, get_rest_uniprot_release
from esmecata import __version__ as esmecata_version

URLLIB_HEADERS = {'User-Agent': 'EsMeCaTa annotation v' + esmecata_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

# Create KEGG instance of bioservices.KEEG.
kegg = KEGG()
uniprot = UniProt(verbose=False)

root = os.path.dirname(__file__)
kegg_archive = os.path.join(*[root, 'data', 'kegg_model.zip'])

def rest_query_uniprot_to_retrieve_kegg_gene(protein_queries):
    """REST query to get KEGG gene ID from proteins.

    Args:
        protein_queries (list): list of proteins

    Returns:
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value
    """
    output_dict = {}

    url = 'http://rest.uniprot.org/idmapping/run'

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


def query_uniprot_kegg_rest(protein_to_search_on_uniprots, output_dict):
    """Query UniProt with REST to find KEGG gene ID

    Args:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value

    Returns:
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value
    """
    # The limit of 15 000 proteins per query comes from the help of Uniprot (inferior to 20 000):
    # https://www.uniprot.org/help/uploadlists
    if len(protein_to_search_on_uniprots) < 15000:
        protein_queries = ','.join(protein_to_search_on_uniprots)
        tmp_output_dict = rest_query_uniprot_to_retrieve_kegg_gene(protein_queries)
        output_dict.update(tmp_output_dict)
        time.sleep(1)
    else:
        protein_chunks = chunks(list(protein_to_search_on_uniprots), 15000)
        for chunk in protein_chunks:
            protein_queries = ','.join(chunk)
            tmp_output_dict = rest_query_uniprot_to_retrieve_kegg_gene(protein_queries)
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
                logger.critical('Stochiometry is not an int for {0} ignore it.'.format(reaction_id))
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


def get_modules(module_file):
    """Using bioservices.KEGG to find the module associated with reacitons.

    Args:
        module_file (str): output file which will contains module ID, name, formula and reactions
    """
    response_text = kegg.link('module', 'reaction')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    module_reactions = {}
    for line in csvreader:
        reaction_id = line[0].split(':')[1]
        module_id = line[1].split(':')[1]
        if module_id not in module_reactions:
            module_reactions[module_id] = [reaction_id]
        else:
            module_reactions[module_id].append(reaction_id)

    response_text = kegg.list('module')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []

    modules = {}
    for line in csvreader:
        module_id = line[0].split(':')[1]
        module_data = line[1].split(', ')
        if len(module_data) > 2:
            module_name = module_data[0]
            module_formula = module_data[-1]
        else:
            module_name = module_data[0]
            module_formula = ''
        modules[module_id] = (module_name, module_formula)

    all_modules = set(list(modules.keys()) + list(module_reactions.keys()))
    with open(module_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['module_id', 'module_name', 'module_formula', 'module_reactions'])
        for module_id in all_modules:
            if module_id in module_reactions:
                module_reaction = ','.join(module_reactions[module_id])
            else:
                module_reaction = ''
            if module_id in modules:
                module_name = modules[module_id][0]
                module_formula = modules[module_id][1]
            else:
                module_name = ''
                module_formula = ''
            csvwriter.writerow([module_id, module_name, module_formula, module_reaction])


def create_sbml_model_from_kegg_file(reaction_folder, compound_file, output_sbml, output_tsv, pathways_tsv):
    """Using the reaction keg files (from retrieve_reactions), the compound file (from get_compound_names),
    create a SBML file (containing all reactions of KEGG) and a tsv file (used to map EC and/or KO to reaction ID)

    Args:
        reaction_folder (str): path to the folder containing keg reaction files
        compound_file (str): path to the tsv file containing compound ID and name
        output_sbml (str): path to the sbml output file
        output_tsv (str): path to an output tsv file mapping reaction ID, with KO and EC
        pathways_tsv (str): path to an output tsv showing the pathway, pathway ID and the associated reactions
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
    pathways = {}
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
        if 'PATHWAY' in reaction_data:
            kegg_pathways = reaction_data['PATHWAY']
            for pathway in kegg_pathways:
                if pathway not in pathways:
                    pathways[pathway] = {}
                    pathways[pathway]['name'] = kegg_pathways[pathway]
                    pathways[pathway]['reaction'] = [reaction_id]
                else:
                    pathways[pathway]['reaction'].append(reaction_id)
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

    with open(pathways_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_reactions'])
        for pathway_id in pathways:
            pathway_name = pathways[pathway_id]['name']
            pathway_reactions = ','.join(set(pathways[pathway_id]['reaction']))
            csvwriter.writerow([pathway_id, pathway_name, pathway_reactions])


def map_protein_to_KO_id(proteins_ids):
    """From a list of UniProt protein IDs, retrieved the corresponding KEGG Orthologs.

    Args:
        proteins_ids (list): List of Uniprot protein IDs

    Returns:
        protein_ko_mapping (dict): mapping between UniProt protein ID as key and KO ID as value
    """
    protein_mapping = {}
    protein_mapping = query_uniprot_kegg_rest(proteins_ids, protein_mapping)
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


def compute_stat_kegg(sbml_folder, stat_file=None):
    """Compute stat associated with the number of reactions and metabolites for each taxonomic affiliations.

    Args:
        sbml_folder (str): pathname to the sbml folder containing draft metabolic networks for each cluster
        stat_file (str): pathname to the tsv stat file

    Returns:
        kegg_numbers (dict): dict containing observation names (as key) associated with reaction and metabolites number (as value)
    """
    kegg_numbers = {}
    for infile in os.listdir(sbml_folder):
        if '.sbml' in infile:
            sbml_input_file_path = os.path.join(sbml_folder, infile)
            kegg_model = read_sbml_model(sbml_input_file_path)
            infile_reactions = kegg_model.reactions
            infile_metabolites = kegg_model.metabolites
            kegg_numbers[infile.replace('.sbml','')] = (len(infile_reactions), len(infile_metabolites))

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_reactions', 'Number_metabolites'])
            for observation_name in kegg_numbers:
                csvwriter.writerow([observation_name, kegg_numbers[observation_name][0], kegg_numbers[observation_name][1]])

    return kegg_numbers


def get_kegg_database_version():
    """Retrieve the KEGG version release

    Returns:
        kegg_version (str): string containing KEGG release version
    """
    response_text = kegg.dbinfo()

    for line in response_text.splitlines():
        if 'Release ' in line:
            kegg_version = line.split('Release ')[1].strip()

    return kegg_version


def create_draft_networks(input_folder, output_folder, mapping_ko=False, recreate_kegg=None):
    """From the output folder of 'esmecata annotation' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        input_folder (str): path to the output folder of esmecata annotation
        output_folder (str): path to the output folder
        mapping_ko (bool): option to use KO ID to retrieve reactions
        recreate_kegg (bool): option to recreate KEGG model by guerying KEGG server
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|kegg| Begin KEGG metabolism mapping.')

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder, 'mapping_ko': mapping_ko}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['cobra'] = cobra_version

    esmecata_kegg_metadata = get_rest_uniprot_release(options)
    esmecata_kegg_metadata['kegg_release_number'] = get_kegg_database_version()

    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(output_folder, 'kegg_model')
    is_valid_dir(kegg_model_path)

    if not recreate_kegg:
        with zipfile.ZipFile(kegg_archive, 'r') as zip_ref:
            zip_ref.extractall(output_folder)

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')

    missing_files = [not os.path.exists(kegg_sbml_model_path), not os.path.exists(kegg_rxn_mapping_path), not os.path.exists(kegg_pathways_path)]
    if any(missing_files):
        if not os.path.exists(kegg_reactions_folder_path):
            logger.info('|EsMeCaTa|kegg| Retrieve reactions from KEGG to create SMBL model.')
            retrieve_reactions(kegg_reactions_folder_path)
        if not os.path.exists(compound_file_path):
            logger.info('|EsMeCaTa|kegg| Retrieve compound IDs and names from KEGG to create SMBL model.')
            get_compound_names(compound_file_path)
        logger.info('|EsMeCaTa|kegg| Create KEGG reference SBML and mapping tsv file.')
        create_sbml_model_from_kegg_file(kegg_reactions_folder_path, compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path, kegg_pathways_path)

    kegg_pathways = {}
    with open(kegg_pathways_path, 'r') as open_kegg_pathways_path:
        csvreader = csv.reader(open_kegg_pathways_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            pathway_id = line[0]
            pathway_name = line[1]
            pathway_reactions = line[2].split(',')
            kegg_pathways[pathway_id] = (pathway_name, pathway_reactions)

    # Create SBML output folder.
    sbml_output_folder_path = os.path.join(output_folder, 'sbml')
    is_valid_dir(sbml_output_folder_path)

    # Create KO output folder.
    ko_output_folder_path = os.path.join(output_folder, 'ko')
    is_valid_dir(ko_output_folder_path)

    # Create pathways annotated output folder.
    pathways_output_folder_path = os.path.join(output_folder, 'pathways')
    is_valid_dir(pathways_output_folder_path)

    ko_to_reactions, ec_to_reactions = retrieve_mapping_dictonaries(kegg_rxn_mapping_path)

    # Retrieve proteins and annotations from esmecata annotation folder.
    #annotation_folder_path = os.path.join(input_folder, 'annotation')
    annotation_reference_folder_path = os.path.join(input_folder, 'annotation_reference')

    for clust_threshold in os.listdir(annotation_reference_folder_path):
        clust_annotation_reference_folder_path = os.path.join(annotation_reference_folder_path, clust_threshold)

        if mapping_ko is True:
            clust_ko_output_folder_path = os.path.join(ko_output_folder_path, clust_threshold)
            is_valid_dir(clust_ko_output_folder_path)

        clust_pathways_output_folder_path = os.path.join(pathways_output_folder_path, clust_threshold)
        is_valid_dir(clust_pathways_output_folder_path)

        clust_sbml_output_folder_path = os.path.join(sbml_output_folder_path, clust_threshold)
        is_valid_dir(clust_sbml_output_folder_path)
        for annot_file in os.listdir(clust_annotation_reference_folder_path):
            annot_file_path = os.path.join(clust_annotation_reference_folder_path, annot_file)

            base_file = os.path.basename(annot_file_path)
            base_filename = os.path.splitext(base_file)[0]

            # Extract protein IDs and EC number from anntotation reference folder.
            protein_ec_numbers = {}
            protein_clusters = {}
            with open(annot_file_path, 'r') as open_annot_file_path:
                csvreader = csv.reader(open_annot_file_path, delimiter='\t')
                next(csvreader)
                for line in csvreader:
                    protein_id = line[0]
                    protein_cluster = line[1].split(',')
                    ec_numbers = line[5].split(',')
                    protein_ec_numbers[protein_id] = ec_numbers
                    protein_clusters[protein_id] = protein_cluster

            taxon_reactions = {}
            # If mapping KO option is used, search for KO terms associated to proteins.
            if mapping_ko is True:
                protein_to_maps = set([protein_id for protein_cluster in protein_clusters for protein_id in protein_clusters[protein_cluster]])
                protein_ko_mapping = map_protein_to_KO_id(protein_to_maps)

                ko_output_file_path = os.path.join(clust_ko_output_folder_path, base_filename+'.tsv')
                with open(ko_output_file_path, 'w') as output_file:
                    csvwriter = csv.writer(output_file, delimiter='\t')
                    csvwriter.writerow(['protein', 'KO'])
                    for protein_cluster in protein_clusters:
                        for protein_id in protein_clusters[protein_cluster]:
                            if protein_id in protein_ko_mapping:
                                ko = protein_ko_mapping[protein_id]
                            else:
                                ko = ''
                            csvwriter.writerow([protein_id, ko])

                # Keep KO and propagate their reaction only if they appear in all protein of the cluster.
                ko_added_reactions = []
                all_kos = []
                for protein_cluster in protein_clusters:
                    ko_ids = [[protein_ko_mapping[protein_id]] if protein_id in protein_ko_mapping else [] for protein_id in protein_clusters[protein_cluster]]
                    protein_all_ko_ids = set([ko_id for subko_ids in ko_ids for ko_id in subko_ids])
                    all_kos.extend(list(protein_all_ko_ids))
                    keep_kos = [ko_id for ko_id in protein_all_ko_ids if sum(subko_ids.count(ko_id) for subko_ids in ko_ids) >= 1 * len(protein_clusters[protein_cluster])]
                    for ko_id in keep_kos:
                        ko_id = ko_id.replace('ko:', '')
                        if ko_id in ko_to_reactions:
                            reaction_ids = ko_to_reactions[ko_id]
                            for reaction_id in reaction_ids:
                                ko_added_reactions.append(reaction_id)
                                if reaction_id not in taxon_reactions:
                                    taxon_reactions[reaction_id] = [protein_cluster]
                                else:
                                    taxon_reactions[reaction_id].append(protein_cluster)
                logger.info('|EsMeCaTa|kegg| Added {0} reactions from {1} KO for taxon {2} (with clustering threshold {3}).'.format(len(set(ko_added_reactions)), len(set(all_kos)), base_filename, clust_threshold))

            # Use EC found to be associated to reference protein to retrieve KEEG reaction.
            ec_added_reactions = []
            all_ecs = []
            for protein in protein_ec_numbers:
                ec_ids = protein_ec_numbers[protein]
                all_ecs.extend(ec_ids)
                for ec_id in ec_ids:
                    if ec_id in ec_to_reactions:
                        reaction_ids = ec_to_reactions[ec_id]
                        for reaction_id in reaction_ids:
                            ec_added_reactions.append(reaction_id)
                            if reaction_id not in taxon_reactions:
                                taxon_reactions[reaction_id] = [protein]
                            else:
                                taxon_reactions[reaction_id].append(protein)
            logger.info('|EsMeCaTa|kegg| Added {0} reactions from {1} EC for taxon {2} (with clustering threshold {3}).'.format(len(set(ec_added_reactions)), len(set(all_ecs)), base_filename, clust_threshold))

            total_added_reactions = list(taxon_reactions.keys())
            if mapping_ko is True:
                logger.info('|EsMeCaTa|kegg| A total of {0} unique reactions are added from EC and KO for taxon {1} (with clustering threshold {2}).'.format(len(total_added_reactions), base_filename, clust_threshold))

            # Create pathway file contening pathway with reacitons in the taxon.
            pathways_output_file_path = os.path.join(clust_pathways_output_folder_path, base_filename+'.tsv')
            with open(pathways_output_file_path, 'w') as open_pathways_output_file_path:
                csvwriter = csv.writer(open_pathways_output_file_path, delimiter='\t')
                csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_compeltion_ratio', 'pathway_reaction_in_taxon', 'pathway_reaction'])
                for pathway in kegg_pathways:
                    pathway_reactions = kegg_pathways[pathway][1]
                    pathway_reaction_in_taxon = set(pathway_reactions).intersection(set(total_added_reactions))
                    if len(pathway_reaction_in_taxon) > 0:
                        pathway_name = kegg_pathways[pathway][0]
                        pathway_completion_ratio = len(pathway_reaction_in_taxon) / len(pathway_reactions)
                        csvwriter.writerow([pathway, pathway_name, pathway_completion_ratio, ','.join(pathway_reaction_in_taxon), ','.join(pathway_reactions)])

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

            # Create file if there si at least 1 reaction.
            if len(species_model.reactions) > 0:
                # Create SBML file.
                sbml_output_file_path = os.path.join(clust_sbml_output_folder_path, base_filename+'.sbml')
                write_sbml_model(species_model, sbml_output_file_path)
            else:
                logger.info('|EsMeCaTa|kegg| No reactions in model for {0}, no SBML file will be created.'.format(base_filename))

    for clust_threshold in os.listdir(sbml_output_folder_path):
        clust_sbml_output_folder_path = os.path.join(sbml_output_folder_path, clust_threshold)
        clust_stat_file = os.path.join(output_folder, 'stat_number_kegg_{0}.tsv'.format(clust_threshold))
        compute_stat_kegg(clust_sbml_output_folder_path, clust_stat_file)

    endtime = time.time()

    duration = endtime - starttime
    esmecata_kegg_metadata['esmecata_kegg_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'esmecata_metadata_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(esmecata_kegg_metadata, ouput_file, indent=4)
    logger.info('|EsMeCaTa|kegg| Draft networks creation complete.')
