# Copyright (C) 2023-2025 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
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
import os
import sys
import pandas as pd
import numpy as np
import datetime
import time
import json
import pronto
import logging
import subprocess
import urllib.request
import gseapy

from bioservices import KEGG
from Bio.ExPASy import Enzyme
from matplotlib import __version__ as matplotlib_version
from esmecata.core.proteomes import get_taxon_obs_name
from esmecata.utils import get_domain_or_superkingdom_from_ncbi_tax_database
from esmecata import __version__ as esmecata_version
from Bio import __version__ as biopython_version

logger = logging.getLogger(__name__)


def extract_annotation(annotation_reference_folder):
    """ From esmecata annotation reference folder, creates a dict indicating the file associated with each annotation.

    Args:
        annotation_reference_folder (str): path to annotation reference folder
        dataset_annotation_file_path (str): path to output dataset annotation file

    Returns:
        annotation_sets (dict): annotation dict: annotation as key and taxon name as value
    """
    annotation_sets = {}
    for annotation_file in os.listdir(annotation_reference_folder):
        annotation_file_name = os.path.splitext(annotation_file)[0]
        annotation_file_path = os.path.join(annotation_reference_folder, annotation_file)

        with open(annotation_file_path, 'r') as open_annotation_input_file_path:
            csvreader = csv.DictReader(open_annotation_input_file_path, delimiter='\t')

            for line in csvreader:
                gos = line['GO'].split(',')
                ecs = line['EC'].split(',')

                if gos != ['']:
                    for go in gos:
                        if go not in annotation_sets:
                            annotation_sets[go] = [annotation_file_name]
                        else:
                            annotation_sets[go].append(annotation_file_name)
                if ecs != ['']:
                    for ec in ecs:
                        if ec not in annotation_sets:
                            annotation_sets[ec] = [annotation_file_name]
                        else:
                            annotation_sets[ec].append(annotation_file_name)

    return annotation_sets

def extract_observation_name_annotation(annotation_reference_folder):
    """ From esmecata annotation reference folder, creates a dict indicating the file associated with each annotation.

    Args:
        annotation_reference_folder (str): path to annotation reference folder
        dataset_annotation_file_path (str): path to output dataset annotation file

    Returns:
        taxa_annotations_sets (dict): taxon name as key and annotation as value
    """
    taxa_annotations_sets = {}
    for annotation_file in os.listdir(annotation_reference_folder):
        annotation_file_name = os.path.splitext(annotation_file)[0]
        annotation_file_path = os.path.join(annotation_reference_folder, annotation_file)

        with open(annotation_file_path, 'r') as open_annotation_input_file_path:
            csvreader = csv.DictReader(open_annotation_input_file_path, delimiter='\t')

            for line in csvreader:
                gos = line['GO'].split(',')
                ecs = line['EC'].split(',')

                if gos != ['']:
                    for go in gos:
                        if annotation_file_name not in taxa_annotations_sets:
                            taxa_annotations_sets[annotation_file_name] = [go]
                        else:
                            taxa_annotations_sets[annotation_file_name].append(go)
                if ecs != ['']:
                    for ec in ecs:
                        if annotation_file_name not in taxa_annotations_sets:
                            taxa_annotations_sets[annotation_file_name] = [ec]
                        else:
                            taxa_annotations_sets[annotation_file_name].append(ec)

    return taxa_annotations_sets


def read_function_table(function_table_path, transpose_table=True):
    """ Read function table and create dictionary compatible with enrichment analysis.
    The table contains function/annotation as columns and organism as rows.

    Args:
        function_table_path (str): path to function table file.
        transpose_table (bool): boolean

    Returns:
        annotation_sets (dict): if transpose_table is True: annotation as key and list of organism possessing them as values, if transpose_table is False: organism name as key and list of functions possessing them as values,
    """
    df_function_table = pd.read_csv(function_table_path, sep='\t|,', index_col=0, engine='python')

    if transpose_table is True:
        df_function_table = df_function_table.T

    annotation_sets = {}
    for index, row in df_function_table.iterrows():
        annotation_sets[index] = [col for col in df_function_table.columns
                    if row[col]>0
                    for index_iter in range(int(row[col]))]

    return annotation_sets


def get_annot_name(enzyme_datfile, gobasic_file):
    """ Extract annotation names from enzyme.dat and go-baisc.obo.

    Args:
        enzyme_datfile (str): path to enzyme.Dat file
        gobasic_file (str): path to go-baisc.obo file

    Returns:
        enzyme_names (dict): EC as key and EC name as value
        go_names (dict): GO term ID as key and GO term name as value
    """
    enzyme_names = {}
    with open(enzyme_datfile) as infile:
        enzymes = Enzyme.parse(infile)
        for enzyme in enzymes:
            enzyme_id = enzyme['ID']
            enzyme_name = enzyme['DE']
            enzyme_names[enzyme_id] = enzyme_name

    go_ontology = pronto.Ontology(gobasic_file)
    go_names = {}
    for go_term in go_ontology:
        if 'GO:' in go_term:
            go = go_ontology[go_term]
            go_names[go_term] = go.name

    return enzyme_names, go_names


def get_ko_names(ko_file):
    """ Extract KEGG Orthologs names from KEGG database.

    Args:
        ko_file (str): path to KEGG orthologs file

    Returns:
        ko_names (dict): KO as key and KO name as value
    """
    if os.path.exists(ko_file):
        ko_df = pd.read_csv(ko_file, sep='\t')
        ko_df.set_index('KO', inplace=True)
        ko_names = ko_df['name'].to_dict()
    else:
        k = KEGG()
        ko_names = {}
        ko_answer = k.list('ko')
        for ko_line in ko_answer.split('\n'):
            ko_dat = ko_line.split('\t')
            ko_entry = 'ko:'+ko_dat[0]
            if ko_entry.startswith('ko:K'):
                if ';' in ko_dat[1]:
                    ko_name = ko_dat[1].split(';')[1]
                else:
                    ko_name = ko_dat[1]
                ko_names[ko_entry] = ko_name
        ko_df = pd.DataFrame(ko_names.items(), columns=['KO', 'name'])
        ko_df.to_csv(ko_file, sep='\t', index=None)

    return ko_names


def run_orsum(annotation_file_gmt_file, enrichr_module_phylum_output, orsum_output_folder, orsum_minterm_size=None):
    """ Run orsum to reduce list of enriched annotations from gseapy.

    Args:
        annotation_file_gmt_file (str): GMT file containing annotations and their labels
        enrichr_module_phylum_output (str): path to output of gseapy containing enriched annotations for each taxon
        orsum_output_folder (str): path to output folder for orsum
        orsum_minterm_size (int): option minTermSize of orsum
    """
    input_files = [os.path.join(enrichr_module_phylum_output, file) for file in os.listdir(enrichr_module_phylum_output)]
    orsum_cmd = ['orsum.py', '--gmt', annotation_file_gmt_file, '--files', *input_files, '--outputFolder', orsum_output_folder]

    if orsum_minterm_size is not None:
        orsum_cmd.append('--minTermSize')
        orsum_cmd.append(str(orsum_minterm_size))

    logger.info('|EsMeCaTa|gseapy_enrichr| orsum command: {0}'.format(' '.join(orsum_cmd)))
    subprocess.call(orsum_cmd)


def get_orsum_version():
    """Returns version of orsum.

    Returns:
        orsum_version (str): version of orsum.
    """
    orsum_version = None
    response = subprocess.Popen(['orsum.py', '-v'], stdout=subprocess.PIPE, start_new_session=True, universal_newlines="")
    for orsum_line in response.stdout:
        orsum_line_decoded = orsum_line.decode('utf-8')
        orsum_version = orsum_line_decoded.strip('\n')

    if orsum_version is None:
        logger.critical('|EsMeCaTa|gseapy_enrichr| esmecata could not find the version of orsum.')
        logger.critical('|EsMeCaTa|gseapy_enrichr| It is possibly an issue with the installation of orsum (maybe it is not in the PATH). Or it can be due to a change in the output of orsum.py -v command.')
        sys.exit()

    return orsum_version


def extract_organisms_selected(taxa_lists_file):
    """ Extract observation names from input file (each column corresponds to a set of observation name to search for enrichment).

    Args:
        taxa_lists_file (str): path to taxa list file

    Returns:
        taxa_lists (dict): dictionary containing as key the column name of the file and as value the list of observation names
    """
    taxa_lists = {}
    df = pd.read_csv(taxa_lists_file, sep='\t')
    df = df.replace(np.nan, '')
    for col in df.columns:
        taxa_lists[col] = [observation_name for observation_name in df[col].tolist() if observation_name != '']
    return taxa_lists


def taxon_rank_annotation_enrichment(annotation_folder_or_file, output_folder, grouping="tax_rank",
                                     taxon_rank='phylum', taxa_lists_file=None, function_lists_file=None,
                                     annot_names_file='download', orsum_minterm_size=None, selected_adjust_pvalue_cutoff=0.05,
                                     taxon_id_file=None):
    """ Run an enrichment analysis on taxon from annotation results of esmecata using gseapy.
    Then filter this list with orsum.

    Args:
        annotation_folder_or_file (str): path to esmecata annotation folder or function table (from esmecata or other methods)
        output_folder (str): path to output folder
        grouping (str): grouping factor, either "tax_rank", "selected" or "selected_function"
        taxon_rank (str): taxon rank to cluster the observation name together (by default, phylum) when selecting grouping "tax_rank"
        taxa_lists_file (str): path to manually created groups of observation names when selecting grouping "selected"
        function_lists_file (str): path to manually created groups of functions (GO or EC) when selecting grouping "selected_function"
        annot_names_file (str): path to json annot names file. If 'download' is given, download EC number, GO Terms and KO names and generate a json file
        orsum_minterm_size (int): option minTermSize of orsum
        selected_adjust_pvalue_cutoff (float): adjust-Pval cutoff for gseapy enrichr, default: 0.05
        taxon_id_file (str): path to taxon ID file if annotation_folder_or_file is a file
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|gseapy_enrichr| Begin enrichment analysis.')

    if grouping == "tax_rank":
        if os.path.isfile(annotation_folder_or_file) and taxon_id_file is None:
            logger.critical('|EsMeCaTa|gseapy_enrichr| If you give an input annotation file for tax_rank analysis, you have to give a taxon_id_file indicating taxon ID associated with organisms.')
            sys.exit()
        domain_superkingdom_tax_rank_name = get_domain_or_superkingdom_from_ncbi_tax_database()
        if taxon_rank is None:
            logger.critical('|EsMeCaTa|gseapy_enrichr| You have to specify a taxon rank for this analysis with --taxon-rank')
            sys.exit()
        taxon_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'domain', 'superkingdom']
        if taxon_rank != 'phylum':
            if taxon_rank not in taxon_ranks:
                logger.critical('|EsMeCaTa|gseapy_enrichr| Incorrect taxon given {0}, possible ranks are: {1}'.format(taxon_rank, ','.join(taxon_ranks)))
                sys.exit()
        if taxon_rank == 'superkingdom':
            if taxon_rank != domain_superkingdom_tax_rank_name and domain_superkingdom_tax_rank_name == 'domain':
                logger.critical('|EsMeCaTa|gseapy_enrichr| superkingdom has been renamed domain in new version of NCBI Taxonomy database, EsMeCaTa will use domain.')
                taxon_rank = domain_superkingdom_tax_rank_name
        if taxon_rank == 'domain':
            if taxon_rank != domain_superkingdom_tax_rank_name and domain_superkingdom_tax_rank_name == 'superkingdom':
                logger.critical('|EsMeCaTa|gseapy_enrichr| domain was named superkingdom in old version of NCBI Taxonomy database (which is the case of your version), EsMeCaTa will use superkingdom.')
                taxon_rank = domain_superkingdom_tax_rank_name

    elif grouping == "selected":
        if taxa_lists_file is None:
            logger.critical('|EsMeCaTa|gseapy_enrichr| You have to specify a taxa lists file for this analysis with --taxa-list.')
            sys.exit()
    elif grouping == "selected_function":
        if function_lists_file is None:
            logger.critical('|EsMeCaTa|gseapy_enrichr| You have to specify a function lists file for this analysis with --function-list.')
            sys.exit()
    elif grouping is None:
        logger.critical('|EsMeCaTa|gseapy_enrichr| You have to choose a grouping factor either "tax_rank", "selected" or "selected_function.')
        sys.exit()

    # Get metadata associated with run.
    options = {'annotation_folder_or_file': annotation_folder_or_file, 'output_folder': output_folder, 'grouping': grouping,
               'taxon_rank': taxon_rank, 'taxa_lists_file': taxa_lists_file,
               'annot_names_file': annot_names_file, 'orsum_minterm_size': orsum_minterm_size,
               'taxon_id_file': taxon_id_file}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['python_package']['pronto'] = pronto.__version__
    options['tool_dependencies']['python_package']['gseapy'] = gseapy.__version__
    options['tool_dependencies']['python_package']['orsum'] = get_orsum_version()
    options['tool_dependencies']['python_package']['pandas'] = pd.__version__
    options['tool_dependencies']['python_package']['biopython'] = biopython_version
    options['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version

    esmecata_metadata = {}
    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')
    esmecata_metadata['access_time'] = date
    esmecata_metadata['tool_options'] = options

    if annot_names_file == 'download':
        annot_names = {}
        # Download enzyme.dat.
        enzyme_data_file = os.path.join(output_folder, 'enzyme.dat')
        urllib.request.urlretrieve('https://ftp.expasy.org/databases/enzyme/enzyme.dat', enzyme_data_file)
        # Download go-basic.obo.
        go_basic_obo_file = os.path.join(output_folder, 'go-basic.obo')
        urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go/go-basic.obo', go_basic_obo_file)
        # Retrieves KO name from KEGG.
        ko_file = os.path.join(output_folder, 'ko_names.tsv')
        ko_names = get_ko_names(ko_file)
        annot_names.update(ko_names)

        enzyme_names, go_names = get_annot_name(enzyme_data_file, go_basic_obo_file)
        annot_names.update(enzyme_names)
        annot_names.update(go_names)
        annotation_names_file = os.path.join(output_folder, 'annotation_names.json')
        with open(annotation_names_file, 'w') as open_annotation_names_file:
            json.dump(annot_names, open_annotation_names_file, indent=4)
    else:
        if os.path.exists(annot_names_file):
            with open(annot_names_file, 'r') as open_annotation_names_file:
                annot_names = json.load(open_annotation_names_file)
        else:
            logger.critical('|EsMeCaTa|gseapy_enrichr| --annot-names/annot_names_file must be a valid file or "download" to either use your own annotation names or EC/GO/KO names.')
            sys.exit()

    if os.path.isdir(annotation_folder_or_file):
        proteome_tax_id_file_path = os.path.join(annotation_folder_or_file, 'proteome_tax_id.tsv')

    # Extract the groups that will be search.
    # Either: - taxonomic groups -> tax_rank.
    #  - manually selected groups of organisms -> selected.
    #  - manually selected groups of functions -> selected_function.
    if grouping == "tax_rank":
        if os.path.isfile(annotation_folder_or_file):
            proteome_tax_id_file_path = taxon_id_file
        identified_groups = get_taxon_obs_name(proteome_tax_id_file_path, taxon_rank)
    elif grouping == "selected":
        identified_groups = extract_organisms_selected(taxa_lists_file)
    elif grouping == "selected_function":
        identified_groups = extract_organisms_selected(function_lists_file)
        obs_names = {}
        if os.path.isdir(annotation_folder_or_file):
            with open(proteome_tax_id_file_path, 'r') as proteome_tax_file:
                csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
                for line in csvreader:
                    observation_name = line['observation_name']
                    tax_name = line['name']
                    obs_names[observation_name] = tax_name

    # Retrieve element that will be searched for over-representation.
    # If input is esmecata annotation folder, iterates on annotation file to generate element sets.
    if os.path.isdir(annotation_folder_or_file):
        annotation_reference_path = os.path.join(annotation_folder_or_file, 'annotation_reference')
        if grouping in ["tax_rank", "selected"]:
            element_sets = extract_annotation(annotation_reference_path)
        elif grouping in ["selected_function"]:
            element_sets = extract_observation_name_annotation(annotation_reference_path)
    # If input is a function table, read it and extract annotations/organisms.
    elif os.path.isfile(annotation_folder_or_file):
        if grouping in ["tax_rank", "selected"]:
            element_sets = read_function_table(annotation_folder_or_file)
        elif grouping in ["selected_function"]:
            element_sets = read_function_table(annotation_folder_or_file, transpose_table=False)

    # Create GMT file for orsum.
    # When using tax_rank or selected parameters, the GMT file contains each annotation, their labels and the observation name in which they were found.
    # When using selected_function, the GMT file contains each taxon, their associated taxon name and the annotations they are linked to. 
    annotation_file_gmt_file = os.path.join(output_folder, 'annotation_file.gmt')
    with open(annotation_file_gmt_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        for element in element_sets:
            element_name = 'na'
            if grouping in ["tax_rank", "selected"]:
                if element in annot_names:
                    element_name = element + ' ' + annot_names[element]
                else:
                    element_name = element
            elif grouping in ["selected_function"]:
                if element in obs_names:
                    element_name = element + ' ' + obs_names[element]
                else:
                    element_name = element
            csvwriter.writerow([element, element_name, *element_sets[element]])

    output_dir = os.path.join(output_folder, 'enrichr_module')
    orsum_input_folder = os.path.join(output_folder, 'orsum_input_folder')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(orsum_input_folder):
        os.mkdir(orsum_input_folder)

    enriched_elements = {}
    for group_name in identified_groups:
        identified_group = identified_groups[group_name]
        # Try to run gseapy enrichr, if no enriched results, continue.
        try:
            gseapy.enrichr(gene_list=identified_group, gene_sets=element_sets, background=None,
                        outdir=os.path.join(output_dir, group_name), cutoff=selected_adjust_pvalue_cutoff)
        except ValueError as error:
            logger.critical('|EsMeCaTa|gseapy_enrichr| No enrichred functions with p-value cutoff < {0} for {1}.'.format(selected_adjust_pvalue_cutoff, group_name))
            continue
        if os.path.exists(os.path.join(output_dir, group_name, 'gs_ind_0.human.enrichr.reports.pdf')):
            # If enriched results, extract the ones with an adjusted p-value inferior to selected_adjust_pvalue_cutoff to output folder (default 0.05).
            df = pd.read_csv(os.path.join(output_dir, group_name, 'gs_ind_0.human.enrichr.reports.txt'), sep='\t')
            df = df[df['Adjusted P-value'] < selected_adjust_pvalue_cutoff]
            enriched_elements[group_name] = df.set_index('Term')['Adjusted P-value'].to_dict()
            df.sort_values('Adjusted P-value', inplace=True)
            df = df['Term']
            df.to_csv(os.path.join(orsum_input_folder, group_name), sep='\t', index=False, header=False)

    all_elements = set([element for org in enriched_elements for element in enriched_elements[org]])
    all_groups = list(identified_groups.keys())

    enrich_matrix_file = os.path.join(output_folder, 'enrich_matrix.tsv')
    with open(enrich_matrix_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['Element', 'Name', *all_groups])
        for element in all_elements:
            if grouping in ["tax_rank", "selected"]:
                if element in annot_names:
                    element_name = element + ' ' + annot_names[element]
                else:
                    element_name = 'na'
            elif grouping in ["selected_function"]:
                if element in obs_names:
                    element_name = element + ' ' + obs_names[element]
                else:
                    element_name = 'na'
            csvwriter.writerow([element, element_name, *[enriched_elements[group_name][element] if element in enriched_elements[group_name] else 'NA'  for group_name in all_groups]])

    # Check that there are files in orsum folder.
    orsum_input_folder_files = os.listdir(orsum_input_folder)
    if len(orsum_input_folder_files) == 0:
        logger.critical('|EsMeCaTa|gseapy_enrichr| No enriched files from gseapy enrichr as input for orsum. It seems that there are no enrichd terms found.')
        sys.exit(1)

    # Run orsum to filter list of enriched annotations.
    logger.info('|EsMeCaTa|gseapy_enrichr| Launch orsum visualisation.')
    orsum_output_folder = os.path.join(output_folder, 'orsum_output_folder')
    run_orsum(annotation_file_gmt_file, orsum_input_folder, orsum_output_folder, orsum_minterm_size)

    endtime = time.time()
    duration = endtime - starttime
    logger.info('|EsMeCaTa|gseapy_enrichr| Enrichment analysis took {0}.'.format(duration))

    esmecata_metadata['esmecata_gseapy_enrichr_duration'] = duration
    gseapy_enrichr_metadata_file = os.path.join(output_folder, 'esmecata_metadata_gseapy_enrichr.json')
    if os.path.exists(gseapy_enrichr_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_gseapy_enrichr' in metadata_file]
        gseapy_enrichr_metadata_file = os.path.join(output_folder, 'esmecata_metadata_gseapy_enrichr{0}.json'.format(len(metadata_files)))
        with open(gseapy_enrichr_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)
    else:
        with open(gseapy_enrichr_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)
