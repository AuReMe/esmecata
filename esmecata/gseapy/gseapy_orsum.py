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

from Bio.ExPASy import Enzyme
from matplotlib import __version__ as matplotlib_version
from esmecata.core.proteomes import get_taxon_obs_name
from esmecata.utils import get_domain_or_superkingdom_from_ncbi_tax_database
from esmecata import __version__ as esmecata_version
from Bio import __version__ as biopython_version

logger = logging.getLogger(__name__)


def extract_annotation(annotation_reference_folder):
    """ From annotation reference folder, creates a dict indicating the file associated with each annotation.

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


def taxon_rank_annotation_enrichment(annotation_folder, output_folder, grouping="tax_rank",
                                     taxon_rank='phylum', taxa_lists_file=None,
                                     enzyme_data_file=None, go_basic_obo_file=None, orsum_minterm_size=None,
                                     selected_adjust_pvalue_cutoff=0.05):
    """ Run an enrichment analysis on taxon from annotation results of esmecata using gseapy.
    Then filter this list with orsum.

    Args:
        annotation_folder (str): path to esmecata annotation folder
        output_folder (str): path to output folder
        grouping (str): grouping factor, either "tax_rank" or "selected"
        taxon_rank (str): taxon rank to cluster the observation name together (by default, phylum) when selecting grouping "tax_rank"
        taxa_lists_file (str): path to manually created groups of observation names when selecting grouping "selected"
        enzyme_data_file (str): path to expasy enzyme.dat file, if not given, download it
        go_basic_obo_file (str): path to Gene Ontology go-basic.obo file, if not given, download it
        orsum_minterm_size (int): option minTermSize of orsum
        selected_adjust_pvalue_cutoff (float): adjust-Pval cutoff for gseapy enrichr, default: 0.05
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|gseapy_enrichr| Begin enrichment analysis.')

    if grouping == "tax_rank":
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
            logger.critical('|EsMeCaTa|gseapy_enrichr| You have to specify a taxa lists file for this analysis with --taxa-list')
            sys.exit()
    elif grouping is None:
        logger.critical('|EsMeCaTa|gseapy_enrichr| You have to choose a grouping factor either "tax_rank" or "selected".')
        sys.exit()

    # Get metadata associated with run.
    options = {'annotation_folder': annotation_folder, 'output_folder': output_folder, 'grouping': grouping,
               'taxon_rank': taxon_rank, 'taxa_lists_file': taxa_lists_file,
               'enzyme_data_file': enzyme_data_file, 'go_basic_obo_file': go_basic_obo_file, 'orsum_minterm_size': orsum_minterm_size}

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

    # If no enzyme.Dat available download it.
    if enzyme_data_file is None:
        enzyme_data_file = os.path.join(output_folder, 'enzyme.dat')
        urllib.request.urlretrieve('https://ftp.expasy.org/databases/enzyme/enzyme.dat', enzyme_data_file)
    # If no go-basic.obo file available, download it.
    if go_basic_obo_file is None:
        go_basic_obo_file = os.path.join(output_folder, 'go-basic.obo')
        urllib.request.urlretrieve('http://purl.obolibrary.org/obo/go/go-basic.obo', go_basic_obo_file)

    enzyme_names, go_names = get_annot_name(enzyme_data_file, go_basic_obo_file)

    proteome_tax_id_file_path = os.path.join(annotation_folder, 'proteome_tax_id.tsv')

    if grouping == "tax_rank":
        taxa_name = get_taxon_obs_name(proteome_tax_id_file_path, taxon_rank)
    elif grouping == "selected":
        taxa_name = extract_organisms_selected(taxa_lists_file)

    annotation_reference_path = os.path.join(annotation_folder, 'annotation_reference')
    annotation_sets = extract_annotation(annotation_reference_path)

    # Create GMT file for orsum.
    # The GMT file contains each annotation, their labels and the observation name in which they were found.
    annotation_file_gmt_file = os.path.join(output_folder, 'annotation_file.gmt')
    with open(annotation_file_gmt_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        for annotation in annotation_sets:
            annot_name = 'na'
            if annotation in enzyme_names:
                annot_name = annotation + ' ' + enzyme_names[annotation]
            if annotation in go_names:
                annot_name = annotation + ' ' + go_names[annotation]
            csvwriter.writerow([annotation, annot_name, *annotation_sets[annotation]])

    output_dir = os.path.join(output_folder, 'enrichr_module')
    orsum_input_folder = os.path.join(output_folder, 'orsum_input_folder')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(orsum_input_folder):
        os.mkdir(orsum_input_folder)

    enriched_elements = {}
    for tax_name in taxa_name:
        organisms = taxa_name[tax_name]
        # Try to run gseapy enrichr, if no enriched results, continue.
        try:
            gseapy.enrichr(gene_list=organisms, gene_sets=annotation_sets, background=None,
                        outdir=os.path.join(output_dir, tax_name), cutoff=selected_adjust_pvalue_cutoff)
        except ValueError as error:
            logger.info('|EsMeCaTa|gseapy_enrichr| No enrichred functions with p-value cutoff < {0} for {1}.'.format(selected_adjust_pvalue_cutoff, tax_name))
            continue
        if os.path.exists(os.path.join(output_dir, tax_name, 'gs_ind_0.human.enrichr.reports.pdf')):
            # If enriched results, extract the ones with an adjusted p-value inferior to selected_adjust_pvalue_cutoff to output folder.
            df = pd.read_csv(os.path.join(output_dir, tax_name, 'gs_ind_0.human.enrichr.reports.txt'), sep='\t')
            df = df[df['Adjusted P-value'] < selected_adjust_pvalue_cutoff]
            enriched_elements[tax_name] = df.set_index('Term')['Adjusted P-value'].to_dict()

            df.sort_values('Adjusted P-value', inplace=True)
            df = df['Term']
            df.to_csv(os.path.join(orsum_input_folder, tax_name), sep='\t', index=False, header=False)

    all_elments = set([element for org in enriched_elements for element in enriched_elements[org]])

    enrich_matrix_file = os.path.join(output_folder, 'enrich_matrix.tsv')
    with open(enrich_matrix_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['Organism', *all_elments])
        for org in enriched_elements:
            csvwriter.writerow([org, *[enriched_elements[org][element] if element in enriched_elements[org] else 'NA' for element in all_elments]])

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