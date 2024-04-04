# Copyright (C) 2023-2024 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
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
import pandas as pd
import pronto
import subprocess
import matplotlib.pyplot as plt
import urllib.request
import gseapy

from Bio.ExPASy import Enzyme
from esmecata.esmecata_analysis.analysis import get_taxon_obs_name

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


def run_orsum(annotation_file_gmt_file, enrichr_module_phylum_output, orsum_output_folder):
    """ Run orsum to reduce list of enriched annotations from gseapy.

    Args:
        annotation_file_gmt_file (str): GMT file containing annotations and their labels
        enrichr_module_phylum_output (str): path to output of gseapy containing enriched annotations for each taxon
        orsum_output_folder (str): path to output folder for orsum
    """
    input_files = [os.path.join(enrichr_module_phylum_output, file) for file in os.listdir(enrichr_module_phylum_output)]
    subprocess.call(['orsum.py', '--gmt', annotation_file_gmt_file, '--files', *input_files, '--outputFolder', orsum_output_folder])


def taxon_rank_annotation_enrichment(annotation_folder, output_folder, enzyme_data_file=None, go_basic_obo_file=None, taxon_rank='phylum'):
    """ Run an enrichment analysis on taxon from annotation results of esmecata using gseapy.
    Then filter this list with orsum.

    Args:
        annotation_folder (str): path to esmecata annotation folder
        output_folder (str): path to output folder
        enzyme_data_file (str): path to expasy enzyme.dat file, if not given, download it
        go_basic_obo_file (str): path to Gene Ontology go-basic.obo file, if not given, download it
        taxon_rank (str): taxon rank to cluster the observation name together (by default, phylum)
    """
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
    taxa_name = get_taxon_obs_name(proteome_tax_id_file_path, taxon_rank)

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
    for tax_name in taxa_name:
        organisms = taxa_name[tax_name]
        # Try to run gseapy enrichr, if no enriched results, continue.
        try:
            gseapy.enrichr(gene_list=organisms, gene_sets=annotation_sets, background=None,
                        outdir=os.path.join(output_dir, tax_name))
        except:
            continue
        if os.path.exists(os.path.join(output_dir, tax_name, 'gs_ind_0.human.enrichr.reports.pdf')):
            # If enriched results, extract the ones with an adjusted p-value inferior to 0.05 to output folder.
            df = pd.read_csv(os.path.join(output_dir, tax_name, 'gs_ind_0.human.enrichr.reports.txt'), sep='\t')
            df = df[df['Adjusted P-value'] < 0.05]
            df.sort_values('Adjusted P-value', inplace=True)
            df = df['Term']
            df.to_csv(os.path.join(orsum_input_folder, tax_name), sep='\t', index=False)

    enriched_elements = {}
    for result_file in os.listdir('enrichr_module_phylum_output'):
        if result_file.endswith('.txt'):
            df = pd.read_csv(os.path.join('enrichr_module_phylum_output', result_file), sep='\t')
            df.set_index('Term', inplace=True)
            base_name = os.path.splitext(os.path.basename(result_file))[0]
            enriched_elements[base_name] = df[df['Adjusted P-value'] < 0.05]['Odds Ratio'].to_dict()

    all_elments = set([element for org in enriched_elements for element in enriched_elements[org]])

    with open('enrich_matrix_phylum.tsv', 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['Organism', *all_elments])
        for org in enriched_elements:
            csvwriter.writerow([org, *[enriched_elements[org][element] if element in enriched_elements[org] else 'NA' for element in all_elments]])
    result_df = pd.read_csv('enrich_matrix_phylum.tsv', sep='\t')
    result_df.set_index('Organism', inplace=True)
    fig, axes = plt.subplots(figsize=(70,25))
    plt.rc('font', size=4)

    # Run orsum to filter list of enriched annotations.
    orsum_output_folder = os.path.join(output_folder, 'orsum_output_folder')
    run_orsum(annotation_file_gmt_file, orsum_input_folder, orsum_output_folder)