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
import datetime
import pandas as pd
import os
import numpy as np
import itertools
import json
import logging
import time
import shutil
import subprocess
import sys

from esmecata.utils import is_valid_dir
from esmecata import __version__ as esmecata_version
from esmecata.annotation import extract_protein_cluster

logger = logging.getLogger(__name__)


def get_eggnog_version():
    """Returns version of eggnog-mapper.

    Returns:
        eggnog_version (str): version of eggnog-mapper.
    """
    eggnog_version = None
    response = subprocess.Popen(['emapper.py', '-v'], stdout=subprocess.PIPE, start_new_session=True, universal_newlines="")
    for eggnog_line in response.stdout:
        eggnog_line_decoded = eggnog_line.decode('utf-8')
        if 'emapper-' in eggnog_line_decoded:
            eggnog_version = eggnog_line_decoded.split(' / ')[0].replace('emapper-', '')

    if eggnog_version is None:
        logger.critical('esmecata could not find the version of eggnog-mapper.')
        logger.critical('It is possibly an issue with the installation of eggnog-mapper (maybe it is not in the PATH). Or it can be due to a change in the output of emapper.py -v command.')
        sys.exit()

    return eggnog_version


def get_proteomes_tax_name(proteomes_tax_id):
    """ Extract tax_id associated with observation name.

    Args:
        proteomes_tax_id (str): pathname to the proteomes_tax_id file.

    Returns:
        proteomes_taxa_names (dict): dict containing observation names (as key) associated with tax name used for proteomes (as value)
    """
    proteomes_taxa_names = {}
    with open(proteomes_tax_id, 'r') as proteome_tax_file:
        csvreader = csv.DictReader(proteome_tax_file, delimiter='\t')
        for line in csvreader:
            observation_name = line['observation_name']
            tax_name = line['name']
            proteomes_taxa_names[observation_name] = tax_name.replace(' ', '_')

    return proteomes_taxa_names


def call_to_emapper(input_path, output_name, output_dir, temporary_dir, eggnog_database_path, nb_cpu):
    """ Calls eggnog-mapper on the protein clusters.

    Args:
        input_path (str): pathname to the fasta file containing consensus sequences for protein clusters.
        output_name (str): observation name associated with the protein sequences.
        output_dir (str)): path to the output folder.
        temporayyr_dir (str): path to temporary folder.
        eggnog_database_path (str): pathname to eggnog database folder.
        nb_cpu (int): number of CPUs to use with eggnog-mapper.
    """
    # Use dbmem for faster run.
    # Use override to avoid issue if already present annotation (if there was a failed run for example).
    eggnog_cmds = ['emapper.py', '--cpu', nb_cpu, '-i', input_path, '--output', output_name,
                   '--data_dir', eggnog_database_path, '--itype', 'proteins', '--output_dir', output_dir,
                   '--dbmem', '--temp_dir', temporary_dir, '--override']

    subprocess.call(eggnog_cmds)


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
                csvreader = csv.DictReader(open_annotation_input_file_path, delimiter='\t')
                for line in csvreader:
                    gos = line['GO'].split(',')
                    ecs = line['EC'].split(',')
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


def read_annotation(eggnog_outfile:str):
    """Read an eggnog-mapper annotation file and retrieve EC numbers and GO terms by genes.
    Args:
        eggnog_outfile (str): path to eggnog-mapper annotation file
    Returns:
        dict: dict of genes and their annotations as {gene1:{EC:'..,..', GOs:'..,..,'}}
    """
    # Look at the twentieth first rows to find the header.
    with open(eggnog_outfile, 'r') as f:
        twentieth_first_rows = list(itertools.islice(f, 20))
        first_row_after_header = min([index for index, str_row in enumerate(twentieth_first_rows) if not str_row.startswith('#')])
        header_row = first_row_after_header - 1
        headers_row = twentieth_first_rows[header_row].lstrip("#").strip().split('\t')

    # Fix issue when header is incomplete (eggnog before version 2.0).
    if len(headers_row) == 17:
        headers_row.extend(['tax_scope', 'eggNOG_OGs', 'bestOG', 'COG_functional_category', 'eggNOG_free_text'])

    to_extract_annotations = ['GOs','EC', 'Preferred_name']
    if 'PFAMs' in headers_row:
        to_extract_annotations.append('PFAMs')
    if 'BiGG_Reaction' in headers_row:
        to_extract_annotations.append('BiGG_Reaction')
    if 'KEGG_Reaction' in headers_row:
        to_extract_annotations.append('KEGG_Reaction')
    if 'CAZy' in headers_row:
        to_extract_annotations.append('CAZy')

    # Use chunk when reading eggnog file to cope with big file.
    chunksize = 10 ** 6
    for annotation_data in pd.read_csv(eggnog_outfile, sep='\t', comment='#', header=None, dtype = str, chunksize = chunksize):
        annotation_data.replace(np.nan, '', inplace=True)
        # Assign the headers
        annotation_data.columns = headers_row
        if 'query_name' in annotation_data.columns:
            # Check if the gene IDs are numeric, if yes add 'gene_' in front of them.
            numeric_row_dataframe = pd.to_numeric(annotation_data['query_name'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query_name'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query_name']
            annotation_dict = annotation_data.set_index('query_name')[to_extract_annotations].to_dict('index')
        # 'query' added for compatibility with eggnog-mapper 2.1.2
        elif 'query' in annotation_data.columns:
            numeric_row_dataframe = pd.to_numeric(annotation_data['query'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query']
            annotation_dict = annotation_data.set_index('query')[to_extract_annotations].to_dict('index')
        for key in annotation_dict:
            yield key, annotation_dict[key]


def write_pathologic(base_filename, annotated_proteins, pathologic_output_file, reference_proteins):
    """Write the annotation associated with a cluster after propagation step into pathologic file for run on Pathway Tools.

    Args:
        base_filename (str): observation name
        annotated_proteins (dict): dict of protein and their annotations as {protein:{EC:'..,..', GOs:'..,..,'}}
        reference_proteins (dict): dict containing representative protein IDs (as key) associated with proteins of the cluster
        pathologic_output_file (str): pathname to output pathologic file
    """
    with open(pathologic_output_file, 'w', encoding='utf-8') as element_file:
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        element_file.write(';; ' + base_filename + '\n')
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        for protein_annots in annotated_proteins:
            protein = protein_annots[0].split('|')[1]

            protein_annot = protein_annots[1]
            element_file.write('ID\t' + protein + '\n')

            protein_name = protein_annot['Preferred_name']
            if protein_name not in ['', '-']:
                element_file.write('NAME\t' + protein_name + '\n')
            else:
                element_file.write('NAME\t' + protein + '\n')
            element_file.write('PRODUCT-TYPE\tP' + '\n')
            element_file.write('PRODUCT-ID\tprot ' + protein + '\n')
            for uniprot_protein in reference_proteins[protein]:
                element_file.write('DBLINK\tUNIPROT:' + uniprot_protein + '\n')
            if 'GOs' in protein_annot:
                gos = [go for go in protein_annot['GOs'].split(',') if go not in ['', '-']]
                for go in gos:
                    element_file.write('GO\t' + go + '\n')
            if 'EC' in protein_annot:
                ecs = [ec for ec in protein_annot['EC'].split(',') if ec not in ['', '-']]
                for ec in ecs:
                    element_file.write('EC\t' + ec + '\n')
            element_file.write('//\n\n')


def write_annotation_reference(protein_annotations, reference_proteins, annotation_reference_file):
    """Write the annotation associated with a cluster after propagation step.

    Args:
        protein_annotations (dict): dict of protein and their annotations as {protein:{EC:'..,..', GOs:'..,..,'}}
        reference_proteins (dict): dict containing representative protein IDs (as key) associated with proteins of the cluster
        annotation_reference_file (str): pathname to output tabulated file
    """
    with open(annotation_reference_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        csvwriter.writerow(['protein_cluster', 'cluster_members', 'gene_name', 'GO', 'EC', 'KEGG_reaction'])
        for protein_annot_tuples in protein_annotations:
            protein = protein_annot_tuples[0].split('|')[1]
            protein_annot = protein_annot_tuples[1]
            gene_name = protein_annot['Preferred_name']
            cluster_members = ','.join(reference_proteins[protein])
            if 'GOs' in protein_annot:
                gos = [go for go in protein_annot['GOs'].split(',') if go not in ['', '-']]
                gos = ','.join(sorted(gos))
            else:
                gos = ''

            if 'EC' in protein_annot:
                ecs = [ec for ec in protein_annot['EC'].split(',') if ec not in ['', '-']]
                ecs = ','.join(sorted(ecs))
            else:
                ecs = ''

            if 'KEGG_Reaction' in protein_annot:
                keggs = [kegg_rxn for kegg_rxn in protein_annot['KEGG_Reaction'].split(',') if kegg_rxn not in ['', '-']]
                keggs = ','.join(sorted(keggs))
            else:
                keggs = ''

            csvwriter.writerow([protein, cluster_members, gene_name, gos, ecs, keggs])


def annotate_with_eggnog(input_folder, output_folder, eggnog_database_path, nb_cpu, eggnog_tmp_dir=None):
    """Write the annotation associated with a cluster after propagation step into pathologic file for run on Pathway Tools.

    Args:
        input_folder (str): pathname to mmseqs clustering output folder as it is the input of annotation.
        output_folder (str): pathname to the output folder.
        eggnog_database_path (str): pathname to eggnog database folder.
        nb_cpu (int): number of CPUs to be used by eggnog-mapper.
        eggnog_tmp_dir (str): pathname to eggnog-mapper temporary folder.
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|annotation-eggnog| Begin annotation.')

    is_valid_dir(output_folder)

    eggnog_output_folder = os.path.join(output_folder, 'eggnog_output')
    is_valid_dir(eggnog_output_folder)

    annotation_reference_folder = os.path.join(output_folder, 'annotation_reference')
    is_valid_dir(annotation_reference_folder)

    pathologic_folder = os.path.join(output_folder, 'pathologic')
    is_valid_dir(pathologic_folder)

    reference_protein_fasta_path = os.path.join(input_folder, 'reference_proteins_consensus_fasta')
    taxa_names = [input_file.replace('.faa', '') for input_file in os.listdir(reference_protein_fasta_path)]

    clustering_taxon_id_file = os.path.join(input_folder, 'proteome_tax_id.tsv')
    annotation_taxon_id_file = os.path.join(output_folder, 'proteome_tax_id.tsv')

    if os.path.exists(annotation_taxon_id_file):
        if not os.path.samefile(clustering_taxon_id_file, annotation_taxon_id_file):
            os.remove(annotation_taxon_id_file)
            shutil.copyfile(clustering_taxon_id_file, annotation_taxon_id_file)
    else:
        shutil.copyfile(clustering_taxon_id_file, annotation_taxon_id_file)
    proteomes_tax_names = get_proteomes_tax_name(annotation_taxon_id_file)

    reference_protein_path = os.path.join(input_folder, 'reference_proteins')

    # Convert CPU int into str.
    nb_cpu = str(nb_cpu)
    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder, 'nb_cpu': nb_cpu,
               'eggnog_database_path': eggnog_database_path, 'eggnog_tmp_dir': eggnog_tmp_dir}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = esmecata_version
    options['tool_dependencies']['eggnog_mapper'] = get_eggnog_version()

    esmecata_metadata = {}
    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')
    esmecata_metadata['access_time'] = date
    esmecata_metadata['tool_options'] = options

    proteome_tax_id_file = os.path.join(input_folder, 'proteome_tax_id.tsv')

    # Run eggnog-mapper on taxon name file.
    for taxa_name in taxa_names:
        # Create temporary folder for eggnog-mapper.
        if eggnog_tmp_dir is None:
            eggnog_temporary_dir = os.path.join(output_folder, 'eggnog_temporary')
        else:
            eggnog_temporary_dir = os.path.join(eggnog_tmp_dir, 'eggnog_temporary')
        is_valid_dir(eggnog_temporary_dir)

        fasta_file_path = os.path.join(reference_protein_fasta_path, taxa_name+'.faa')

        # Check if the annotation by eggnog-mapper has been performed.
        eggnog_mapper_annotation_file = os.path.join(eggnog_output_folder, taxa_name+'.emapper.annotations')
        if not os.path.exists(eggnog_mapper_annotation_file):
            logger.info('|EsMeCaTa|annotation| Launch eggnog-mapper on %s.',  taxa_name)
            call_to_emapper(fasta_file_path, taxa_name, eggnog_output_folder, eggnog_temporary_dir, eggnog_database_path, nb_cpu)
        else:
            logger.info('|EsMeCaTa|annotation| Results of eggnog-mapper already present for %s.',  taxa_name)

        # If temporary folder exists, remove it.
        if os.path.exists(eggnog_temporary_dir):
            shutil.rmtree(eggnog_temporary_dir)

    # Create annotation reference and pathologic files for each observation name.
    for observation_name in proteomes_tax_names:
        proteomes_tax_name = proteomes_tax_names[observation_name]
        reference_protein_pathname = os.path.join(reference_protein_path, proteomes_tax_name+'.tsv')
        reference_proteins, set_proteins = extract_protein_cluster(reference_protein_pathname)

        # Read eggnog output.
        eggnog_mapper_annotation_file = os.path.join(eggnog_output_folder, proteomes_tax_name+'.emapper.annotations')
        annotated_proteins = read_annotation(eggnog_mapper_annotation_file)
        annotated_proteins = list(annotated_proteins)

        gos = [go for protein_id, protein_annot in annotated_proteins for go in protein_annot['GOs'].split(',') if go not in ['', '-']]
        unique_gos = set(gos)
        ecs = [ec for protein_id, protein_annot in annotated_proteins for ec in protein_annot['EC'].split(',') if ec not in ['', '-']]
        unique_ecs = set(ecs)
        logger.info('|EsMeCaTa|annotation| %d Go Terms (with %d unique GO Terms) and %d EC numbers (with %d unique EC) associated with %s.', len(gos),
                                                                                                len(unique_gos), len(ecs), len(unique_ecs), observation_name)

        # Create annotation reference file.
        annotation_reference_file = os.path.join(annotation_reference_folder, observation_name+'.tsv')
        if not os.path.exists(annotation_reference_file):
            write_annotation_reference(annotated_proteins, reference_proteins, annotation_reference_file)

        # Create pathologic files.
        pathologic_organism_folder = os.path.join(pathologic_folder, observation_name)
        # Add _1 to pathologic file as genetic element cannot have the same name as the organism.
        pathologic_file = os.path.join(pathologic_organism_folder, observation_name+'_1.pf')

        if not os.path.exists(pathologic_file):
            if len(annotated_proteins) > 0:
                is_valid_dir(pathologic_organism_folder)
                write_pathologic(observation_name, annotated_proteins, pathologic_file, reference_proteins)
            elif len(annotated_proteins) == 0:
                logger.critical('|EsMeCaTa|annotation-eggnog| No reference proteins for %s, esmecata will not create a pathologic folder for it.', observation_name)

        else:
            logger.critical('|EsMeCaTa|annotation-eggnog| Annotation already performed for %s.', observation_name)

    # Create mpwt taxon ID file.
    clustering_taxon_id = {}
    with open(proteome_tax_id_file, 'r') as input_taxon_id_file:
        taxon_id_csvreader = csv.DictReader(input_taxon_id_file, delimiter='\t')
        for line in taxon_id_csvreader:
            clustering_taxon_id[line['observation_name']] = line['tax_id']

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
    esmecata_metadata['esmecata_annotation_method'] = 'eggnog-mapper'
    esmecata_metadata['esmecata_annotation_duration'] = duration
    eggnog_metadata_file = os.path.join(output_folder, 'esmecata_metadata_annotation.json')
    if os.path.exists(eggnog_metadata_file):
        metadata_files = [metadata_file for metadata_file in os.listdir(output_folder) if 'esmecata_metadata_annotation' in metadata_file]
        eggnog_metadata_file = os.path.join(output_folder, 'esmecata_metadata_annotation_{0}.json'.format(len(metadata_files)))
        with open(eggnog_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)
    else:
        with open(eggnog_metadata_file, 'w') as ouput_file:
            json.dump(esmecata_metadata, ouput_file, indent=4)

    logger.info('|EsMeCaTa|annotation-eggnog| Annotation complete in {0}s.'.format(duration))
