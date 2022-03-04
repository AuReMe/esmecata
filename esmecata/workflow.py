import csv
import json
import os
import time

from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins

def compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file):
    result_folder = os.path.join(proteomes_output_folder, 'result')
    proteome_numbers = {}
    for folder in os.listdir(result_folder):
        result_folder_path = os.path.join(result_folder, folder)
        number_proteome = len([proteome for proteome in os.listdir(result_folder_path)])
        proteome_numbers[folder] = number_proteome

    clustering_folder = os.path.join(clustering_output_folder, 'reference_proteins')
    clustering_numbers = {}
    for clustering_file in os.listdir(clustering_folder):
        clustering_file_path = os.path.join(clustering_folder, clustering_file)
        num_lines = sum(1 for line in open(clustering_file_path))
        clustering_numbers[clustering_file.replace('.tsv', '')] = num_lines

    annotation_reference_folder = os.path.join(annotation_output_folder, 'annotation_reference')
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

    all_observation_names = {*proteome_numbers.keys(), *clustering_numbers.keys(), *annotation_numbers.keys()}
    with open(stat_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_proteomes', 'Number_shared_proteins', 'Number_go_terms', 'Number_ecs'])
        for observation_name in all_observation_names:
            if observation_name in proteome_numbers:
                nb_proteomes = proteome_numbers[observation_name]
            else:
                nb_proteomes = 'NA'
            if observation_name in clustering_numbers:
                nb_shared_proteins = clustering_numbers[observation_name]
            else:
                nb_shared_proteins = 'NA'
            if observation_name in annotation_numbers:
                nb_gos = annotation_numbers[observation_name][0]
                nb_ecs = annotation_numbers[observation_name][1]
            else:
                nb_proteomes = 'NA'
            csvwriter.writerow([observation_name, nb_proteomes, nb_shared_proteins, nb_gos, nb_ecs])


def perform_workflow(input_file, output_folder, busco_percentage_keep=80, ignore_taxadb_update=None,
                        all_proteomes=None, uniprot_sparql_endpoint=None, remove_tmp=None,
                        limit_maximal_number_proteomes=99, rank_limit=None,
                        nb_cpu=1, clust_threshold=1, mmseqs_options=None,
                        linclust=None, propagate_annotation=None, uniref_annotation=None,
                        expression_annotation=None, beta=None):
    starttime = time.time()
    workflow_metadata = {}

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteomes_output_folder = os.path.join(output_folder, '0_proteomes')
    retrieve_proteomes(input_file, proteomes_output_folder, busco_percentage_keep,
                        ignore_taxadb_update, all_proteomes, uniprot_sparql_endpoint,
                        remove_tmp, limit_maximal_number_proteomes, beta, rank_limit)

    clustering_output_folder = os.path.join(output_folder, '1_clustering')
    make_clustering(proteomes_output_folder, clustering_output_folder, nb_cpu, clust_threshold, mmseqs_options, linclust, remove_tmp)

    annotation_output_folder = os.path.join(output_folder, '2_annotation')
    annotate_proteins(clustering_output_folder, annotation_output_folder, uniprot_sparql_endpoint, propagate_annotation, uniref_annotation, expression_annotation, beta)

    stat_file = os.path.join(output_folder, 'stat_number_workflow.tsv')
    compute_stat_workflow(proteomes_output_folder, clustering_output_folder, annotation_output_folder, stat_file)

    proteomes_metadata_file = os.path.join(proteomes_output_folder, 'esmecata_metadata_proteomes.json')
    with open(proteomes_metadata_file, 'r') as json_idata:
        proteomes_metadata = json.load(json_idata)
    clustering_metadata_file = os.path.join(clustering_output_folder, 'esmecata_metadata_clustering.json')
    with open(clustering_metadata_file, 'r') as json_idata:
        clustering_metadata = json.load(json_idata)
    annotation_metadata_file = os.path.join(annotation_output_folder, 'esmecata_metadata_annotation.json')
    with open(annotation_metadata_file, 'r') as json_idata:
        annotation_metadata = json.load(json_idata)

    workflow_metadata['proteomes_metadata'] = proteomes_metadata
    workflow_metadata['clustering_metadata'] = clustering_metadata
    workflow_metadata['annotation_metadata'] = annotation_metadata

    endtime = time.time()
    duration = endtime - starttime
    workflow_metadata['esmecata_workflow_duration'] = duration
    workflow_metadata_file = os.path.join(output_folder, 'esmecata_metadata_workflow.json')
    with open(workflow_metadata_file, 'w') as ouput_file:
        json.dump(workflow_metadata, ouput_file, indent=4)
