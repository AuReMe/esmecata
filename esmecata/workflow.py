import os

from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins

def perform_workflow(input_file, output_folder, busco_percentage_keep=80, ignore_taxadb_update=None,
                        all_proteomes=None, uniprot_sparql_endpoint=None, remove_tmp=None,
                        limit_maximal_number_proteomes=99, rank_limit=None,
                        nb_cpu=1, clust_threshold=1, mmseqs_options=None,
                        linclust=None, propagate_annotation=None, uniref_annotation=None,
                        expression_annotation=None, beta=None):

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