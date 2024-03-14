import os
import json
import logging
import pandas as pd

from esmecata.esmecata_analysis.analysis import create_dataset_annotation_file
from esmecata.esmecata_analysis.ec_sunburst_per_model import create_sunburst_ec
from esmecata.esmecata_analysis.viz_esmecata_datapane import create_datapane
from esmecata.esmecata_analysis.esmecata2taxontology import get_fig_parameters, generate_sunburst_fig
from esmecata.esmecata_analysis.stats_workflow_figures import distributions_by_ranks, compare_ranks_in_out, _format_taxo, create_proteome_representativeness_lineplot_px, \
                                                    data_proteome_representativeness, proteomes_representativeness_details
from esmecata.esmecata_analysis.esmecata_compression import esmecata_compression_taxonomy_file

logger = logging.getLogger(__name__)


def run_analysis_pipeline(input_file, input_folder, output_folder):
    logger.info("--- Launch creation of report for workflow ---")
    annotation_referene_folder_path = os.path.join(input_folder, '2_annotation', 'annotation_reference')
    dataset_annotation_ec_file_path = os.path.join(output_folder, 'dataset_annotation_ec.tsv')
    create_dataset_annotation_file(annotation_referene_folder_path, dataset_annotation_ec_file_path, content='EC')

    dataset_annotation_go_file_path = os.path.join(output_folder, 'dataset_annotation_go.tsv')
    create_dataset_annotation_file(annotation_referene_folder_path, dataset_annotation_go_file_path, content='GO')
    #create_sunburst_ec(input_folder, output_folder)
    create_datapane(input_file, input_folder, output_folder)


def run_proteomes_report_creation(input_folder, output_folder):
    logger.info("--- Launch creation of report for proteomes ---")
    df_stat = pd.read_csv(os.path.join(input_folder, 'stat_number_proteome.tsv'), header=0, index_col='observation_name', sep='\t')
    df_proteome_tax_id = pd.read_csv(os.path.join(input_folder, 'proteome_tax_id.tsv'), header=0, index_col='observation_name', sep='\t')
    df_taxon_proteome_tax_id = _format_taxo(df_proteome_tax_id)
    df_stat_join = df_stat.join(df_proteome_tax_id)
    df_stat_join = df_stat_join.join(df_taxon_proteome_tax_id)
    fig1 = distributions_by_ranks(df_stat_join, output_folder, 15, 'phylum')

    with open(os.path.join(input_folder, 'association_taxon_taxID.json')) as f:
        association_taxon_tax_id = json.load(f)
    fig4 = compare_ranks_in_out(df_proteome_tax_id, association_taxon_tax_id, output_folder)

    tax_comp_file = os.path.join(input_folder, 'taxonomy_diff.tsv')
    output_file = os.path.join(output_folder, 'esmecata2taxonomy')
    data = get_fig_parameters(tax_comp_file)
    fig = generate_sunburst_fig(data, output_file)

    output_file = os.path.join(output_folder, 'esmecata_compression')
    fig10 = esmecata_compression_taxonomy_file(tax_comp_file, output_file, False)


def run_clustering_report_creation(input_folder, output_folder):
    logger.info("--- Launch creation of report for clustering ---")

    df_proteome_tax_id = pd.read_csv(os.path.join(input_folder, 'proteome_tax_id.tsv'), header=0, index_col='observation_name', sep='\t')
    computed_threshold_folder = os.path.join(input_folder, 'computed_threshold')
    df_clustering = data_proteome_representativeness(df_proteome_tax_id, computed_threshold_folder)

    metadata_clustering = json.load(open(os.path.join(input_folder, 'esmecata_metadata_clustering.json'), 'r'))

    fig12 = create_proteome_representativeness_lineplot_px(df_clustering,
        metadata_clustering["tool_options"]["clust_threshold"],
        output_folder)

    fig12_details = proteomes_representativeness_details(df_clustering,
        metadata_clustering["tool_options"]["clust_threshold"],
        output_folder)
