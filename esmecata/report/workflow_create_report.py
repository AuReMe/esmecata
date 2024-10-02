# Copyright (C) 2023-2024  Victor Mataigne and Pauline Hamon-Giraud - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
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

import os
import json
import logging
import pandas as pd

from esmecata.core.annotation import create_dataset_annotation_file
from esmecata.report.report_creation import create_esmecata_report
from esmecata.report.esmecata2taxontology import get_fig_parameters, generate_sunburst_fig
from esmecata.report.stats_workflow_figures import distributions_by_ranks, compare_ranks_in_out, _format_taxo, create_proteome_representativeness_lineplot_px, \
                                                    data_proteome_representativeness, proteomes_representativeness_details, create_annot_obs_df,annot_frequencies_in_obs, \
                                                    annot_frequencies_in_obs_hist, ec_sunburst
from esmecata.report.esmecata_compression import esmecata_compression_taxonomy_file

logger = logging.getLogger(__name__)


def run_create_workflow_report(input_file, input_folder, output_folder, create_svg=False):
    """ Create report from esmecata workflow output.

    Args:
        input_file (str): pathname to the input file of esmecata.
        input_folder (str): pathname to the output folder of esmecata worklow.
        output_folder (str): output folder path
        create_svg (bool): boolean to create svg files of the figure
    """
    logger.info("--- Launch creation of report for workflow ---")
    create_esmecata_report(input_file, input_folder, output_folder, create_svg)


def run_proteomes_report_creation(input_folder, output_folder, create_svg=False):
    """ Create report from esmecata proteomes output.

    Args:
        input_folder (str): pathname to the output folder of esmecata proteomes.
        output_folder (str): output folder path
        create_svg (bool): boolean to create svg files of the figure
    """
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


def run_clustering_report_creation(input_folder, output_folder, create_svg=False):
    """ Create report from esmecata clustering output.

    Args:
        input_folder (str): pathname to the output folder of esmecata clustering.
        output_folder (str): output folder path
        create_svg (bool): boolean to create svg files of the figure
    """
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


def run_annotation_report_creation(input_folder, output_folder, create_svg=False):
    """ Create report from esmecata annotation output.

    Args:
        input_folder (str): pathname to the output folder of esmecata annotation.
        output_folder (str): output folder path
        create_svg (bool): boolean to create svg files of the figure
    """
    annotation_reference_folder = os.path.join(input_folder, 'annotation_reference')
    dataset_annotation_ec_file = os.path.join(output_folder, 'dataset_annotation_ec.tsv')
    dataset_annotation_go_file = os.path.join(output_folder, 'dataset_annotation_go.tsv')

    _ = create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_ec_file, 'EC')

    _ = create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_go_file, 'GO')

    data_ec_number = create_annot_obs_df(dataset_annotation_ec_file, output_folder, 'EC numbers')

    data_go_term = create_annot_obs_df(dataset_annotation_go_file, output_folder, "GO terms")

    fig5 = annot_frequencies_in_obs(data_ec_number["df_annot_frequencies"], output_folder, "EC numbers")
    #fig6 = fraction_of_all_annot_in_obs(data_ec_number["df_fractionin_obs"], DATA["DF_STATS"], output_folder, "EC numbers")
    fig7 = annot_frequencies_in_obs_hist(data_ec_number["df_annot_frequencies"], output_folder, "EC numbers")
    #fig8 = fraction_of_all_annot_in_obs_hist(data_ec_number["df_fractionin_obs"], DATA["DF_STATS"], output_folder, "EC numbers")

    fig5b = annot_frequencies_in_obs(data_go_term["df_annot_frequencies"], output_folder, "GO terms")
    #fig6b = fraction_of_all_annot_in_obs(data_go_term["df_fractionin_obs"], DATA["DF_STATS"], output_folder, "GO terms")
    fig7b = annot_frequencies_in_obs_hist(data_go_term["df_annot_frequencies"], output_folder, "GO terms")
    #fig8b = fraction_of_all_annot_in_obs_hist(data_go_term["df_fractionin_obs"], DATA["DF_STATS"], output_folder, "GO terms")

    fig11 = ec_sunburst(data_ec_number["df_annot_frequencies"].index, output_folder)