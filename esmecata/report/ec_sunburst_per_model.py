# Copyright (C) 2023-2024 Victor Mataigne - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
""" 
Uses esmecata sunburst on lists of EC, one per taxonomic rank.
The output of the function `create_dataset_annotation_file()` is 
needed (in its default path)
"""

import os
import pandas as pd
from esmecata.report.stats_workflow_figures import ec_sunburst_per_model


def create_sunburst_ec(esmecata_output_folder, output_folder):
    dataset_annotation_ec_file = os.path.join(output_folder, 'dataset_annotation_ec.tsv')
    proteome_tax_id_filepath = os.path.join(esmecata_output_folder, '0_proteomes', 'proteome_tax_id.tsv')
    # Get EC data of biogas dataset and taxonomic lineages
    df_ec = pd.read_table(dataset_annotation_ec_file,
        header=0,
        index_col=0,
        sep='\t')

    # Used to group inputs with the same esmecata output (same rank = same annotation)
    taxo = pd.read_table(proteome_tax_id_filepath,
        header=0,
        index_col=0,
        sep='\t')

    if not os.path.exists(os.path.join(output_folder, 'ec_sunburst_per_model')):
        os.mkdir(os.path.join(output_folder, 'ec_sunburst_per_model'))

    # observations_and_groups = {}

    for taxid in taxo["name"].unique():
        group = taxo[taxo["name"] == taxid].index
        # observations_and_groups[taxid] = group
        df_ec_reduced = df_ec.loc[group,]
        ecs = df_ec_reduced.loc[:, (df_ec_reduced != 0).all()].columns.tolist()
        fig = ec_sunburst_per_model(ecs, output_folder, taxid, True)
