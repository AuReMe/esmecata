#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Uses esmecata sunburst on lists of EC, one per taxonomic rank.
The output of the function `create_dataset_annotation_file()` is 
needed (in its default path)
"""

import os
import json
import argparse
import pandas as pd
from stats_workflow_figures import ec_sunburst_per_model

__author__ = "Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

parser = argparse.ArgumentParser()
parser.add_argument("-d","--directory",
    required=True, 
    help="EsMeCaTa output folder (must contains the output of create_dataset_annotation_file() in the '3_analysis' directory)")
args = parser.parse_args()

# Get EC data of biogas dataset and taxonomic lineages
df_ec = pd.read_table(os.path.join(args.directory, "3_analysis/dataset_annotation_ec.tsv"),
    header=0,
    index_col=0,
    sep='\t')

# Used to group inputs with the same esmecata output (same rank = same annotation)
taxo = pd.read_table(os.path.join(args.directory, "2_annotation/proteome_tax_id.tsv"),
    header=0,
    index_col=0,
    sep='\t')

if not os.path.exists(os.path.join(args.directory, "3_analysis/annotation_figures/ec_sunburst_per_model")):
    os.mkdir(os.path.join(args.directory, "3_analysis/annotation_figures/ec_sunburst_per_model"))

# observations_and_groups = {}

for taxid in taxo["name"].unique():
    group = taxo[taxo["name"] == taxid].index
    # observations_and_groups[taxid] = group
    df_ec_reduced = df_ec.loc[group,]
    ecs = df_ec_reduced.loc[:, (df_ec_reduced != 0).all()].columns.tolist()
    fig = ec_sunburst_per_model(ecs, args.directory, taxid, True)
