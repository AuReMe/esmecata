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
import csv

from typing import Dict, List, Tuple, Any
import plotly.graph_objects as go
import pandas as pd

# CONSTANTS
# ======================================================================================================================

# TAXONOMIC RANK RELATED CONSTANTS
# --------------------------------

RANK2COL = {  # KINGDOM (RED)
    'superkingdom': '#6f1a12', 'kingdom': '#872015', 'subkingdom': '#a2261a',
    # PHYLUM (ORANGE)
    'superphylum': '#b26f15', 'phylum': '#c87c14', 'subphylum': '#de8710', 'infraphylum': '#fa960e',
    # CLASS (YELLOW)
    'superclass': '#c3b21f', 'class': '#d1be14', 'subclass': '#e2cd10', 'infraclass': '#f0db1a',
    # COHORT (KAKHI)
    'cohort': '#a2af21', 'subcohort': '#c1d02e',
    # ORDER (GREEN)
    'superorder': '#699924', 'order': '#7ab12a', 'suborder': '#8ecc33', 'infraorder': '#9dde3e', 'parvorder': '#a9e94c',
    # FAMILY (TURQUOISE)
    'superfamily': '#29b28d', 'family': '#29c198', 'subfamily': '#36d5aa',
    # TRIBE (TURQUOISE - BLUE)
    'tribe': '#3097a4', 'subtribe': '#3eb5c4',
    # GENUS (BLUE)
    'genus': '#2a51ba', 'subgenus': '#3863d7',
    # SECTION (PURPLE)
    'section': '#3a31c6', 'subsection': '#554ce5',
    # SERIES (MAGENTA)
    'series': '#7838c4', 'subseries': '#904ae3',
    # SPECIES (PINK)
    'species group': '#883c91', 'species subgroup': '#9c3ea7', 'species': '#af41bb', 'forma specialis': '#c356cf',
    'subspecies': '#d66be2',
    # (PINK - RED)
    'varietas': '#7e2d48', 'subvariety': '#903855', 'forma': '#a44564', 'serogroup': '#bd5979', 'serotype': '#c96685',
    'strain': '#d57392', 'isolate': '#e085a2',
    # NO RANK (GREY)
    'clade': '#585858', 'environmental samples': '#686868', 'incertae sedis': '#828181', 'unclassified': '#939393',
    'no rank': '#acacac'}

RANK_SORTED = ['isolate', 'strain', 'serotype', 'serogroup', 'forma', 'subvariety', 'varietas',
               'subspecies', 'forma specialis', 'species', 'species subgroup', 'species group',
               'subseries', 'series',
               'subsection', 'section',
               'subgenus', 'genus',
               'subtribe', 'tribe',
               'subfamily', 'family', 'superfamily',
               'parvorder', 'infraorder', 'suborder', 'order', 'superorder',
               'subcohort', 'cohort',
               'infraclass', 'subclass', 'class', 'superclass',
               'infraphylum', 'subphylum', 'phylum', 'superphylum',
               'subkingdom', 'kingdom', 'superkingdom',
               'clade', 'environmental samples', 'incertae sedis', 'unclassified', 'no rank']

# FIGURE COMPOSITION CONSTANTS
# ----------------------------

DEFAULT_NODE_COLOR = 'rgba(191, 191, 191, 0.3)'
DEFAULT_LINK_COLOR = 'rgba(191, 191, 191, 0.3)'
CHANGE_LINK_COLOR = 'rgba(47, 133, 129, 0.3)'

# DATA CONSTANT
# -------------

# INPUT_DATA KEYS
IN_LST = 'input_lst'
OUT_LST = 'output_lst'
VAL_D = 'val_d'
RANK_D = 'rank_d'

# FIG_DATA KEYS
LABELS = 'labels'
NODE_COLORS = 'node_colors'
RANKS = 'ranks'
SOURCES = 'sources'
TARGETS = 'targets'
VALUES = 'values'
LINK_COLORS = 'link_colors'
TEXT_ANNOT = 'text_annot'
X_POS = 'x_pos'
Y_POS = 'y_pos'

# SUFFIXES
IN_OTU_SUFFIX = ' Organisms'
OUT_TAXA_SUFFIX = ' '


# FUNCTIONS
# ======================================================================================================================
def sort_tax_by_ranks(input_data: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """ Sort the 'input_lst' and 'output_list' lists of 'input_data' according to taxonomic rank : more specific to
    less specific (~ species --> kingdom)

    Parameters
    ----------
    input_data: Dict[str, Any]
        Dictionary containing input data information

    Returns
    -------
    List[str]
        List of input taxa sorted
    List[str]
        List of output taxa sorted
    """
    input_rvalue = [RANK_SORTED.index(input_data[RANK_D][s]) for s in input_data[IN_LST]]
    output_rvalue = [RANK_SORTED.index(input_data[RANK_D][s]) for s in input_data[OUT_LST]]
    df = pd.DataFrame([input_data[IN_LST], input_data[OUT_LST], input_rvalue, output_rvalue],
                      index=['i_tax', 'o_tax', 'i_val', 'o_val']).T
    df = df.sort_values(by=['i_val', 'o_val', 'o_tax'])
    return list(df['i_tax']), list(df['o_tax'])


def get_input_data(input_file: str) -> Dict[str, Any]:
    """ Read taxonomy_diff.tsv file to extract input and output taxa from a run.

    Parameters
    ----------
    input_file: str
        taxonomy_diff.tsv file from proteome Esmecata step.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing input data information
    """
    input_data = {IN_LST: [], OUT_LST: [], VAL_D: {}, RANK_D: {}}
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for line in reader:
            input_data[IN_LST].append(line['Input Name'])
            input_data[OUT_LST].append(line['Esmecata Name'] + OUT_TAXA_SUFFIX)
            input_data[VAL_D][line['Input Name']] = len(line['Observation Names'].split(';'))
            input_data[RANK_D][line['Input Name']] = line['Input Rank']
            input_data[RANK_D][line['Input Name'] + IN_OTU_SUFFIX] = line['Input Rank']
            input_data[RANK_D][line['Esmecata Name'] + OUT_TAXA_SUFFIX] = line['Esmecata Rank']
    input_lst, output_lst = sort_tax_by_ranks(input_data)
    input_data[IN_LST] = input_lst
    input_data[OUT_LST] = output_lst
    return input_data


def get_fig_parameters(input_data: Dict[str, Any]) -> Dict[str, List[Any]]:
    """ Generate parameters for the sankey figure from input data.

    Parameters
    ----------
    input_data: Dict[str, Any]
        Dictionary containing input data information

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information
    """
    # Init data
    output_set_lst = list(set(input_data[OUT_LST]))
    input_otu = [x + IN_OTU_SUFFIX for x in input_data[IN_LST]]
    output_str = 'Output Diversity'
    input_str = '  '
    labels = input_otu + input_data[IN_LST] + output_set_lst + [str(output_str), str(input_str)]
    nb_otu_input = sum(input_data[VAL_D].values())
    nb_input_div = len(set(input_data[IN_LST]))
    nb_output_div = len(output_set_lst)
    fig_data = {LABELS: labels,
                NODE_COLORS: [RANK2COL[input_data[RANK_D][sp]] if sp in input_data[RANK_D] else DEFAULT_NODE_COLOR
                              for sp in labels],
                RANKS: [input_data[RANK_D][sp] if sp in input_data[RANK_D] else '' for sp in labels],
                SOURCES: [],
                TARGETS: [],
                VALUES: [],
                LINK_COLORS: [],
                TEXT_ANNOT: [f'Input : {nb_otu_input} Organisms',
                             f'Input : {nb_input_div} taxonomic diversity',
                             f'Output : {nb_output_div} taxonomic diversity'],
                X_POS: [0.01] * len(labels),
                Y_POS: [0.0] * len(labels)}
    ii_y = 0.001
    i_y = 0.001
    o_y = 0.001

    # Set Output node position
    fig_data[X_POS][labels.index(output_str)] = 0.99
    fig_data[Y_POS][labels.index(output_str)] = 0.5
    # Make links
    for sp_i in input_data[IN_LST]:
        sp_o = input_data[OUT_LST][input_data[IN_LST].index(sp_i)]
        sp_ii = input_otu[input_data[IN_LST].index(sp_i)]
        # Link Input node (input_str) to Input OTU (sp_ii)
        fig_data, ii_y = link_input_node_to_input_otu(input_data, fig_data, labels, input_str, sp_i, sp_ii,
                                                      nb_otu_input, ii_y)
        # Link Input OTU (sp_ii) to Input Taxa (sp_i)
        fig_data, i_y = link_input_otu_to_input_taxa(fig_data, labels, sp_i, sp_ii, nb_input_div, i_y)
        # Link Input taxa (sp_i) to Output taxa (sp_o) (+ Link Output taxa (sp_o) to Output node (output_str))
        if sp_i != sp_o:
            fig_data, o_y = link_input_taxa_to_output_taxa(input_data, fig_data, labels, sp_i, sp_o, nb_input_div,
                                                           output_str, o_y)
    return fig_data


def set_node_position(node: str, step: float, x_pos: float, increment: float, fig_data: Dict[str, List[Any]],
                      labels: List[str]) -> Tuple[Dict[str, List[Any]], float]:
    """ Set the x and y position of a node in the figure.

    Parameters
    ----------
    node: str
        Label of the node
    step: float
        Proportion of the node height for incrementation of y-axis position
    x_pos: float
        The x-axis position
    increment: float
        Current y-axis position of the iteration
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    labels: List[str]
        List of nodes labels

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current node position)
    float
        New y-axis position of the iteration
    """
    fig_data[X_POS][labels.index(node)] = x_pos
    fig_data[Y_POS][labels.index(node)] = increment + (step / 2)
    increment += step
    return fig_data, increment


def link_input_node_to_input_otu(input_data: Dict[str, Any], fig_data: Dict[str, List[Any]], labels: List[str],
                                 input_str: str, sp_i: str, sp_ii: str, nb_otu_input: int, ii_y: float)\
        -> Tuple[Dict[str, List[Any]], float]:
    """ Link the input node to an input OTU node (col 1 to col 2).

    Parameters
    ----------
    input_data: Dict[str, Any]
        Dictionary containing input data information
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    labels: List[str]
        List of nodes labels
    input_str: str
        Label of the input node
    sp_i: str
        Label of the input taxa node
    sp_ii: str
        Label of the input OTU node
    nb_otu_input: int
        Total number of OTUs in the input
    ii_y: float
        Previous y-axis position of the input OTU (sp_ii) node

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current link information)
    float
        New input OTUs y-axis position
    """
    source = labels.index(input_str)
    target = labels.index(sp_ii)
    value = input_data[VAL_D][sp_i]
    fig_data = add_link_data(fig_data, source, target, value, DEFAULT_LINK_COLOR)
    # Set sp_ii position
    step = input_data[VAL_D][sp_i] / nb_otu_input
    return set_node_position(sp_ii, step, 0.2, ii_y, fig_data, labels)


def link_input_otu_to_input_taxa(fig_data: Dict[str, List[Any]], labels: List[str], sp_i: str, sp_ii: str,
                                 nb_input_div: int, i_y: float) -> Tuple[Dict[str, List[Any]], float]:
    """ Link the input OTU nodes do the input taxa nodes (col 2 to col 3).

    Parameters
    ----------
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    labels: List[str]
        List of nodes labels
    sp_i: str
        Label of the input taxa node
    sp_ii: str
        Label of the input OTU node
    nb_input_div: int
        Total number of taxa in the input
    i_y: float
        Previous y-axis position of the input taxa (sp_i) node

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current link information)
    float
        New input taxa y-axis position
    """
    source = labels.index(sp_ii)
    target = labels.index(sp_i)
    fig_data = add_link_data(fig_data, source, target, 1, DEFAULT_LINK_COLOR)
    # Set sp_i position
    step = 1 / nb_input_div
    return set_node_position(sp_i, step, 0.5, i_y, fig_data, labels)


def link_input_taxa_to_output_taxa(input_data: Dict[str, Any], fig_data: Dict[str, List[Any]], labels: List[str],
                                   sp_i: str, sp_o: str, nb_input_div: int, output_str: str, o_y: float)\
        -> Tuple[Dict[str, List[Any]], float]:
    """ Link the input taxa nodes to output taxa nodes (col 3 to col 4).

    Parameters
    ----------
    input_data: Dict[str, Any]
        Dictionary containing input data information
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    labels: List[str]
        List of nodes labels
    sp_i: str
        Label of the input taxa node
    sp_o: str
        Label of the output taxa node
    nb_input_div: int
        Total number of input taxonomic diversity
    output_str: str
        Label of the output node
    o_y: float
        Previous y-axis position of the output taxa (sp_o) node

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current link information)
    float
        New output taxa y-axis position

    """
    source = labels.index(sp_i)
    target = labels.index(sp_o)
    value = 1
    if sp_i == sp_o[:-len(OUT_TAXA_SUFFIX)]:
        fig_data = add_link_data(fig_data, source, target, value, DEFAULT_LINK_COLOR)
    else:
        fig_data = add_link_data(fig_data, source, target, value, CHANGE_LINK_COLOR)

    if fig_data[Y_POS][labels.index(sp_o)] == 0.0:
        # Set sp_o position
        step = input_data[OUT_LST].count(sp_o) / nb_input_div
        fig_data, o_y = set_node_position(sp_o, step, 0.8, o_y, fig_data, labels)
        # Link Output taxa (sp_o) to Output node (output_str)
        fig_data = link_output_taxa_to_output_node(fig_data, labels, sp_o, output_str)
    return fig_data, o_y


def link_output_taxa_to_output_node(fig_data: Dict[str, List[Any]], labels: List[str], sp_o: str, output_str: str)\
        -> Dict[str, List[Any]]:
    """ Link the output taxa nodes to the output node (col 4 to col 5).

    Parameters
    ----------
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    labels: List[str]
        List of nodes labels
    sp_o: str
        Label of the output taxa node
    output_str: str
        Label of the output node

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current link information)
    """
    source = labels.index(sp_o)
    target = labels.index(output_str)
    return add_link_data(fig_data, source, target, 1, DEFAULT_LINK_COLOR)


def add_link_data(fig_data: Dict[str, List[Any]], source: int, target: int, value: int, color: str):
    """ Add a link parameters to fig_data information.

    Parameters
    ----------
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    source: int
        Index of source label
    target: int
        Index of target label
    value: int
        Value of link (flow / size / height)
    color: str
        Color of the link

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing figure parameters information (+ current link information)
    """
    fig_data[SOURCES].append(source)
    fig_data[TARGETS].append(target)
    fig_data[VALUES].append(value)
    fig_data[LINK_COLORS].append(color)
    return fig_data


def draw_sankey_fig(fig_data: Dict[str, List[Any]], output: str, show_fig: bool):
    """ Draw the Sankey graph showing compression between input and output taxa of an Esmecata run
    from the fig parameters.

    Parameters
    ----------
    fig_data: Dict[str, List[Any]]
        Dictionary containing figure parameters information
    output: str
        Name of the output file
    show_fig: bool
        True to show figure, False to only save in html
    """
    line_width = 1
    if len(fig_data[LABELS]) > 1000:
        line_width = 0    

    # opacity = 0.4
    # fig_data[""][0]["link"] = [data['data'][0]['node']['color'][src].replace("0.8", str(opacity))
    #                                 for src in data['data'][0]['link']['source']]

    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node={'line': dict(color='black', width=line_width),
              'label': fig_data[LABELS],
              'thickness': 40,
              'pad': 2,
              'color': fig_data[NODE_COLORS],
              'customdata': fig_data[RANKS],
              'hovertemplate': '<b>Name :</b> %{label}<br>'
                               '<b>Rank :</b> %{customdata}<br>'
                               '<extra><b>%{value:,.0f}</b><br></extra>',
              'x': fig_data[X_POS],
              'y': fig_data[Y_POS]},

        link={'source': fig_data[SOURCES],
              'target': fig_data[TARGETS],
              'value': fig_data[VALUES],
              'color': fig_data[LINK_COLORS],
              'line': dict(color='white', width=line_width)}),
    )

    text_pos = [0.2, 0.5, 0.8]
    for i in range(0, 3):
        fig.add_annotation(dict(font=dict(color='darkslategrey', size=18), x=text_pos[i], y=1.05,
                                showarrow=False, text=fig_data[TEXT_ANNOT][i], textangle=0,
                                xanchor='center', xref='paper', yref='paper'))
    # fig.update_layout(height=fig_height)
    # fig.write_html(output)
    if show_fig:
        fig.show()
    return fig


# MAIN FUNCTION
# ======================================================================================================================
def esmecata_compression_taxonomy_file(input_file, output, show_fig):
    """ Draw a Sankey graph showing compression and taxonomic ranks rearrangement between input and output taxa of an
    Esmecata run.

    Parameters
    ----------
    input_file: str
        Path to the taxonomy_diff file in proteomes folder
    output: str
        Name of the output file
    show_fig: bool (default=False)
        True to show figure, False to only save in html
    """
    input_data = get_input_data(input_file)
    fig_data = get_fig_parameters(input_data)
    
    fig_height = 10*len(fig_data[LABELS]) # bug report : Datapane does not uses height correctly if coded inside draw_sankey_fig. Code here instead works (why ?)
    if fig_height < 600:
        fig_height = 600
    
    fig = draw_sankey_fig(fig_data, f'{output}.html', show_fig)
    fig.update_layout(height=fig_height, 
                      font=dict(family="Courier New, monospace",
                                size=14,
                                color="black"))
    fig.write_html(f'{output}.html')

    return fig


def esmecata_compression(run_name: str, output: str, show_fig: bool = False):
    """ Draw a Sankey graph showing compression and taxonomic ranks rearrangement between input and output taxa of an
    Esmecata run.

    Parameters
    ----------
    run_name: str
        Name of the folder of Esmecata results
    output: str
        Name of the output file
    show_fig: bool (default=False)
        True to show figure, False to only save in html
    """
    tax_diff_file = os.path.join(run_name, '0_proteomes', 'taxonomy_diff.tsv')
    fig = esmecata_compression_taxonomy_file(tax_diff_file, output, show_fig)
    return fig