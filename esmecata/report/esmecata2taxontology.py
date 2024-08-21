# Copyright (C) 2023-2024 Victor Mataigne and Pauline Hamon-Giraud - Inria, Univ Rennes, CNRS, IRISA Dyliss
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

import csv
import os.path

import plotly.graph_objects as go
from ete3 import NCBITaxa
from typing import Dict, List, Set, Tuple

from esmecata.report.esmecata_compression import RANK2COL

# CONSTANTS ==========================================================================================================

NCBI = NCBITaxa()


# FUNCTIONS ==========================================================================================================
def get_fig_parameters(tax_diff_file: str) -> Dict[str, List]:
    """ Fill a dictionary with input parameters for the sunburst figure creation.

    Parameters
    ----------
    tax_diff_file: str
        Path of "taxonomy_diff.tsv" file from proteome step.

    Returns
    -------
    Dict[str, List]
        Dictionary of parameters for sunburst figure. Multiple lists of size N with following information for each
         taxon i in [1:n]:
            - ID: tax ID
            - Name: taxon name
            - Rank: taxonomic rank
            - Parent: tax ID of the parent taxon
            - Count Input: Number of OTU corresponding to the taxon in input data
            - Count Esmecata: Number of OTU corresponding to the taxon after esmecata proteome selection
            - Color: Hexadecimal color value corresponding of the taxonomic rank
            - Shape: Shape of the sunburst cell, hatched if the taxon is not kept after esmecata filter, plain otherwise
    """
    count_input = dict()
    count_esmecata = dict()
    tax_names = dict()
    tax_ranks = dict()
    parent = dict()
    with open(tax_diff_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        reader.__next__()
        for line in reader:
            count, input_tax_id, esmecata_tax_id = len(line[6].split(';')), int(line[4]), int(line[5])
            lineage_ids = NCBI.get_lineage(input_tax_id)
            count_input = fill_count_input(count_input, lineage_ids, count)
            tax_names = fill_tax_names(tax_names, lineage_ids)
            parent = fill_parent(parent, lineage_ids)
            tax_ranks = fill_tax_ranks(tax_ranks, lineage_ids)
            count_esmecata = fill_count_esmecata(count_esmecata, count, esmecata_tax_id)

    to_exclude, parent = select_root(count_input, parent)
    return fill_parameters(tax_names, tax_ranks, parent, count_input, count_esmecata, to_exclude)


def fill_count_input(count_input: Dict[int, int], lineage_ids: List[int], count: int) -> Dict[int, int]:
    """ Fill the count_input dictionary for Count input fig parameter for one taxon.

    Parameters
    ----------
    count_input: Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon.
    lineage_ids: List[int]
        List of tax id representing the lineage of the taxon (root to taxon order)
    count: int
        Number of OTU corresponding to the taxon.

    Returns
    -------
    Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon.
    """
    for tax_id in lineage_ids:
        if tax_id not in count_input.keys():
            count_input[tax_id] = count
        else:
            count_input[tax_id] += count
    return count_input


def fill_tax_names(tax_names: Dict[int, str], lineage_ids: List[int]) -> Dict[int, str]:
    """ Fill the tax_names dictionary for Name fig parameter for one taxon.

    Parameters
    ----------
    tax_names: Dict[int, str]
        Dictionary associating to each tax id, its taxon name.
    lineage_ids: List[int]
        List of tax id representing the lineage of the taxon (root to taxon order)

    Returns
    -------
    Dict[int, str]
        Dictionary associating to each tax id, its taxon name.
    """
    lineage_names = NCBI.translate_to_names(lineage_ids)
    for i in range(len(lineage_ids)):
        tax_names[lineage_ids[i]] = lineage_names[i]
    return tax_names


def fill_parent(parent: Dict[int, int or str], lineage_ids: List[int]) -> Dict[int, int or str]:
    """ Fill the parent dictionary for Parent fig parameter for one taxon.

    Parameters
    ----------
    parent: Dict[int, int or str]
        Dictionary associating to each tax id, its parent taxon ID.
    lineage_ids: List[int]
        List of tax id representing the lineage of the taxon (root to taxon order)

    Returns
    -------
    Dict[int, int or str]
        Dictionary associating to each tax id, its parent taxon ID ('' if root).
    """
    for i in range(len(lineage_ids)):
        if lineage_ids[i] != 1:
            parent[lineage_ids[i]] = lineage_ids[i - 1]
        else:
            parent[lineage_ids[i]] = ''
    return parent


def fill_tax_ranks(tax_ranks: Dict[int, str], lineage_ids: List[int]) -> Dict[int, str]:
    """ Fill the tax_ranks dictionary for Rank fig parameter for one taxon.

    Parameters
    ----------
    tax_ranks: Dict[int, str]
        Dictionary associating to each tax id, its taxonomic rank.
    lineage_ids: List[int]
        List of tax id representing the lineage of the taxon (root to taxon order)

    Returns
    -------
    Dict[int, str]
        Dictionary associating to each tax id, its taxonomic rank.
    """
    lineage_ranks = NCBI.get_rank(lineage_ids)
    return {**tax_ranks, **lineage_ranks}


def fill_count_esmecata(count_esmecata: Dict[int, int], count: int, esmecata_id: int) -> Dict[int, int]:
    """ Fill the count_esmecata dictionary for "Count esmecata" fig parameter for one esmecata taxon kept.

    Parameters
    ----------
    count_esmecata: Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon for taxon kept.
    count: int
        Number of OTU corresponding to the taxon.
    esmecata_id
        Tax ID of the esmecata taxon kept.
    Returns
    -------
    Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon for taxon kept.

    """
    lineage_ids = NCBI.get_lineage(esmecata_id)
    for tax_id in lineage_ids:
        if tax_id not in count_esmecata.keys():
            count_esmecata[tax_id] = count
        else:
            count_esmecata[tax_id] += count
    return count_esmecata


def select_root(count_input: Dict[int, int], parent: Dict[int, int or str]) -> Tuple[Set[int], Dict[int, int]]:
    """ Select the root to use and return tax IDs to exclude from the figure. Keep as root the tax ID shared by all
    taxa having the lowest taxonomic rank.

    Parameters
    ----------
    count_input: Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon.
    parent: Dict[int, int or str]
        Dictionary associating to each tax id, its parent taxon ID.

    Returns
    -------
    Tuple[Set[int], Dict[int, int]]
        Set of tax IDs to exclude and Parent dictionary modified for the root ID.
    """
    max_count = max(count_input.values())
    roots = [x for x, y in count_input.items() if y == max_count]
    root = 1
    for ch_id in parent.values():
        if ch_id != '':
            if parent[ch_id] in roots and ch_id not in roots:
                root = parent[ch_id]
    parent[root] = ''

    return set(roots).difference({root}), parent


def fill_parameters(tax_names: Dict[int, str], tax_ranks: Dict[int, str], parent: Dict[int, int or str],
                    count_input: Dict[int, int], count_esmecata: Dict[int, int], to_exclude: Set[int]) \
        -> Dict[str, List]:
    """ Fill a dictionary with parameters for sunburst figure.

    Parameters
    ----------
    tax_names: Dict[int, str]
        Dictionary associating to each tax id, its taxon name.
    tax_ranks: Dict[int, str]
        Dictionary associating to each tax id, its taxonomic rank.
    parent: Dict[int, int or str]
        Dictionary associating to each tax id, its parent taxon ID.
    count_input: Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon.
    count_esmecata: Dict[int, int]
        Dictionary associating to each tax id, the number of OTU associated to the taxon for taxon kept.
    to_exclude: Set[int]
        Set of tax IDs to exclude for root selection.

    Returns
    -------
    Dict[str, List]:
        Dictionary of parameters for sunburst figure. Multiple lists of size N with following information for each
        taxon i in [1:n]:
            - ID: tax ID
            - Name: taxon name
            - Rank: taxonomic rank
            - Parent: tax ID of the parent taxon
            - Count Input: Number of OTU corresponding to the taxon in input data
            - Count Esmecata: Number of OTU corresponding to the taxon after esmecata proteome selection
            - Color: Hexadecimal color value corresponding of the taxonomic rank
            - Shape: Shape of the sunburst cell, hatched if the taxon is not kept after esmecata filter, plain otherwise
    """
    data = {'ID': list(),
            'Name': list(),
            'Rank': list(),
            'Parent': list(),
            'Count Input': list(),
            'Count Esmecata': list(),
            'Color': list(),
            'Shape': list()}
    for tax_id in tax_names.keys():
        data['ID'].append(tax_id)
        data['Name'].append(tax_names[tax_id])
        data['Rank'].append(tax_ranks[tax_id])
        data['Parent'].append(parent[tax_id])
        data['Count Input'].append(count_input[tax_id])
        # data['Color'].append(RANK2COL[tax_ranks[tax_id]])
        if tax_id in count_esmecata.keys():
            data['Count Esmecata'].append(count_esmecata[tax_id])
            data['Color'].append(RANK2COL[tax_ranks[tax_id]])
            #data['Shape'].append('')
        else:
            data['Count Esmecata'].append(0)
            data['Color'].append('#ffffff')
            #data['Shape'].append('/')

    return data


def generate_sunburst_fig(data: Dict[str, List[str or int]], output: str):
    """ Generate a Sunburst figure and save it to output path.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary of parameters for sunburst figure. Multiple lists of size N with following information for each
         taxon i in [1:n]:
            - ID: tax ID
            - Name: taxon name
            - Rank: taxonomic rank
            - Parent: tax ID of the parent taxon
            - Count Input: Number of OTU corresponding to the taxon in input data
            - Count Esmecata: Number of OTU corresponding to the taxon after esmecata proteome selection
            - Color: Hexadecimal color value corresponding of the taxonomic rank
            - Shape: Shape of the sunburst cell, hatched if the taxon is not kept after esmecata filter, plain otherwise
    output: str
        Path to output to save the figure without extension
    """
    # Create Sunburst
    fig = go.Figure(go.Sunburst(labels=data['Name'], parents=data['Parent'], values=data['Count Input'], ids=data['ID'],
                                branchvalues='total', hovertext=[f'Rank: {x}' for x in data['Rank']], opacity = 0.8,
                                hoverinfo='label+value+percent entry+text', maxdepth=10,
                                insidetextfont=dict(color='#000000', size=16), #cccccc
                                marker=dict(colors=data['Color'],
                                            line=dict(color='#000000'),                                            
                                            pattern=dict(shape=data['Shape'])
                                            )
                                )
                    )

    # Add ranks color annotation
    m = len(set(data['Rank']))
    i = 0
    for rank, color in RANK2COL.items():
        if rank in data['Rank']:
            i += 1
            fig.add_annotation(dict(font=dict(color='#ffffff', size=16), x=0, y=1 - (i / m), showarrow=False, text=rank,
                                    textangle=0, xanchor='left', xref="paper", yref="paper",
                                    bgcolor=color,
                                    opacity=0.8))

    # Update html page layout
    fig.update_layout(paper_bgcolor="#ffffff", #1b1b1b
                      font_color='#ffffff',
                      # width=2000,
                      height=1300,
                      font=dict(family="Courier New, monospace",size=14,color="black"),
                      margin=go.layout.Margin(
                        l=0, #left margin
                        r=0, #right margin
                        b=0, #bottom margin
                        t=0) #top margin
    )
    fig.write_html(f'{output}.html')

    return fig


# MAIN FUNCTION =====================================================================================================

def esmecata2taxonomy(run_name: str, output: str):
    """ Generate a Sunburst representing the taxonomic variety of OTU of a run. Compare ranks from input data to ranks
    kept by the proteome selection.

    Parameters
    ----------
    run_name: str
    output: str
    """
    tax_comp_file = os.path.join(run_name, '0_proteomes', 'taxonomy_diff.tsv')
    data = get_fig_parameters(tax_comp_file)
    fig = generate_sunburst_fig(data, output)
    
    return fig
