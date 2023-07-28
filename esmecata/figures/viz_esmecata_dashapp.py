#!/usr/bin/env python

from os import path
import json
import argparse

import dash
import dash_bootstrap_components as dbc

import stats_workflow_figures as swf
from esmecata2taxontology import esmecata2taxonomy

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', 
                    action='store', 
                    type=argparse.FileType('r'), 
                    required=True, 
                    help='EsMeCaTa input table')
parser.add_argument('-o','--outdir', 
                    action='store', 
                    required=True, 
                    help='Path to EsMeCaTa output directory')
args = parser.parse_args()


INPUT_DATA, DISCARDED, N_DISCARDED, DF_STATS, N_IN, N_OUT, PROTEOME_TAX_ID, ASSOCIATION_TAXON_TAX_ID = swf.post_analysis_config(args.input, args.outdir)
RANK = 'Phylum'

CSS = {'main-title': {'textAlign': 'center', 'color': '#7FDBFF'},
       'headers': {'textAlign': 'left', 'color': '#7FDBFF'},
       'div-blocks': {'width': '49%', 'display': 'inline-block'}}

CONFIG = {'modeBarButtonsToRemove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

def reproducibility_tokens(outdir):
    '''
    Loads json metadata of all the EsmeCata workflow and converts it to a string.

        Parameters:
            outdir (str) : path to the results of an EsMeCaTa workflow

        Returns:
            metadata (str) : all metadata of the workflow
    '''
    reprometa_proteomes = json.load(open(path.join(outdir, '0_proteomes', 'esmecata_metadata_proteomes.json'), 'r'))
    reprometa_clustering = json.load(open(path.join(outdir, '1_clustering', 'esmecata_metadata_clustering.json'), 'r'))
    reprometa_annotation = json.load(open(path.join(outdir, '2_annotation', 'esmecata_metadata_annotation.json'), 'r'))

    metadata = {"proteomes": reprometa_proteomes,
        "clustering": reprometa_clustering,
        "annotation": reprometa_annotation}
    
    metadata = json.dumps(metadata, indent=4)

    return metadata

if __name__ == '__main__':
    # app = dash.Dash(__name__)
    app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

    # TODO : small texts explaining figures

    fig1 = swf.distributions_by_ranks(DF_STATS, args.outdir, RANK)
    fig2 = swf.n_prot_ec_go_correlations(DF_STATS, args.outdir)
    fig3 = swf.taxo_ranks_contribution(PROTEOME_TAX_ID, args.outdir)
    fig4 = swf.compare_ranks_in_out(PROTEOME_TAX_ID, ASSOCIATION_TAXON_TAX_ID, args.outdir)
    fig5 = esmecata2taxonomy(args.outdir, path.join(args.outdir, 'esmecata2taxonomy'))

    metadata = reproducibility_tokens(args.outdir)

    app.layout = dash.html.Div(children=[
        dash.html.H1(
            children=f'EsMeCaTa - Graphical summary',
            style=CSS['main-title']
        ),

        dash.html.Div(children=[
            dash.html.H2(
                children='Distribution of the number of proteomes found',
                style=CSS['headers']),
            dash.dcc.Graph(
                id='graph1', 
                figure=fig1,
                config=CONFIG)],
            style=CSS['div-blocks']),

        dash.html.Div(children=[
            dash.html.H2(
                children='Correlation between the number of EC, GO terms, and proteomes',
                style=CSS['headers']),
            dash.dcc.Graph(
                id='graph2', 
                figure=fig2,
                config=CONFIG)],
            style=CSS['div-blocks']),

        dash.html.Div(children=[
            dash.html.H2(
                children=f'Proteomes per taxonomic rank (up to the {RANK} limit).',
                style=CSS["headers"]),
            dash.dcc.Graph(
                id='graph3', 
                figure=fig3,
                config=CONFIG)],
            style=CSS['div-blocks']),

        dash.html.Div(children=[
            dash.html.H2(
                children='Input and output ranks difference',
                style=CSS["headers"]),
            dash.dcc.Graph(
                id='graph4', 
                figure=fig4, 
                config={'staticPlot': True})],
            style=CSS['div-blocks']),

        dash.html.Div(children=[
            dash.html.H2(
                children='Taxonomic diversity',
                style=CSS["headers"]),
            dash.dcc.Graph(id='graph5', figure=fig5)]),

        dash.html.H2(
            children="Reproducibility metadata",
            style=CSS["headers"]
        ),

        # nb : some really weird behavior here, clean indentation messes up the results on screen
        dash.dcc.Markdown('''Full workflow : 
```json{}
'''.format(metadata))
     ])

app.run_server(debug=True)
