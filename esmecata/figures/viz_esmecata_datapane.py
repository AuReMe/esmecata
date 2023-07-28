#!/usr/bin/env python

from os import path
import json
import argparse
import datapane as dp
import stats_workflow_figures as swf
from esmecata2taxontology import esmecata2taxonomy

__author__ = "Pauline Hamon-Giraud, Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

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

INPUT_DATA, DISCARDED, N_DISCARDED, DF_STATS, N_IN, N_OUT, PROTEOME_TAX_ID, ASSOCIATION_TAXON_TAX_ID = swf.post_analysis_config(args.input, args.outdir)
RANK = 'phylum'

CSS = {'main-title': {'textAlign': 'center', 'color': '#7FDBFF'},
       'headers': {'textAlign': 'left', 'color': '#7FDBFF'},
       'div-blocks': {'width': '49%', 'display': 'inline-block'}}

CONFIG = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

fig1 = swf.distributions_by_ranks(DF_STATS, args.outdir, RANK)
fig2 = swf.n_prot_ec_go_correlations(DF_STATS, args.outdir, RANK)
fig3 = swf.taxo_ranks_contribution(PROTEOME_TAX_ID, args.outdir)
fig4 = swf.compare_ranks_in_out(PROTEOME_TAX_ID, ASSOCIATION_TAXON_TAX_ID, args.outdir)
fig5 = esmecata2taxonomy(args.outdir, path.join(args.outdir, 'esmecata2taxonomy'))

metadata = reproducibility_tokens(args.outdir)

report = dp.Blocks(
    dp.Page(
        title="EsMeCaTa - Graphical summary",
        blocks=[
            dp.Group(
                dp.Plot(fig1.update_layout(modebar=CONFIG), label='Distribution of the number of proteomes found'),
                dp.Plot(fig2.update_layout(modebar=CONFIG), label='Correlation between the number of EC, GO terms, and proteomes'),
                columns=2,
            ),
            dp.Group(
                dp.Plot(fig3.update_layout(modebar=CONFIG), label=f'Proteomes per taxonomic rank (up to the {RANK} limit).'),
                dp.Plot(fig4.update_layout(modebar=CONFIG), label='Input and output ranks difference'),
                columns=2,
            ),
            dp.Plot(fig5, label='Taxonomic diversity')
        ],
    ),
    dp.Page(
        title="Data", blocks=[
            dp.Group(
                dp.BigNumber(heading="Number of inputs", value=N_IN),
                dp.BigNumber(heading="Kept by EsMeCaTa", value=N_OUT),
                dp.BigNumber(heading="Discarded", value=N_DISCARDED),
                columns=3,
            ),
            dp.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete3</p>"),
            dp.DataTable(DF_STATS, 
                        label="Data"),
            dp.HTML("<h2>Discarded</h2><p>Taxonomic ranks were not inferred; only names are displayed</p>"),
            dp.DataTable(DISCARDED, 
                        label="Data")
        ]
    ),
    dp.Page(title="Metadata", blocks=[dp.Code(code=metadata, language="json", label="Metadata")]),
)

dp.save_report(report, 
    path=path.join(args.outdir, 'esmecata_summary.html'),
    formatting=dp.Formatting(
        accent_color="rgb(204, 255, 204)",
        font=dp.FontChoice.MONOSPACE,
        width=dp.Width.FULL
    )
)

