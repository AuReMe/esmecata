#!/usr/bin/env python

from os import path
import json
import argparse
import datapane as dp
import stats_workflow_figures as swf
from esmecata2taxontology import esmecata2taxonomy
from analysis import create_dataset_annotation_file

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

# INPUT_DATA, DISCARDED, N_DISCARDED, DF_STATS, N_IN, N_OUT, PROTEOME_TAX_ID, ASSOCIATION_TAXON_TAX_ID, ECS_MATRIX = swf.post_analysis_config(args.input, args.outdir)

DATA = swf.post_analysis_config(args.input, args.outdir)
_ = create_dataset_annotation_file(path.join(args.outdir, "2_annotation/annotation_reference"), path.join(args.outdir, "3_analysis/dataset_annotation.tsv"))
DATA2 = swf.create_ec_obs_df(path.join(args.outdir, "3_analysis/dataset_annotation.tsv"), args.outdir)
RANK = 'phylum'

CSS = {'main-title': {'textAlign': 'center', 'color': '#7FDBFF'},
       'headers': {'textAlign': 'left', 'color': '#7FDBFF'},
       'div-blocks': {'width': '49%', 'display': 'inline-block'}}

CONFIG = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

fig1 = swf.distributions_by_ranks(DATA["DF_STATS"], args.outdir, RANK)
fig2 = swf.n_prot_ec_go_correlations(DATA["DF_STATS"], args.outdir, RANK)
fig3 = swf.taxo_ranks_contribution(DATA["PROTEOME_TAX_ID"], args.outdir)
fig4 = swf.compare_ranks_in_out(DATA["PROTEOME_TAX_ID"], DATA["ASSOCIATION_PROTEOME_TAX_ID"], args.outdir)
fig5 = swf.ecs_frequencies_in_obs(DATA2["df_ec_frequencies"], args.outdir)
fig6 = swf.fraction_of_all_ec_in_obs(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir)
fig7 = swf.ecs_frequencies_in_obs_hist(DATA2["df_ec_frequencies"], args.outdir)
fig8 = swf.fraction_of_all_ec_in_obs_hist(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir)
fig9 = esmecata2taxonomy(args.outdir, path.join(args.outdir, '3_analysis/esmecata2taxonomy'))

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
            dp.Group(
                dp.BigNumber(heading="EC with a frequency > 0.9", value=DATA2["n_9_ec"]),
                dp.BigNumber(heading="EC with a frequency < 0.1", value=DATA2["n_1_ec"]),
                dp.BigNumber(heading="EC with a frequency between 0.1 and 0.9", value=DATA2["n91_ec"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig5.update_layout(modebar=CONFIG), label='EC frequencies in observations'),
                dp.Plot(fig7.update_layout(modebar=CONFIG), label='EC frequencies in observations (barplot)'),
                columns=2,
                #widths=[0.75, 0.25],
            ),
            dp.Group(
                dp.BigNumber(heading="Observations with > 0.9 of ECs", value=DATA2["n_9_ob"]),
                dp.BigNumber(heading="Observations with < 0.1 of ECs", value=DATA2["n_1_ob"]),
                dp.BigNumber(heading="Observations having between 0.1 and O.9 of ECs", value=DATA2["n91_ob"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig6.update_layout(modebar=CONFIG), label="Observations' content in EC numbers"),
                dp.Plot(fig8.update_layout(modebar=CONFIG), label="Observations' content in EC numbers (barplot)"),
                columns=2,
                #widths=[0.75, 0.25],                
            ),
            dp.Plot(fig9, label='Taxonomic diversity')
        ],
    ),
    dp.Page(
        title="Data", blocks=[
            dp.Group(
                dp.BigNumber(heading="Number of inputs", value=DATA["N_IN"]),
                dp.BigNumber(heading="Kept by EsMeCaTa", value=DATA["N_OUT"]),
                dp.BigNumber(heading="Discarded", value=DATA["N_DISCARDED"]),
                columns=3,
            ),
            dp.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete3</p>"),
            dp.DataTable(DATA["DF_STATS"], 
                        label="Data"),
            dp.HTML("<h2>Discarded</h2><p>Taxonomic ranks were not inferred; only names are displayed</p>"),
            dp.DataTable(DATA["DISCARDED"], 
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

