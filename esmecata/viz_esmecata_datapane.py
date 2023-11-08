#!/usr/bin/env python

from os import path
import json
import argparse
import datapane as dp
import stats_workflow_figures as swf
from esmecata2taxontology import esmecata2taxonomy
from analysis import create_dataset_annotation_file
from esmecata_compression import esmecata_compression

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

# ======
# CONFIG
# ======

DATA = swf.post_analysis_config(args.input, args.outdir)
_ = create_dataset_annotation_file(path.join(args.outdir, "2_annotation/annotation_reference"), path.join(args.outdir, "3_analysis/dataset_annotation_ec.tsv"), "EC")
_ = create_dataset_annotation_file(path.join(args.outdir, "2_annotation/annotation_reference"), path.join(args.outdir, "3_analysis/dataset_annotation_go.tsv"), "GO")
DATA2 = swf.create_annot_obs_df(path.join(args.outdir, "3_analysis/dataset_annotation_ec.tsv"), args.outdir, "EC numbers")
DATA3 = swf.create_annot_obs_df(path.join(args.outdir, "3_analysis/dataset_annotation_go.tsv"), args.outdir, "GO terms")

RANK = 'phylum'

CSS = {'main-title': {'textAlign': 'center', 
                      'color': '#7FDBFF'},
       'headers': {'textAlign': 'left', 
                   'color': '#7FDBFF'},
       'div-blocks': {'width': '49%', 
                      'display': 'inline-block'}}
CONFIG = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

# =======
# Figures
# =======

fig1 = swf.distributions_by_ranks(DATA["DF_STATS"], args.outdir, RANK)
fig2 = swf.n_prot_ec_go_correlations(DATA["DF_STATS"], args.outdir, RANK)
fig3 = swf.taxo_ranks_contribution(DATA["PROTEOME_TAX_ID"], args.outdir)
fig4 = swf.compare_ranks_in_out(DATA["PROTEOME_TAX_ID"], DATA["ASSOCIATION_PROTEOME_TAX_ID"], args.outdir)

fig5 = swf.annot_frequencies_in_obs(DATA2["df_annot_frequencies"], args.outdir, "EC numbers")
fig6 = swf.fraction_of_all_annot_in_obs(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "EC numbers")
fig7 = swf.annot_frequencies_in_obs_hist(DATA2["df_annot_frequencies"], args.outdir, "EC numbers")
fig8 = swf.fraction_of_all_annot_in_obs_hist(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "EC numbers")

fig5b = swf.annot_frequencies_in_obs(DATA3["df_annot_frequencies"], args.outdir, "GO terms")
fig6b = swf.fraction_of_all_annot_in_obs(DATA3["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "GO terms")
fig7b = swf.annot_frequencies_in_obs_hist(DATA3["df_annot_frequencies"], args.outdir, "GO terms")
fig8b = swf.fraction_of_all_annot_in_obs_hist(DATA3["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "GO terms")

fig9 = esmecata2taxonomy(args.outdir, path.join(args.outdir, '3_analysis/esmecata2taxonomy'))
fig10 = esmecata_compression(args.outdir, path.join(args.outdir, '3_analysis/esmecata_compression'), False)

fig11 = swf.ec_sunburst(DATA2["df_annot_frequencies"].index, args.outdir)
fig12 = swf.create_proteome_representativeness_lineplot_px(DATA["PROTEOME_TAX_ID"], path.join(args.outdir, '1_clustering/computed_threshold'), args.outdir)


if not DATA["DISCARDED"].empty:    
    panel_content = dp.DataTable(DATA["DISCARDED"],label="Data")
else:
    panel_content = dp.HTML("<p>None of input the taxonomic observations were discarded</p>")

metadata = swf.reproducibility_tokens(args.outdir)

# ======
# Report
# ======

report = dp.Blocks(

    # Page 1 : EsMeCaTa's inputs and outputs taxonomic ranks
    dp.Page(
        title="Inputs and outputs ranks comparison",
        blocks=[
            dp.HTML("<h2>Reduction of taxonomic diversity between EsMeCaTa's inputs and outputs</h2>"),            
            dp.Plot(fig9), #.update_traces(marker=dict(pattern=dict(shape=tmp_data['Shape'])))),
            dp.HTML("<h2>\"Compression\" of taxonomic diversity between inputs and outputs</h2>"),
            dp.Plot(fig10, label='Taxonomic compression')
        ],
    ),

    # Step 1 (Page 2) : EsMeCaTa proteomes
    dp.Page(
        title="Proteomes summary",
        blocks=[
            dp.Group(
                dp.HTML("<p>When proteomes of the lowest taxonomic rank of a taxonomic observation are not found, EsMeCaTa goes up in the taxonomic levels to found suited proteomes. In addition to figure 2, the heatmap belows details the input taxonomic rank of all taxonomic observations (y-axis), and the rank given by EsMeCaTa to the output (x-axis).</p>"),
                dp.HTML("<p>The barplot belows details how many proteomes were assigned to each taxonomic observation (colored by phylum). The number of proteomes found for each taxonomic observation can vary from only a few to dozens or hundreds, depending on the data available on UniProt.</p>"),
                dp.Plot(fig4.update_layout(modebar=CONFIG), label='Input and output ranks difference'),
                dp.Plot(fig1.update_layout(modebar=CONFIG), label='Distribution of the number of proteomes found'),                
                columns=2,
            )
        ],
    ),

    # Step 2 (Page 3): EsMeCaTa clustering
    dp.Page(
        title="Clustering summary",
        blocks=[
            dp.HTML("<p>In the clustering step, for each set of proteomes associated to an input, EsMeCaTa searchs for clusters of similar proteic sequences among these proteomes. This step is supervised by a threshold between 0 and 1, 1 meaning that a cluster contains sequences from each proteomes, 0.5 meaning that a cluster contains sequences from at least a half of the proteomes (etc). This threshold can then be an approximation of the pan-proteome (0) and the core-proteome (1). Clusters below the threshold are discarded. The figure below displays how many clusters were retained at each threshold value.</p>"),
            dp.Plot(fig12.update_layout(modebar=CONFIG), label='proteomes representativeness')
        ]
    ),

    # Step 3 (Page 4): EsMeCaTa annotation
    dp.Page(
        title="Annotation summary", blocks=[

            # EC numbers figures
            dp.HTML("<h2>EC numbers frequencies among taxa</h2>"),
            dp.HTML("<p>Numbers and figures below detail how EC numbers vary among taxa, i.e. if an EC number is found in many taxa (a generic EC number) or only in very few taxa (a specific EC number).</p>"),
            
            dp.Group(
                dp.BigNumber(heading="EC found in more than 90% of taxa", value=DATA2["n_9_an"]),
                dp.BigNumber(heading="EC found in less than 10% of taxa", value=DATA2["n_1_an"]),
                dp.BigNumber(heading="EC found between 10% and 90% of taxa", value=DATA2["n91_an"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig5.update_layout(modebar=CONFIG), label='EC numbers frequencies in observations'),
                dp.Plot(fig7.update_layout(modebar=CONFIG), label='EC numbers frequencies in observations (barplot version)'),
                columns=2,
            ),

            dp.HTML("<h2>Taxa content in EC numbers</h2>"),
            dp.HTML("<p>Numbers and figures below detail which taxa are annotated with most of the EC numbers of the whole dataset or only a few EC numbers. For instance, taxa with many EC numbers could have bigger genomes or be generalist species.</p>"),
            
            dp.Group(
                dp.BigNumber(heading="Taxa with more than 90% of all the dataset's EC numbers", value=DATA2["n_9_ob"]),
                dp.BigNumber(heading="Taxa with less than 10% of all the dataset's EC numbers", value=DATA2["n_1_ob"]),
                dp.BigNumber(heading="Taxa with between 10% and 90% of all the dataset's EC numbers", value=DATA2["n91_ob"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig6.update_layout(modebar=CONFIG), label="Taxa content in EC numbers"),
                dp.Plot(fig8.update_layout(modebar=CONFIG), label="Taxa content in EC numbers (barplot version)"),
                columns=2,
            ),

            dp.HTML("<h2>EC numbers classes, counts and proportions</h2>"),
            dp.Plot(fig11.update_layout(modebar=CONFIG), label="EC numbers classes, counts and proportions"),

            # GO terms figures
            dp.HTML("<h2>GO terms frequencies among taxa</h2>"),
            dp.HTML("<p>Numbers and figures below detail how GO terms vary among taxa, i.e. if a GO term is found in many taxa (a generic GO term) or only in very few taxa (a specific GO term).</p>"),
            
            dp.Group(
                dp.BigNumber(heading="GO terms found in more than 90% of taxa", value=DATA3["n_9_an"]),
                dp.BigNumber(heading="GO terms found in less than 10% of taxa", value=DATA3["n_1_an"]),
                dp.BigNumber(heading="GO terms found between 10% and 90% of taxa", value=DATA3["n91_an"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig5b.update_layout(modebar=CONFIG), label='Go terms frequencies in taxa'),
                dp.Plot(fig7b.update_layout(modebar=CONFIG), label='Go terms frequencies in taxa (barplot version)'),
                columns=2,
            ),

            dp.HTML("<h2>Taxa content in GO terms</h2>"),
            dp.HTML("<p>Numbers and figures below detail which taxa are annotated with most of the GO terms of the whole dataset or only a few Go terms. For instance, taxa with many Go terms could have bigger genomes or be generalist species.</p>"),
            
            dp.Group(
                dp.BigNumber(heading="Taxa with more than 90% of all the dataset's GO terms", value=DATA3["n_9_ob"]),
                dp.BigNumber(heading="Taxa with less than 10% of all the dataset's GO terms", value=DATA3["n_1_ob"]),
                dp.BigNumber(heading="Taxa with between 10% and 90% of all the dataset's GO terms", value=DATA3["n91_ob"]),
                columns=3,
            ),
            dp.Group(
                dp.Plot(fig6b.update_layout(modebar=CONFIG), label="Observations' content in GO terms"),
                dp.Plot(fig8b.update_layout(modebar=CONFIG), label="Observations' content in GO terms (barplot version)"),
                columns=2,
            ),
        ]
    ),

    # Page 4 : dataframes of default summary statistics
    dp.Page(
        title="Data summary", blocks=[
            dp.Group(
                dp.BigNumber(heading="Number of inputs", value=DATA["N_IN"]),
                dp.BigNumber(heading="Kept by EsMeCaTa", value=DATA["N_OUT"]),
                dp.BigNumber(heading="Discarded", value=DATA["N_DISCARDED"]),
                columns=3,
            ),

            dp.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete3</p>"),
            dp.DataTable(DATA["DF_STATS"], label="Data"),
            dp.HTML("<h2>Discarded</h2><p>Taxonomic ranks were not inferred; only names are displayed</p>"),
            panel_content
        ]
    ),

    # Page 5 : metadata with tools versions, files paths (...) for reproducibility
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
