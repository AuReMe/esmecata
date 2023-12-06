#!/usr/bin/env python

from os import path, mkdir, makedirs
import json
import argparse
import datapane as dp
import stats_workflow_figures as swf
from esmecata2taxontology import esmecata2taxonomy
from analysis import create_dataset_annotation_file
from esmecata_compression import esmecata_compression
from plotly.io import write_json

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

# if not path.exists(path.join(args.outdir, "3_analysis")):
#     mkdir(path.join(args.outdir, "3_analysis"))

dirs = ["inputs_outputs_figures", "proteomes_figures", "clustering_figures", "annotation_figures"]
for dir in dirs:
    try: makedirs(path.join(args.outdir, "3_analysis", dir))
    except OSError: pass

print("Getting data")
DATA = swf.post_analysis_config(args.input, args.outdir)
_ = create_dataset_annotation_file(path.join(args.outdir, "2_annotation/annotation_reference"), path.join(args.outdir, "3_analysis/dataset_annotation_ec.tsv"), "EC")
_ = create_dataset_annotation_file(path.join(args.outdir, "2_annotation/annotation_reference"), path.join(args.outdir, "3_analysis/dataset_annotation_go.tsv"), "GO")
DATA2 = swf.create_annot_obs_df(path.join(args.outdir, "3_analysis/dataset_annotation_ec.tsv"), args.outdir, "EC numbers")
DATA3 = swf.create_annot_obs_df(path.join(args.outdir, "3_analysis/dataset_annotation_go.tsv"), args.outdir, "GO terms")
DF_CLUSTERING = swf.data_proteome_representativeness(DATA["PROTEOME_TAX_ID"], path.join(args.outdir, '1_clustering/computed_threshold'))

RANK = 'phylum'

CSS = {'main-title': {'textAlign': 'center', 
                      'color': '#7FDBFF'},
       'headers': {'textAlign': 'left', 
                   'color': '#7FDBFF'},
       'div-blocks': {'width': '49%', 
                      'display': 'inline-block'}}
CONFIG = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

metadata = swf.reproducibility_tokens(args.outdir)

# =======
# Figures
# =======

print("Building Inputs summary figures")
fig9 = esmecata2taxonomy(args.outdir, path.join(args.outdir, '3_analysis/inputs_outputs_figures/esmecata2taxonomy'))
fig10 = esmecata_compression(args.outdir, path.join(args.outdir, '3_analysis/inputs_outputs_figures/esmecata_compression'), False)

write_json(fig9, 
    path.join(args.outdir, "3_analysis/inputs_outputs_figures/esmecata2taxonomy.json"), 
    pretty=True)
write_json(fig10, 
    path.join(args.outdir, "3_analysis/inputs_outputs_figures/esmecata_compression.json"), 
    pretty=True)

print("Building proteomes summary figures")
fig1 = swf.distributions_by_ranks(DATA["DF_STATS"], args.outdir, 15, RANK)
# fig2 = swf.n_prot_ec_go_correlations(DATA["DF_STATS"], args.outdir, RANK)
# fig3 = swf.taxo_ranks_contribution(DATA["PROTEOME_TAX_ID"], args.outdir)
fig4 = swf.compare_ranks_in_out(DATA["PROTEOME_TAX_ID"], DATA["ASSOCIATION_PROTEOME_TAX_ID"], args.outdir)

print("Building clustering summary figures")
fig12 = swf.create_proteome_representativeness_lineplot_px(DF_CLUSTERING,
    metadata["clustering"]["tool_options"]["clust_threshold"],
    args.outdir)

fig12_details = swf.proteomes_representativeness_details(DF_CLUSTERING,
    metadata["clustering"]["tool_options"]["clust_threshold"],
    args.outdir)

print("Building annotation summary figures")
fig5 = swf.annot_frequencies_in_obs(DATA2["df_annot_frequencies"], args.outdir, "EC numbers")
fig6 = swf.fraction_of_all_annot_in_obs(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "EC numbers")
fig7 = swf.annot_frequencies_in_obs_hist(DATA2["df_annot_frequencies"], args.outdir, "EC numbers")
fig8 = swf.fraction_of_all_annot_in_obs_hist(DATA2["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "EC numbers")

fig5b = swf.annot_frequencies_in_obs(DATA3["df_annot_frequencies"], args.outdir, "GO terms")
fig6b = swf.fraction_of_all_annot_in_obs(DATA3["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "GO terms")
fig7b = swf.annot_frequencies_in_obs_hist(DATA3["df_annot_frequencies"], args.outdir, "GO terms")
fig8b = swf.fraction_of_all_annot_in_obs_hist(DATA3["df_fractionin_obs"], DATA["DF_STATS"], args.outdir, "GO terms")

fig11 = swf.ec_sunburst(DATA2["df_annot_frequencies"].index, args.outdir)

print("Formatting summary dataframes")
if not DATA["DISCARDED"].empty:    
    df_discarded_panel_content = dp.DataTable(DATA["DISCARDED"],label="Data")
else:
    df_discarded_panel_content = dp.HTML("<p>None of input the taxonomic observations were discarded</p>")

print("Formatting metadata")
metadata = json.dumps(metadata, indent=4)

# ======
# Report
# ======

report = dp.Blocks(

    # Page 1 : EsMeCaTa's inputs and outputs taxonomic ranks
    dp.Page(
        title="Inputs and outputs ranks comparison",
        blocks=[
            dp.HTML("<h2>Reduction of taxonomic diversity between EsMeCaTa's inputs and outputs</h2>"),
            dp.HTML("<p>The sunburst below displays the taxonomic diversity of inputs, i.e. all the taxonomic affiliations listed in the tsv input file. White cells indicates inputs whose ranks were not retained by EsMeCaTa, which went up to the superior taxonomic rank to find proteomes.</p>"),           
            dp.Plot(fig9), #.update_traces(marker=dict(pattern=dict(shape=tmp_data['Shape'])))),
            dp.HTML("<h2>\"Compression\" of taxonomic diversity between inputs and outputs</h2>"),
            dp.HTML("<p>The sankey diagram below also displays the initial taxonomic diversity of taxonomic affiliations given to EsMeCaTa as inputs. The first column details the content of the input file : one block is equivalent to one taxonomic affiliation, and the block height indicates the number of times this affiliation appears in the input file. The second column details which of these taxonomic affiliations were merged together (i.e. 'compressed') because their taxonomyies were identical in the input file. Finally, The first column displays the ranks EsMeCaTa attributed to each of them : at this step, some taxonomic affiliations can be merged because EsMeCaTa went up the taxonomic ranks to find proteomes.</p>"),
            dp.Plot(fig10, label='Taxonomic compression')
        ],
    ),

    # Step 1 (Page 2) : EsMeCaTa proteomes
    dp.Page(
        title="Proteomes summary",
        blocks=[
            dp.Group(
                dp.HTML("<p>When proteomes of the lowest taxonomic rank of an input taxa are not found, EsMeCaTa goes up in the taxonomic levels to found suited proteomes. In complement to the figures in the first panel, the heatmap below details the difference between the lowest taxonomic rank of all input taxonomic observations (y-axis) and the lowest rank atributed by EsMeCaTa (x-axis).</p>"),
                dp.HTML("<p>The barplot belows details how many proteomes were assigned to inputs taxonomic affiliations (grouped and colored by phylum). The number of proteomes found can vary from only a few to dozens or hundreds, depending on the data available on UniProt for each taxonomic group.</p>"),
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
            dp.HTML("<p>In the clustering step, for the set of proteomes associated to each input taxa, EsMeCaTa searches for clusters of similar proteic sequences among these proteomes. This step is supervised by a threshold between 0 and 1, 1 meaning that a cluster contains sequences from each proteomes, 0.5 meaning that a cluster contains sequences from at least a half of the proteomes (etc). This threshold can then be an approximation of the pan-proteome (0) and the core-proteome (1). Clusters below the threshold are discarded. The figure below displays how many clusters were retained at each threshold value. The vertical red line displays the threshold set dor this run of EsMeCata.</p>"),
            dp.Plot(fig12.update_layout(modebar=CONFIG), label='proteomes representativeness'),
            dp.Plot(fig12_details.update_layout(modebar=CONFIG), label="proteomes representativeness details")
        ]
    ),

    # Step 3 (Page 4): EsMeCaTa annotation
    dp.Page(
        title="Annotation summary", 
        blocks=[

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
        title="Data summary", 
        blocks=[
            dp.Group(
                dp.BigNumber(heading="Number of inputs", value=DATA["N_IN"]),
                dp.BigNumber(heading="Kept by EsMeCaTa", value=DATA["N_OUT"]),
                dp.BigNumber(heading="Discarded", value=DATA["N_DISCARDED"]),
                columns=3,
            ),

            dp.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete3</p>"),
            dp.DataTable(DATA["DF_STATS"], label="Data"),
            dp.HTML("<h2>Discarded</h2><p>Taxonomic ranks were not inferred; only names are displayed</p>"),
            df_discarded_panel_content
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
