#!/usr/bin/env python

# Copyright (C) 2023-2024 Pauline Hamon-Giraud and Victor Mataigne - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import datapane as dp

from plotly.io import write_json

from esmecata.core.annotation import create_dataset_annotation_file
from esmecata.report import stats_workflow_figures as swf
from esmecata.report.esmecata2taxontology import esmecata2taxonomy
from esmecata.report.esmecata_compression import esmecata_compression


def create_panel_report(CONFIG, DATA, DATA2, DATA3, fig1, fig4, fig5, fig5b, fig6, fig6b,
                        fig7, fig7b, fig8, fig8b, fig9, fig10, fig11,
                        fig12, fig12_details, metadata, output_file, reduce=None):
    from bokeh.resources import INLINE, CDN
    import panel as pn

    styles = {
        'background-color': '#F6F6F6', 'border': '2px solid black',
        'border-radius': '5px', 'padding': '10px'
    }

    page_one = pn.Column(
        pn.pane.Markdown("""# Reduction of taxonomic diversity between EsMeCaTa's inputs and outputs
                    The sunburst below displays the taxonomic diversity of inputs, i.e. all the taxonomic affiliations listed in the tsv input file. White cells indicates inputs whose ranks were not retained by EsMeCaTa, which went up to the superior taxonomic rank to find proteomes.""",
                    styles=styles, renderer='markdown'),
        pn.Column('## Sunburst of taxon selected by EsMeCaTa',
            pn.pane.Plotly(fig9.update_layout(autosize=True)),
            sizing_mode='scale_both'),

        pn.pane.Markdown("""# \"Compression\" of taxonomic diversity between inputs and outputs
                    The sankey diagram below also displays the initial taxonomic diversity of taxonomic affiliations given to EsMeCaTa as inputs. The first column details the content of the input file : one block is equivalent to one taxonomic affiliation, and the block height indicates the number of times this affiliation appears in the input file. The second column details which of these taxonomic affiliations were merged together (i.e. 'compressed') because their taxonomyies were identical in the input file. Finally, The first column displays the ranks EsMeCaTa attributed to each of them : at this step, some taxonomic affiliations can be merged because EsMeCaTa went up the taxonomic ranks to find proteomes.""", styles=styles),
        pn.Column('## Sankey diagram of taxonomic compression',
            pn.pane.Plotly(fig10.update_layout(autosize=True)),
            sizing_mode='scale_both'),
    )
    page_two = pn.Column(
        pn.pane.Markdown("""# Proteomes summary
                     When proteomes of the lowest taxonomic rank of an input taxa are not found, EsMeCaTa goes up in the taxonomic levels to found suited proteomes. In complement to the figures in the first panel, the heatmap below details the difference between the lowest taxonomic rank of all input taxonomic observations (y-axis) and the lowest rank atributed by EsMeCaTa (x-axis).""", styles=styles),
        pn.Column('## Input and output ranks difference',
            pn.pane.Plotly(fig4.update_layout(modebar=CONFIG, autosize=True)),
            sizing_mode='scale_both'),
        pn.pane.HTML("""<p>The barplot belows details how many proteomes were assigned to inputs taxonomic affiliations (grouped and colored by phylum). The number of proteomes found can vary from only a few to dozens or hundreds, depending on the data available on UniProt for each taxonomic group.</p>""", styles=styles),
        pn.Column('## Distribution of the number of proteomes found',
            pn.pane.Plotly(fig1.update_layout(modebar=CONFIG, autosize=True)),
            sizing_mode='scale_both'),
    )
    page_three = pn.Column(
            pn.pane.HTML("""<h2>Clustering summary</h2>
                         <p>In the clustering step, for the set of proteomes associated to each input taxa, EsMeCaTa searches for clusters of similar proteic sequences among these proteomes. This step is supervised by a threshold between 0 and 1, 1 meaning that a cluster contains sequences from each proteomes, 0.5 meaning that a cluster contains sequences from at least a half of the proteomes (etc). This threshold can then be an approximation of the pan-proteome (0) and the core-proteome (1). Clusters below the threshold are discarded. The figure below displays how many clusters were retained at each threshold value. The vertical red line displays the threshold set dor this run of EsMeCata.</p>""", styles=styles),
        pn.Column('## Proteomes representativeness',
            pn.pane.Plotly(fig12.update_layout(modebar=CONFIG, autosize=True)),
            sizing_mode='scale_both'),
        pn.Column('## Proteomes representativeness details',
            pn.pane.Plotly(fig12_details.update_layout(modebar=CONFIG, autosize=True)),
            sizing_mode='scale_both'),
    )
    page_four = pn.Column(
        # EC numbers figures
        pn.FlexBox(
            pn.pane.HTML("""<h2>EC numbers frequencies among taxa</h2>
                        <p>Numbers and figures below detail how EC numbers vary among taxa, i.e. if an EC number is found in many taxa (a generic EC number) or only in very few taxa (a specific EC number).</p>""", styles=styles),
            pn.indicators.Number(name='EC found in more than 90% of taxa', value=DATA2["n_9_an"], format='{value}'),
            pn.indicators.Number(name='EC found in less than 10% of taxa', value=DATA2["n_1_an"], format='{value}'),
            pn.indicators.Number(name='EC found between 10% and 90% of taxa', value=DATA2["n91_an"], format='{value}'),
        ),
        pn.Column(
            pn.Column('## EC numbers frequencies in observations',
                pn.pane.Plotly(fig5.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
            pn.Column('## EC numbers frequencies in observations (barplot version)',
                pn.pane.Plotly(fig7.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
        ),

        pn.FlexBox(
            pn.pane.HTML("""<h2>Taxa content in EC numbers</h2>
                        <p>Numbers and figures below detail which taxa are annotated with most of the EC numbers of the whole dataset or only a few EC numbers. For instance, taxa with many EC numbers could have bigger genomes or be generalist species.</p>""", styles=styles),
            pn.indicators.Number(name="Taxa with more than 90% of all the dataset's EC numbers", value=DATA2["n_9_ob"], format='{value}'),
            pn.indicators.Number(name="Taxa with less than 10% of all the dataset's EC numbers", value=DATA2["n_1_ob"], format='{value}'),
            pn.indicators.Number(name="Taxa with between 10% and 90% of all the dataset's EC numbers", value=DATA2["n91_ob"], format='{value}'),
        ),
        pn.Column(
            pn.Column('## Taxa content in EC numbers',
                pn.pane.Plotly(fig6.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
            pn.Column('## Taxa content in EC numbers (barplot version)',
                pn.pane.Plotly(fig8.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
        ),

        pn.pane.HTML("<h2>EC numbers classes, counts and proportions</h2>", styles=styles),
        pn.Column('## EC numbers classes, counts and proportions',
            pn.pane.Plotly(fig11.update_layout(modebar=CONFIG, autosize=True)),
            sizing_mode='scale_both'),

        # GO terms figures
        pn.FlexBox(
            pn.pane.HTML("""<h2>GO terms frequencies among taxa</h2>
                        <p>Numbers and figures below detail how GO terms vary among taxa, i.e. if a GO term is found in many taxa (a generic GO term) or only in very few taxa (a specific GO term).</p>""", styles=styles),
            pn.indicators.Number(name="GO terms found in more than 90% of taxa", value=DATA3["n_9_an"], format='{value}'),
            pn.indicators.Number(name="GO terms found in less than 10% of taxa", value=DATA3["n_1_an"], format='{value}'),
            pn.indicators.Number(name="GO terms found between 10% and 90% of taxa", value=DATA3["n91_an"], format='{value}'),
        ),
        pn.Column(
            pn.Column('## GO terms frequencies in taxa',
                pn.pane.Plotly(fig5b.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
            pn.Column('## GO terms frequencies in taxa (barplot version)',
                pn.pane.Plotly(fig7b.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
        ),

        pn.FlexBox(
        pn.pane.HTML("""<h2>Taxa content in GO terms</h2>
                     <p>Numbers and figures below detail which taxa are annotated with most of the GO terms of the whole dataset or only a few Go terms. For instance, taxa with many Go terms could have bigger genomes or be generalist species.</p>""", styles=styles),
            pn.indicators.Number(name="Taxa with more than 90% of all the dataset's GO terms", value=DATA3["n_9_ob"], format='{value}'),
            pn.indicators.Number(name="Taxa with less than 10% of all the dataset's GO terms", value=DATA3["n_1_ob"], format='{value}'),
            pn.indicators.Number(name="Taxa with between 10% and 90% of all the dataset's GO terms", value=DATA3["n91_ob"], format='{value}'),
        ),
        pn.Column(
            pn.Column("## Observations' content in GO terms",
                pn.pane.Plotly(fig6b.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
            pn.Column("## Observations' content in GO terms (barplot version)",
                pn.pane.Plotly(fig8b.update_layout(modebar=CONFIG, autosize=True)),
                sizing_mode='scale_both'),
        ),
        height=5000
    )

    if not DATA["DISCARDED"].empty:    
        discarded_panel_content_panel = pn.pane.DataFrame(DATA["DISCARDED"])
    else:
        discarded_panel_content_panel =  pn.pane.HTML("<p>None of input the taxonomic observations were discarded</p>", styles=styles)

    page_five = pn.Column(
                pn.FlexBox(
                    pn.indicators.Number(name="Number of inputs", value=DATA["N_IN"], format='{value}'),
                    pn.indicators.Number(name="Kept by EsMeCaTa", value=DATA["N_OUT"], format='{value}'),
                    pn.indicators.Number(name="Discarded", value=DATA["N_DISCARDED"], format='{value}'),
                ),
                pn.FlexBox(
                    pn.pane.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete3</p>", styles=styles),
                    pn.pane.DataFrame(DATA["DF_STATS"]),
                ),
                pn.FlexBox(
                    pn.pane.HTML("<h2>Discarded</h2><p>Taxonomic ranks were not inferred; only names are displayed</p>", styles=styles),
                    discarded_panel_content_panel
                ),
    )

    page_six = pn.Column(
        pn.pane.JSON(metadata, name="Metadata", depth=-1),
    )

    tabs = pn.Tabs(("Inputs and outputs ranks comparison", page_one),
                   ("Proteomes summary", page_two),
                   ("Clustering summary", page_three),
                   ("Annotation summary", page_four),
                   ("Data summary", page_five),
                   ("Metadata", page_six),
                   height_policy="fit", sizing_mode="fixed")

    panel = pn.template.VanillaTemplate(
        title="EsMeCaTa report",
        main=[tabs],
        header_color='black',
        header_background='white'
    )
 
    if reduce is True:
        panel.save(output_file, resources=CDN)
    else:
        panel.save(output_file, resources=INLINE)


def create_datapane_report(CONFIG, DATA, DATA2, DATA3, fig1, fig4, fig5, fig5b, fig6, fig6b,
                        fig7, fig7b, fig8, fig8b, fig9, fig10, fig11,
                        fig12, fig12_details, metadata, output_file):
    print("Formatting summary dataframes")
    if not DATA["DISCARDED"].empty:    
        df_discarded_panel_content = dp.DataTable(DATA["DISCARDED"],label="Data")
    else:
        df_discarded_panel_content = dp.HTML("<p>None of input the taxonomic observations were discarded</p>")

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
        path=output_file,
        formatting=dp.Formatting(
            accent_color="rgb(204, 255, 204)",
            font=dp.FontChoice.MONOSPACE,
            width=dp.Width.FULL
        )
    )

def create_esmecata_report(esmecata_input_file, esmecata_core_output_folder, output_folder, create_svg=False):

    # ======
    # CONFIG
    # ======
    output_dirs = ["inputs_outputs_figures", "proteomes_figures", "clustering_figures", "annotation_figures"]
    for output_dir in output_dirs:
        try: os.makedirs(os.path.join(output_folder, output_dir))
        except OSError: pass

    print("Getting data")
    computed_threshold_folder = os.path.join(esmecata_core_output_folder, '1_clustering', 'computed_threshold')
    annotation_reference_folder = os.path.join(esmecata_core_output_folder, '2_annotation', 'annotation_reference')
    dataset_annotation_ec_file = os.path.join(output_folder, 'dataset_annotation_ec.tsv')
    dataset_annotation_go_file = os.path.join(output_folder, 'dataset_annotation_go.tsv')

    DATA = swf.post_analysis_config(esmecata_input_file, esmecata_core_output_folder)

    data_summary_file_path = os.path.join(output_folder, 'data_summary.tsv')
    DATA["DF_STATS"].to_csv(data_summary_file_path, sep='\t')

    _ = create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_ec_file, 'EC')

    _ = create_dataset_annotation_file(annotation_reference_folder, dataset_annotation_go_file, 'GO')

    print("Getting stats from workflow")
    DATA2 = swf.create_annot_obs_df(dataset_annotation_ec_file, output_folder, 'EC numbers')

    DATA3 = swf.create_annot_obs_df(dataset_annotation_go_file, output_folder, "GO terms")

    print("Getting data for proteome representativeness")
    DF_CLUSTERING = swf.data_proteome_representativeness(DATA['PROTEOME_TAX_ID'], computed_threshold_folder)

    RANK = 'phylum'

    CSS = {'main-title': {'textAlign': 'center', 
                        'color': '#7FDBFF'},
        'headers': {'textAlign': 'left', 
                    'color': '#7FDBFF'},
        'div-blocks': {'width': '49%', 
                        'display': 'inline-block'}}
    CONFIG = {'remove': ['select', 'zoomIn', 'zoomOut', 'autoScale', 'lasso2d']}

    metadata = swf.reproducibility_tokens(esmecata_core_output_folder)

    # =======
    # Figures
    # =======

    print("Building Inputs summary figures")
    esmecata2taxonomy_output_folder = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata2taxonomy')
    fig9 = esmecata2taxonomy(esmecata_core_output_folder, esmecata2taxonomy_output_folder)
    esmecata_compression_output_folder = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata_compression')
    fig10 = esmecata_compression(esmecata_core_output_folder, esmecata_compression_output_folder, False)

    esmecata2taxonomy_json_file_path = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata2taxonomy.json')
    write_json(fig9, esmecata2taxonomy_json_file_path, pretty=True)
    if create_svg is True:
        esmecata2taxonomy_svg_file_path = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata2taxonomy.svg')
        fig9.write_image(esmecata2taxonomy_svg_file_path)

    esmecata_compression_json_file_path = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata_compression.json')
    write_json(fig10, esmecata_compression_json_file_path, pretty=True)
    if create_svg is True:
        esmecata_compression_svg_file_path = os.path.join(output_folder, 'inputs_outputs_figures', 'esmecata_compression.svg')
        fig10.write_image(esmecata_compression_svg_file_path)

    print("Building proteomes summary figures")
    output_proteomes_figure_folder = os.path.join(output_folder, 'proteomes_figures')
    fig1 = swf.distributions_by_ranks(DATA["DF_STATS"], output_proteomes_figure_folder, 15, RANK)
    if create_svg is True:
        proteomes_distribution_by_ranks_svg_file = os.path.join(output_proteomes_figure_folder, 'proteomes_distribution_by_ranks.svg')
        fig1.write_image(proteomes_distribution_by_ranks_svg_file)

    # fig2 = swf.n_prot_ec_go_correlations(DATA["DF_STATS"], args.outdir, RANK)
    # fig3 = swf.taxo_ranks_contribution(DATA["PROTEOME_TAX_ID"], args.outdir)
    fig4 = swf.compare_ranks_in_out(DATA["PROTEOME_TAX_ID"], DATA["ASSOCIATION_PROTEOME_TAX_ID"], output_proteomes_figure_folder)
    if create_svg is True:
        input_and_output_ranks_svg_file = os.path.join(output_proteomes_figure_folder, 'input_and_output_ranks.svg')
        fig4.write_image(input_and_output_ranks_svg_file)

    print("Building clustering summary figures")
    output_clustering_figures_folder = os.path.join(output_folder, 'clustering_figures')
    fig12 = swf.create_proteome_representativeness_lineplot_px(DF_CLUSTERING,
        metadata["clustering"]["tool_options"]["clust_threshold"],
        output_clustering_figures_folder)
    if create_svg is True:
        proteome_representativeness_svg_file = os.path.join(output_clustering_figures_folder, 'proteome_representativeness.svg')
        fig12.write_image(proteome_representativeness_svg_file)

    fig12_details = swf.proteomes_representativeness_details(DF_CLUSTERING,
        metadata["clustering"]["tool_options"]["clust_threshold"],
        output_clustering_figures_folder)
    if create_svg is True:
        proteome_representativeness_details_svg_file = os.path.join(output_clustering_figures_folder, 'proteome_representativeness_details.svg')
        fig12_details.write_image(proteome_representativeness_details_svg_file)

    print("Building annotation summary figures")
    output_annotation_figures_folder = os.path.join(output_folder, 'annotation_figures')
    fig5 = swf.annot_frequencies_in_obs(DATA2["df_annot_frequencies"], output_annotation_figures_folder, "EC numbers")
    fig6 = swf.fraction_of_all_annot_in_obs(DATA2["df_fractionin_obs"], DATA["DF_STATS"], output_annotation_figures_folder, "EC numbers")
    fig7 = swf.annot_frequencies_in_obs_hist(DATA2["df_annot_frequencies"], output_annotation_figures_folder, "EC numbers")
    fig8 = swf.fraction_of_all_annot_in_obs_hist(DATA2["df_fractionin_obs"], DATA["DF_STATS"], output_annotation_figures_folder, "EC numbers")
    if create_svg is True:
        EC_numbers_frequencies_in_taxa_svg_file = os.path.join(output_annotation_figures_folder, 'EC_numbers_frequencies_in_taxa.svg')
        fig5.write_image(EC_numbers_frequencies_in_taxa_svg_file)
        EC_numbers_fraction_per_taxa_svg_file = os.path.join(output_annotation_figures_folder, 'EC_numbers_fraction_per_taxa.svg')
        fig6.write_image(EC_numbers_fraction_per_taxa_svg_file)
        EC_numbers_frequencies_in_taxa_hist_svg_file = os.path.join(output_annotation_figures_folder, 'EC_numbers_frequencies_in_taxa_hist.svg')
        fig7.write_image(EC_numbers_frequencies_in_taxa_hist_svg_file)
        EC_numbers_fraction_per_taxa_hist_svg_file = os.path.join(output_annotation_figures_folder, 'EC_numbers_fraction_per_taxa_hist.svg')
        fig8.write_image(EC_numbers_fraction_per_taxa_hist_svg_file)

    fig5b = swf.annot_frequencies_in_obs(DATA3["df_annot_frequencies"], output_annotation_figures_folder, "GO terms")
    fig6b = swf.fraction_of_all_annot_in_obs(DATA3["df_fractionin_obs"], DATA["DF_STATS"], output_annotation_figures_folder, "GO terms")
    fig7b = swf.annot_frequencies_in_obs_hist(DATA3["df_annot_frequencies"], output_annotation_figures_folder, "GO terms")
    fig8b = swf.fraction_of_all_annot_in_obs_hist(DATA3["df_fractionin_obs"], DATA["DF_STATS"], output_annotation_figures_folder, "GO terms")

    fig11 = swf.ec_sunburst(DATA2["df_annot_frequencies"].index, output_annotation_figures_folder)

    if create_svg is True:
        GO_terms_frequencies_in_taxa_svg_file = os.path.join(output_annotation_figures_folder, 'GO_terms_frequencies_in_taxa.svg')
        fig5b.write_image(GO_terms_frequencies_in_taxa_svg_file)
        GO_terms_fraction_per_taxa_svg_file = os.path.join(output_annotation_figures_folder, 'GO_terms_fraction_per_taxa.svg')
        fig6b.write_image(GO_terms_fraction_per_taxa_svg_file)
        GO_terms_frequencies_in_taxa_hist_svg_file = os.path.join(output_annotation_figures_folder, 'GO_terms_frequencies_in_taxa_hist.svg')
        fig7b.write_image(GO_terms_frequencies_in_taxa_hist_svg_file)
        GO_terms_fraction_per_taxa_hist_svg_file = os.path.join(output_annotation_figures_folder, 'GO_terms_fraction_per_taxa_hist.svg')
        fig8b.write_image(GO_terms_fraction_per_taxa_hist_svg_file)
        ec_classes_sunburst_svg_file = os.path.join(output_annotation_figures_folder, 'ec_classes_sunburst.svg')
        fig11.write_image(ec_classes_sunburst_svg_file)

    print("Formatting metadata")
    metadata = json.dumps(metadata, indent=4)

    # ======
    # Report
    # ======
    """
    esmecata_summary_panel = os.path.join(output_folder, 'esmecata_summary_panel.html')
    create_panel_report(CONFIG, DATA, DATA2, DATA3, fig1, fig4, fig5, fig5b, fig6, fig6b,
                            fig7, fig7b, fig8, fig8b, fig9, fig10, fig11,
                            fig12, fig12_details, metadata, esmecata_summary_panel, True)
    """
    esmecata_summary = os.path.join(output_folder, 'esmecata_summary.html')
    create_datapane_report(CONFIG, DATA, DATA2, DATA3, fig1, fig4, fig5, fig5b, fig6, fig6b,
                            fig7, fig7b, fig8, fig8b, fig9, fig10, fig11,
                            fig12, fig12_details, metadata, esmecata_summary)
