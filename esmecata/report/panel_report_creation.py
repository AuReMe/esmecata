#!/usr/bin/env python

# Copyright (C) 2023-2025 Pauline Hamon-Giraud and Alice Mataigne - Inria, Univ Rennes, CNRS, IRISA Dyliss
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

from bokeh.resources import INLINE, CDN
import panel as pn


def create_panel_report(CONFIG, DATA, DATA2, DATA3, fig1, fig4, fig5, fig5b, fig6, fig6b,
                        fig7, fig7b, fig8, fig8b, fig9, fig10, fig10a, fig11,
                        fig12, fig12_details, metadata, output_file, reduce=None):
    """ Create HTML report using panel from plotly figures.
    This function was created as a temporary fix after datapane was no longer maintained.
    Now it is replaced by the function using arakawa report. 

    Args:
        CONFIG (dict): metadata dictionary for plotly.
        DATA (dict): gathers all useful data for the analysis of an esmecata run (inputs, summary tables, taxonomy ...)
        DATA2 (dict): dictionary storing dataframes corresponding to EC number with summary numbers
        DATA3 (dict): dictionary storing dataframes corresponding to GO terms with summary numbers
        fig1 (plotly figure): histogram showing the distribution of proteomes according to taxonomic rank
        fig4 (plotly figure): heatmap showing how many input taxonomic ranks were assigned to the same rank by esmecata and how many were assigned to higher ranks
        fig5 (plotly figure): scatterplot showing the frequency of each EC in taxa
        fig5b (plotly figure): scatterplot showing the frequency of each GO Term in taxa
        fig6 (plotly figure): scatterplot showing the fraction of all EC of the dataset in each taxa
        fig6b (plotly figure): scatterplot showing the fraction of all GO Term of the dataset in each taxa
        fig7 (plotly figure): histogram showing how many EC numbers belongs to each taxa of this rank
        fig7b (plotly figure): histogram showing how many GO Terms belongs to each taxa of this rank
        fig8 (plotly figure): histogram of the fraction of all ECs of the dataset in each taxa
        fig8b (plotly figure): histogram of the fraction of all GO Terms of the dataset in each taxa
        fig9 (plotly figure): sunburst representing the taxonomic variety of taxa in the dataset.
        fig10 (plotly figure): sankey graph showing compression and taxonomic ranks rearrangement between input and output taxa of an Esmecata run
        fig10a (plotly figure): barplot showing the number proteomes according to the taxonomic rank
        fig11 (plotly figure): sunburst summarizing EC numbers categories, counts and proportions for the whole dataset
        fig12 (plotly figure): scatterplot showing the representativeness ratio and the number of associated protein clusters according to the taxonomic rank
        fig12_details (plotly figure): scatterplot showing details of the representativeness ratio and the number of associated protein clusters according to each cluster
        metadata (dict): all metadata of the workflow (made by merging esmecata metadata json from proteomes, clustering and annotation)
        output_file (str): path to HTML output file
    """
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
        pn.Column('## Distribution of proteomes according to taxonomic ranks',
            pn.pane.Plotly(fig10a.update_layout(modebar=CONFIG, autosize=True)),
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
                    pn.pane.HTML("<h2>Output stats</h2><p>Taxonomic ranks are NCBI ranks returned by ete4</p>", styles=styles),
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

