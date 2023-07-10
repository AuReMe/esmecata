#!/usr/bin/env python
# -*- coding: utf-8 -*-

from esmecata.figures.stats_workflow_figures import *

__author__ = "Victor Mataigne"
__email__ = "victor.mataigne@irisa.fr"

INPUT_TABLE = "figures_input/input_table.tsv"
RESULTS_PATH = "figures_input"
OUTPATH = "figures_input"
RANK = "Phylum"

def main():
    input_data, df_stats, n_in, n_out, proteome_tax_id, association_taxon_tax_id = post_analysis_config(INPUT_TABLE, RESULTS_PATH)
    distributions_by_ranks(df_stats, n_in, n_out, OUTPATH, RANK)
    n_prot_ec_go_correlations(df_stats, n_in, n_out, OUTPATH)
    taxo_ranks_contribution(proteome_tax_id, n_in, n_out, OUTPATH)
    compare_ranks_in_out(input_data, proteome_tax_id, association_taxon_tax_id, OUTPATH)

if __name__ == "__main__":
    main()
