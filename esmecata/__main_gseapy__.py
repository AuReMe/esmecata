# Copyright (C) 2021-2025 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
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

import argparse
import logging
import os
import sys
import time

from esmecata.gseapy.gseapy_orsum import taxon_rank_annotation_enrichment
from esmecata.utils import is_valid_dir
from esmecata import __version__ as VERSION

MESSAGE = '''
Create enrichment analysis files from esmecata results.
'''
REQUIRES = '''
Requires: bioservices, pronto, gseapy and orsum
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'esmecata_gseapy',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION + '\n')

    parent_parser_f = argparse.ArgumentParser(add_help=False)
    parent_parser_f.add_argument(
        '-f',
        '--folder',
        dest='input_file_folder',
        required=True,
        help='Annotation input folder or file. For folder, it corresponds to EsMeCaTa annotation output folder. For file, it ocrresponds to a function table containing annotation as column and organism as row.',
        metavar='INPUT_FOLDER')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_grouping = argparse.ArgumentParser(add_help=False)
    parent_parser_grouping.add_argument(
        '--grouping',
        dest='grouping',
        required=True,
        help='Grouping factor used for enrichment analysis (either "tax_rank", "selected" or "selected_function").',
        metavar='STR',
        default='tax_rank')

    parent_parser_t = argparse.ArgumentParser(add_help=False)
    parent_parser_t.add_argument(
        '--taxon-rank',
        dest='taxon_rank',
        required=False,
        help='Taxon rank to cluster observation names of EsMeCaTa together (default "phylum").',
        metavar='INPUT_TAXON',
        default='phylum')

    parent_parser_taxa_list = argparse.ArgumentParser(add_help=False)
    parent_parser_taxa_list.add_argument(
        '--taxa-list',
        dest='taxa_list',
        required=False,
        help='When using value "selected" for option "grouping", you have to give this file indicating the different groups of taxon to compare.',
        metavar='STR',
        default='tax_rank')

    parent_parser_function_list = argparse.ArgumentParser(add_help=False)
    parent_parser_function_list.add_argument(
        '--function-list',
        dest='function_list',
        required=False,
        help='When using value "selected_function" for option "grouping", you have to give this file indicating the different groups of functions to compare.',
        metavar='STR',
        default='function_list')

    parent_parser_e = argparse.ArgumentParser(add_help=False)
    parent_parser_e.add_argument(
        '--enzyme',
        dest='enzyme_file',
        required=False,
        help='Pathanme to enzyme.dat file of expasy. If no path is given, it will be downloaded and put in the output folder.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        '--go',
        dest='go_file',
        required=False,
        help='Pathname to go-basic.obo file. If no path is given, it will be downloaded and put in the output folder.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_taxon_id = argparse.ArgumentParser(add_help=False)
    parent_parser_taxon_id.add_argument(
        '--taxon-id',
        dest='taxon_id',
        required=False,
        help='Pathname to taxon ID file indicating taxonomic ID associated with organisms name. Required for tax_rank analysis if the annotation input is a file and not an esmecata annotation folder.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_ko = argparse.ArgumentParser(add_help=False)
    parent_parser_ko.add_argument(
        '--ko',
        dest='ko',
        help='If KEGG Orthologs are in the function table, add this parameter to get their names (require bioservices).',
        required=False,
        action='store_true',
        default=False)

    parent_parser_orsumMinTermSize = argparse.ArgumentParser(add_help=False)
    parent_parser_orsumMinTermSize.add_argument(
        '--orsumMinTermSize',
        dest='orsumMinTermSize',
        required=False,
        help='MinTermSize of orsum.',
        metavar='INT',
        default=None)

    parent_parser_gseapyCutOff = argparse.ArgumentParser(add_help=False)
    parent_parser_gseapyCutOff.add_argument(
        '--gseapyCutOff',
        dest='gseapyCutOff',
        required=False,
        type=float,
        help='Adjust-Pval cutoff for gseapy enrichr, default: 0.05 (--cut-off argument of gseapy).',
        metavar='FLOAT',
        default=0.05)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    gseapy_taxon_parser = subparsers.add_parser(
        'gseapy_enrichr',
        help='Extract enriched functions from groups (either chosen from tax_rank or manually selected) using gseapy enrichr and orsum.',
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_grouping, parent_parser_t, parent_parser_taxa_list, parent_parser_function_list,
            parent_parser_e, parent_parser_g, parent_parser_orsumMinTermSize, parent_parser_gseapyCutOff, parent_parser_taxon_id, parent_parser_ko
            ],
        allow_abbrev=False)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    is_valid_dir(args.output)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.output, f'esmecata_gseapy_{args.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args.cmd == 'gseapy_enrichr':
        taxon_rank_annotation_enrichment(args.input_file_folder, args.output, args.grouping, taxon_rank=args.taxon_rank, taxa_lists_file=args.taxa_list,
                                         function_lists_file=args.function_list,
                                         enzyme_data_file=args.enzyme_file, go_basic_obo_file=args.go_file, orsum_minterm_size=args.orsumMinTermSize,
                                         selected_adjust_pvalue_cutoff=args.gseapyCutOff, taxon_id_file=args.taxon_id, ko_annotation=args.ko)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
