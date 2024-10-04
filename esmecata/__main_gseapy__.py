# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
Requires: pronto, gseapy and orsum
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
        dest='input_folder',
        required=True,
        help='EsMeCaTa annotation output folder.',
        metavar='INPUT_FOLDER')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_e = argparse.ArgumentParser(add_help=False)
    parent_parser_e.add_argument(
        '-e',
        '--enzyme',
        dest='enzyme_file',
        required=False,
        help='Pathanme to enzyme.dat file of expasy. If no path is given, it will be downloaded and put in the output folder.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_g = argparse.ArgumentParser(add_help=False)
    parent_parser_g.add_argument(
        '-g',
        '--go',
        dest='go_file',
        required=False,
        help='Pathname to go-basic.obo file. If no path is given, it will be downloaded and put in the output folder.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_t = argparse.ArgumentParser(add_help=False)
    parent_parser_t.add_argument(
        '-t',
        '--taxon',
        dest='taxon_rank',
        required=False,
        help='Taxon rank to cluster observation names of EsMeCaTa together (default "phylum").',
        metavar='INPUT_TAXON',
        default='phylum')

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    gseapy_taxon_parser = subparsers.add_parser(
        'gseapy_taxon',
        help='Extract enriched functions from taxon using gseapy and orsum.',
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_e, parent_parser_g, parent_parser_t
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

    if args.cmd == 'gseapy_taxon':
        taxon_rank_annotation_enrichment(args.input_folder, args.output, enzyme_data_file=args.enzyme_file,
                                         go_basic_obo_file=args.go_file, taxon_rank=args.taxon_rank)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
