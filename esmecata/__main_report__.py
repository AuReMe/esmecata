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

from esmecata.report.workflow_create_report import run_create_workflow_report, run_proteomes_report_creation, run_clustering_report_creation, \
                                                                    run_annotation_report_creation
from esmecata.utils import is_valid_dir
from esmecata import __version__ as VERSION

MESSAGE = '''
Create report files from esmecata output folder.
'''
REQUIRES = '''
Requires: datapane, plotly, kaleido, ontosunburst.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'esmecata_report',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION + '\n')

    parent_parser_i = argparse.ArgumentParser(add_help=False)
    parent_parser_i.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='EsMeCaTa input file.',
        metavar='INPUT_FILE')

    parent_parser_f = argparse.ArgumentParser(add_help=False)
    parent_parser_f.add_argument(
        '-f',
        '--folder',
        dest='input_folder',
        required=True,
        help='EsMeCaTa output folder.',
        metavar='INPUT_FOLDER')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_svg = argparse.ArgumentParser(add_help=False)
    parent_parser_svg.add_argument(
        '--svg',
        dest='create_svg',
        action='store_true',
        required=False,
        help='Create SVG files of picture.')

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    create_report_parser = subparsers.add_parser(
        'create_report',
        help='Create report from esmecata output folder of workflow or precomputed subcommands.',
        parents=[
            parent_parser_i, parent_parser_f, parent_parser_o, parent_parser_svg
            ],
        allow_abbrev=False)

    create_report_proteomes_parser = subparsers.add_parser(
        'create_report_proteomes',
        help='Create report from esmecata output folder of proteomes subcommand.',
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_svg
            ],
        allow_abbrev=False)

    create_report_clustering_parser = subparsers.add_parser(
        'create_report_clustering',
        help='Create report from esmecata output folder of clustering subcommand.',
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_svg
            ],
        allow_abbrev=False)

    create_report_annotation_parser = subparsers.add_parser(
        'create_report_annotation',
        help='Create report from esmecata output folder of annotation subcommand.',
        parents=[
            parent_parser_f, parent_parser_o, parent_parser_svg
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
    log_file_path = os.path.join(args.output, f'esmecata_analysis_{args.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args.cmd == 'create_report':
        run_create_workflow_report(args.input, args.input_folder, args.output, args.create_svg)
    elif args.cmd == 'create_report_proteomes':
        run_proteomes_report_creation(args.input_folder, args.output, args.create_svg)
    elif args.cmd == 'create_report_clustering':
        run_clustering_report_creation(args.input_folder, args.output, args.create_svg)
    elif args.cmd == 'create_report_annotation':
        run_annotation_report_creation(args.input_folder, args.output, args.create_svg)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
