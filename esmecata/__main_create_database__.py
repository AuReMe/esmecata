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

from esmecata.precomputed.create_database import create_database_from_esmecata_run, merge_db_files
from esmecata.utils import is_valid_dir
from esmecata import __version__ as VERSION

MESSAGE = '''
Create database file from esmecata run.
'''
REQUIRES = '''
Requires: esmecata, pandas.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'esmecata_create_db',
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
        help='EsMeCaTa output folder.',
        metavar='INPUT_FOLDER')

    parent_parser_i_zip = argparse.ArgumentParser(add_help=False)
    parent_parser_i_zip.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Multiple path of esmecata database zip files separated by ",".',
        metavar='INPUT_FILES')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        '-c',
        '--core',
        dest='core',
        help='Number of CPU-cores for multiprocessing.',
        required=False,
        type=int,
        default=1)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    from_workflow = subparsers.add_parser(
        'from_workflow',
        help='Create database from esmecata workflow output.',
        parents=[
            parent_parser_i, parent_parser_o, parent_parser_c
            ],
        allow_abbrev=False)

    merge_db = subparsers.add_parser(
        'merge_db',
        help='Merge multiple zip files corresponding to EsMeCaTa databases.',
        parents=[
            parent_parser_i_zip, parent_parser_o
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
    log_file_path = os.path.join(args.output, f'esmecata_create_db_{args.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args.cmd == 'from_workflow':
        esmecata_proteomes_folder = os.path.join(args.input, '0_proteomes')
        esmecata_clustering_folder = os.path.join(args.input, '1_clustering')
        esmecata_annotation_folder = os.path.join(args.input, '2_annotation')
        create_database_from_esmecata_run(esmecata_proteomes_folder, esmecata_clustering_folder, esmecata_annotation_folder, args.output, args.core)
    elif args.cmd == 'merge_db':
        if ',' in args.input:
            list_db_files = args.input.split(',')
        if ' ' in args.input:
            list_db_files = args.input.split(' ')
        merge_db_files(list_db_files, args.output)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
