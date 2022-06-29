# Copyright (C) 2021-2022 Arnaud Belcour - Inria Dyliss
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

from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins
from esmecata.workflow import perform_workflow
from esmecata.utils import limited_integer_type, range_limited_float_type, is_valid_dir
from esmecata import __version__ as VERSION

MESSAGE = '''
From taxonomic affiliation to metabolism using Uniprot.
'''
REQUIRES = '''
Requires: mmseqs2 and an internet connection (for REST and SPARQL queries, except if you have a local Uniprot SPARQL endpoint).
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'esmecata',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION + '\n')

    parent_parser_i_taxon = argparse.ArgumentParser(add_help=False)
    parent_parser_i_taxon.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input taxon file (excel, tsv or csv) containing a column associating ID to a taxonomic affiliation (separated by ;).',
        metavar='INPUT_FILE')

    parent_parser_i_clustering_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_i_clustering_folder.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='This input folder of clustering is the output folder of proteomes command.',
        metavar='INPUT_DIR')

    parent_parser_i_annotation_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_i_annotation_folder.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='This input folder of annotation is the output folder of clustering command.',
        metavar='INPUT_DIR')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_b = argparse.ArgumentParser(add_help=False)
    parent_parser_b.add_argument(
        '-b',
        '--busco',
        dest='busco',
        required=False,
        type=range_limited_float_type,
        help='BUSCO percentage between 0 and 1. This will remove all the proteomes without BUSCO score and the score before the selected ratio of completion.',
        metavar='BUSCO',
        default=0.8)
    parent_parser_taxadb = argparse.ArgumentParser(add_help=False)
    parent_parser_taxadb.add_argument(
        '--ignore-taxadb-update',
        dest='ignore_taxadb_update',
        help='If you have a not up-to-date version of the NCBI taxonomy database with ete3, use this option to bypass the warning message and use the old version.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_all_proteomes = argparse.ArgumentParser(add_help=False)
    parent_parser_all_proteomes.add_argument(
        '--all-proteomes',
        dest='all_proteomes',
        help='Download all proteomes associated with a taxon even if they are no reference proteomes.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_limit_maximal_number_proteomes = argparse.ArgumentParser(add_help=False)
    parent_parser_limit_maximal_number_proteomes.add_argument(
        '-l',
        '--limit-proteomes',
        dest='limit_maximal_number_proteomes',
        required=False,
        type=limited_integer_type,
        help='Choose the maximal number of proteomes after which the tool will select a subset of proteomes instead of using all the available proteomes (default is 99).',
        default=99)
    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        '-c',
        '--cpu',
        dest='cpu',
        help='CPU number for multiprocessing.',
        required=False,
        type=int,
        default=1)
    parent_parser_thr = argparse.ArgumentParser(add_help=False)
    parent_parser_thr.add_argument(
        '-t',
        '--threshold',
        dest='threshold_clustering',
        help='Proportion [0 to 1] of proteomes required to occur in a proteins cluster for that cluster to be kept in core proteome assembly.',
        required=False,
        type=range_limited_float_type,
        default=0.95)
    parent_parser_mmseqs_options = argparse.ArgumentParser(add_help=False)
    parent_parser_mmseqs_options.add_argument(
        '-m',
        '--mmseqs',
        dest='mmseqs_options',
        help='String containing mmseqs options for cluster command (except --threads which is already set by --cpu command and -v). If nothing is given, esmecata will used the option "--min-seq-id 0.3 -c 0.8"',
        required=False,
        type=str,
        default=None)
    parent_parser_linclust = argparse.ArgumentParser(add_help=False)
    parent_parser_linclust.add_argument(
        '--linclust',
        dest='linclust',
        help='Use mmseqs linclust (clustering in lienar time) to cluster proteins sequences. It is faster than mmseqs cluster (default behaviour) but less sensitive.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_remove_tmp = argparse.ArgumentParser(add_help=False)
    parent_parser_remove_tmp.add_argument(
        '--remove-tmp',
        dest='remove_tmp',
        help='Delete tmp files to limit the disk space used: files in tmp_proteome for esmecata proteomes and files created by mmseqs (in mmseqs_tmp).',
        required=False,
        action='store_true',
        default=None)
    parent_parser_propagate = argparse.ArgumentParser(add_help=False)
    parent_parser_propagate.add_argument(
        '-p',
        '--propagate',
        dest='propagate_annotation',
        help='Proportion [0 to 1] of the occurrence of an annotation to be propagated from the protein of a cluster to the reference protein of the cluster. 0 mean the annotations from all proteins are propagated to the reference and 1 only the annotation occurring in all the proteins of the cluster (default).',
        required=False,
        type=range_limited_float_type,
        default=1)
    parent_parser_uniref = argparse.ArgumentParser(add_help=False)
    parent_parser_uniref.add_argument(
        '--uniref',
        dest='uniref',
        help='Use uniref cluster to extract more annotations from the representative member of the cluster associated with the proteins. Needs the --sparql option.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_expression = argparse.ArgumentParser(add_help=False)
    parent_parser_expression.add_argument(
        '--expression',
        dest='expression',
        help='Extract expression information associated with the proteins. Needs the --sparql option.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_sparql = argparse.ArgumentParser(add_help=False)
    parent_parser_sparql.add_argument(
        '-s',
        '--sparql',
        dest='sparql',
        help='Use sparql endpoint instead of REST queries on Uniprot.',
        required=False)
    parent_parser_rank_limit = argparse.ArgumentParser(add_help=False)
    parent_parser_rank_limit.add_argument(
        '-r',
        '--rank-limit',
        dest='rank_limit',
        required=False,
        help='This option limit the rank used by the tool for searching for proteomes. The given rank and all the superior ranks will be ignored. Look at the readme for more information (and a list of possible rank).',
        default=None)
    parent_parser_minimal_number_proteomes = argparse.ArgumentParser(add_help=False)
    parent_parser_minimal_number_proteomes.add_argument(
        '--minimal-nb-proteomes',
        dest='minimal_number_proteomes',
        required=False,
        type=limited_integer_type,
        help='Choose the minimal number of proteomes to be selected by EsMeCaTa. If a taxon has less proteomes, it will be ignored and a higher taxonomic rank will be used. Default is 1.',
        default=1)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    proteomes_parser = subparsers.add_parser(
        'proteomes',
        help='Download proteomes associated with taxon from Uniprot Proteomes.',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b,
            parent_parser_taxadb, parent_parser_all_proteomes, parent_parser_sparql,
            parent_parser_remove_tmp, parent_parser_limit_maximal_number_proteomes,
            parent_parser_rank_limit, parent_parser_minimal_number_proteomes
            ],
        allow_abbrev=False)
    clustering_parser = subparsers.add_parser(
        'clustering',
        help='Cluster the proteins of the different proteomes of a taxon into a single set of representative shared proteins.',
        parents=[
            parent_parser_i_clustering_folder, parent_parser_o, parent_parser_c,
            parent_parser_thr, parent_parser_mmseqs_options, parent_parser_linclust,
            parent_parser_remove_tmp
            ],
        allow_abbrev=False)
    annotation_parser = subparsers.add_parser(
        'annotation',
        help='Retrieve protein annotations from Uniprot.',
        parents=[
            parent_parser_i_annotation_folder, parent_parser_o, parent_parser_sparql,
            parent_parser_propagate, parent_parser_uniref, parent_parser_expression
            ],
        allow_abbrev=False)
    workflow_parser = subparsers.add_parser(
        'workflow',
        help='Run all esmecata steps (proteomes, clustering and annotation).',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b, parent_parser_c,
            parent_parser_taxadb, parent_parser_all_proteomes, parent_parser_sparql,
            parent_parser_remove_tmp, parent_parser_limit_maximal_number_proteomes,
            parent_parser_thr, parent_parser_mmseqs_options, parent_parser_linclust,
            parent_parser_propagate, parent_parser_uniref, parent_parser_expression,
            parent_parser_rank_limit, parent_parser_minimal_number_proteomes
            ],
        allow_abbrev=False)

    args = parser.parse_args()

    is_valid_dir(args.output)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.output, f'esmecata_{args.cmd}.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.cmd in ['proteomes', 'annotation', 'workflow']:
        if args.sparql is None:
            uniprot_sparql_endpoint = None
        elif args.sparql == 'uniprot':
            uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
        else:
            uniprot_sparql_endpoint = args.sparql

    if args.cmd in ['proteomes', 'workflow']:
        if args.busco is not None:
            busco_score = 100*args.busco

    if args.cmd == 'proteomes':
        retrieve_proteomes(args.input, args.output, busco_score, args.ignore_taxadb_update,
                            args.all_proteomes, uniprot_sparql_endpoint, args.remove_tmp,
                            args.limit_maximal_number_proteomes, args.rank_limit, args.minimal_number_proteomes)
    elif args.cmd == 'clustering':
        make_clustering(args.input, args.output, args.cpu, args.threshold_clustering, args.mmseqs_options, args.linclust, args.remove_tmp)
    elif args.cmd == 'annotation':
        annotate_proteins(args.input, args.output, uniprot_sparql_endpoint, args.propagate_annotation, args.uniref, args.expression)
    elif args.cmd == 'workflow':
        perform_workflow(args.input, args.output, busco_score, args.ignore_taxadb_update,
                            args.all_proteomes, uniprot_sparql_endpoint, args.remove_tmp,
                            args.limit_maximal_number_proteomes, args.rank_limit,
                            args.cpu, args.threshold_clustering, args.mmseqs_options,
                            args.linclust, args.propagate_annotation, args.uniref,
                            args.expression, args.minimal_number_proteomes)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
