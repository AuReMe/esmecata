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

from esmecata.proteomes import check_proteomes, retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins
from esmecata.workflow import perform_workflow, perform_workflow_eggnog
from esmecata.eggnog import annotate_with_eggnog
from esmecata.utils import limited_integer_type, range_limited_float_type, is_valid_dir
from esmecata.analysis import perform_analysis
from esmecata import __version__ as VERSION

MESSAGE = '''
From taxonomic affiliation to metabolism using Uniprot.
'''
REQUIRES = '''
Requires: mmseqs2 and an internet connection (for REST and SPARQL queries, except if you have a local Uniprot SPARQL endpoint).
Annotation can be performed with UniProt or eggnog-mapper (which is then a requirement if the option is selected).
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

    parent_parser_i_analysis_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_i_analysis_folder.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='This input folder of analysis is the output folder of annotation command.',
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
    parent_parser_update_affiliation = argparse.ArgumentParser(add_help=False)
    parent_parser_update_affiliation.add_argument(
        '--update-affiliations',
        dest='update_affiliations',
        help='''If the taxonomic affiliations were assigned from an outdated taxonomic database, this can lead to taxon not be found in ete3 database. \
            This option tries to udpate the taxonomic affiliations using the lowest taxon name.''',
        required=False,
        action='store_true',
        default=None)
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
        help='Proportion [0 to 1] of proteomes required to occur in a proteins cluster for that cluster to be kept in core proteome assembly. Default is 0.5.',
        required=False,
        type=range_limited_float_type,
        default=0.5)
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
        help='Use mmseqs linclust (clustering in linear time) to cluster proteins sequences. It is faster than mmseqs cluster (default behaviour) but less sensitive.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_remove_tmp = argparse.ArgumentParser(add_help=False)
    parent_parser_remove_tmp.add_argument(
        '--remove-tmp',
        dest='remove_tmp',
        help='Delete tmp files to limit the disk space used: files created by mmseqs (in mmseqs_tmp).',
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
        help='''This option limits the rank used when searching for proteomes. All the ranks superior to the given rank will be ignored. \
            For example, if 'family' is given, only taxon ranks inferior or equal to family will be kept. \
            Look at the readme for more information (and a list of rank names).''',
        default=None)
    parent_parser_minimal_number_proteomes = argparse.ArgumentParser(add_help=False)
    parent_parser_minimal_number_proteomes.add_argument(
        '--minimal-nb-proteomes',
        dest='minimal_number_proteomes',
        required=False,
        type=limited_integer_type,
        help='Choose the minimal number of proteomes to be selected by EsMeCaTa. If a taxon has less proteomes, it will be ignored and a higher taxonomic rank will be used. Default is 5.',
        default=5)
    parent_parser_annotation_file = argparse.ArgumentParser(add_help=False)
    parent_parser_annotation_file.add_argument(
        '--annotation-files',
        dest='annotation_files',
        required=False,
        help='Use UniProt annotation files (uniprot_trembl.txt and uniprot_sprot.txt) to avoid querying UniProt REST API. Need both paths to these files separated by a ",".',
        default=None)
    parent_parser_bioservices = argparse.ArgumentParser(add_help=False)
    parent_parser_bioservices.add_argument(
        '--bioservices',
        dest='option_bioservices',
        help='Use bioservices instead of esmecata functions for protein annotation.',
        required=False,
        action='store_true',
        default=None)
    parent_parser_eggnog_database = argparse.ArgumentParser(add_help=False)
    parent_parser_eggnog_database.add_argument(
        '-e',
        '--eggnog',
        dest='eggnog_database',
        help='Path to eggnog database.',
        required=True)
    parent_parser_nb_digit = argparse.ArgumentParser(add_help=False)
    parent_parser_nb_digit.add_argument(
        '--nb-digit',
        dest='nb_digit',
        required=False,
        type=limited_integer_type,
        help='Number of the digit to keep on the clustermap (1, 2, 3 or 4). Defualt 3.',
        default=3)
    parent_parser_taxon_rank = argparse.ArgumentParser(add_help=False)
    parent_parser_taxon_rank.add_argument(
        '--taxon-rank',
        dest='taxon_rank',
        required=False,
        help='Taxon rank to merge organisms (default family).',
        default='family')
    parent_parser_eggnog_tmp_dir = argparse.ArgumentParser(add_help=False)
    parent_parser_eggnog_tmp_dir.add_argument(
        '--eggnog-tmp',
        dest='eggnog_tmp_dir',
        help='Path to eggnog tmp dir.',
        required=False,
        default=None)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    check_parser = subparsers.add_parser(
        'check',
        help='Check proteomes associated with taxon in Uniprot Proteomes database.',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b,
            parent_parser_taxadb, parent_parser_all_proteomes, parent_parser_sparql,
            parent_parser_limit_maximal_number_proteomes, parent_parser_rank_limit,
            parent_parser_minimal_number_proteomes, parent_parser_update_affiliation,
            parent_parser_bioservices
            ],
        allow_abbrev=False)
    proteomes_parser = subparsers.add_parser(
        'proteomes',
        help='Download proteomes associated with taxon from Uniprot Proteomes.',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b,
            parent_parser_taxadb, parent_parser_all_proteomes, parent_parser_sparql,
            parent_parser_limit_maximal_number_proteomes, parent_parser_rank_limit,
            parent_parser_minimal_number_proteomes, parent_parser_update_affiliation,
            parent_parser_bioservices
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
    annotation_uniprot_parser = subparsers.add_parser(
        'annotation_uniprot',
        help='Retrieve protein annotations from Uniprot.',
        parents=[
            parent_parser_i_annotation_folder, parent_parser_o, parent_parser_sparql,
            parent_parser_propagate, parent_parser_uniref, parent_parser_expression,
            parent_parser_annotation_file, parent_parser_bioservices
            ],
        allow_abbrev=False)
    annotation_eggnog_parser = subparsers.add_parser(
        'annotation',
        help='Annotate protein clusters using eggnog-mapper.',
        parents=[
            parent_parser_i_annotation_folder, parent_parser_o, parent_parser_eggnog_database,
            parent_parser_c, parent_parser_eggnog_tmp_dir
            ],
        allow_abbrev=False)
    workflow_uniprot_parser = subparsers.add_parser(
        'workflow_uniprot',
        help='Run all esmecata steps (proteomes, clustering and annotation).',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b, parent_parser_c,
            parent_parser_taxadb, parent_parser_all_proteomes, parent_parser_sparql,
            parent_parser_remove_tmp, parent_parser_limit_maximal_number_proteomes,
            parent_parser_thr, parent_parser_mmseqs_options, parent_parser_linclust,
            parent_parser_propagate, parent_parser_uniref, parent_parser_expression,
            parent_parser_rank_limit, parent_parser_minimal_number_proteomes,
            parent_parser_annotation_file, parent_parser_update_affiliation,
            parent_parser_bioservices
            ],
        allow_abbrev=False)
    workflow_eggnog_parser = subparsers.add_parser(
        'workflow',
        help='Run all esmecata steps (proteomes, clustering and annotation with eggnog-mapper).',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_eggnog_database,
            parent_parser_b, parent_parser_c, parent_parser_taxadb,
            parent_parser_all_proteomes, parent_parser_sparql, parent_parser_remove_tmp,
            parent_parser_limit_maximal_number_proteomes, parent_parser_thr, parent_parser_mmseqs_options,
            parent_parser_linclust, parent_parser_rank_limit, parent_parser_minimal_number_proteomes,
            parent_parser_update_affiliation, parent_parser_bioservices, parent_parser_eggnog_tmp_dir
            ],
        allow_abbrev=False)
    analysis_parser = subparsers.add_parser(
        'analysis',
        help='Create clustermap for EC.',
        parents=[
            parent_parser_i_analysis_folder, parent_parser_o, parent_parser_taxon_rank,
            parent_parser_nb_digit
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

    if args.cmd in ['proteomes', 'annotation_uniprot', 'workflow_uniprot', 'workflow', 'check']:
        if args.sparql is None:
            uniprot_sparql_endpoint = None
        elif args.sparql == 'uniprot':
            uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
        else:
            uniprot_sparql_endpoint = args.sparql

    if args.cmd in ['proteomes', 'workflow_uniprot', 'workflow', 'check']:
        if args.busco is not None:
            busco_score = 100*args.busco

    if args.cmd == 'proteomes':
        retrieve_proteomes(args.input, args.output, busco_score, args.ignore_taxadb_update,
                            args.all_proteomes, uniprot_sparql_endpoint, args.limit_maximal_number_proteomes,
                            args.rank_limit, args.minimal_number_proteomes, args.update_affiliations,
                            args.option_bioservices)
    if args.cmd == 'check':
        check_proteomes(args.input, args.output, busco_score, args.ignore_taxadb_update,
                            args.all_proteomes, uniprot_sparql_endpoint, args.limit_maximal_number_proteomes,
                            args.rank_limit, args.minimal_number_proteomes, args.update_affiliations,
                            args.option_bioservices)
    elif args.cmd == 'clustering':
        make_clustering(args.input, args.output, args.cpu, args.threshold_clustering, args.mmseqs_options, args.linclust, args.remove_tmp)
    elif args.cmd == 'annotation_uniprot':
        annotate_proteins(args.input, args.output, uniprot_sparql_endpoint,
                        args.propagate_annotation, args.uniref, args.expression,
                        args.annotation_files, args.option_bioservices)
    elif args.cmd == 'workflow_uniprot':
        perform_workflow(args.input, args.output, busco_score, args.ignore_taxadb_update,
                            args.all_proteomes, uniprot_sparql_endpoint, args.remove_tmp,
                            args.limit_maximal_number_proteomes, args.rank_limit,
                            args.cpu, args.threshold_clustering, args.mmseqs_options,
                            args.linclust, args.propagate_annotation, args.uniref,
                            args.expression, args.minimal_number_proteomes, args.annotation_files,
                            args.update_affiliations, args.option_bioservices)
    elif args.cmd == 'annotation':
        annotate_with_eggnog(args.input, args.output, args.eggnog_database, args.cpu,
                             args.eggnog_tmp_dir)
    elif args.cmd == 'workflow':
        perform_workflow_eggnog(args.input, args.output, args.eggnog_database, busco_score,
                                args.ignore_taxadb_update, args.all_proteomes, uniprot_sparql_endpoint,
                                args.remove_tmp, args.limit_maximal_number_proteomes, args.rank_limit,
                                args.cpu, args.threshold_clustering, args.mmseqs_options,
                                args.linclust, args.minimal_number_proteomes, args.update_affiliations,
                                args.option_bioservices, args.eggnog_tmp_dir)
    elif args.cmd == 'analysis':
        perform_analysis(args.input, args.output, args.taxon_rank, args.nb_digit)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
