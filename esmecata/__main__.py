import argparse
import sys

from esmecata.proteomes import retrieve_proteomes
from esmecata.clustering import make_clustering
from esmecata.annotation import annotate_proteins
from esmecata.utils import range_limited_float_type
from esmecata import __version__ as VERSION

MESSAGE = '''
From taxonomy to metabolism using Uniprot.
'''
REQUIRES = '''
Requires: mmseqs2 and an internet connection (for REST and SPARQL queries, except if you have a local Uniprot SPARQL endpoint).
'''

def main():
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
        help='Input taxon file (excel, tsv or csv) containing a column associating ID to a taxonomy (separated by ;).',
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
        metavar='BUSCO')
    parent_parser_taxadb = argparse.ArgumentParser(add_help=False)
    parent_parser_taxadb.add_argument(
        '--ignore-taxadb-update',
        dest='ignore_taxadb_update',
        help='If you have a not up-to-date version of the NCBI taxonomy database with ete3, use this option to bypass the warning message and use the old version.',
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
        help='Proportion [0 to 1] of proteomes required to occur in a proteins cluster for that cluster to be kept in core proteome assembly.',
        required=False,
        type=range_limited_float_type,
        default=1)
    parent_parser_propagate = argparse.ArgumentParser(add_help=False)
    parent_parser_propagate.add_argument(
        '-p',
        '--propagate',
        dest='propagate_annotation',
        help='Proportion [0 to 1] of the reccurence of an annotation to be propagated from the protein of a cluster to the reference protein of the cluster. 0 mean the annotaitons from all proteins are propagated to the reference and 1 only the annotation occurring in all the proteins of the cluster.',
        required=False,
        type=range_limited_float_type,
        default=1)
    parent_parser_sparql = argparse.ArgumentParser(add_help=False)
    parent_parser_sparql.add_argument(
        '-s',
        '--sparql',
        dest='sparql',
        help='Use sparql endpoint instead of REST queries on Uniprot.',
        required=False)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    proteomes_parser = subparsers.add_parser(
        'proteomes',
        help='Download proteomes associated to taxon from Uniprot Proteomes.',
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b,
            parent_parser_taxadb, parent_parser_sparql
        ])
    clustering_parser = subparsers.add_parser(
        'clustering',
        help='Cluster the proteins of the different proteomes of a taxon into a single set of representative shared proteins.',
        parents=[
            parent_parser_i_clustering_folder, parent_parser_o, parent_parser_c, parent_parser_thr
        ])
    annotation_parser = subparsers.add_parser(
        'annotation',
        help='Retrieve protein annotations from Uniprot.',
        parents=[
            parent_parser_i_annotation_folder, parent_parser_o, parent_parser_sparql,
            parent_parser_propagate
        ])

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.cmd in ['proteomes', 'annotation']:
        if args.sparql is None:
            uniprot_sparql_endpoint = None
        elif args.sparql == 'uniprot':
            uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
        else:
            uniprot_sparql_endpoint = args.sparql

    if args.cmd == 'proteomes':
        if args.busco:
            busco_score = 100*args.busco
        else:
            busco_score = None

    if args.cmd == 'proteomes':
        retrieve_proteomes(args.input, args.output, busco_score, args.ignore_taxadb_update, uniprot_sparql_endpoint)
    elif args.cmd == 'clustering':
        make_clustering(args.input, args.output, args.cpu, args.threshold_clustering)
    elif args.cmd == 'annotation':
        annotate_proteins(args.input, args.output, uniprot_sparql_endpoint, args.propagate_annotation)

if __name__ == '__main__':
    main()
