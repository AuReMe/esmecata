
import argparse

import pkg_resources

from esmecata.retrieve_proteome import retrieve_proteome
from esmecata.coreproteome import create_coreproteome

VERSION = pkg_resources.get_distribution("metage2metabo").version

MESSAGE = """
From taxonomy to metabolism
"""
REQUIRES = """
Requires: mmseqs2
"""

def main():
    parser = argparse.ArgumentParser(
        "esmecata",
        description=MESSAGE + " For specific help on each subcommand use: esmecata {cmd} --help",
        epilog=REQUIRES
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n")

    parent_parser_i_taxon = argparse.ArgumentParser(add_help=False)
    parent_parser_i_taxon.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        help="input taxon file (excel, tsv or csv) containing a column associating ID to a taxonomy (separated by ;)",
        metavar="INPUT_FILE")

    parent_parser_i_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_i_folder.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        help="This input folder of clustering is the output folder of proteomes command",
        metavar="INPUT_DIR")

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="output directory path",
        metavar="OUPUT_DIR")

    parent_parser_b = argparse.ArgumentParser(add_help=False)
    parent_parser_b.add_argument(
        "-b",
        "--busco",
        dest="busco",
        required=False,
        help="busco percentage between 0 and 100. This will remove all the proteomes without BSUCO score and the score before the selected percentage.",
        metavar="BUSCO")
    parent_parser_taxadb = argparse.ArgumentParser(add_help=False)
    parent_parser_taxadb.add_argument(
        "--ignore-taxadb-update",
        dest="ignore_taxadb_update",
        help="If you have a not up-to-date version of NCBI taxonomy database with ete3, use this option to use this version and bypass the warning message.",
        required=False,
        action="store_true",
        default=None)
    parent_parser_c = argparse.ArgumentParser(add_help=False)
    parent_parser_c.add_argument(
        "-c",
        "--cpu",
        dest="cpu",
        help="cpu number for multiprocessing",
        required=False,
        type=int,
        default=1)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest="cmd")

    proteomes_parser = subparsers.add_parser(
        "proteomes",
        help="Download proteomes associated to taxon from Uniprot Proteomes.",
        parents=[
            parent_parser_i_taxon, parent_parser_o, parent_parser_b,
            parent_parser_taxadb
        ])
    clustering_parser = subparsers.add_parser(
        "clustering",
        help="Cluster proteins proteomes for a taxon into a single set of representative shared proteins.",
        parents=[
            parent_parser_i_folder, parent_parser_o, parent_parser_c
        ])

    args = parser.parse_args()
    print(args)
    if args.cmd == "proteomes":
        retrieve_proteome(args.input, args.output, args.busco, args.ignore_taxadb_update)
    elif args.cmd == "clustering":
        create_coreproteome(args.input, args.output, args.cpu)


if __name__ == "__main__":
    main()