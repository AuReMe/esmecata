import pandas as pd
import os
import shutil
import subprocess

from esmecata.gseapy.gseapy_orsum import taxon_rank_annotation_enrichment


def test_taxon_rank_annotation_enrichment_selected():
    input_folder = 'annotation_output'
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = 'taxa_list.tsv'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, taxa_lists_file=taxa_list_file, orsum_minterm_size=4)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_cli():
    input_folder = 'annotation_output'
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = 'taxa_list.tsv'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--taxa-list', taxa_list_file, '--orsumMinTermSize', '4'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank():
    input_folder = 'annotation_output'
    output_folder = 'output_folder'
    grouping = 'tax_rank'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, orsum_minterm_size=4)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank_cli():
    input_folder = 'annotation_output'
    output_folder = 'output_folder'
    grouping = 'tax_rank'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--orsumMinTermSize', '4'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)
    