import pandas as pd
import os
import shutil
import subprocess

from esmecata.gseapy.gseapy_orsum import taxon_rank_annotation_enrichment

ANNOTATION_NAMES_JSON = os.path.join('test_data', 'annotation_names.json')


def test_taxon_rank_annotation_enrichment_selected():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, taxa_lists_file=taxa_list_file, annot_names_file=ANNOTATION_NAMES_JSON, orsum_minterm_size=4)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_download():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, taxa_lists_file=taxa_list_file, annot_names_file='download', orsum_minterm_size=4)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_cli():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--taxa-list', taxa_list_file, '--orsumMinTermSize', '4',
                     '--annot-names', ANNOTATION_NAMES_JSON])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_cutoff():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, taxa_lists_file=taxa_list_file, annot_names_file=ANNOTATION_NAMES_JSON, orsum_minterm_size=4, selected_adjust_pvalue_cutoff=0.1)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_cutoff_cli():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--taxa-list', taxa_list_file, '--annot-names', ANNOTATION_NAMES_JSON,
                     '--orsumMinTermSize', '4', '--gseapyCutOff', '0.1'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'tax_rank'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, annot_names_file=ANNOTATION_NAMES_JSON, orsum_minterm_size=4)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank_cli():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'tax_rank'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping,
                     '--annot-names', ANNOTATION_NAMES_JSON, '--orsumMinTermSize', '4'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank_cli_input_file():
    input_folder = os.path.join('test_data', 'annotation_output', 'dataset_annotation.tsv')
    taxon_id_file = os.path.join('test_data', 'annotation_output', 'proteome_tax_id.tsv')

    output_folder = 'output_folder'
    grouping = 'tax_rank'
    expected_terms = ['1.1.1.86', '2.8.4.1']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping,
                     '--annot-names', ANNOTATION_NAMES_JSON, '--orsumMinTermSize', '4', '--taxon-id', taxon_id_file])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_selected_cli_input_file_picrust2():
    input_folder = os.path.join('test_data', 'picrust_test_predicted_func.tsv')

    output_folder = 'output_folder'
    grouping = 'selected'
    taxa_list_file = os.path.join('test_data', 'taxa_list.tsv')
    expected_terms = ['3.4.21.19', 'ko:K07785']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--taxa-list', taxa_list_file,
                     '--annot-names', ANNOTATION_NAMES_JSON, '--orsumMinTermSize', '4'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_tax_rank_cli_input_file_picrust2():
    input_annot_file = os.path.join('test_data', 'picrust_test_predicted_func.tsv')
    taxon_id_file = os.path.join('test_data', 'picrust_taxon_id.tsv')

    output_folder = 'output_folder'
    grouping = 'tax_rank'

    expected_terms = ['3.4.21.19', 'ko:K07785']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_annot_file, '-o', output_folder, '--grouping', grouping, '--taxon-id', taxon_id_file,
                     '--annot-names', ANNOTATION_NAMES_JSON, '--orsumMinTermSize', '4'])

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_function_selected():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected_function'
    function_list_file = os.path.join('test_data', 'function_list.tsv')
    expected_terms = ['org_06', 'org_01']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, function_lists_file=function_list_file, annot_names_file=ANNOTATION_NAMES_JSON,
                                     orsum_minterm_size=4, selected_adjust_pvalue_cutoff=1)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_function_selected_cli():
    input_folder = os.path.join('test_data', 'annotation_output')
    output_folder = 'output_folder'
    grouping = 'selected_function'
    function_list_file = os.path.join('test_data', 'function_list.tsv')
    expected_terms = ['org_06', 'org_01']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    subprocess.call(['esmecata_gseapy', 'gseapy_enrichr', '-f', input_folder, '-o', output_folder, '--grouping', grouping, '--function-list', function_list_file,
                     '--annot-names', ANNOTATION_NAMES_JSON, '--orsumMinTermSize', '4'])

    taxon_rank_annotation_enrichment(input_folder, output_folder, grouping, function_lists_file=function_list_file, annot_names_file=ANNOTATION_NAMES_JSON,
                                     orsum_minterm_size=4, selected_adjust_pvalue_cutoff=1)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


def test_taxon_rank_annotation_enrichment_function_selected_picrust2():
    input_annot_file = os.path.join('test_data', 'picrust_test_predicted_func.tsv')
    output_folder = 'output_folder'
    grouping = 'selected_function'
    function_list_file = os.path.join('test_data', 'picrust_function_list.tsv')
    expected_terms = ['org_06', 'org_01']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    taxon_rank_annotation_enrichment(input_annot_file, output_folder, grouping, function_lists_file=function_list_file, annot_names_file=ANNOTATION_NAMES_JSON,
                                     orsum_minterm_size=4, selected_adjust_pvalue_cutoff=0.9)

    result_file = os.path.join(output_folder, 'orsum_output_folder', 'filteredResult-Summary.tsv')
    df = pd.read_csv(result_file, sep='\t')
    representing_terms = df['Representing term id']

    assert sorted(representing_terms) == sorted(expected_terms)
    shutil.rmtree(output_folder)


if __name__ == "__main__":
    test_taxon_rank_annotation_enrichment_selected_cli_input_file_picrust2()
