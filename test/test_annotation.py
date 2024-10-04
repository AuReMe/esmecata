import os

from esmecata.core.annotation import extract_protein_cluster, search_already_annotated_protein, query_uniprot_annotation_rest, \
                                query_uniprot_annotation_sparql, propagate_annotation_in_cluster, extract_protein_annotation_from_files

from Bio import SeqIO

ANOTATIONS = {'Q7CGB6': ['Protein translocase subunit SecA', 'UniProtKB reviewed (Swiss-Prot)',
            ['GO:0008564', 'GO:0006605', 'GO:0005886', 'GO:0031522', 'GO:0043952', 'GO:0046872', 'GO:0005737', 'GO:0065002', 'GO:0017038', 'GO:0005524'],
            ['7.4.2.8'],
            ['IPR027417','IPR004027','IPR000185','IPR020937','IPR011115','IPR014018','IPR011130','IPR011116','IPR036266','IPR036670','IPR014001','IPR044722'],
            [],
            'secA']}

UP000119554_ANOTATIONS = {'O91464': ['Genome polyprotein', 'UniProtKB reviewed (Swiss-Prot)',
                            ['GO:0039618', 'GO:0044162', 'GO:0006508', 'GO:0003723', 'GO:0039694', 'GO:0016020', 'GO:0019062',
                             'GO:0003968', 'GO:0016887', 'GO:0005524', 'GO:0044178', 'GO:0046718', 'GO:0003724', 'GO:0005198',
                             'GO:0039522', 'GO:0006351', 'GO:0004197', 'GO:0034220', 'GO:0015267'],
                            ['2.7.7.48', '3.6.4.13'],
                            ['IPR000199', 'IPR000605', 'IPR001205', 'IPR001676', 'IPR004004',
                            'IPR007053', 'IPR007094', 'IPR009003', 'IPR014759', 'IPR027417',
                            'IPR029053', 'IPR033703', 'IPR043128', 'IPR043502', 'IPR043504'],
                            ['13065', '21248'],
                            '']}

TREMBL_ANNOTATIONS = {'A0A1B2H8S9': ['Siroheme synthase', False,
        ['GO:0009236', 'GO:0043115', 'GO:0032259', 'GO:0019354', 'GO:0004851', 'GO:0051266', 'GO:0051287'],
        ['2.1.1.107', '4.99.1.4', '1.3.1.76'],
        ['IPR019478', 'IPR000878', 'IPR037115', 'IPR014776', 'IPR036291', 'IPR012409', 'IPR006367', 'IPR006366', 'IPR014777', 'IPR035996', 'IPR003043'],
        ['RHEA:32459', 'RHEA:15613', 'RHEA:24360'], 'cysG'],
        'A0A0M4HE72': ['Dual-specificity RNA methyltransferase RlmN', False,
        ['GO:0070040', 'GO:0070475', 'GO:0005737', 'GO:0000049', 'GO:0002935', 'GO:0051539', 'GO:0046872', 'GO:0019843'],
        ['2.1.1.192'], ['IPR004383', 'IPR027492', 'IPR007197', 'IPR013785', 'IPR040072'], ['RHEA:43332', 'RHEA:42916'], 'rlmN'],
        'A0A5A7R956': ['', False, ['GO:0003682'], [], ['IPR039276', 'IPR032001', 'IPR009057'], [], '']}

SWISSPROT_ANNOTATIONS = {'P57136': ['Ribosomal RNA small subunit methyltransferase D', True,
    ['GO:0052913', 'GO:0003676'], ['2.1.1.171'], ['IPR004398', 'IPR029063', 'IPR002052'], ['RHEA:23548'], 'rsmD'],
    'P57634': ['Mannitol-1-phosphate 5-dehydrogenase', True, ['GO:0019594', 'GO:0008926'], ['1.1.1.17'],
    ['IPR036291', 'IPR013328', 'IPR008927', 'IPR000669', 'IPR023028', 'IPR013118', 'IPR013131'], ['RHEA:19661'], 'mtlD'],
    'P57406': ['Protease HtpX', True, ['GO:0008270', 'GO:0004222', 'GO:0006508', 'GO:0005886'],
    ['3.4.24.-'], ['IPR022919', 'IPR001915'], [], 'htpX']
    }


def compare_annotation_dict(expected_dict, result_dict, propagate_test=None):
    for protein in result_dict:
        assert expected_dict[protein][0] == result_dict[protein][0]
        # For propagation 1 is GO term where for other dict it is reviewed
        if propagate_test is None:
            assert expected_dict[protein][1] == result_dict[protein][1]
        else:
            assert set(expected_dict[protein][1]) == set(result_dict[protein][1])
        assert set(expected_dict[protein][2]) == set(result_dict[protein][2])
        assert set(expected_dict[protein][3]) == set(result_dict[protein][3])
        # No more annotations for protein_annotation dict from propagation
        if propagate_test is None:
            assert set(expected_dict[protein][4]) == set(result_dict[protein][4])
            assert set(expected_dict[protein][5]) == set(result_dict[protein][5])
            assert expected_dict[protein][6] == result_dict[protein][6]


def test_extract_protein_cluster_offline():
    reference_proteins, set_proteins = extract_protein_cluster(os.path.join('annotation_input', 'reference_proteins', 'Buchnera_c32_aphidicola__taxid__9.tsv'))
    assert len(reference_proteins) == 573
    assert len(set_proteins) == 9854


def test_search_already_annotated_protein_offline():
    output_dict = {}
    already_annotated_proteins = {'P59488': ['Shikimate kinase', True, ['GO:0000287', 'GO:0004765', 'GO:0005524', 'GO:0005737', 'GO:0008652', 'GO:0009073', 'GO:0009423'], ['2.7.1.71'],
                                    ['IPR000623', 'IPR023000', 'IPR027417', 'IPR031322'], ['RHEA:13121'], ['aroK']]}
    reference_proteins, set_proteins = extract_protein_cluster(os.path.join('annotation_input', 'reference_proteins', 'Buchnera_c32_aphidicola__taxid__9.tsv')) 
    protein_to_search_on_uniprots, output_dict = search_already_annotated_protein(set_proteins, already_annotated_proteins, output_dict)

    assert len(set_proteins)-1 == len(protein_to_search_on_uniprots)

def test_query_uniprot_annotation_rest_online():
    protein_to_search_on_uniprots = ['Q7CGB6', 'NOTAGOODIDEAOFPROTEIN']
    output_dict = {}
    output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict)
    compare_annotation_dict(ANOTATIONS, output_dict)


def test_query_uniprot_annotation_rest_mocked(mocker):
    mock_data = ANOTATIONS

    # Create a mock response object with a .json() method that returns the mock data
    mock_response = mocker.MagicMock()
    mock_response.json.return_value = mock_data

    # Patch 'requests.get' to return the mock response
    mocker.patch("requests.get", return_value=mock_response)
    protein_to_search_on_uniprots = ['Q7CGB6', 'NOTAGOODIDEAOFPROTEIN']
    output_dict = {}
    output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict)
    compare_annotation_dict(ANOTATIONS, output_dict)


def test_query_uniprot_annotation_rest_bioservices_online():
    protein_to_search_on_uniprots = ['Q7CGB6', 'NOTAGOODIDEAOFPROTEIN']
    output_dict = {}
    output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict, option_bioservices=True)
    compare_annotation_dict(ANOTATIONS, output_dict)


def test_query_uniprot_annotation_rest_bioservices_mocked(mocker):
    mock_data = ANOTATIONS

    # Create a mock response object with a .json() method that returns the mock data
    mock_response = mocker.MagicMock()
    mock_response.json.return_value = mock_data

    # Patch 'requests.get' to return the mock response
    mocker.patch("requests.get", return_value=mock_response)
    protein_to_search_on_uniprots = ['Q7CGB6', 'NOTAGOODIDEAOFPROTEIN']
    output_dict = {}
    output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict, option_bioservices=True)
    compare_annotation_dict(ANOTATIONS, output_dict)


def test_query_uniprot_annotation_sparql_online():
    proteomes = ['UP000119554']
    uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
    output_dict = {}
    output_dict = query_uniprot_annotation_sparql(proteomes, uniprot_sparql_endpoint, output_dict)
    compare_annotation_dict(UP000119554_ANOTATIONS, output_dict)


def test_propagate_annotation_in_cluster_offline():
    output_dict = {'prot_1': ['function_1', True, ['GO:0031522', 'GO:0004765'], ['7.4.2.8'], ['IPR027417'], [], 'gene_1'],
                    'prot_2': ['function_2', True, ['GO:0031522', 'GO:0005737'], ['7.4.2.8', '2.7.1.71'], ['IPR027417'], [], 'gene_1'],
                    'prot_3': ['function_1', True, ['GO:0031522', 'GO:0005737'], ['7.4.2.8'], [], [], 'gene_3']}
    reference_proteins = {'prot_1': ['prot_1', 'prot_2', 'prot_3']}
    uniref_output_dict = None

    # Default behaviour no propagation only annotation from representative proteins.
    propagate_annotation = None
    protein_annotations = propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict)
    expected_result_no_propagation = {'prot_1': ['function_1', ['GO:0031522', 'GO:0004765'], ['7.4.2.8'], 'gene_1']}
    compare_annotation_dict(expected_result_no_propagation, protein_annotations, True)

    # Propagation with threshold at 0, meaning all the annotations from all the protein in the cluster will be used (union of annotations).
    propagate_annotation = 0
    protein_annotations = propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict)
    expected_result_propagation_0_union = {'prot_1': ['function_1', ['GO:0004765', 'GO:0005737', 'GO:0031522'], ['2.7.1.71', '7.4.2.8'], 'gene_1']}
    compare_annotation_dict(expected_result_propagation_0_union, protein_annotations, True)

    # Propagation with threshold at 0.66, an annotation is kept if it appears at least in 2 proteins on the 3 of the test.
    propagate_annotation = 0.66
    protein_annotations = propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict)
    expected_result_propagation_0_66 = {'prot_1': ['function_1', ['GO:0005737', 'GO:0031522'], ['7.4.2.8'], 'gene_1']}
    compare_annotation_dict(expected_result_propagation_0_66, protein_annotations, True)

    # Propagation with threshold at 1, an annotation is kept only if it occurs in all protein of the clsuter (intersection of annotation).
    propagate_annotation = 1
    protein_annotations = propagate_annotation_in_cluster(output_dict, reference_proteins, propagate_annotation, uniref_output_dict)
    expected_result_propagation_1_intersection = {'prot_1': ['function_1', ['GO:0031522'], ['7.4.2.8'], 'gene_1']}
    compare_annotation_dict(expected_result_propagation_1_intersection, protein_annotations, True)


def test_annotation_from_files_offline():
    uniprot_trembl_index = SeqIO.index('uniprot_trembl.txt', 'swiss')

    uniprot_sprot_index = SeqIO.index('uniprot_sprot.txt', 'swiss')

    protein_to_search_on_uniprots = ['A0A5A7R956', 'A0A0M4HE72', 'A0A1B2H8S9', 'P57406', 'P57634', 'P57136']
    output_dict = {}
    output_dict = extract_protein_annotation_from_files(protein_to_search_on_uniprots, uniprot_trembl_index, uniprot_sprot_index, output_dict)

    expected_dict = {}
    expected_dict.update(TREMBL_ANNOTATIONS)
    expected_dict.update(SWISSPROT_ANNOTATIONS)

    compare_annotation_dict(expected_dict, output_dict)


if __name__ == "__main__":
    test_extract_protein_cluster_offline()
    test_search_already_annotated_protein_offline()
    test_query_uniprot_annotation_rest_online()
    test_query_uniprot_annotation_rest_bioservices_online()
    test_query_uniprot_annotation_sparql_online()
    test_propagate_annotation_in_cluster_offline()
    test_annotation_from_files_offline()
