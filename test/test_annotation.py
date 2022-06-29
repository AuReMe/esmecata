from esmecata.annotation import rest_query_uniprot_to_retrieve_function, sparql_query_uniprot_to_retrieve_function, extract_protein_cluster, \
                                search_already_annotated_protein, query_uniprot_annotation_rest, query_uniprot_annotation_sparql, \
                                propagate_annotation_in_cluster

ANOTATIONS = {'Q7CGB6': ['Protein translocase subunit SecA', True,
            ['GO:0031522','GO:0005829','GO:0005887','GO:0015462','GO:0005524','GO:0046872','GO:0065002','GO:0017038','GO:0006605','GO:0043952','GO:0008564'],
            ['7.4.2.8'],
            ['IPR027417','IPR004027','IPR000185','IPR020937','IPR011115','IPR014018','IPR011130','IPR011116','IPR036266','IPR036670','IPR014001','IPR044722'],
            [],
            'secA']}

UP000119554_ANOTATIONS = {'O91464': ['Genome polyprotein', True,
                            ['GO:0003723', 'GO:0003724', 'GO:0003968', 'GO:0004197', 'GO:0005198',
                                'GO:0005216', 'GO:0005524', 'GO:0006351', 'GO:0016021', 'GO:0016887',
                                'GO:0018144', 'GO:0019028', 'GO:0019062', 'GO:0039522', 'GO:0039618',
                                'GO:0039657', 'GO:0039690', 'GO:0039694', 'GO:0039707', 'GO:0044162',
                                'GO:0044178', 'GO:0044385', 'GO:0046718', 'GO:0051259'],
                            ['2.7.7.48', '3.6.4.13'],
                            ['IPR000199', 'IPR000605', 'IPR001205', 'IPR001676', 'IPR004004',
                            'IPR007053', 'IPR007094', 'IPR009003', 'IPR014759', 'IPR027417',
                            'IPR029053', 'IPR033703', 'IPR043128', 'IPR043502', 'IPR043504'],
                            ['13065', '21248'],
                            '']}

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


def test_extract_protein_cluster():
    reference_proteins, set_proteins = extract_protein_cluster('annotation_input/reference_proteins/Cluster_1.tsv')
    assert len(reference_proteins) == 460
    assert len(set_proteins) == 936


def test_search_already_annotated_protein():
    output_dict = {}
    already_annotated_proteins = {'P59488': ['Shikimate kinase', True, ['GO:0000287', 'GO:0004765', 'GO:0005524', 'GO:0005737', 'GO:0008652', 'GO:0009073', 'GO:0009423'], ['2.7.1.71'],
                                    ['IPR000623', 'IPR023000', 'IPR027417', 'IPR031322'], ['RHEA:13121'], ['aroK']]}
    reference_proteins, set_proteins = extract_protein_cluster('annotation_input/reference_proteins/Cluster_1.tsv') 
    protein_to_search_on_uniprots, output_dict = search_already_annotated_protein(set_proteins, already_annotated_proteins, output_dict)

    assert len(set_proteins)-1 == len(protein_to_search_on_uniprots)

def test_query_uniprot_annotation_rest():
    protein_to_search_on_uniprots = ['Q7CGB6']
    output_dict = {}
    output_dict = query_uniprot_annotation_rest(protein_to_search_on_uniprots, output_dict)
    compare_annotation_dict(ANOTATIONS, output_dict)


def test_query_uniprot_annotation_sparql():
    proteomes = ['UP000119554']
    uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
    output_dict = {}
    output_dict = query_uniprot_annotation_sparql(proteomes, uniprot_sparql_endpoint, output_dict)
    compare_annotation_dict(UP000119554_ANOTATIONS, output_dict)


def test_propagate_annotation_in_cluster():
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


if __name__ == "__main__":
    test_extract_protein_cluster()
    test_search_already_annotated_protein()
    test_propagate_annotation_in_cluster()
    test_query_uniprot_annotation_rest()
    test_query_uniprot_annotation_sparql()
