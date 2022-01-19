from esmecata.annotation import rest_query_uniprot_to_retrieve_function, sparql_query_uniprot_to_retrieve_function

ANOTATIONS = {'Q7CGB6': ['Protein translocase subunit SecA (EC 7.4.2.8)', True,
            ['GO:0031522','GO:0005829','GO:0005887','GO:0015462','GO:0005524','GO:0046872','GO:0065002','GO:0017038','GO:0006605','GO:0043952','GO:0008564'],
            ['7.4.2.8'],
            ['IPR027417','IPR004027','IPR000185','IPR020937','IPR011115','IPR014018','IPR011130','IPR011116','IPR036266','IPR036670','IPR014001','IPR044722'],
            [],
            'secA']}

def test_rest_query_uniprot_to_retrieve_function():
    output_dict = {}
    output_dict = rest_query_uniprot_to_retrieve_function('Q7CGB6')

    for protein in output_dict:
        assert ANOTATIONS[protein][0] == output_dict[protein][0]
        assert ANOTATIONS[protein][1] == output_dict[protein][1]
        assert set(ANOTATIONS[protein][2]) == set(output_dict[protein][2])
        assert set(ANOTATIONS[protein][3]) == set(output_dict[protein][3])
        assert set(ANOTATIONS[protein][4]) == set(output_dict[protein][4])
        assert set(ANOTATIONS[protein][5]) == set(output_dict[protein][5])
        assert ANOTATIONS[protein][6] == output_dict[protein][6]


def test_sparql_query_uniprot_to_retrieve_function():
    output_dict = {}
    output_dict = sparql_query_uniprot_to_retrieve_function('Q7CGB6', 'https://sparql.uniprot.org/sparql')

    for protein in output_dict:
        assert ANOTATIONS[protein][0] == output_dict[protein][0]
        assert ANOTATIONS[protein][1] == output_dict[protein][1]
        assert set(ANOTATIONS[protein][2]) == set(output_dict[protein][2])
        assert set(ANOTATIONS[protein][3]) == set(output_dict[protein][3])
        assert set(ANOTATIONS[protein][4]) == set(output_dict[protein][4])
        assert set(ANOTATIONS[protein][5]) == set(output_dict[protein][5])
        assert ANOTATIONS[protein][6] == output_dict[protein][6]


if __name__ == "__main__":
    test_rest_query_uniprot_to_retrieve_function()
    test_sparql_query_uniprot_to_retrieve_function()
