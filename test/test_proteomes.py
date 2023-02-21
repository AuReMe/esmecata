import os
import shutil
import time

from collections import OrderedDict, Counter
from ete3 import NCBITaxa

from esmecata.proteomes import taxonomic_affiliation_to_taxon_id, associate_taxon_to_taxon_id, \
                                disambiguate_taxon, find_proteomes_tax_ids, filter_rank_limit, \
                                rest_query_proteomes, sparql_query_proteomes, subsampling_proteomes, \
                                update_taxonomy

TAXONOMIES = {'id_1': 'cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia;species not found'}


def test_taxonomic_affiliation_to_taxon_id():
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_taxonomic_affiliations = OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])

    tax_ids_to_names, taxon_ids = taxonomic_affiliation_to_taxon_id('test', 'cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia;species not found')
    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_ids_to_names[tax_id]
    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == taxon_ids[taxon]


def test_associate_taxon_to_taxon_id():
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])}

    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_id_names[tax_id]
    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_disambiguate_taxon():
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_filter_rank_limit():
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, input_json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    input_json_taxonomic_affiliations = disambiguate_taxon(input_json_taxonomic_affiliations, ncbi)
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    # Keep genus and inferior (here Yersinia).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'genus')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]

    # Keep family and inferior (here Yersiniaceae).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'family')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]

    # Keep order and inferior (here Enterobacterales).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'order')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]),
                                                                ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]

    # Keep class and inferior (here Gammaproteobacteria).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'class')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]),
                                                                ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]

    # Keep infraphylum and inferior (here Gammaproteobacteria, as the phylum is superior to infraphylum).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'infraphylum')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]),
                                                                ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]

    # Keep superkingdom and inferior (here Bacteria).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'superkingdom')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Bacteria', [2]), ('Proteobacteria', [1224]),
                                        ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]),
                                        ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]


def test_organism_ids():
    ncbi = NCBITaxa()
    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', '9', 'Buchnera aphidicola', 80, None, None)

    selected_proteomes = subsampling_proteomes(organism_ids, 10, ncbi)
    expected_organism_ids = {'224915': ['UP000000601'], '107806': ['UP000001806']}
    for org_id in expected_organism_ids:
        assert organism_ids[org_id] == expected_organism_ids[org_id]


def test_subsampling_proteomes():
    ncbi = NCBITaxa()
    organism_ids = {'2562891': ['UP000477739'], '208962': ['UP000292187', 'UP000193938', 'UP000650750', 'UP000407502', 'UP000003042'],
                    '562': ['UP000000625', 'UP000000558', 'UP000464341', 'UP000219757', 'UP000092491',
                            'UP000567387', 'UP000234906', 'UP000016096', 'UP000017268', 'UP000017618'],
                    '564': ['UP000510927', 'UP000392711', 'UP000000745', 'UP000033773']}
    revert_organism_ids ={}
    for org_id, proteomes in organism_ids.items():
        revert_organism_ids.update({proteome_id: org_id for proteome_id in proteomes})
    limit_maximal_number_proteomes = 10

    selected_proteomes = subsampling_proteomes(organism_ids, limit_maximal_number_proteomes, ncbi)
    selected_organisms = [revert_organism_ids[proteome] for proteome in selected_proteomes]

    # 562 has 10 proteomes among a total of 20 proteomes, so 50% of proteomes, subsampling will select 50% out of 10: 5 proteomes.
    # 208962 has 5 proteomes among the 20, so we expect at least 25% of 10 proteomes: ~3 proteomes.
    # 564 has 4 proteoemes among the 20, so we expect at least 20% of 10 proteomes: 2 proteomes.
    # 2562891 has 1 proteoùe with 20 proteomes the subsampling will ensure that there is this proteome among the 10 selected.
    # This leads to 11 proteomes.
    expected_proteomes_representation = {'562': 5, '2562891': 1, '208962': 3, '564': 2}
    for org_id in Counter(selected_organisms):
        assert Counter(selected_organisms)[org_id] == expected_proteomes_representation[org_id]


def test_update_taxonomy():
    outdated_taxonomic_affiliation = 'Firmicutes;Bacilli;Bacillales;Bacilliaceae;Bacillus'
    new_taxonomic_affiliation = update_taxonomy('test', outdated_taxonomic_affiliation)

    expected_taconomic_affiliation = 'root;cellular organisms;Bacteria;Terrabacteria group;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus'

    assert new_taxonomic_affiliation == expected_taconomic_affiliation


def test_rest_query_proteomes():
    expected_proteoems = ['UP000255169', 'UP000000815']
    expected_organism_ids = {'632': ['UP000000815'], '29486': ['UP000255169']}
    expected_proteome_data = [['UP000255169', 98.86363636363636, 'full', '29486', True], ['UP000000815', 99.77272727272727, 'full', '632', True]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_rest_query_proteomes_bioservices():
    expected_proteoems = ['UP000255169', 'UP000000815']
    expected_organism_ids = {'632': ['UP000000815'], '29486': ['UP000255169']}
    expected_proteome_data = [['UP000255169', 98.86363636363636, 'full', '29486', True], ['UP000000815', 99.77272727272727, 'full', '632', True]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None, option_bioservices=True)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_non_reference_proteome_rest():
    expected_proteoems = ['UP000829720']
    expected_organism_ids = {'1534307': ['UP000829720']}
    expected_proteome_data = [['UP000824540', None, 'full', '121402', False], ['UP000829720', 85.71428571428571, 'full', '1534307', False]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 54906, 'Albuliformes', 0.8, all_proteomes=None)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_proteome_rest_all_proteomes():
    expected_proteoems = {'UP000036680', 'UP000267096', 'UP000036681', 'UP000267007', 'UP000050794', 'UP000031036'}
    expected_organism_ids = {'6252': ['UP000036681'], '6265': ['UP000050794', 'UP000031036', 'UP000031036', 'UP000267007'], '6269': ['UP000267096', 'UP000267096', 'UP000036680']}
    expected_proteome_data =  [['UP000031036', 87.70360907058448, 'full', '6265', True], ['UP000267096', 64.58000638773555, 'full', '6269', True],
                            ['UP000036680', 64.58000638773555, 'full', '6269', False], ['UP000036681', 49.66464388374322, 'full', '6252', False],
                            ['UP000050794', 77.96231236026829, 'full', '6265', False], ['UP000267007', 77.8984350047908, 'full', '6265', False]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 33256, 'Ascaridoidea', 0.8, all_proteomes=True)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_non_reference_proteome_sparql():
    expected_proteoems = ['UP000829720']
    expected_organism_ids = {'1534307': ['UP000829720']}
    expected_proteome_data = [['UP000824540', None, 'full', '121402', False], ['UP000829720', 85.71428571428571, 'full', '1534307', False]]

    proteomes, organism_ids, proteomes_data = sparql_query_proteomes('test', 54906, 'Albuliformes', 0.8, all_proteomes=None)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_proteome_sparql_all_proteomes():
    expected_proteoems = {'UP000036680', 'UP000267096', 'UP000036681', 'UP000267007', 'UP000050794', 'UP000031036'}
    expected_organism_ids = {'6252': ['UP000036681'], '6265': ['UP000050794', 'UP000031036', 'UP000031036', 'UP000267007'], '6269': ['UP000267096', 'UP000267096', 'UP000036680']}
    expected_proteome_data = [['UP000036681', 49.66464388374322, 'full', '6252', False], ['UP000050794', 77.96231236026829, 'full', '6265', False],
                            ['UP000031036', 87.70360907058448, 'full', '6265', False], ['UP000031036', 87.70360907058448, 'full', '6265', True],
                            ['UP000267007', 77.8984350047908, 'full', '6265', False], ['UP000267096', 64.58000638773555, 'full', '6269', False],
                            ['UP000267096', 64.58000638773555, 'full', '6269', True], ['UP000036680', 64.58000638773555, 'full', '6269', False]]

    proteomes, organism_ids, proteomes_data = sparql_query_proteomes('test', 33256, 'Ascaridoidea', 0.8, all_proteomes=True)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000000815', 'UP000255169'])}
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    proteomes_description_folder = 'proteomes_description'
    os.mkdir(proteomes_description_folder)
    proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_taxonomic_affiliations=json_taxonomic_affiliations, ncbi=ncbi, proteomes_description_folder=proteomes_description_folder,
                                                                        busco_percentage_keep=90, all_proteomes=None)
    time.sleep(1)
    shutil.rmtree(proteomes_description_folder)
    for taxon in expected_proteomes_ids:
        assert expected_proteomes_ids[taxon][0] == proteomes_ids[taxon][0]
        assert set(expected_proteomes_ids[taxon][1]) == set(proteomes_ids[taxon][1])


def test_sparql_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000000815', 'UP000255169'])}
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    proteomes_description_folder = 'proteomes_description'
    os.mkdir(proteomes_description_folder)
    proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_taxonomic_affiliations=json_taxonomic_affiliations, ncbi=ncbi, proteomes_description_folder=proteomes_description_folder,
                                                                            busco_percentage_keep=90, all_proteomes=None, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql')
    shutil.rmtree(proteomes_description_folder)
    for taxon in expected_proteomes_ids:
        assert expected_proteomes_ids[taxon][0] == proteomes_ids[taxon][0]
        assert set(expected_proteomes_ids[taxon][1]) == set(proteomes_ids[taxon][1])


if __name__ == "__main__":
    #test_find_proteomes_tax_ids()
    #test_disambiguate_taxon()
    #test_subsampling_proteomes()
    #test_organism_ids()
    #test_sparql_find_proteomes_tax_ids()
    #test_rest_query_proteomes()
    #test_find_non_reference_proteome_rest()
    #test_find_proteome_rest_all_proteomes()
    #test_find_proteome_sparql_all_proteomes()
    #test_filter_rank_limit()
    test_update_taxonomy()