import os
import shutil
import time

from collections import OrderedDict, Counter
from ete3 import NCBITaxa

from esmecata.proteomes import taxonomic_affiliation_to_taxon_id, associate_taxon_to_taxon_id, \
                                disambiguate_taxon, find_proteomes_tax_ids, filter_rank_limit, \
                                rest_query_proteomes, subsampling_proteomes

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
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])}

    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_id_names[tax_id]
    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_disambiguate_taxon():
    ncbi = NCBITaxa()
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_filter_rank_limit():
    ncbi = NCBITaxa()
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    json_taxonomic_affiliations = filter_rank_limit(json_taxonomic_affiliations, ncbi, 'superkingdom')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    for taxon in expected_json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]


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


def test_rest_query_proteomes():
    expected_proteoems = ['UP000255169', 'UP000000815']
    expected_organism_ids = {632: ['UP000000815'], 29486: ['UP000255169']}
    expected_proteome_data = [['UP000255169', 98.86363636363636, 'full', 29486, True], ['UP000000815', 99.77272727272727, 'full', 632, True]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None)
    time.sleep(1)
    print(organism_ids)
    assert set(expected_proteoems) == set(proteomes)
    for organism in organism_ids:
        assert set(organism_ids[organism]) == set(expected_organism_ids[organism])
    for index, data in enumerate(sorted(proteomes_data)):
        assert sorted(expected_proteome_data)[index] == data


def test_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000000815', 'UP000255169'])}
    ncbi = NCBITaxa()
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
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
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
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
    test_find_proteomes_tax_ids()
    test_disambiguate_taxon()
    test_find_proteomes_tax_ids()
    test_subsampling_proteomes()
    #test_sparql_find_proteomes_tax_ids()