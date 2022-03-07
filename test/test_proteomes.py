import os
import shutil
import time

from collections import OrderedDict
from ete3 import NCBITaxa

from esmecata.proteomes import taxonomic_affiliation_to_taxon_id, associate_taxon_to_taxon_id, \
                                disambiguate_taxon, find_proteomes_tax_ids, filter_rank_limit, \
                                rest_query_proteomes

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


def test_rest_query_proteomes():
    expected_proteoems = ['UP000000815']
    expected_organism_ids = {'632': ['UP000000815']}
    expected_proteome_data = [['UP000000815', 99.8, 'full', '632', True]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None, beta=None)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert expected_organism_ids[organism] == organism_ids[organism]
    for index, data in enumerate(sorted(expected_proteome_data)):
        assert data == sorted(proteomes_data)[index]


def test_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000000815'])}
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

"""
def test_sparql_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000255169', 'UP000000815'])}
    ncbi = NCBITaxa()
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    json_taxonomic_affiliations = filter_taxon(json_taxonomic_affiliations, ncbi)
    proteomes_description_folder = 'proteomes_description'
    os.mkdir(proteomes_description_folder)
    proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_taxonomic_affiliations=json_taxonomic_affiliations, ncbi=ncbi, proteomes_description_folder=proteomes_description_folder,
                                                                            busco_percentage_keep=90, all_proteomes=None, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql')
    shutil.rmtree(proteomes_description_folder)
    for taxon in expected_proteomes_ids:
        assert expected_proteomes_ids[taxon][0] == proteomes_ids[taxon][0]
        assert set(expected_proteomes_ids[taxon][1]) == set(proteomes_ids[taxon][1])
"""

if __name__ == "__main__":
    test_find_proteomes_tax_ids()
    test_disambiguate_taxon()
    test_find_proteomes_tax_ids()
    #test_sparql_find_proteomes_tax_ids()