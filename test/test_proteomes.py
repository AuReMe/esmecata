from collections import OrderedDict
from ete3 import NCBITaxa

from esmecata.retrieve_proteome import associate_taxon_to_taxon_id, filter_taxon, find_proteomes_tax_ids

TAXONOMIES = {'id_1': 'cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia;species not found'}

def test_associate_taxon_to_taxon_id():
    ncbi = NCBITaxa()
    tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_cluster_taxons = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])}

    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_id_names[tax_id]
    for taxon in expected_json_cluster_taxons:
        assert expected_json_cluster_taxons[taxon] == json_cluster_taxons[taxon]

def test_filter_taxon():
    ncbi = NCBITaxa()
    tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    json_cluster_taxons = filter_taxon(json_cluster_taxons, ncbi)
    expected_json_cluster_taxons = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    for taxon in expected_json_cluster_taxons:
        assert expected_json_cluster_taxons[taxon] == json_cluster_taxons[taxon]

def test_find_proteomes_tax_ids():
    expected_proteomes_ids = {'id_1': (629, ['UP000255169', 'UP000000815'])}
    ncbi = NCBITaxa()
    tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(TAXONOMIES, ncbi)
    json_cluster_taxons = filter_taxon(json_cluster_taxons, ncbi)
    proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_cluster_taxons, 90)
    for taxon in expected_proteomes_ids:
        assert expected_proteomes_ids[taxon] == proteomes_ids[taxon]


if __name__ == "__main__":
    test_find_proteomes_tax_ids()
    test_filter_taxon()
    test_find_proteomes_tax_ids()