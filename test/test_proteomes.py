import csv
import os
import shutil
import subprocess
import time

from collections import OrderedDict, Counter
from ete4 import NCBITaxa

from esmecata.core.proteomes import taxonomic_affiliation_to_taxon_id, associate_taxon_to_taxon_id, \
                                disambiguate_taxon, find_proteomes_tax_ids, filter_rank_limit, \
                                rest_query_proteomes, sparql_query_proteomes, subsampling_proteomes, \
                                update_taxonomy

TAXONOMIES = {'id_1': 'cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia;species not found'}
EXPECTED_PROTEOMES_DATA = [['UP000000815', 94.54545454545455, 'full', '632', True, [['Chromosome', 'Yersinia pestis CO92 complete genome'], ['Plasmid pCD1', 'Yersinia pestis CO92 plasmid pCD1'], ['Plasmid pMT1', 'Yersinia pestis CO92 plasmid pMT1'], ['Plasmid pPCP1', 'Yersinia pestis CO92 plasmid pPCP1']]], ['UP000255169', 98.86363636363636, 'full', '29486', True, [['Unassembled WGS sequence', 'Yersinia ruckeri strain NCTC10476 genome assembly']]], ['UP000000642', 100.0, 'full', '393305', False, [['Chromosome', 'Yersinia enterocolitica subsp. enterocolitica 8081 complete genome'], ['Plasmid pYVe8081', 'Yersinia enterocolitica subsp. enterocolitica 8081 plasmid pYVe8081 complete genome']]], ['UP000001011', 100.0, 'full', '273123', False, [['Chromosome', 'Yersinia pseudotuberculosis IP32953 genome'], ['Plasmid pYV', 'Yersinia pseudotuberculosis IP32953 pYV plasmid'], ['Plasmid pYptb32953', 'Yersinia pseudotuberculosis IP32953 cryptic plasmid']]], ['UP000001019', 99.31818181818181, 'full', '632', False, [['Chromosome', 'Yersinia pestis biovar Microtus str. 91001, complete genome.'], ['Plasmid pCD1', 'Yersinia pestis biovar Microtus str. 91001 plasmid pCD1, complete sequence.'], ['Plasmid pCRY', 'Yersinia pestis biovar Microtus str. 91001 plasmid pCRY, complete sequence.'], ['Plasmid pMT1', 'Yersinia pestis biovar Microtus str. 91001 plasmid pMT1, complete sequence.'], ['Plasmid pPCP1', 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.']]], ['UP000001971', 99.0909090909091, 'full', '360102', False, [['Chromosome', 'Yersinia pestis Antiqua'], ['Plasmid pMT', 'Yersinia pestis Antiqua plasmid pMT'], ['Plasmid pPCP', 'Yersinia pestis Antiqua plasmid pPCP'], ['Plasmid pCD', 'Yersinia pestis Antiqua plasmid pCD']]], ['UP000002412', 99.54545454545455, 'full', '349747', False, [['Chromosome', 'Yersinia pseudotuberculosis IP 31758'], ['Plasmid plasmid_59kb', 'Yersinia pseudotuberculosis IP 31758 plasmid plasmid_59kb'], ['Plasmid plasmid_153kb', 'Yersinia pseudotuberculosis IP 31758 plasmid plasmid_153kb']]], ['UP000002490', 98.18181818181819, 'full', '632', False, [['Chromosome', 'Yersinia pestis KIM10+'], ['Plasmid pMT-1', 'Yersinia pestis KIM10+ plasmid pMT-1']]], ['UP000004430', 94.0909090909091, 'full', '373665', False, [['Plasmid pIP1202', 'Yersinia pestis biovar Orientalis str. IP275 plasmid pIP1202'], ['Unassembled WGS sequence', 'Yersinia pestis biovar Orientalis str. IP275']]], ['UP000008084', 98.18181818181819, 'full', '930944', False, [['Chromosome', 'Yersinia enterocolitica subsp. palearctica Y11'], ['Plasmid pYV03', 'Yersinia enterocolitica subsp. palearctica Y11 plasmid pYVO3 complete sequence']]], ['UP000008936', 99.54545454545455, 'full', '377628', False, [['Chromosome', 'Yersinia pestis Nepal516'], ['Plasmid pMT', 'Yersinia pestis Nepal516 plasmid pMT'], ['Plasmid pPCP', 'Yersinia pestis Nepal516 plasmid pPCP']]], ['UP000038204', 99.0909090909091, 'full', '367190', False, [['Unassembled WGS sequence', 'Yersinia similis']]], ['UP000038750', 99.0909090909091, 'full', '631', False, [['Unassembled WGS sequence', 'Yersinia intermedia']]], ['UP000040088', 98.63636363636363, 'full', '263819', False, [['Unassembled WGS sequence', 'Yersinia aleksiciae']]], ['UP000040841', 98.86363636363636, 'full', '33060', False, [['Unassembled WGS sequence', 'Yersinia mollaretii']]], ['UP000041356', 99.0909090909091, 'full', '630', False, [['Unassembled WGS sequence', 'Yersinia enterocolitica']]], ['UP000041595', 98.63636363636363, 'full', '29483', False, [['Unassembled WGS sequence', 'Yersinia aldovae']]], ['UP000041882', 99.0909090909091, 'full', '2890319', False, [['Unassembled WGS sequence', 'Yersinia thracica']]], ['UP000042054', 98.86363636363636, 'full', '29485', False, [['Unassembled WGS sequence', 'Yersinia rohdei']]], ['UP000043316', 98.63636363636363, 'full', '631', False, [['Unassembled WGS sequence', 'Yersinia intermedia']]], ['UP000045824', 99.0909090909091, 'full', '28152', False, [['Unassembled WGS sequence', 'Yersinia kristensenii']]], ['UP000045840', 98.18181818181819, 'full', '1288385', False, [['Unassembled WGS sequence', 'Yersinia pekkanenii']]], ['UP000046784', 98.86363636363636, 'full', '29484', False, [['Unassembled WGS sequence', 'Yersinia frederiksenii']]], ['UP000048841', 99.0909090909091, 'full', '630', False, [['Unassembled WGS sequence', 'Yersinia enterocolitica']]], ['UP000196440', 99.54545454545455, 'full', '631', False, [['Unassembled WGS sequence', 'Yersinia intermedia']]], ['UP000220513', 100.0, 'full', '28152', False, [['Unassembled WGS sequence', 'Yersinia kristensenii']]], ['UP000229378', 99.54545454545455, 'full', '634', False, [['Unassembled WGS sequence', 'Yersinia bercovieri strain SCPM-O-B-7607, whole genome shotgun sequencing project.']]], ['UP000230961', 99.0909090909091, 'full', '2339259', False, [['Chromosome', 'Yersinia enterocolitica LC20, complete genome.'], ['Plasmid p2_50K', 'Yersinia enterocolitica LC20 plasmid plasmid2_50K, complete sequence.'], ['Plasmid p1_80K', 'Yersinia enterocolitica LC20 plasmid plasmid1_80K, complete sequence.']]], ['UP000251879', 98.4090909090909, 'full', '632', False, [['Unassembled WGS sequence', 'Yersinia pestis']]], ['UP000254835', 99.0909090909091, 'full', '29484', False, [['Unassembled WGS sequence', 'Yersinia frederiksenii']]], ['UP000255087', 97.72727272727273, 'full', '633', False, [['Unassembled WGS sequence', 'Yersinia pseudotuberculosis']]], ['UP000265864', 100.0, 'full', '1604335', False, [['Chromosome', 'Yersinia rochesterensis strain ATCC BAA-2637 chromosome.'], ['Plasmid pUnnamed1', 'Yersinia rochesterensis strain ATCC BAA-2637 plasmid pUnnamed1.'], ['Plasmid pUnnamed2', 'Yersinia rochesterensis strain ATCC BAA-2637 plasmid pUnnamed2.']]], ['UP000464402', 93.86363636363636, 'full', '2607663', False, [['Chromosome', 'Yersinia canariae strain NCTC 14382 chromosome']]], ['UP000595309', 100.0, 'full', '630', False, [['Chromosome', 'Yersinia enterocolitica strain FDAARGOS_1082 chromosome']]], ['UP000698240', 100.0, 'full', '419257', False, [['Unassembled WGS sequence', 'Yersinia massiliensis']]], ['UP000712947', 100.0, 'full', '33060', False, [['Unassembled WGS sequence', 'Yersinia mollaretii']]], ['UP001167864', 99.77272727272727, 'full', '685706', False, [['Unassembled WGS sequence', 'Yersinia nurmii']]], ['UP001182355', 99.77272727272727, 'full', '630', False, [['Unassembled WGS sequence', 'Yersinia enterocolitica']]], ['UP000005382', 95.9090909090909, 'full', '992138', False, [['Unassembled WGS sequence', 'Yersinia pestis PY-12, whole genome shotgun sequencing project.']]]]


def test_taxonomic_affiliation_to_taxon_id_offline():
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_taxonomic_affiliations = OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])

    tax_ids_to_names, taxon_ids = taxonomic_affiliation_to_taxon_id('test', 'cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia;species not found')
    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_ids_to_names[tax_id]
    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == taxon_ids[taxon]


def test_associate_taxon_to_taxon_id_offline():
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    expected_tax_id_names = {2: 'Bacteria', 131567: 'cellular organisms', 91347: 'Enterobacterales', 1236: 'Gammaproteobacteria', 1224: 'Proteobacteria', 629: 'Yersinia', 444888: 'Yersinia', 1903411: 'Yersiniaceae'}
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629, 444888]), ('species not found', ['not_found'])])}

    for tax_id in expected_tax_id_names:
        assert expected_tax_id_names[tax_id] == tax_id_names[tax_id]
    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_disambiguate_taxon_offline():
    ncbi = NCBITaxa()
    update_affiliations = None
    tax_id_names, json_taxonomic_affiliations = associate_taxon_to_taxon_id(TAXONOMIES, update_affiliations, ncbi)
    json_taxonomic_affiliations = disambiguate_taxon(json_taxonomic_affiliations, ncbi)
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('cellular organisms', [131567]), ('Bacteria', [2]), ('Proteobacteria', [1224]), ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]), ('Yersinia', [629]), ('species not found', ['not_found'])])}

    for taxon in expected_json_taxonomic_affiliations:
        assert expected_json_taxonomic_affiliations[taxon] == json_taxonomic_affiliations[taxon]


def test_filter_rank_limit_offline():
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

    # Keep domain and inferior (here Bacteria).
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'domain')
    expected_json_taxonomic_affiliations = {'id_1': OrderedDict([('Bacteria', [2]), ('Proteobacteria', [1224]),
                                        ('Gammaproteobacteria', [1236]), ('Enterobacterales', [91347]), ('Yersiniaceae', [1903411]),
                                        ('Yersinia', [629]), ('species not found', ['not_found'])])}
    for taxon in json_taxonomic_affiliations['id_1']:
        assert expected_json_taxonomic_affiliations['id_1'][taxon] == json_taxonomic_affiliations['id_1'][taxon]


def test_filter_rank_limit_family_synonym():
    # In this case, two taxa have the same tax IDs but different names (synonyms): Actinobacteriota/Actinobacteria.
    # Check that both are deleted.
    input_json_taxonomic_affiliations = {'Cluster_234': OrderedDict({'Bacteria': [2], 'Actinobacteriota': [201174], 'Actinobacteria': [201174], 'Multi-affiliation': ['not_found']})}

    ncbi = NCBITaxa()
    json_taxonomic_affiliations = filter_rank_limit(input_json_taxonomic_affiliations, ncbi, 'family')

    expected_json_taxonomic_affiliations = {'Cluster_234': OrderedDict({'Multi-affiliation': ['not_found']})}
    for taxon in json_taxonomic_affiliations['Cluster_234']:
        assert expected_json_taxonomic_affiliations['Cluster_234'][taxon] == json_taxonomic_affiliations['Cluster_234'][taxon]


def test_organism_ids_online():
    ncbi = NCBITaxa()
    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', '9', 'Buchnera aphidicola', 80, None, None)

    selected_proteomes = subsampling_proteomes(organism_ids, 10, ncbi)
    expected_organism_ids = {'224915': ['UP000000601'], '107806': ['UP000001806']}
    for org_id in expected_organism_ids:
        assert organism_ids[org_id] == expected_organism_ids[org_id]


def test_subsampling_proteomes_offline():
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


def test_update_taxonomy_bacillus_offline():
    outdated_taxonomic_affiliation = 'Firmicutes;Bacilli;Bacillales;Bacilliaceae;Bacillus'
    new_taxonomic_affiliation = update_taxonomy('test', outdated_taxonomic_affiliation)

    expected_taconomic_affiliation = 'root;cellular organisms;Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus'

    assert new_taxonomic_affiliation == expected_taconomic_affiliation


def test_update_taxonomy_yersinia_offline():
    outdated_taxonomic_affiliation = 'Bacteria;Yersinia'
    new_taxonomic_affiliation = update_taxonomy('test', outdated_taxonomic_affiliation)

    expected_taconomic_affiliation = 'root;cellular organisms;Bacteria;Pseudomonadati;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Yersiniaceae;Yersinia'

    assert new_taxonomic_affiliation == expected_taconomic_affiliation


def test_rest_query_proteomes_online():
    expected_proteoems = ['UP000255169', 'UP000000815']
    expected_organism_ids = {'632': ['UP000000815'], '29486': ['UP000255169']}

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in EXPECTED_PROTEOMES_DATA:
        assert data in proteomes_data


def test_rest_query_proteomes_bioservices_online():
    expected_proteoems = ['UP000255169', 'UP000000815']
    expected_organism_ids = {'632': ['UP000000815'], '29486': ['UP000255169']}

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 629, 'Yersinia', 0.8, all_proteomes=None, option_bioservices=True)
    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in EXPECTED_PROTEOMES_DATA:
        assert data in proteomes_data


def test_find_non_reference_proteome_rest_online():
    expected_proteoems = ['UP000829720']
    expected_organism_ids = {'1534307': ['UP000829720']}
    expected_proteome_data = ['UP000824540', 67.91208791208791, 'full', '121402', True, [['Unassembled WGS sequence', 'Albula glossodonta']]], ['UP000829720', 85.71428571428571, 'full', '1534307', True, [['Chromosome 1', 'Albula goreensis ecotype Florida chromosome 1, whole genome shotgun sequence.'], ['Chromosome 2', 'Albula goreensis ecotype Florida chromosome 2, whole genome shotgun sequence.'], ['Chromosome 3', 'Albula goreensis ecotype Florida chromosome 3, whole genome shotgun sequence.'], ['Chromosome 4', 'Albula goreensis ecotype Florida chromosome 4, whole genome shotgun sequence.'], ['Chromosome 5', 'Albula goreensis ecotype Florida chromosome 5, whole genome shotgun sequence.'], ['Chromosome 6', 'Albula goreensis ecotype Florida chromosome 6, whole genome shotgun sequence.'], ['Chromosome 7', 'Albula goreensis ecotype Florida chromosome 7, whole genome shotgun sequence.'], ['Chromosome 8', 'Albula goreensis ecotype Florida chromosome 8, whole genome shotgun sequence.'], ['Chromosome 9', 'Albula goreensis ecotype Florida chromosome 9, whole genome shotgun sequence.'], ['Chromosome 10', 'Albula goreensis ecotype Florida chromosome 10, whole genome shotgun sequence.'], ['Chromosome 11', 'Albula goreensis ecotype Florida chromosome 11, whole genome shotgun sequence.'], ['Chromosome 12', 'Albula goreensis ecotype Florida chromosome 12, whole genome shotgun sequence.'], ['Chromosome 13', 'Albula goreensis ecotype Florida chromosome 13, whole genome shotgun sequence.'], ['Chromosome 14', 'Albula goreensis ecotype Florida chromosome 14, whole genome shotgun sequence.'], ['Chromosome 15', 'Albula goreensis ecotype Florida chromosome 15, whole genome shotgun sequence.'], ['Chromosome 16', 'Albula goreensis ecotype Florida chromosome 16, whole genome shotgun sequence.'], ['Chromosome 17', 'Albula goreensis ecotype Florida chromosome 17, whole genome shotgun sequence.'], ['Chromosome 18', 'Albula goreensis ecotype Florida chromosome 18, whole genome shotgun sequence.'], ['Chromosome 19', 'Albula goreensis ecotype Florida chromosome 19, whole genome shotgun sequence.'], ['Chromosome 20', 'Albula goreensis ecotype Florida chromosome 20, whole genome shotgun sequence.'], ['Chromosome 21', 'Albula goreensis ecotype Florida chromosome 21, whole genome shotgun sequence.'], ['Chromosome 22', 'Albula goreensis ecotype Florida chromosome 22, whole genome shotgun sequence.'], ['Chromosome 23', 'Albula goreensis ecotype Florida chromosome 23, whole genome shotgun sequence.'], ['Chromosome 24', 'Albula goreensis ecotype Florida chromosome 24, whole genome shotgun sequence.'], ['Chromosome 25', 'Albula goreensis ecotype Florida chromosome 25, whole genome shotgun sequence.'], ['Unassembled WGS sequence', 'Albula goreensis']]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 54906, 'Albuliformes', 80, all_proteomes=None)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_proteome_rest_all_proteomes_online():
    expected_proteoems = {'UP000050794', 'UP000036681', 'UP000031036', 'UP000267096', 'UP000887564', 'UP000887569'}
    expected_organism_ids = {'6265': ['UP000031036', 'UP000050794'], '6269': ['UP000267096'], '6256': ['UP000887564'], '6257': ['UP000887569'], '6252': ['UP000036681']}
    expected_proteome_data =  [['UP000031036', 87.70360907058448, 'full', '6265', True, [['Unassembled WGS sequence', 'Toxocara canis']]], ['UP000267096', 64.58000638773555, 'full', '6269', True, [['Unassembled WGS sequence', 'Anisakis simplex genome assembly']]], ['UP000887564', 9.581603321622485, 'full', '6256', True, [['Unplaced', 'Unplaced_6256']]], ['UP000887569', 91.05717023315235, 'full', '6257', True, [['Unplaced', 'Unplaced_6257']]], ['UP000050794', 77.96231236026829, 'full', '6265', False, [['Unassembled WGS sequence', 'Toxocara canis genome assembly']]], ['UP000036681', 49.63270520600447, 'full', '6252', False, [['Genome', 'Genome'], ['Unplaced', 'ASCLU_Unplaced']]]]

    proteomes, organism_ids, proteomes_data = rest_query_proteomes('test', 33256, 'Ascaridoidea', 0.8, all_proteomes=True)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        assert data in proteomes_data


def test_find_non_reference_proteome_sparql_online():
    expected_proteoems = ['UP000829720']
    expected_organism_ids = {'1534307': ['UP000829720']}
    expected_proteome_data = [['UP000829720', 85.71428571428571, 'full', '1534307', True, [['Chromosome 20', 'Albula goreensis ecotype Florida chromosome 20, whole genome shotgun sequence.'], ['Chromosome 20', 'Albula goreensis ecotype Florida chromosome 20, whole genome shotgun sequence.'], ['Chromosome 14', 'Albula goreensis ecotype Florida chromosome 14, whole genome shotgun sequence.'], ['Chromosome 14', 'Albula goreensis ecotype Florida chromosome 14, whole genome shotgun sequence.'], ['Chromosome 9', 'Albula goreensis ecotype Florida chromosome 9, whole genome shotgun sequence.'], ['Chromosome 9', 'Albula goreensis ecotype Florida chromosome 9, whole genome shotgun sequence.'], ['Chromosome 11', 'Albula goreensis ecotype Florida chromosome 11, whole genome shotgun sequence.'], ['Chromosome 11', 'Albula goreensis ecotype Florida chromosome 11, whole genome shotgun sequence.'], ['Chromosome 24', 'Albula goreensis ecotype Florida chromosome 24, whole genome shotgun sequence.'], ['Chromosome 24', 'Albula goreensis ecotype Florida chromosome 24, whole genome shotgun sequence.'], ['Chromosome 13', 'Albula goreensis ecotype Florida chromosome 13, whole genome shotgun sequence.'], ['Chromosome 13', 'Albula goreensis ecotype Florida chromosome 13, whole genome shotgun sequence.'], ['Chromosome 25', 'Albula goreensis ecotype Florida chromosome 25, whole genome shotgun sequence.'], ['Chromosome 25', 'Albula goreensis ecotype Florida chromosome 25, whole genome shotgun sequence.'], ['Chromosome 1', 'Albula goreensis ecotype Florida chromosome 1, whole genome shotgun sequence.'], ['Chromosome 1', 'Albula goreensis ecotype Florida chromosome 1, whole genome shotgun sequence.'], ['Chromosome 16', 'Albula goreensis ecotype Florida chromosome 16, whole genome shotgun sequence.'], ['Chromosome 16', 'Albula goreensis ecotype Florida chromosome 16, whole genome shotgun sequence.'], ['Chromosome 17', 'Albula goreensis ecotype Florida chromosome 17, whole genome shotgun sequence.'], ['Chromosome 17', 'Albula goreensis ecotype Florida chromosome 17, whole genome shotgun sequence.'], ['Chromosome 10', 'Albula goreensis ecotype Florida chromosome 10, whole genome shotgun sequence.'], ['Chromosome 10', 'Albula goreensis ecotype Florida chromosome 10, whole genome shotgun sequence.'], ['Chromosome 19', 'Albula goreensis ecotype Florida chromosome 19, whole genome shotgun sequence.'], ['Chromosome 19', 'Albula goreensis ecotype Florida chromosome 19, whole genome shotgun sequence.'], ['Chromosome 12', 'Albula goreensis ecotype Florida chromosome 12, whole genome shotgun sequence.'], ['Chromosome 12', 'Albula goreensis ecotype Florida chromosome 12, whole genome shotgun sequence.'], ['Chromosome 15', 'Albula goreensis ecotype Florida chromosome 15, whole genome shotgun sequence.'], ['Chromosome 15', 'Albula goreensis ecotype Florida chromosome 15, whole genome shotgun sequence.'], ['Chromosome 23', 'Albula goreensis ecotype Florida chromosome 23, whole genome shotgun sequence.'], ['Chromosome 23', 'Albula goreensis ecotype Florida chromosome 23, whole genome shotgun sequence.'], ['Unassembled WGS sequence', 'Albula goreensis'], ['Unassembled WGS sequence', 'Albula goreensis'], ['Chromosome 18', 'Albula goreensis ecotype Florida chromosome 18, whole genome shotgun sequence.'], ['Chromosome 18', 'Albula goreensis ecotype Florida chromosome 18, whole genome shotgun sequence.'], ['Chromosome 3', 'Albula goreensis ecotype Florida chromosome 3, whole genome shotgun sequence.'], ['Chromosome 3', 'Albula goreensis ecotype Florida chromosome 3, whole genome shotgun sequence.'], ['Chromosome 22', 'Albula goreensis ecotype Florida chromosome 22, whole genome shotgun sequence.'], ['Chromosome 22', 'Albula goreensis ecotype Florida chromosome 22, whole genome shotgun sequence.'], ['Chromosome 6', 'Albula goreensis ecotype Florida chromosome 6, whole genome shotgun sequence.'], ['Chromosome 6', 'Albula goreensis ecotype Florida chromosome 6, whole genome shotgun sequence.'], ['Chromosome 4', 'Albula goreensis ecotype Florida chromosome 4, whole genome shotgun sequence.'], ['Chromosome 4', 'Albula goreensis ecotype Florida chromosome 4, whole genome shotgun sequence.'], ['Chromosome 5', 'Albula goreensis ecotype Florida chromosome 5, whole genome shotgun sequence.'], ['Chromosome 5', 'Albula goreensis ecotype Florida chromosome 5, whole genome shotgun sequence.'], ['Chromosome 21', 'Albula goreensis ecotype Florida chromosome 21, whole genome shotgun sequence.'], ['Chromosome 21', 'Albula goreensis ecotype Florida chromosome 21, whole genome shotgun sequence.'], ['Chromosome 2', 'Albula goreensis ecotype Florida chromosome 2, whole genome shotgun sequence.'], ['Chromosome 2', 'Albula goreensis ecotype Florida chromosome 2, whole genome shotgun sequence.'], ['Chromosome 8', 'Albula goreensis ecotype Florida chromosome 8, whole genome shotgun sequence.'], ['Chromosome 8', 'Albula goreensis ecotype Florida chromosome 8, whole genome shotgun sequence.'], ['Chromosome 7', 'Albula goreensis ecotype Florida chromosome 7, whole genome shotgun sequence.'], ['Chromosome 7', 'Albula goreensis ecotype Florida chromosome 7, whole genome shotgun sequence.']]], ['UP000824540', 67.91208791208791, 'full', '121402', True, [['Unassembled WGS sequence', 'Albula glossodonta'], ['Unassembled WGS sequence', 'Albula glossodonta']]]]

    proteomes, organism_ids, proteomes_data = sparql_query_proteomes('test', 54906, 'Albuliformes', 80, all_proteomes=None)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in proteomes_data:
        if data[0] == 'UP000829720':
            assert data[0:5] == expected_proteome_data[0][0:5]
            assert sorted(data[5]) == sorted(expected_proteome_data[0][5])


def test_find_proteome_sparql_all_proteomes_online():
    expected_proteoems = {'UP000031036', 'UP000887569'}
    expected_organism_ids = {'6265': ['UP000031036']}
    expected_proteome_data = [['UP000036680', 64.58000638773555, 'full', '6269', False, [['Genome', 'Genome'], ['Unplaced', 'ANISI_Unplaced']]], ['UP000267096', 64.58000638773555, 'full', '6269', True, [['Unassembled WGS sequence', 'Anisakis simplex'], ['Unassembled WGS sequence', 'Anisakis simplex']]], ['UP000887569', 91.05717023315235, 'full', '6257', True, [['Unplaced', 'Unplaced_6257'], ['Unplaced', 'Unplaced_6257']]], ['UP000031036', 87.70360907058448, 'full', '6265', True, [['Unassembled WGS sequence', 'Toxocara canis'], ['Unassembled WGS sequence', 'Toxocara canis']]], ['UP000050794', 77.96231236026829, 'full', '6265', False, [['Genome assembly', 'Genome assembly'], ['Unplaced', 'TOXCA_Unplaced']]], ['UP000267007', 77.8984350047908, 'full', '6265', False, [['Unassembled WGS sequence', 'Toxocara canis']]], ['UP000887564', 9.581603321622485, 'full', '6256', True, [['Unplaced', 'Unplaced_6256'], ['Unplaced', 'Unplaced_6256']]], ['UP000036681', 49.63270520600447, 'full', '6252', False, [['Genome', 'Genome'], ['Unplaced', 'ASCLU_Unplaced']]]]

    proteomes, organism_ids, proteomes_data = sparql_query_proteomes('test', 33256, 'Ascaridoidea', 80, all_proteomes=True)

    time.sleep(1)

    assert set(expected_proteoems) == set(proteomes)
    for organism in expected_organism_ids:
        assert set(expected_organism_ids[organism]).issubset(set(organism_ids[organism]))
    for data in expected_proteome_data:
        if data[0] == 'UP000036680':
            assert data[0:5] == expected_proteome_data[0][0:5]
            assert sorted(data[5]) == sorted(expected_proteome_data[0][5])


def test_find_proteomes_tax_ids_online():
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


def test_sparql_find_proteomes_tax_ids_online():
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


def test_check_cli_online():
    output_folder = 'proteomes_output'
    subprocess.call(['esmecata', 'check', '-i', 'buchnera_workflow.tsv', '-o', output_folder])
    expected_results = []
    output_stat_file = os.path.join(output_folder, 'proteome_tax_id.tsv')
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.DictReader(stat_file_read, delimiter='\t')
        for line in csvreader:
            expected_results.append(line['tax_rank'])

    assert expected_results == ['species']

if __name__ == "__main__":
    #test_find_proteomes_tax_ids_online()
    #test_disambiguate_taxon_offline()
    #test_subsampling_proteomes_offline()
    #test_organism_ids_online()
    #test_sparql_find_proteomes_tax_ids_online()
    #test_rest_query_proteomes_online()
    #test_find_non_reference_proteome_rest_online()
    #test_find_proteome_rest_all_proteomes_online()
    #test_find_proteome_sparql_all_proteomes_online()
    #test_filter_rank_limit_offline()
    test_update_taxonomy_bacillus_offline()
    test_update_taxonomy_yersinia_offline()