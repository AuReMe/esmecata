import csv
import json
import os.path
import sys
from ete3 import NCBITaxa


ncbi = NCBITaxa()


def create_krona_input_in(association_taxon_id_json, name, output):
    krona_input = os.path.join(output, f'{name}_krona_input.tsv')
    with open(krona_input, 'w') as f:
        f.write('#queryID\t#taxID')

    with open(association_taxon_id_json, 'r') as f:
        input_taxon_data = json.load(f)

    with open(krona_input, 'a') as f:
        for observation_name in input_taxon_data:
            reversed_affiliation_taxa = list(reversed(input_taxon_data[observation_name]))
            for tax_name in reversed_affiliation_taxa:
                if tax_name != 'unknown' and input_taxon_data[observation_name][tax_name] != ['not_found']:
                    tax_id = input_taxon_data[observation_name][tax_name][0]
                    f.write(f'\n{observation_name}\t{tax_id}')
                    break

    krona_output = os.path.join(output, f'{name}_krona_input.html')
    os.system(f'ktImportTaxonomy {krona_input} -o {krona_output}')


def create_krona_input_esmecata(proteome_tax_id_tsv, name, output):
    krona_input = os.path.join(output, f'{name}_krona_esmecata.tsv')
    with open(krona_input, 'w') as o_f:
        o_f.write('#queryID\t#taxID')

        with open(proteome_tax_id_tsv, 'r') as i_f:
            lines = list(csv.reader(i_f, delimiter='\t'))
            for l in lines[1:]:
                o_f.write(f'\n{l[0]}\t{l[2]}')

    krona_output = os.path.join(output, f'{name}_krona_esmecata.html')
    os.system(f'ktImportTaxonomy {krona_input} -o {krona_output}')


def create_comp_taxonomy_file(association_taxon_id_json, proteome_tax_id_tsv, name, output):
    with open(association_taxon_id_json, 'r') as f:
        input_taxon_data = json.load(f)

    d_tax = dict()
    with open(proteome_tax_id_tsv, 'r') as f:
        lines = list(csv.reader(f, delimiter='\t'))
        for l in lines[1:]:
            observation_name = l[0]
            reversed_affiliation_taxa = list(reversed(input_taxon_data[observation_name]))
            for tax_name in reversed_affiliation_taxa:
                if tax_name != 'unknown' and input_taxon_data[observation_name][tax_name] != ['not_found']:
                    if tax_name not in d_tax.keys():
                        tax_id = input_taxon_data[observation_name][tax_name][0]
                        tax_rank = ncbi.get_rank([tax_id])
                        if tax_id in tax_rank.keys():
                            tax_rank = tax_rank[tax_id]
                        else:
                            tax_rank = ''
                        d_tax[tax_name] = {'Esmecata ID': l[2], 'Esmecata Rank': l[3], 'Esmecata Name': l[1],
                                           'Input ID': tax_id, 'Input Rank': tax_rank}
                    if 'obs' not in d_tax[tax_name].keys():
                        d_tax[tax_name]['obs'] = list()
                    d_tax[tax_name]['obs'].append(observation_name)
                    break

    output = os.path.join(output, f'{name}_taxonomy_diff.tsv')
    with open(output, 'w') as f:
        f.write('\t'.join(['Input Name', 'Esmecata Name', 'Input Rank', 'Esmecata Rank',
                           'Input ID', 'Esmecata ID', 'Observation Names']))
        for name, ref in d_tax.items():
            f.write('\n' + '\t'.join([name, ref['Esmecata Name'], ref['Input Rank'], ref['Esmecata Rank'],
                                      str(ref['Input ID']), str(ref['Esmecata ID']), ';'.join(ref['obs'])]))


# MAIN
if __name__ == "__main__":
    RUN_NAME = sys.argv[1]
    INPUT_RANK_FILE = os.path.join(RUN_NAME, '0_proteomes', 'association_taxon_taxID.json')
    ESMECATA_RANK_FILE = os.path.join(RUN_NAME, '0_proteomes', 'proteome_tax_id.tsv')
    OUTPUT = os.path.join(RUN_NAME, 'Krona')
    if not os.path.exists(OUTPUT):
        os.mkdir(OUTPUT)
    create_krona_input_in(INPUT_RANK_FILE, RUN_NAME, OUTPUT)
    create_krona_input_esmecata(ESMECATA_RANK_FILE, RUN_NAME, OUTPUT)
    create_comp_taxonomy_file(INPUT_RANK_FILE, ESMECATA_RANK_FILE, RUN_NAME, OUTPUT)
