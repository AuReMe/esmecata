# Copyright (C) 2024-2025 Arnaud Belcour -  Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import csv
from esmecata.utils import send_uniprot_sparql_query

from ete3 import NCBITaxa
def get_taxon_proteomes(output_file, taxon_rank='Genus'):
    uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql'

    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT DISTINCT ?taxon_name ?taxon

    WHERE
    {{
        ?taxon up:rank up:{0} .
        ?taxon up:scientificName ?taxon_name .
        ?proteome rdf:type up:Proteome .
        ?proteome up:organism ?organism .
        ?organism rdfs:subClassOf ?taxon .
    }}
    """.format(taxon_rank)
    csvreader = send_uniprot_sparql_query(uniprot_sparql_query, uniprot_sparql_endpoint)
    taxon_proteomes = []
    for line in csvreader:
        taxon = line[0].split('/')[-1]
        taxon_id = line[1].split('/')[-1]
        taxon_proteomes.append([taxon, taxon_id])
    ncbi = NCBITaxa()
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['observation_name', 'taxonomic_affiliation', 'taxon_id'])
        for taxon, taxon_id in taxon_proteomes:
            try:
                lineage = ncbi.get_lineage(taxon_id)
                convert_taxon = ''.join(c if c.isalnum() or c in '-_' else ('_c' + str(ord(c)) + '_') for c in taxon)
                observation_name = convert_taxon + '__taxid__' + str(taxon_id)
                lineage_names = [ncbi.get_taxid_translator([lineage_tax_id])[lineage_tax_id] for lineage_tax_id in lineage]
                csvwriter.writerow([observation_name, ';'.join(lineage_names), taxon_id])
            except:
                continue
get_taxon_proteomes('esmecata_precomputed_species.tsv', 'Species')
get_taxon_proteomes('esmecata_precomputed_genus.tsv')
get_taxon_proteomes('esmecata_precomputed_family.tsv', 'Family')
get_taxon_proteomes('esmecata_precomputed_order.tsv', 'Order')
get_taxon_proteomes('esmecata_precomputed_class.tsv', 'Class')
get_taxon_proteomes('esmecata_precomputed_phylum.tsv', 'Phylum')