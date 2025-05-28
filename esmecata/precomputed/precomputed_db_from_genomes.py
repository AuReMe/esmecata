# Copyright (C) 2024-2025 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Microcosme
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
import os
import gzip
import pandas as pd
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete4 import NCBITaxa

from esmecata.core.clustering import make_clustering

def create_tree(taxonomic_file, proteome_tax_id_file):
    ncbi = NCBITaxa()
    df_taxonomic = pd.read_csv(taxonomic_file, sep='\t', index_col=0)
    all_taxa = df_taxonomic['taxon_id'].unique()
    all_taxa_expanded = []
    for taxon in all_taxa:
        all_taxa_expanded.extend(ncbi.get_lineage(taxon))
    
    taxon_id_to_names = ncbi.get_taxid_translator(all_taxa_expanded)
    taxon_names = {str(taxon): taxon_id_to_names[taxon] for taxon in taxon_id_to_names}
    tree = ncbi.get_topology([org_tax_id for org_tax_id in all_taxa_expanded])

    taxon_genomes = {str(tax_id): df_taxonomic[df_taxonomic['taxon_id']==tax_id].index.tolist() for tax_id in df_taxonomic['taxon_id']}

    tax_id_proteomes = {}
    for node in tree.traverse():
        node_to_checks = [desc.name for desc in node.descendants()]
        node_to_checks.append(node.name)
        tax_id_proteome = [proteome for tax_id in node_to_checks if tax_id in taxon_genomes
                                        for proteome in taxon_genomes[tax_id]]
        tax_id_proteomes[node.name] = tax_id_proteome

    # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
    with open(proteome_tax_id_file, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        csvwriter.writerow(['observation_name', 'name', 'tax_id', 'tax_id_name', 'tax_rank', 'proteome'])
        for tax_id in tax_id_proteomes:
            tax_name = taxon_names[tax_id]
            proteomes = tax_id_proteomes[tax_id]
            # Convert characters in tax_name into unicode code.
            convert_tax_name = ''.join(c if c.isalnum() or c in '-_' else ('_c' + str(ord(c)) + '_') for c in tax_name)
            tax_id_name = convert_tax_name + '__taxid__' + str(tax_id)
            tax_rank = ncbi.get_rank([int(tax_id)])[int(tax_id)]
            if tax_rank in ['species', 'genus', 'family', 'order', 'class', 'phylum']:
                csvwriter.writerow([tax_id_name, tax_name, tax_id, tax_id_name, tax_rank, ','.join(proteomes)])


def read_gbk_prots(gbk_file, qualifier='locus_tag'):
    protein_records = []
    with open(gbk_file, "r") as gbk:
        for seq_record in SeqIO.parse(gbk, "genbank"):
            seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
            for seq_feature in seq_feature_cds:
                fasta_id = seq_feature.qualifiers[qualifier][0]
                if 'translation' in seq_feature.qualifiers:
                    fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0]), id=fasta_id, description=fasta_id)
                    protein_records.append(fasta_record)

    return protein_records

def run_proteomes_step(input_taxon_id, input_proteome_folder, output_folder, nb_core=1, clust_threshold=0.5, mmseqs_options=None, linclust=None, remove_tmp=None):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteome_step_folder = os.path.join(output_folder, '0_proteomes')
    if not os.path.exists(proteome_step_folder):
        os.mkdir(proteome_step_folder)

    proteome_tax_id_file = os.path.join(proteome_step_folder, 'proteome_tax_id.tsv')
    create_tree(input_taxon_id, proteome_tax_id_file)
    output_proteome_folder = os.path.join(proteome_step_folder, 'proteomes')
    if not os.path.exists(output_proteome_folder):
        os.mkdir(output_proteome_folder)

    for proteome_file in os.listdir(input_proteome_folder):
        proteome_filepath = os.path.join(input_proteome_folder, proteome_file)
        output_proteome_filepath = os.path.join(output_proteome_folder, proteome_file.replace('gbff', 'faa'))
        output_proteome_gz_filepath = os.path.join(output_proteome_folder, proteome_file.replace('gbff', 'faa')+'.gz')

        protein_records = read_gbk_prots(proteome_filepath, qualifier='locus_tag')
        SeqIO.write(protein_records, output_proteome_filepath, 'fasta')
        with open(output_proteome_filepath, 'rb') as input_file:
            with gzip.open(output_proteome_gz_filepath, 'wb') as compressed_file:
                shutil.copyfileobj(input_file, compressed_file)
        os.remove(output_proteome_filepath)

    clustering_step_folder = os.path.join(output_folder, '1_clustering')
    if not os.path.exists(clustering_step_folder):
        os.mkdir(clustering_step_folder)

    make_clustering(proteome_step_folder, clustering_step_folder, nb_core, clust_threshold, mmseqs_options, linclust, remove_tmp)

run_proteomes_step('taxon_id.tsv', 'test_prot', 'output_folder', 35)