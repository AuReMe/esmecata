import csv
import gzip
import json
import math
import os
import pandas as pd
import random
import requests
import shutil
import sys
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from ete3 import NCBITaxa, is_taxadb_up_to_date
from io import StringIO
from SPARQLWrapper import SPARQLWrapper, TSV

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release, is_valid_file, is_valid_dir

def associate_taxon_to_taxon_id(taxonomies, ncbi):
    tax_id_names = {}
    clusters_dicts = {}
    json_cluster_taxons = {}

    if taxonomies == {}:
        print('Empty taxonomy dictionary.')
        return

    # For each taxonomy find the taxonomy ID corresponding to each taxon.
    for cluster in taxonomies:
        if isinstance(taxonomies[cluster], str):
            taxons = [taxon for taxon in taxonomies[cluster].split(';')]
            names = ncbi.get_name_translator(taxons)
            taxon_ids = OrderedDict()
            for taxon in taxons:
                taxon_translations = ncbi.get_name_translator([taxon])
                if taxon_translations == {}:
                    print('For {0}, no taxon ID has been found associated to the taxon "{1}" in the NCBI taxonomy of ete3.'.format(cluster, taxon))
                    taxon_translations = {taxon: ['not_found']}
                taxon_ids.update(taxon_translations)
            invert_names = {v: k for k, vs in names.items() for v in vs }
            tax_id_names.update(invert_names)
            clusters_dicts[cluster] = [names[tax_name] for tax_name in taxonomies[cluster].split(';') if tax_name in names]
            json_cluster_taxons[cluster] = taxon_ids

    return tax_id_names, json_cluster_taxons

def filter_taxon(json_cluster_taxons, ncbi):
    # If there is multiple taxon ID for a taxon, use the taxonomy to find the most relevant taxon.
    # The most relevant taxon is the one with the most overlapping lineage with the taxonomy.
    taxon_to_modify = {}
    for cluster in json_cluster_taxons:
        cluster_taxons = json_cluster_taxons[cluster].values()
        for taxon in json_cluster_taxons[cluster]:
            taxon_ids = json_cluster_taxons[cluster][taxon]
            if len(taxon_ids) > 1:
                relevant_taxon = None
                previous_nb_shared_ids = 0
                for taxon_id in taxon_ids:
                    taxons_ids = list(next(zip(*cluster_taxons)))
                    lineage = ncbi.get_lineage(taxon_id)
                    nb_shared_ids = len(set(taxons_ids).intersection(set(lineage)))
                    if nb_shared_ids > previous_nb_shared_ids:
                        relevant_taxon = taxon_id
                    previous_nb_shared_ids = nb_shared_ids
                taxon_to_modify[cluster] = {}
                taxon_to_modify[cluster][taxon] = [relevant_taxon]

    for cluster in taxon_to_modify:
        for taxon in taxon_to_modify[cluster]:
            json_cluster_taxons[cluster][taxon] = taxon_to_modify[cluster][taxon]

    return json_cluster_taxons

def rest_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes):
    proteomes = []
    organism_ids = {}
    # Find proteomes associated with taxon.
    # Take reference proteomes with "reference:yes".
    # Avoid redundant and excluded proteomes with "redundant%3Ano+excluded%3Ano".
    # Use "format=tab" to easily handle the ouput.
    http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+reference:yes+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)

    if all_proteomes:
        http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)

    # If esmecata does not find proteomes with only reference, search for all poroteomes even if they are not reference.
    all_http_str = 'https://www.uniprot.org/proteomes/?query=taxonomy:{0}+redundant%3Ano+excluded%3Ano&format=tab'.format(tax_id)

    response = requests.get(http_str)
    # Raise error if we have a bad request.
    response.raise_for_status()

    response_proteome_status = False

    with requests.get(http_str) as proteome_response:
        proteome_response_text = proteome_response.text
        if proteome_response_text != '':
            csvreader = csv.reader(proteome_response_text.splitlines(), delimiter='\t')
            response_proteome_status = True
            # Avoid header.
            next(csvreader)
        else:
            csvreader = []

    if response_proteome_status is False:
        time.sleep(1)
        print('{0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(taxon, tax_id, tax_name))
        with requests.get(all_http_str) as proteome_response:
            proteome_response_text = proteome_response.text
            if proteome_response_text != '':
                csvreader = csv.reader(proteome_response_text.splitlines(), delimiter='\t')
                response_proteome_status = True
                # Avoid header.
                next(csvreader)
            else:
                csvreader = []

    for line in csvreader:
        proteome = line[0]
        completness = line[6]
        org_tax_id = line[2]

        # Check that proteome has busco score.
        if busco_percentage_keep:
            if line[4] != '':
                busco_percentage = float(line[4].split(':')[1].split('%')[0])
                if busco_percentage >= busco_percentage_keep and completness == 'full':
                    proteomes.append(proteome)
                    if org_tax_id not in organism_ids:
                        organism_ids[org_tax_id] = [proteome]
                    else:
                        organism_ids[org_tax_id].append(proteome)
        else:
            if completness == 'full':
                proteomes.append(proteome)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome]
                else:
                    organism_ids[org_tax_id].append(proteome)

    return proteomes, organism_ids


def sparql_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint='https://sparql.uniprot.org/sparql'):
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)

    # Test with SPARQL query
    # First FILTER NOT EXISTS to avoid redundant proteomes.
    # Second FILTER NOT EXISTS to avoid excluded proteomes.
    # OPTIONAL allows to retrive reference and non reference proteomes and split them.
    # From: uniprot graph location https://sparql.uniprot.org/.well-known/void%3E
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX busco: <http://busco.ezlab.org/schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

    SELECT DISTINCT ?proteome ?score ?fragmented ?missing ?organism ?completion ?type
    FROM <http://sparql.uniprot.org/proteomes>
    FROM <http://sparql.uniprot.org/taxonomy>
    WHERE
    {{
        ?proteome rdf:type up:Proteome .
        ?proteome up:organism ?organism .
        ?organism rdfs:subClassOf* taxon:{0} .
        OPTIONAL {{
            ?proteome busco:has_score ?busco .
            ?busco busco:complete ?score .
            ?busco busco:fragmented ?fragmented .
            ?busco busco:missing ?missing .
        }}
        OPTIONAL {{
            ?proteome rdfs:seeAlso ?metadataID .
            ?metadataID rdfs:label ?completion .
        }}
        FILTER NOT EXISTS {{
        ?proteome up:redundantTo ?proteome2 .
            }}
        FILTER NOT EXISTS {{
        ?proteome up:exclusionReason ?exclusion .
                }}
        OPTIONAL {{
            ?proteome rdf:type ?type .
            VALUES ?type {{up:Representative_Proteome}}
        }}
    }}""".format(tax_id)

    sparql.setQuery(uniprot_sparql_query)
    # Parse output.
    sparql.setReturnFormat(TSV)

    results = sparql.query().convert().decode('utf-8')
    if results != '':
        csvreader = csv.reader(StringIO(results), delimiter='\t')
        # Avoid header.
        next(csvreader)
    else:
        csvreader = []

    proteomes = []
    organism_ids = {}
    reference_proteomes = []
    other_proteomes = []
    for line in csvreader:
        proteome_id = line[0].split('/')[-1]
        if line[1] != '':
            score = float(line[1].split('^^')[0])
        else:
            score = None
        if line[2] != '':
            fragmented = float(line[2].split('^^')[0])
        else:
            fragmented = None
        if line[3] != '':
            missing = float(line[3].split('^^')[0])
        else:
            missing = None
        org_tax_id = line[4].split('/')[-1]
        completness = line[5]
        if line[6] == '':
            reference = False
        else:
            reference = True

        if all([score, fragmented, missing]) is True:
            busco_percentage = (score / (score+fragmented+missing)) * 100
        else:
            busco_percentage = None
        # Check that proteome has busco score.
        if busco_percentage_keep and busco_percentage:
            if busco_percentage >= busco_percentage_keep and completness == 'full':
                if reference == True:
                    reference_proteomes.append(proteome_id)
                else:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)
        else:
            if completness == 'full':
                if reference == True:
                    reference_proteomes.append(proteome_id)
                else:
                    other_proteomes.append(proteome_id)
                if org_tax_id not in organism_ids:
                    organism_ids[org_tax_id] = [proteome_id]
                else:
                    organism_ids[org_tax_id].append(proteome_id)

    if all_proteomes:
        proteomes = other_proteomes
    else:
        if len(reference_proteomes) == 0:
            print('{0}: No reference proteomes found for {1} ({2}) try non-reference proteomes.'.format(taxon, tax_id, tax_name))
            proteomes = other_proteomes
        else:
            proteomes = reference_proteomes

    return proteomes, organism_ids


def find_proteomes_tax_ids(json_cluster_taxons, ncbi, busco_percentage_keep=None, all_proteomes=None, uniprot_sparql_endpoint=None):
    # Query the Uniprot proteomes to find all the proteome IDs associated to taxonomy.
    # If there is more thant 100 proteomes we do not keep it because there is too many proteome.
    print('Find proteome ID associated to taxonomy')
    proteomes_ids = {}
    single_proteomes = {}
    tax_id_not_founds = {}
    tax_id_founds = {}
    for taxon in json_cluster_taxons:
        for tax_name in reversed(json_cluster_taxons[taxon]):
            tax_id = json_cluster_taxons[taxon][tax_name][0]

            # If tax_id has already been found use the corresponding proteomes without new requests.
            if tax_id in tax_id_founds:
                proteomes_ids[taxon] = (tax_id, tax_id_founds[tax_id])
                if len(tax_id_founds[tax_id]) == 1:
                    single_proteomes[taxon] = (tax_id, tax_id_founds[tax_id])
                break

            # If tax_id has not been found with a request do not try a new request with the same tax_id.
            if tax_id in tax_id_not_founds:
                continue

            if uniprot_sparql_endpoint:
                proteomes, organism_ids = sparql_query_proteomes(taxon, tax_id, tax_name, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint)
            else:
                proteomes, organism_ids = rest_query_proteomes(taxon, tax_id, tax_name, all_proteomes, busco_percentage_keep)

            # Answer is empty no corresponding proteomes to the tax_id.
            if len(proteomes) == 0:
                if tax_id not in tax_id_not_founds:
                    tax_id_not_founds[tax_id] = [tax_name]
                else:
                    tax_id_not_founds[tax_id].append(tax_name)
                continue

            elif len(proteomes) > 0 and len(proteomes) < 100:
                print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(proteomes)))
                proteomes_ids[taxon] = (tax_id, proteomes)
                tax_id_founds[tax_id] = proteomes
                if len(proteomes) == 1:
                    single_proteomes[taxon] = (tax_id, proteomes)
                break

            elif len(proteomes) >= 100:
                print('More than 99 proteomes are associated to the taxa {0} associated to {1}, esmecata will randomly select around 100 proteomes with respect to the taxonomy proportion.'.format(taxon, tax_name))
                tree = ncbi.get_topology([org_tax_id for org_tax_id in organism_ids])

                # For each direct descendant taxon of the tree root (our tax_id), we will look for the proteomes inside these subtaxons.
                childs = {}
                for parent_node in tree.get_children():
                    parent_tax_id = str(parent_node.taxid)
                    if parent_node.get_descendants() != []:
                        for child_node in parent_node.get_descendants():
                            child_tax_id = str(child_node.taxid)
                            if parent_tax_id not in childs:
                                if child_tax_id in organism_ids:
                                    childs[parent_tax_id] = organism_ids[child_tax_id]
                            else:
                                if child_tax_id in organism_ids:
                                    childs[parent_tax_id].extend(organism_ids[child_tax_id])
                    else:
                        if parent_tax_id in organism_ids:
                            childs[parent_tax_id] = organism_ids[parent_tax_id]
                        else:
                            childs[parent_tax_id].extend(organism_ids[parent_tax_id])

                # For each direct descendant taxon, compute the number of proteomes in their childrens.
                elements = [v for k,v in childs.items()]
                elements_counts = [len(element) for element in elements]
                # Compute the proportion of proteomes for each direct descendant taxon compare to the total number of proteomes in all the descendant taxons.
                percentages = [(i/sum(elements_counts))*100  for i in elements_counts]
                # Can be superior to 100 if there is a lot of data with around 0.xxx percentage.
                percentages_round = [math.ceil(percentage) if percentage < 1 else math.floor(percentage) for percentage in percentages]

                # Choose randomly a number of proteomes corresponding to the computed percentage.
                all_proteomes = []
                for index, element in enumerate(elements):
                    percentage_to_keep = percentages_round[index]
                    proteomes_to_keep = random.sample(element, percentage_to_keep)
                    all_proteomes.extend(proteomes_to_keep)
                proteomes_ids[taxon] = (tax_id, all_proteomes)
                tax_id_founds[tax_id] = all_proteomes
                print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(taxon, tax_name, len(all_proteomes)))
                break

            time.sleep(1)

    return proteomes_ids, single_proteomes, tax_id_not_founds

def sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint):
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)

    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX proteome: <http://purl.uniprot.org/proteomes/>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>

    SELECT ?protein ?name ?sequence ?sequenceaa ?review
    FROM <http://sparql.uniprot.org/uniprot>
    FROM <http://sparql.uniprot.org/proteomes>
    WHERE
    {{
    ?protein a up:Protein ;
            up:mnemonic ?name ;
            up:reviewed ?review ;
            up:proteome ?genomicComponent .
    ?proteome skos:narrower ?genomicComponent .
    ?protein up:sequence ?sequence .
    ?sequence rdf:value ?sequenceaa .
    VALUES (?proteome) {{ (proteome:{0}) }}
    }}""".format(proteome)

    sparql.setQuery(uniprot_sparql_query)
    # Parse output.
    sparql.setReturnFormat(TSV)

    results = sparql.query().convert().decode('utf-8')
    if results != '':
        csvreader = csv.reader(StringIO(results), delimiter='\t')
        # Avoid header.
        next(csvreader)
    else:
        csvreader = []

    records = []
    already_added_proteins = []
    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        protein_name = line[1]
        protein_isoform = line[2].split('/')[-1]
        protein_seq = line[3]
        prot_review = line[4].split('^^')[0]
        # Check protein review.
        # If review == Swissprot else Trembl.
        if prot_review == 'true':
            prefix_record = 'sp'
        elif prot_review == 'false':
            prefix_record = 'tr'
        else:
            print(protein_id, prot_review)

        # Take only the canonical sequences (ending with -1):
        if protein_isoform.endswith('-1'):
            record_id = f'{prefix_record}|{protein_id}|{protein_name}'
            records.append(SeqRecord(Seq(protein_seq), id=record_id, description=''))
            already_added_proteins.append(protein_id)

    intermediary_file = output_proteome_file[:-3]
    SeqIO.write(records, intermediary_file, 'fasta')

    with open(intermediary_file, 'rb') as input_file:
        with gzip.open(output_proteome_file, 'wb') as output_file:
            shutil.copyfileobj(input_file, output_file)
    os.remove(intermediary_file)


def retrieve_proteomes(input_file, output_folder, busco_percentage_keep=None, ignore_taxadb_update=None, all_proteomes=None, uniprot_sparql_endpoint=None):
    if is_valid_file(input_file) == False:
        print('The input {0} is not a valid file pathname.'.format(input_file))
        sys.exit()

    try:
        ete_taxadb_up_to_date = is_taxadb_up_to_date()
    except:
        ncbi = NCBITaxa()
        ete_taxadb_up_to_date = is_taxadb_up_to_date()

    if ete_taxadb_up_to_date is False:
        print('''WARNING: ncbi taxonomy database is not up to date with the last NCBI Taxonomy. Update it using:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()''')
        if ignore_taxadb_update is None:
            print('If you want to stil use esmecata with the old taxonomy database use the option --ignore-taxadb-update/ignore_taxadb_update.')
            return
        else:
            print('--ignore-taxadb-update/ignore_taxadb_update option detected, esmecata will continue with this version.')

    is_valid_dir(output_folder)

    if '.xlsx' in input_file:
        df = pd.read_excel(input_file)
    if '.tsv' in input_file:
        df = pd.read_csv(input_file, sep='\t')
    if '.csv' in input_file:
        df = pd.read_csv(input_file, sep=',')

    # Set index on the column containing the OTU name.
    df.set_index('observation_name', inplace=True)

    # taxonomy is the column containing the taxonomy separated by ';': phylum;class;order;family;genus;genus + species
    taxonomies = df.to_dict()['taxonomy']

    proteome_cluster_tax_id_file = os.path.join(output_folder, 'proteome_cluster_tax_id.tsv')

    if not os.path.exists(proteome_cluster_tax_id_file):
        ncbi = NCBITaxa()

        tax_id_names, json_cluster_taxons = associate_taxon_to_taxon_id(taxonomies, ncbi)

        json_cluster_taxons = filter_taxon(json_cluster_taxons, ncbi)

        json_log = os.path.join(output_folder, 'log.json')
        with open(json_log, 'w') as ouput_file:
            json.dump(json_cluster_taxons, ouput_file, indent=4)


    if not os.path.exists(proteome_cluster_tax_id_file):
        proteomes_ids, single_proteomes, tax_id_not_founds = find_proteomes_tax_ids(json_cluster_taxons, ncbi, busco_percentage_keep, all_proteomes, uniprot_sparql_endpoint)

        proteome_to_download = []
        for proteomes_id in proteomes_ids:
            proteome_to_download.extend(proteomes_ids[proteomes_id][1])
        proteome_to_download = set(proteome_to_download)

        # Write for each taxon the corresponding tax ID, the name of the taxon and the proteome associated with them.
        with open(proteome_cluster_tax_id_file, 'w') as out_file:
            csvwriter = csv.writer(out_file, delimiter='\t')
            csvwriter.writerow(['cluster', 'name', 'tax_id', 'proteome'])
            for cluster in proteomes_ids:
                csvwriter.writerow([cluster, tax_id_names[int(proteomes_ids[cluster][0])], proteomes_ids[cluster][0], ','.join(proteomes_ids[cluster][1])])
    else:
        proteome_to_download = []
        proteomes_ids = {}
        single_proteomes = {}
        with open(proteome_cluster_tax_id_file, 'r') as input_file:
            csvreader = csv.reader(input_file, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                observation_name = line[0]
                name = line[1]
                tax_id = line[2]
                proteomes = line[3].split(',')
                proteome_to_download.extend(proteomes)
                if len(proteomes) == 1:
                    single_proteomes[observation_name] = (tax_id, proteomes)
                proteomes_ids[observation_name] = (tax_id, proteomes)
                print('{0} will be associated to the taxon "{1}" with {2} proteomes.'.format(observation_name, name, len(proteomes)))

        proteome_to_download = set(proteome_to_download)

    # Download all the proteomes in tmp folder.
    print('Download proteome')
    tmp_folder = os.path.join(output_folder, 'tmp_proteome')
    is_valid_dir(tmp_folder)

    for proteome in proteome_to_download:
        output_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
        if not os.path.exists(output_proteome_file):
            if uniprot_sparql_endpoint is not None:
                sparql_get_protein_seq(proteome, output_proteome_file, uniprot_sparql_endpoint)
            else:
                with requests.get('https://www.uniprot.org/uniprot/?query=proteome:{0}&format=fasta&compress=yes'.format(proteome)) as proteome_response:
                    with open(output_proteome_file, 'wb') as f:
                        f.write(proteome_response.content)
        time.sleep(1)

    # Download Uniprot metadata and create a json file containing them.
    if uniprot_sparql_endpoint:
        uniprot_releases = get_sparql_uniprot_release(uniprot_sparql_endpoint)
    else:
        uniprot_releases = get_rest_uniprot_release()

    uniprot_metadata_file = os.path.join(output_folder, 'uniprot_release_metadata_proteomes.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    print('Creating result folder')
    # Create a result folder which contains one sub-folder per OTU.
    # Each OTU sub-folder will contain the proteome found.
    result_folder = os.path.join(output_folder, 'result')
    is_valid_dir(result_folder)

    for cluster in proteomes_ids:
        output_cluster = os.path.join(result_folder, cluster)
        is_valid_dir(output_cluster)
        for proteome in proteomes_ids[cluster][1]:
            input_proteome_file = os.path.join(tmp_folder, proteome+'.faa.gz')
            output_proteome = os.path.join(output_cluster, proteome+'.faa')
            with gzip.open(input_proteome_file, 'rb') as f_in:
                with open(output_proteome, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
