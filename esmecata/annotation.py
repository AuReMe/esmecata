import csv
import json
import os
import re
import urllib.parse
import urllib.request

from esmecata.utils import get_rest_uniprot_release, get_sparql_uniprot_release

def rest_query_uniprot_to_retrieve_function(protein_queries, output_dict):
    url = 'https://www.uniprot.org/uploadlists/'

    # Column names can be found at: https://www.uniprot.org/help/uniprotkb_column_names
    # from and to are linked to mapping ID: https://www.uniprot.org/help/api%5Fidmapping
    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
        'query': protein_queries,
        'columns': 'id,protein names,reviewed,go,ec,interpro(db_abbrev),rhea-id,genes(PREFERRED)'
    }
    go_pattern = re.compile(r'\[GO:(?P<go>\d{7})\]')
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        csvreader = csv.reader(f.read().decode('utf-8').splitlines(), delimiter='\t')
    # Avoid header
    next(csvreader)

    results = {}
    for line in csvreader:
        protein_name = line[1]
        if line[2] == 'reviewed':
            review = True
        else:
            review = False
        go_result = go_pattern.findall(line[3])
        gos = ['GO:'+ go_term for go_term in go_result]
        ec_numbers = [ec for ec in line[4].split('; ') if ec != '']
        interpros =  [interpro for interpro in line[5].split('; ') if interpro != '']
        rhea_ids = [rhea_id for rhea_id in line[6].split('; ') if rhea_id != '']
        gene_name = line[7]
        results[line[0]] = [protein_name, review, gos, ec_numbers, interpros, rhea_ids, gene_name]
        output_dict.update(results)

    return output_dict


def sparql_query_uniprot_to_retrieve_function(proteomes, output_dict, uniprot_sparql_endpoint):
    from SPARQLWrapper import SPARQLWrapper, TSV
    from io import StringIO

    proteomes = ' '.join(['( proteome:'+proteome+' )' for proteome in proteomes])

    sparql = SPARQLWrapper(uniprot_sparql_endpoint)

    # uniprotkb:ID
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
    PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX proteome: <http://purl.uniprot.org/proteomes/>

    SELECT ?protein
    (GROUP_CONCAT(DISTINCT ?fullName; separator=";") AS ?name)
    (GROUP_CONCAT(DISTINCT ?goTerm; separator=";") AS ?go)
    (GROUP_CONCAT(DISTINCT ?ecNumber; separator=";") AS ?ec)
    (GROUP_CONCAT(DISTINCT ?ecNumber2; separator=";") AS ?ec2)
    (GROUP_CONCAT(DISTINCT ?interpro; separator=";") AS ?ipr)
    (GROUP_CONCAT(DISTINCT ?rheaReaction; separator=";") AS ?rhea)
    (GROUP_CONCAT(DISTINCT ?reviewed; separator=";") AS ?review)
    (GROUP_CONCAT(DISTINCT ?geneLabel; separator=";") AS ?geneName)
    (GROUP_CONCAT(DISTINCT ?subName; separator=";") AS ?submitName)
    WHERE {{
        ?protein a up:Protein ;
            up:proteome ?genomicComponent .
        ?proteome skos:narrower ?genomicComponent .
        OPTIONAL {{
            ?protein up:reviewed ?reviewed  .
        }}
        OPTIONAL {{
            ?protein up:annotation ?annot .
            ?annot a up:Catalytic_Activity_Annotation ;
                up:catalyticActivity ?catalyticAct .
            ?catalyticAct up:catalyzedReaction ?rheaReaction .
            ?catalyticAct up:enzymeClass ?ecNumber2 .
        }}
        OPTIONAL {{
            ?protein up:classifiedWith ?goTerm .
            FILTER (regex(str(?goTerm), "GO")) .
        }}
        OPTIONAL {{
            ?protein rdfs:seeAlso ?interpro .
            FILTER (regex(str(?interpro), "interpro")) .
        }}
        OPTIONAL {{
            ?protein up:enzyme ?ecNumber .
        }}
        OPTIONAL {{
            ?protein up:recommendedName ?recommendName .
            ?recommendName up:fullName ?fullName .
        }}
        OPTIONAL {{
            ?protein up:submittedName ?submittedName .
            ?submittedName up:fullName ?subName .
        }}
        OPTIONAL {{
            ?protein up:encodedBy ?gene .
            ?gene skos:prefLabel ?geneLabel .
        }}
    VALUES (?proteome) {{ {0} }}
    }}
    GROUP BY ?protein
    """.format(proteomes)

    sparql.setQuery(uniprot_sparql_query)
    # Parse output.
    sparql.setReturnFormat(TSV)
    results = sparql.query().convert().decode('utf-8')
    csvreader = csv.reader(StringIO(results), delimiter='\t')
    # Avoid header.
    next(csvreader)
    results = {}
    for line in csvreader:
        protein_id = line[0].split('/')[-1]
        protein_name = line[1].split('^^')[0]
        go_terms = [go_uri.split('obo/')[1].replace('_', ':') for go_uri in line[2].split('^^')[0].split(';') if go_uri != '']
        ec_numbers = [ec_uri.split('enzyme/')[1] for ec_uri in line[3].split('^^')[0].split(';') if ec_uri != '']
        ec_numbers.extend([ec_uri.split('enzyme/')[1] for ec_uri in line[4].split('^^')[0].split(';') if ec_uri != ''])
        ec_numbers = list(set(ec_numbers))
        interpros = [ipr_uri.split('interpro/')[1] for ipr_uri in line[5].split('^^')[0].split(';') if ipr_uri != '']
        rhea_ids = [rhea_uri.split('rdf.rhea-db.org/')[1] for rhea_uri in line[6].split('^^')[0].split(';') if rhea_uri != '']
        review = line[7].split('^^')[0]
        if review == '0':
            review = False
        elif review == '1':
            review = True
        gene_name = line[8].split('^^')[0]

        submitted_name = [submitted_name for submitted_name in line[9].split('^^')[0].split(';') if submitted_name != '']
        if submitted_name != []:
            submitted_name = submitted_name[0]
        else:
            submitted_name = ''

        # If there is no recommendedName, use submittedName instead.
        if protein_name == '':
            if submitted_name != '':
                protein_name = submitted_name

        results[protein_id] = [protein_name, review, go_terms, ec_numbers, interpros, rhea_ids, gene_name]
        output_dict.update(results)

    return output_dict


def chunks(lst, n):
    """Yield successive n-sized chunks from list.
    Form: https://stackoverflow.com/a/312464
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def create_pathologic(base_filename, protein_annotations, protein_set, pathologic_output_file):
    with open(pathologic_output_file, 'w', encoding='utf-8') as element_file:
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        element_file.write(';; ' + base_filename + '\n')
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        for protein in protein_annotations:
            if protein in protein_set:
                element_file.write('ID\t' + protein + '\n')
                if protein_annotations[protein][3] != '':
                    element_file.write('NAME\t' + protein_annotations[protein][3] + '\n')
                else:
                    element_file.write('NAME\t' + protein + '\n')
                if protein_annotations[protein][0] != '':
                    element_file.write('FUNCTION\t' + protein_annotations[protein][0] + '\n')
                element_file.write('PRODUCT-TYPE\tP' + '\n')
                element_file.write('PRODUCT-ID\tprot ' + protein + '\n')
                element_file.write('DBLINK\tUNIPROT:' + protein + '\n')
                for go in protein_annotations[protein][1]:
                    element_file.write('GO\t' + go + '\n')
                for ec in protein_annotations[protein][2]:
                    element_file.write('EC\t' + ec + '\n')
                element_file.write('//\n\n')


def annotate_proteins(input_folder, output_folder, uniprot_sparql_endpoint):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Download Uniprot metadata and create a json file containing them.
    if uniprot_sparql_endpoint:
        uniprot_releases = get_sparql_uniprot_release(uniprot_sparql_endpoint)
    else:
        uniprot_releases = get_rest_uniprot_release()

    uniprot_metadata_file = os.path.join(output_folder, 'uniprot_release_metadata.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(uniprot_releases, ouput_file, indent=4)

    if uniprot_sparql_endpoint:
        input_proteomes = {}
        proteome_cluster_tax_id_file = os.path.join(input_folder, 'proteome_cluster_tax_id.tsv')
        with open(proteome_cluster_tax_id_file, 'r') as tax_id_file:
            csvreader = csv.reader(tax_id_file, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                input_proteomes[line[0]] = line[3].split(',')

    reference_protein_path = os.path.join(input_folder, 'reference_proteins')
    for input_file in os.listdir(reference_protein_path):
        base_file = os.path.basename(input_file)
        base_filename = os.path.splitext(base_file)[0]
        proteins = []
        reference_protein_list_file = os.path.join(reference_protein_path, input_file)

        reference_proteins = {}
        with open(reference_protein_list_file, 'r') as input_file:
            csvreader = csv.reader(input_file, delimiter='\t')
            for line in csvreader:
                proteins.extend(line)
                reference_proteins[line[0]] = line[1:]

        set_proteins = list(set(proteins))

        # Query Uniprot to get the annotation of each proteins.
        # The limit of 20 000 proteins per query comes from the help of Uniprot:
        # https://www.uniprot.org/help/uploadlists
        output_dict = {}
        if len(set_proteins) < 20000:
            protein_queries = ' '.join(set_proteins)
            if uniprot_sparql_endpoint:
                proteomes = input_proteomes[base_filename]
                sparql_query_uniprot_to_retrieve_function(proteomes, output_dict, uniprot_sparql_endpoint)
            else:
                rest_query_uniprot_to_retrieve_function(protein_queries, output_dict)
        else:
            protein_chunks = chunks(proteins, 20000)
            for chunk in protein_chunks:
                protein_queries = ' '.join(chunk)
                if uniprot_sparql_endpoint:
                    proteomes = input_proteomes[base_filename]
                    output_dict = sparql_query_uniprot_to_retrieve_function(proteomes, output_dict, uniprot_sparql_endpoint)
                else:
                    output_dict = rest_query_uniprot_to_retrieve_function(protein_queries, output_dict)

        annotation_folder = os.path.join(output_folder, 'annotation')
        if not os.path.exists(annotation_folder):
            os.mkdir(annotation_folder)
        annotation_file = os.path.join(annotation_folder, base_filename+'.tsv')
        with open(annotation_file, 'w') as output_tsv:
            csvwriter = csv.writer(output_tsv, delimiter='\t')
            csvwriter.writerow(['protein_id', 'protein_name', 'review', 'gos', 'ecs', 'interpros', 'rhea_ids', 'gene_name'])
            for protein in output_dict:
                protein_name = output_dict[protein][0]
                protein_review_satus = str(output_dict[protein][1])
                go_tersm = ','.join(output_dict[protein][2])
                ec_numbers = ','.join(output_dict[protein][3])
                interpros = ','.join(output_dict[protein][4])
                rhea_ids = ','.join(output_dict[protein][5])
                gene_name = output_dict[protein][6]
                csvwriter.writerow([protein, protein_name, protein_review_satus, go_tersm, ec_numbers, interpros, rhea_ids, gene_name])

        # For each reference protein, get the annotation of the other proteins clustered with it and add to its annotation.
        protein_annotations = {}
        for reference_protein in reference_proteins:
            gos = set()
            ecs = set()
            gene_name = ''
            for protein in reference_proteins[reference_protein]:
                if protein in output_dict:
                    if output_dict[protein][2] != []:
                        gos.update(set(output_dict[protein][2]))
                    if output_dict[protein][3] != []:
                        ecs.update(set(output_dict[protein][3]))
                    if output_dict[protein][6] != '':
                        gene_name = output_dict[protein][6]
            protein_annotations[reference_protein] = [output_dict[protein][0], gos, ecs, gene_name]
  
        annotation_reference_folder = os.path.join(output_folder, 'annotation_reference')
        if not os.path.exists(annotation_reference_folder):
            os.mkdir(annotation_reference_folder)
        annotation_reference_file = os.path.join(annotation_reference_folder, base_filename+'.tsv')
        with open(annotation_reference_file, 'w') as output_tsv:
            csvwriter = csv.writer(output_tsv, delimiter='\t')
            csvwriter.writerow(['protein', 'GO', 'EC'])
            for protein in protein_annotations:
                csvwriter.writerow([protein, ','.join(list(protein_annotations[protein][1])), ','.join(list(protein_annotations[protein][2]))])

        # Create PathoLogic file and folder for each input.
        pathologic_folder = os.path.join(output_folder, 'pathologic')
        if not os.path.exists(pathologic_folder):
            os.mkdir(pathologic_folder)

        pathologic_organism_folder = os.path.join(pathologic_folder, base_filename)
        if not os.path.exists(pathologic_organism_folder):
            os.mkdir(pathologic_organism_folder)
        pathologic_file = os.path.join(pathologic_organism_folder, base_filename+'.pf')
        create_pathologic(base_filename, protein_annotations, set_proteins, pathologic_file)

    # Create mpwt taxon ID file.
    clustering_taxon_id_file = os.path.join(input_folder, 'proteome_cluster_tax_id.tsv')
    clustering_taxon_id = {}
    with open(clustering_taxon_id_file, 'r') as input_taxon_id_file:
        taxon_id_csvreader = csv.reader(input_taxon_id_file, delimiter='\t')
        next(taxon_id_csvreader)
        for line in taxon_id_csvreader:
            clustering_taxon_id[line[0]] = line[2]

    pathologic_taxon_id_file = os.path.join(pathologic_folder, 'taxon_id.tsv')
    with open(pathologic_taxon_id_file, 'w') as taxon_id_file:
        taxon_id_csvwriter = csv.writer(taxon_id_file, delimiter='\t')
        taxon_id_csvwriter.writerow(['species', 'taxon_id'])
        for species in clustering_taxon_id:
            taxon_id_csvwriter.writerow([species, clustering_taxon_id[species]])
