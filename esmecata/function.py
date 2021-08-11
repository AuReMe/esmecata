import csv
import re
import os
import urllib.parse
import urllib.request

def annotation_sparql_query(missing_proteins):
    from SPARQLWrapper import SPARQLWrapper, JSON

    uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)
    missing_proteins_query = ['uniprotkb:'+prot_id for prot_id in missing_proteins]

    """
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rh: <http://rdf.rhea-db.org/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
                SELECT DISTINCT ?protein ?enzyme WHERE {{
                    # ECO 269 is experimental evidence
                    BIND (<http://purl.obolibrary.org/obo/ECO_0000269> as ?evidence)
                    ?protein a up:Protein .
                    ?protein up:reviewed true.
                    ?protein up:enzyme ?enzyme .
                    [] a rdf:Statement ;
                            rdf:subject ?protein ;
                            rdf:predicate up:enzyme ;
                            rdf:object ?enzyme ;
                            up:attribution ?attribution .
                    ?attribution up:evidence ?evidence .
    VALUES ?protein {{ {0} }} 
    }}""".format()

    sparql.setQuery(missing_proteins_query)

    # Parse output.
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

def query_uniprot_to_retrieve_function(protein_queries, output_file, output_dict):
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': 'ACC+ID',
        'to': 'ACC',
        'format': 'tab',
        'query': protein_queries,
        'columns': 'id,protein names,reviewed,go,ec'
    }
    p = re.compile(r'\[GO:(?P<go>\d{7})\]')
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
        result = p.findall(line[3])
        gos = ['GO:'+ go_term for go_term in result]
        ec_numbers = line[4].split('; ')
        results[line[0]] = [protein_name, review, gos, ec_numbers]
        output_dict.update(results)

    with open(output_file, 'w') as output_tsv:
        csvwriter = csv.writer(output_tsv, delimiter='\t')
        csvwriter.writerow(['protein_id', 'protein_name', 'review', 'gos', 'ecs'])
        for protein in results:
            csvwriter.writerow([protein, results[protein][0], str(results[protein][1]), ','.join(results[protein][2]), ','.join(results[protein][3])])

    return output_dict

def chunks(lst, n):
    """Yield successive n-sized chunks from list.
    Form: https://stackoverflow.com/a/312464
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def annotate_coreproteome(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    proteins = []
    with open(input_folder+'/coreproteome.tsv', 'r') as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        for line in csvreader:
            proteins.extend(line)

    set_proteins = list(set(proteins))

    output_dict = {}
    if len(set_proteins) < 20000:
        protein_queries = ' '.join(set_proteins)

        query_uniprot_to_retrieve_function(protein_queries, output_folder+'/annotation_protein.tsv', output_dict)
    else:
        protein_chunks = chunks(proteins, 20000)
        for chunk in protein_chunks:
            output_dict = query_uniprot_to_retrieve_function(chunk, output_folder+'/annotation_protein.tsv', output_dict)

def create_pathologic(protein_annotations, protein_set):
    with open('test.pf', 'w', encoding='utf-8') as element_file:

        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        element_file.write(';; ' + 'test' + '\n')
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        for protein in protein_annotations:
            if protein in protein_set:
                element_file.write('ID\t' + protein + '\n')
                element_file.write('NAME\t' + protein_annotations[protein][0] + '\n')
                element_file.write('PRODUCT-TYPE\tP' + '\n')
                element_file.write('PRODUCT-ID\tprot ' + protein + '\n')
                for go in protein_annotations[protein][2].split(','):
                    element_file.write('GO\t' + go + '\n')
                for ec in protein_annotations[protein][3].split(','):
                    element_file.write('EC\t' + ec + '\n')
                element_file.write('//\n\n')

def annotate_proteins(input_folder, output_folder, cpu):
    annotate_coreproteome(input_folder, output_folder)