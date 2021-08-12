import csv
import re
import os
import urllib.parse
import urllib.request

def old_annotation_sparql_query(missing_proteins):
    """
    Not used and deprecated.
    """
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

def query_uniprot_to_retrieve_function(protein_queries, output_dict):
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
    reference_protein = os.path.join(input_folder, 'reference_proteins')
    for input_file in os.listdir(reference_protein):
            base_file = os.path.basename(input_file)
            base_filename = os.path.splitext(base_file)[0]
            proteins = []
            reference_proteins = {}
            reference_protein_list_file = os.path.join(reference_protein, input_file)
            with open(reference_protein_list_file, 'r') as input_file:
                csvreader = csv.reader(input_file, delimiter='\t')
                for line in csvreader:
                    proteins.extend(line)
                    reference_proteins[line[0]] = line[1:]

            set_proteins = list(set(proteins))

            output_dict = {}
            if len(set_proteins) < 20000:
                protein_queries = ' '.join(set_proteins)

                query_uniprot_to_retrieve_function(protein_queries, output_dict)
            else:
                protein_chunks = chunks(proteins, 20000)
                for chunk in protein_chunks:
                    output_dict = query_uniprot_to_retrieve_function(chunk, output_dict)

            annotation_folder = os.path.join(output_folder, 'annotation')
            if not os.path.exists(annotation_folder):
                os.mkdir(annotation_folder)
            annotation_file = os.path.join(annotation_folder, base_filename+'.tsv')
            with open(annotation_file, 'w') as output_tsv:
                csvwriter = csv.writer(output_tsv, delimiter='\t')
                csvwriter.writerow(['protein_id', 'protein_name', 'review', 'gos', 'ecs'])
                for protein in output_dict:
                    csvwriter.writerow([protein, output_dict[protein][0], str(output_dict[protein][1]), ','.join(output_dict[protein][2]), ','.join(output_dict[protein][3])])

            protein_annotations = {}
            for reference_protein in reference_proteins:
                gos = set()
                ecs = set()
                for protein in reference_proteins[reference_protein]:
                    if protein in output_dict:
                        if output_dict[protein][2] != ['']:
                            gos.update(set(output_dict[protein][2]))
                        if output_dict[protein][3] != ['']:
                            ecs.update(set(output_dict[protein][3]))
                protein_annotations[reference_protein] = [output_dict[protein][0], gos, ecs]

            annotation_reference_folder = os.path.join(output_folder, 'annotation_reference')
            if not os.path.exists(annotation_reference_folder):
                os.mkdir(annotation_reference_folder)
            annotation_reference_file = os.path.join(annotation_reference_folder, base_filename+'.tsv')
            with open(annotation_reference_file, 'w') as output_tsv:
                csvwriter = csv.writer(output_tsv, delimiter='\t')
                csvwriter.writerow(['protein', 'GO', 'EC'])
                for protein in protein_annotations:
                    csvwriter.writerow([protein, ','.join(list(protein_annotations[protein][0])), ','.join(list(protein_annotations[protein][1]))])

            """
            with open(output_folder+'/diff_annotation_reference_protein.tsv', 'w') as output_tsv:
                csvwriter = csv.writer(output_tsv, delimiter='\t')
                csvwriter.writerow(['protein', 'GO', 'EC'])
                for protein in protein_annotations:
                    csvwriter.writerow([protein, ','.join(list(protein_annotations[protein][0])), ','.join(list(output_dict[protein][3])), ','.join(list(protein_annotations[protein][1])), ','.join(list(output_dict[protein][4]))])
            """

            pathologic_folder = os.path.join(output_folder, 'pathologic')
            if not os.path.exists(pathologic_folder):
                os.mkdir(pathologic_folder)
            pathologic_file = os.path.join(pathologic_folder, base_filename+'.pf')
            create_pathologic(base_filename, protein_annotations, set_proteins, pathologic_file)

def create_pathologic(base_filename, protein_annotations, protein_set, pathologic_output_file):
    with open(pathologic_output_file, 'w', encoding='utf-8') as element_file:
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        element_file.write(';; ' + base_filename + '\n')
        element_file.write(';;;;;;;;;;;;;;;;;;;;;;;;;\n')
        for protein in protein_annotations:
            if protein in protein_set:
                element_file.write('ID\t' + protein + '\n')
                element_file.write('NAME\t' + protein_annotations[protein][0] + '\n')
                element_file.write('PRODUCT-TYPE\tP' + '\n')
                element_file.write('PRODUCT-ID\tprot ' + protein + '\n')
                for go in protein_annotations[protein][1]:
                    element_file.write('GO\t' + go + '\n')
                for ec in protein_annotations[protein][2]:
                    element_file.write('EC\t' + ec + '\n')
                element_file.write('//\n\n')

def annotate_proteins(input_folder, output_folder, cpu):
    annotate_coreproteome(input_folder, output_folder)