import csv

from SPARQLWrapper import SPARQLWrapper, TSV
from io import StringIO

def old_retrieve_proteome(tax_id):
    proteomes = []
    organism_ids = {}
    uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)

    # Test with SPARQL query
    # First FILTER NOT EXISTS to avoid redundant proteomes.
    # Second FILTER NOT EXISTS to avoid excluded proteomes.
    # OPTIONAL allows to retrive reference and non reference proteomes and split them.
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX busco: <http://busco.ezlab.org/schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

    SELECT DISTINCT ?proteome ?score ?fragmented ?missing ?organism ?completion ?type
    WHERE
    {{
        ?proteome rdf:type up:Proteome .
        ?proteome up:organism ?organism .
        ?organism rdfs:subClassOf taxon:{0} .
        ?proteome busco:has_score ?busco .
        ?busco busco:complete ?score .
        ?busco busco:fragmented ?fragmented .
        ?busco busco:missing ?missing .
        ?proteome rdfs:seeAlso ?metadataID .
        ?metadataID rdfs:label ?completion .
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
    csvreader = csv.reader(StringIO(results), delimiter='\t')
    # Avoid header.
    next(csvreader)

    for line in csvreader:
        proteome_id = line[0].split('/')[-1]
        score = float(line[1].split('^^')[0])
        fragmented = float(line[2].split('^^')[0])
        missing = float(line[3].split('^^')[0])
        org_tax_id = line[4].split('/')[-1]
        completness = line[5]
        if line[6] == '':
            reference = False
        else:
            reference = True

        busco_percentage = (score / (score+fragmented+missing)) * 100
        print(proteome_id, busco_percentage, org_tax_id, completness, reference)

old_retrieve_proteome(629)
def old_annotation_sparql_query(missing_proteins):
    """
    Not used and deprecated.
    """

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