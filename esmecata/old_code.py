from SPARQLWrapper import SPARQLWrapper, JSON

def old_retrieve_proteome():
    tax_ids = []
    uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
    sparql = SPARQLWrapper(uniprot_sparql_endpoint)
    for tax_id in tax_ids:
        # Test with SPARQL query
        uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX busco: <http://busco.ezlab.org/schema#>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

        SELECT DISTINCT ?proteome ?score ?fragmented ?missing ?organism ?completion
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
        }}""".format(tax_id)

        sparql.setQuery(uniprot_sparql_query)
        # Parse output.
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

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