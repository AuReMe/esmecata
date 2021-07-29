from SPARQLWrapper import SPARQLWrapper, JSON
tax_ids = []
uniprot_sparql_endpoint = 'https://sparql.uniprot.org/sparql'
sparql = SPARQLWrapper(uniprot_sparql_endpoint)
for tax_id in tax_ids:
    # Test with SPARQL query
    uniprot_sparql_query = """PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX busco: <http://busco.ezlab.org/schema#>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>

    SELECT DISTINCT ?proteome ?score ?fragmented ?missing
    WHERE
    {{
        ?proteome rdf:type up:Proteome.
        ?protoeme up:organism taxon:{0}.
        ?proteome busco:has_score ?busco .
        ?busco busco:complete ?score .
        ?busco busco:fragmented ?fragmented .
        ?busco busco:missing ?missing .

    }}""".format(tax_id)

    sparql.setQuery(uniprot_sparql_query)
    # Parse output.
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()