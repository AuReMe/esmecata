# EsMeCaTa: *Es*timating *Me*tabolic *Ca*pabilties from *Ta*xonomy

EsMeCaTa is a tool to estimate metabolic capabilities from a taxonomy (for example from 16S RNA sequencing). This is useful when we do not have sequencged genomes or proteome associated to our data.

## Requirements

EsMeCaTa needs the following python packages:
 
- [biopython](https://pypi.org/project/biopython/)
- [pandas](https://pypi.org/project/pandas/)
- [requests](https://pypi.org/project/requests/)
- [ete3](https://pypi.org/project/ete3/)

And requires mmseqs2 for protein clustering:

- [mmseqs2](https://github.com/soedinglab/MMseqs2)

Optionnaly, you can use SPARQL queries instead of REST queries. This can be done either with the [Uniprot SPARQL Endpoint](https://sparql.uniprot.org/) (with the option `--sparql uniprot`) or with a Uniprot SPARQL Endpoint that you created locally (it is supposed to work but not tested, only SPARQL queries on the Uniprot SPARQL endpoint have been tested). This option needs a new package:

- [SPARQLwrapper](https://pypi.org/project/SPARQLWrapper/)

**Warning**: using SPARQL queries will lead to less functional annotations and metabolic reactions due to difference in retrieving the annotations between REST query and SPARQL query.

## Input

EsMeCaTa takes as input a tabulated or an excel file with two columns one with the ID corresponding to the taxonomy (for example the OTU ID for 16S RNA sequencing) and a second column with taxonomy separated by ';'. For example:

| observation_name | taxonomy                                                                                                     |
|------------------|--------------------------------------------------------------------------------------------------------------|
| Cluster_1        | Bacteria;Spirochaetes;Spirochaetia;Spirochaetales;Spirochaetaceae;Sphaerochaeta;unknown species              |
| Cluster_2        | Bacteria;Chloroflexi;Anaerolineae;Anaerolineales;Anaerolineaceae;ADurb.Bin120;unknown species                |
| Cluster_3        | Bacteria;Cloacimonetes;Cloacimonadia;Cloacimonadales;Cloacimonadaceae;Candidatus Cloacimonas;unknown species |
| Cluster_4        | Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Rikenellaceae;Rikenellaceae RC9 gut group;unknown species   |
| Cluster_5        | Bacteria;Cloacimonetes;Cloacimonadia;Cloacimonadales;Cloacimonadaceae;W5;unknown species                     |
| Cluster_6        | Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Dysgonomonadaceae;unknown genus;unknown species             |
| Cluster_7        | Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;unknown species                      |

## EsMeCaTa functions

### Retrieve proteomes associated to taxonomy

For each taxon in each taxonomy EsMeCaTa will use ete3 to find the corresponding taxon ID. Then it will search for proteomes assocaited to these taxon ID in the Uniprot Proteomes database.

If there is only 1 proteomes, it will be put aside.

If there is more than 100 proteomes, esmecata will apply a specific method: (1) use the taxon ID associated to each proteomes to create a taxonomy tree with ete3, (2) from the root of the tree (the input taxon), esmecata will find the direct deescendant (sub-taxons), (3) then esmecata will compute the number of proteomes associated to each sub-taxon, (4) the corresponding proportions will be used to select randomly a number of proteomes corresponding to the propotion.
For example: for the taxon Clostridiales, 645 proteomes are found. Using the organism taxon ID associated to the 645 proteomes we found that there is 17 direct sub-taxons. Then for each sub-taxon we compute the percentage of proportion of proteomes given by the sub-taxon to the taxon Clostridiales.
There is 198 proteomes associated to the sub-taxon Clostridiaceae, the percentage will be computed as follow: 198 / 645 = 30% (if a percentage is superior to 1 it will be round down and if the percentage is lower than 1 it will be round up to keep all the low proportion sub-taxons). We will use this 30% to select randomly 30 proteomes amongst the 198 proteomes of Clostridiaceae. This is done for all the other sub-taxons, so we get a number of proteomes around 100 (here it will be 102). Due to the different rounds (up or down) the total number of proteomes will not be equal to exactly 100 but it will be around it.

Then the proteomes found will be downloaded.

### Proteins clustering

For each taxon (a row in the table) EsMeCaTa will use mmseqs2 to cluster the proteins. Then if a cluster contains at least one protein from each proteomes, it will be kept (this threshold can be change using the --threshold option). The representative proteins from the cluster will be used. A fasta file of all the representative proteins will be created for each taxon.

### Retrive proteins annotations

For each of the representative proteins conserved, esmecata will look for the annotation (GO terms, EC number, function, gene name, Interpro) in Uniprot. it will also look for the annotation of the protein of the same cluster than the representative one. And all the anntoations found for a cluster will be propagated to the representative proteins.

Then esmecata will create a tabulated file for each row of the input file and also a folde rcontaing PathoLogic file that can be used as input for Pathway Tools.