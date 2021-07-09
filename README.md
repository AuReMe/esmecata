# EsMeCaTa: *Es*timating *Me*tabolic *Ca*pabilties from *Ta*xonomy

EsMeCaTa is a tool to estimate metabolic capabilities from a taxonomy (for example from 16S RNA sequencing). This is useful when we do not have sequencged genomes or proteome associated to our data.

## Requirements

EsMeCaTa needs the following ptyhon packages:
 
- biopython
- pandas
- requests
- ete3

And requires mmseqs2 for protein clustering:

- mmseqs2

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

### Retriece proteomes associated to taxonomy

For each taxon in each taxonomy EsMeCaTa will use ete3 to find the corresponding taxon ID. Then it will search for proteomes assocaited to these taxon ID in the Uniprot Proteomes database. The proteomes foudn will be downloaded.

### Pseudo-coreproteome creation

For each taxon (a row in the table) EsMeCaTa will use mmseqs2 to cluster the proteins. Then if a cluster contains at least one protein from each proteomes, it will be kept. The representative proteins from the cluster will be used. A fasta file of all the representative proteins will be created for each taxon.
