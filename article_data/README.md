# Input files from the article

## Table of contents
- [Input files from the article](#input-files-from-the-article)
  - [Table of contents](#table-of-contents)
  - [Experiments](#experiments)
    - [Taxonomic affiliations from Gammaproteobacteria and Alveolata](#taxonomic-affiliations-from-gammaproteobacteria-and-alveolata)
    - [Taxonomic affiliations from 16S and rpoB](#taxonomic-affiliations-from-16s-and-rpob)

## Experiments

In the article two experiments were performed:

- one with 13 taxonomic affiliations from Gammaproteobacteria and Alveolata.
- one with taxonomic affiliations from the article of [Ogier et al (2019)](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-019-1546-z).

Both experiments were performed on a cluster with 20 CPUs (especially for the MMseqs part).

### Taxonomic affiliations from Gammaproteobacteria and Alveolata

Two input files were used (Gammaproteobacteria.tsv and Alveolata.tsv).

The following commands were used:

```bash

esmecata proteomes -i input_file.tsv -o 0_proteomes --remove-tmp

esmecata clustering -i 0_proteomes -o 1_clustering_t05 -t 0.5 --cpu 20 --remove-tmp

esmecata annotation -i 1_clustering_t05 -o 2_annotation_t05

esmecata clustering -i 0_proteomes -o 1_clustering_t095 -t 0.95 --cpu 20 --remove-tmp

esmecata annotation -i 1_clustering_t095 -o 2_annotation_t095

esmecata clustering -i 0_proteomes -o 1_clustering_t0 -t 0 --cpu 20 --remove-tmp

esmecata annotation -i 1_clustering_t0 -o 2_annotation_t0
```

### Taxonomic affiliations from 16S and rpoB

Two input files were used (ogier_16s_frogs.tsv and ogier_rpoB_frogs.tsv).

```bash
esmecata proteomes -i input_file.tsv -o 00_proteomes --remove-tmp

esmecata clustering -i 00_proteomes -o 1_clustering_t095 -t 0.95 --cpu 20 --remove-tmp

esmecata annotation -i 1_clustering_t095 -o 2_annotation_t095
```