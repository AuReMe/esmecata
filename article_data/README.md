# Input files from the article

## Table of contents
- [Input files from the article](#input-files-from-the-article)
  - [Table of contents](#table-of-contents)
  - [Experiments](#experiments)
    - [Manually selected taxa](#manually-selected-taxa)
    - [MGnify validation](#mgnify-validation)
    - [Biogas reactor](#biogas-reactor)
    - [Algae symbionts](#algae-symbionts)
  - [Old Experiments](#old-experiments)
    - [Taxonomic affiliations from Gammaproteobacteria and Alveolata](#taxonomic-affiliations-from-gammaproteobacteria-and-alveolata)
    - [Taxonomic affiliations from 16S and rpoB](#taxonomic-affiliations-from-16s-and-rpob)

## Experiments

### Manually selected taxa

A folder containing an input file for esmecata containing 13 manually selected taxonomic affiliations from Gammaproteobacteria and Alveolata.

### MGnify validation

4 input files associated with dataset from MGnify:
- [honeybee gut v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/honeybee-gut/v1.0/)
- [marine v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v1.0/)
- [human oral v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/)
- [pig gut v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/pig-gut/v1.0/)

It contains the input files for esmecata (`*_esmecata.tsv` files) and metadata files (`*_genomes-all_metadata.tsv` files) indicating for each MAGs different measures (such as the completness and contamination).

### Biogas reactor

An OTU table from a biogas reactor experiment containing the taxonomic assignment and the abundance of each OTU in different samples (corresponding to time of measurements of the biogas reactor).

### Algae symbionts

Microbiotes from Metagenomes experiment ([Burgunter et al. 2020](https://doi.org/10.3389/fmars.2020.00085) and [KleinJan et al. 2023](https://doi.org/10.1111/mec.16766)). 35 MAGs were selected from KleinJan et al. 2023 with more than 90% completion.

## Old Experiments

In the [first preprint](https://www.biorxiv.org/content/10.1101/2022.03.16.484574v1) two experiments were performed:

- one with 13 taxonomic affiliations from Gammaproteobacteria and Alveolata.
- one with taxonomic affiliations from the article of [Ogier et al (2019)](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-019-1546-z) using the supplementary [Additional File 10](https://static-content.springer.com/esm/art%3A10.1186%2Fs12866-019-1546-z/MediaObjects/12866_2019_1546_MOESM10_ESM.xlsx) and the FROGs taxonomic affiliations (sheets `rpoB_FROGS` and `16S_FROGs` and column `blast_taxonomy`).

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