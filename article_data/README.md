# Input files from the article

## Table of contents
- [Input files from the article](#input-files-from-the-article)
  - [Table of contents](#table-of-contents)
  - [Experiments (preprint of 2025)](#experiments-preprint-of-2025)
    - [Datasets](#datasets)
      - [Manually selected taxa](#manually-selected-taxa)
      - [MGnify validation](#mgnify-validation)
      - [Methanogenic reactor](#methanogenic-reactor)
      - [Algae symbionts](#algae-symbionts)
    - [Reproduce experiments](#reproduce-experiments)
  - [Old Experiments (preprint of 2022)](#old-experiments-preprint-of-2022)
    - [Taxonomic affiliations from Gammaproteobacteria and Alveolata](#taxonomic-affiliations-from-gammaproteobacteria-and-alveolata)
    - [Taxonomic affiliations from 16S and rpoB](#taxonomic-affiliations-from-16s-and-rpob)

## Experiments (preprint of 2025)

### Datasets

#### Manually selected taxa

A folder containing an input file for esmecata containing 13 manually selected taxonomic affiliations from Gammaproteobacteria and Alveolata.

#### MGnify validation

4 input files associated with dataset from MGnify:
- [honeybee gut v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/honeybee-gut/v1.0/)
- [marine v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v1.0/)
- [human oral v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/)
- [pig gut v1.0](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/pig-gut/v1.0/)

For each dataset, it contains all the metadata associated with the genomes of the dataset. The columns `observation_name` and `taxonomic_affiliation` are added so the file can be used by EsMeCaTa.
The column `Completness` is used to filter the genomes (for the article by using the threshold of at least 90% of Completness).

#### Methanogenic reactor

An OTU table from a methanogenic reactor experiment containing the taxonomic assignment and the abundance of each OTU in different samples (corresponding to time of measurements of the biogas reactor).

#### Algae symbionts

Microbiotes from Metagenomes experiment ([Burgunter et al. 2020](https://doi.org/10.3389/fmars.2020.00085) and [KleinJan et al. 2023](https://doi.org/10.1111/mec.16766)). 35 MAGs were selected from KleinJan et al. 2023 with more than 90% completion.

### Reproduce experiments

The experimetns made in the article of EsMeCaTa (a preprint is available at [biorXiv](https://doi.org/10.1101/2022.03.16.484574
)) are available in this [Zenodo archive](https://zenodo.org/records/14502342).

Furthermore, in this archive, there are several precomptued database present (to try) to reproduce these experiments.

To run these experiments, you will have to install esmecata with: `pip install esmecata`

To do so, download one of the precomputed database associated with the dataset you want to run (such as `precomputed_db_honeybee.zip`). Then you can use `esmecata precomputed` command to use this database on an input file. For example:

```
esmecata precomputed -i honeybee_esmecata_metdata.tsv -d precomputed_db_honeybee.zip -o output_folder_honeybee
```

For a better change to reproduce the results, it is recommended to use the NCBI Taxonomy database associacted with the dataset. The NCBI Taxonomy database can be dwonload from the Zenodo archive (file `ncbi_taxonomy_database.zip`). To know which version use, refer to the following table:

| dataset                   | UniProt | NCBI Taxonomy |
|---------------------------|---------|---------------|
| Toy example               | 2023_04 | 09-2023       |
| E. siliculosus microbiota | 2023_02 | 04-2023       |
| Honeybee gut              | 2023_05 | 12-2023       |
| Human Oral                | 2023_05 | 12-2023       |
| Marine                    | 2023_05 | 12-2023       |
| Pig Gut                   | 2023_05 | 12-2023       |
| Methanogenic reactor      | 2024_01 | 01-2024       |

You can update the NCBI Taxonomy database used by ete3 with the following command (here we use as example `taxdmp_2024-01.tar.gz` associated with the methanogenic reactor):

```
python3 -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(); ncbi.update_taxonomy_database('taxdmp_2024-01.tar.gz')"
```

This will create output folder containing the predictions for the associated organisms of the community.

## Old Experiments (preprint of 2022)

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