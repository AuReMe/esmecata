# Changelog

# EsMeCaTa v0.4.0 (2022-12-10)

WARNING:
* change in intermediary files of `clustering` and `annotation` in order to reduce disk space used by EsMeCaTa and the number of operations performed by the methods.
* annotation with eggnog-mapper is now the default workflow methods of EsMeCaTa. The previous annotation methods with UniProt has been moved to `annotation_uniprot` and `workflow_uniprot`.

## Add

* Add a new step called `check`, which corresponds to the `proteomes` step of EsMeCaTa without downloading the proteomes. This allows to look at the knowledge associated with the dataset in the UniProt databases.


## Fix


## Modify

* Modify the intermerdiary files of `clustering` and `annotation` (made with an idea propsoed by @PaulineGHG). This change replaces tsv files, that were created for each observation names. Now they will be created for each taxon instead. This means that observation names with the same taxon will be associated with the same file. This reduces the redundancy of the file and decreases the number of operations made by EsMeCaTa.
* Modify how the log json files are created so if a run failed, a new log json file is created instead of erasing the previous ones.

# EsMeCaTa v0.3.0 (2023-12-10)

Add a new way to annotate protein clusters using eggnog-mapper. From test on metagenomcis data, it is more accurate than the methods with UniProt.
Also modify the default option of EsMeCaTa for option with better results on tested data (minimal number of proteomes from 1 to 5 and clustering threshold from 0.95 to 0.5).

## Add

* Add a new method to annotate protein clusters using eggnog-mapper: new script `eggnog.py`, new commands `annotation_eggnog` and `workflow_eggnog`.
* Add option to query uniprot dat files during annotation (`--annotation-files`, needs `biopython`>=`1.81`).
* Add an option to use bioservices for annotation queries (`--bioservices`, requires `bioservices`>=`1.11.2`).
* Add more tests for proteomes selection.
* Add an option to update taxonomic affiliations (`--update-affiliations`).
* Show the failedIDs during mapping for annotation.
* Add an option to specify eggnog-mapper tmp fodler (`--eggnog-tmp`). By default, it is in esmecata output folder.
* Add KEGG reaction in annotation_reference file when using eggnog-mapper.
* Add a function to compare Input taxa information to esmecata taxa information (taxa name, taxa ID, taxa rank) + precise OTUs associated. Thanks to @PaulineGHG.

## Fix

* Do not use already annotated proteins when using annotation files.
* Fix issue in esmecata proteomes, not using non-reference proteome.
* Fix issue with missing reference proteome when parsing SPARQL results.
* Fix issue in main with cli.
* Fix an issue with minimal-nb-proteomes and non-reference proteomes.

## Modify

* Modify default options to `--minimal-nb-proteomes` of 5 (from 1 in previous version) and `-t` (clustering threshold from 0.95 to 0.5).
* Modify rank_limit option to make it more understandable. Only taxon ranks inferior or equal to the one given will be kept.
* Remove several output folders (proteomes `result`, clustering `fasta_consensus` and clustering `fasta_representative`) to reduce size of EsMeCaTa results.
* Rename `tmp_proteome` into `proteomes`.
* Remove FROM in SPARQL queries to speed up the queries (could speed up SPARQL queries).
* Change header for annotation files, especially: 'gos', 'ecs', 'interpros', 'rhea_ids' into 'GO', 'EC', 'InterPro', 'Rhea'.
* Add column `cluster_members` in annotation reference file and renamed column `protein` into `protein_cluster`.
* Do not create fasta file when there are no protein clusters.
* Update license year.
* Update esmecata worfklow picture.
* Update the doc of esmecata.

# EsMeCaTa v0.2.12 (2022-12-01)

## Fix

* esmecata clustering should now avoid reclustering already clustered data.

# EsMeCaTa v0.2.11 (2022-11-25)

## Modify

* remove already downloaded proteomes from the list of proteomes to download.

## Fix

* esmecata annotation should now avoid reannotating already annotated data (as it was not correctly done in previous version).

# EsMeCaTa v0.2.10 (2022-10-19)

## Modify

* esmecata annotation should now avoid reannotating already annotated data.

## Fix

* KeyError issues when retrieving annotation for proteins.
* issue with uncorrect review information for proteins.
* issue in tutorial.

# EsMeCaTa v0.2.9 (2022-09-28)

## Fix

* issue with subsampling.
* issue when checking job status.

# EsMeCaTa v0.2.8 (2022-09-22)

## Fix

* issue in Uniprot mapping requests.

# EsMeCaTa v0.2.7 (2022-09-19)

## Fix

* issue in requirements.txt (that induces the failure of the 0.2.6 release).

# EsMeCaTa v0.2.6 (2022-09-19)

## Modify

* use UniProt functions to wait for the query.
* update ete3 taxaonomy database it it is possible.

## Fix

* multiple KeyError issues in annotation.
* issue in SPARQL queriy for version number.
* typos.

# EsMeCaTa v0.2.5 (2022-08-30)

## Add

* more logs for proteomes download and annotation.

## Fix

* issue with UniProt API, using UniProt functions.

# EsMeCaTa v0.2.4 (2022-06-29)

## Modify

* use new UniProt API.
* remove beta option and old UniProt queries.

# EsMeCaTa v0.2.3 (2022-05-25)

## Fix

* issue with new URL for UniProt mapping.
* issue in test with new UniProt release.

# EsMeCaTa v0.2.2 (2022-04-20)

## Fix

* issue with time out in `proteomes.py`.

# EsMeCaTa v0.2.1 (2022-03-17)

## Add

* `--minimal-nb-proteomes`: Choose the minimal number of proteomes to be selected by EsMeCaTa, if a taxon has less proteomes, it will be ignored (issue #11).
* add a jupyter notebook explaining EsMeCaTa method.

## Fix

* issue when parsing SPARQL query result. The filtering of proteomes with BUSCO score in SPARQL query was incorrect.

# EsMeCaTa v0.2.0 (2022-03-15)

Refactor the code in EsMeCaTa to make it clearer.

## Add

* New functions with the refactoring.
* Comments have been added to make the code more understandable.
* Singularity recipe.
* Timeout in query, this should help for issue #7.
* Tests for the new functions.
* Logging to show log and creates a log file at each run with the command line.

## Fix

* Decrease the number of proteins given to Uniprot in a query, to avoid issue #7 (from 20,000 to 15,000).
* Fix issue in sub sampling with incorrect threshold.
* Fix numerous typos.

## Modify

* Update esmecata workflow picture and readme.

# EsMeCaTa v0.1.0 (2022-03-04)

## Add

* EsMeCaTa step `workflow` (performing steps `proteomes`, `clustering` and `annotation`).
* Time measures in metadata jsons.
* Article data.
* Workflow test using conda (for MMseqs2) and GitHub Actions.
* `--beta` should work to be used with new UniProt REST API (still in development).

## Fix

* Issue in statistics computation (wrong number return).

# EsMeCaTa v0.0.2 (2022-01-19)

## Fix

* Fix readme in PyPI.

# EsMeCaTa v0.0.1 (2022-01-19)

First test release for EsMeCaTa.

## Add

* EsMeCaTa steps: proteomes, clustering and annotation.
* Readme.
* Licence.
* Release on PyPI.
* Test (on GitHub Actions).

