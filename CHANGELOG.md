# Changelog

# EsMeCaTa v0.6.5 (2025-03-31)

## Fix

* Issue with `cellular root` in `esmecata_compression.py`.

## Modify

* Complete the move to ete4.
* Move test data files.

# EsMeCaTa v0.6.4 (2025-03-07)

## Add

* Function `get_domain_or_superkingdom_from_ncbi_tax_database` in `esmecata.utils` to check if domain is called domain or superkingdom in NCBI Taxonomy database.

## Fix

* Taxonomic rank `superkingdom` has been renamed to `domain` in recent version of the NCBI Taxonomy database. Fix several related issues in different parts of EsMeCaTa.
* `esmecata_gseapy` was not working with results from `esmecata precomputed` due to missing `proteome_tax_id.tsv` in `2_annotation` folder.

# EsMeCaTa v0.6.3 (2025-03-07)

## Fix

* KeyError when using precomputed database from EsMeCaTa article (as there is no `KEGG_reaction` in these databases).

# EsMeCaTa v0.6.2 (2025-03-07)

## Modify

* Replace `datapane` by `arakawa` for HTML report creation in `esmecata_report` (issue #20). Make HTMLs standalone (they are a lot bigger but bo not require internet to display).
* Replace `ete3` by `ete4` as ete3 is no longer maintained (issue #18).

# EsMeCaTa v0.6.1 (2025-03-03)

## Fix

* Issue with tax synonyms on the same tax rank when filtering rank in `proteomes`.

## Modify

* Use taxon ID to search for taxon presence in precomputed database.
* Use Trusted publisher to publish on PyPI.

# EsMeCaTa v0.6.0 (2025-01-27)

## Add

* New command `esmecata_create_db` to create database from different output folders of esmecata (`from_runs`).
* Full release of `esmecata precomputed` associated with the first version of [esmecata precomputed database](https://doi.org/10.5281/zenodo.13354073).
* Option threshold (`-t`) to precomputed.
* Add `--gseapyCutOff` option to `gseapy_enrichr`.
* A check after database creation to detect taxon with few predicted proteins compared to higher affiliated taxon. 
* Check the good format of the gzip file.
* Header `KEGG_reaction` in annotation_reference from `annotation_uniprot` to avoid issues with `esmecata_create_db`.

## Fix

* Issue with protein IDs from UniParc during annotation (incorrect split on '|').
* Fix issue in `get_taxon_obs_name` function.
* Issues in test.

## Modify

* Add database version in log.
* Rename `test_workflow.py` into `test_workflow_uniprot.py`, to better reflect what is done.
* Update workflow figure.
* Update readme.
* Update article_data folder and the associated readme.

# EsMeCaTa v0.5.4 (2024-11-06)

## Fix

* Issue with proteomes from UniParc during clustering (incorrect split on '|').
* Issue in test with updated taxonomic group.

# EsMeCaTa v0.5.3 (2024-10-31)

## Fix

* Handle an issue with requests.exceptions.ChunkedEncodingError.

## Modify

* Remove unused header in output file of gseapy.

# EsMeCaTa v0.5.2 (2024-10-21)

## Add

* New way to search proteomes by using UniParc. Some proteomes, when downloaded directly from UniProt are empty. A solution is to search for them in UniParc and retrieved the associated protein sequences.
* New plot in report showing proteomes according to tax_rank.
* Database number version when creating precomputed database.
* The possibility to give a file containing manually selected groups of observation names for `esmecata_gseapy gseapy_enrichr`.
* Tests for `esmecata_gseapy gseapy_enrichr`.

## Fix

* Issue in creating heatmap of proteomes (missing taxon rank) in report creation.
* Issue when creating database: there was a possibility that a taxon without consensus proteomes and associated annotations was kept.

## Modify

* Update parameter description for SPARQL option to indicate the value to query SPARQL UniProt Endpoint.
* Rename `esmecata_gseapy gseapy_taxon` into `esmecata_gseapy gseapy_enrichr` to reflect the changes in the command.
* Modify how `esmecata_gseapy gseapy_enrichr` works by adding a grouping parameters allowing to choose either groups according to taxon_rank or with a file created by the user and containing manually selected groups of observation names.
* Update readme according to the different changes made in this release.

## TODO

* Investigate and solve memory leak when mapping UniParc IDs to UniProt with bioservices.
* Add handling of UniParc IDs with SPARQL queries.

# EsMeCaTa v0.5.1 (2024-10-04)

## Add

* Metadata file for `esmecata_gseapy gseapy_taxon`.

## Fix

* Several issues in `esmecata_gseapy gseapy_taxon`.
* Issues in tests due to UniProt updates.

## Modify

* Update metadata files for annotation and eggnog by adding missing dependencies.

# EsMeCaTa v0.5.0 (2024-10-02)

**WARNING**: Changes in the structure of the python package of EsMeCaTa.
If you have been importing the package in Python, you will need to modify your import.

## Add

* New command `esmecata_report` to create a report from the output folder of EsMeCaTa. Scripts of `esmecata_report` allow to create `html`, `pdf` and `tsv` reports from EsMeCaTa (work of @alimatai and @PaulineGHG). This command has several subcommands:
    * (1) `create_report` to create a report from the output folder of the `esmecata workflow` subcommand (only way to have the complete HTML report).
    * (2) `create_report_proteomes` to create report files from output of `esmecata proteomes` subcommand.
    * (3) `create_report_clustering` to create report files from output of `esmecata clustering` subcommand.
    * (4) `create_report_annotation` to create report files from output of `esmecata annotation` subcommand.
* New command `esmecata_gseapy` to create enrichment analysis of functions predicted by EsMeCaTa according to taxon rank.
* New optional dependencies required for `esmecata_report`: [datapane](https://github.com/datapane/datapane), [plotly](https://github.com/plotly/plotly.py), [kaleido](https://github.com/plotly/Kaleido), [ontosunburst](https://github.com/AuReMe/Ontology_sunburst). As `datapane` is no more maintained, an alternative with [panel](https://github.com/holoviz/panel) is currently developed.
* New optional dependencies required for `esmecata_gseapy`: [gseapy](https://github.com/zqfang/GSEApy) and [orsum](https://github.com/ozanozisik/orsum).
* New file indicating the EC numbers and GO Terms for the different observation name of the dataset (file `function_table.tsv`).
* New subcommand `esmecata precomputed`. This subcommand uses a precomputed database to make predictions from the input file (using EsMeCaTa default parameters). It has been added to avoid creating the same prediction every run and to have a fast way to make predictions with EsMeCaTa. It is necessary to download the precompiled database before using it. At the moment of this release, the database is not available, these scripts are present for testing purposes.
* Prototype for precomputed database creation: several scripts are added in `esmecata/precomputed` folder to create the input and the precomputed database.
* Check that the proteome files are not completely empty, which could cause problems with mmseqs2.
* Tests for precomputed database, report creation, database creation and eggnog annotation. Add mock on sevral functions to perform the test. Required `pytest-mock`.
* Add readme in test folder.

## Fix

* Issue in proteomes SPARQL query (missing `PREFIX`).

## Modify

* Modification of the structure of the EsMeCaTa package, now divided into 4 main folders: (1) `esmecata/core` (for the script previously contained in the EsMeCaTa folder) and used for the workflow, (2) `esmecata/report` to generate a report from the esmecata output folder, (3) `esmecata/gseapy` to perform enrichment analysis on the esmecata output, and (4) `esmecata/precomputed` to create precomputed database (in development).
* Change the name of intermediary files in `clustering` and `annotation` to avoid issues with ambiguous taxon names.
* Modify test according to changes of packaging structure.
* Modify the behaviour of annotation by eggnog-mapper. Now it merges protein sequences from clustering into bigger files (associated with superkingdom). This increases the performance of eggnog-mapper. Modification made with @megyl. Use `--tax_scope` with eggnog-mapper.
* Update article_data folder.
* Update CI tests of github workflow according to the new tests and the new dependencies. 

## Remove

* Remove `esmecata analysis` subcommand as it was not used and not very useful.

# EsMeCaTa v0.4.2 (2024-02-26)

## Fix

* Issue in `pyproject.toml` and how the package is installed with `pip`.

# EsMeCaTa v0.4.1 (2024-02-16)

## Fix

* Some issues on the PyPI page of the package (due to change from `setup.py` to `pyproject.toml`).

# EsMeCaTa v0.4.0 (2024-02-16)

WARNING:
* change in intermediary files of `clustering` and `annotation` in order to reduce disk space used by EsMeCaTa and the number of operations performed by the methods. Instead of assigning one file per observation name, this new version assigns one file per taxon used by EsMeCaTa. This removes a lot of redundant work that slowed EsMeCaTa and could lead to issue.
* annotation with eggnog-mapper is now the default workflow methods of EsMeCaTa. The previous annotation methods with UniProt has been moved to `annotation_uniprot` and `workflow_uniprot`.


## Add

* Add sub-commands `annotation_uniprot` and `workflow_uniprot` to use the old method of protein annotation.
* Add `check` subcommand that performs the first step of EsMeCaTa without downloading the proteomes. This is helpful when you want to have a glimpse on the available knowledge for your dataset.
* Error message if incorrect extension is given as input to esmecata.

## Fix

* Missing import in proteomes.
* Github Actions.

## Modify

* Modify intermediary files to associate them with taxon name selected by EsMeCaTa instead of the observation name (based on an idea of @PaulineGHG). This change replaces tsv files, that were created for each observation names. Now they will be created for each taxon instead. This means that observation names with the same taxon will be associated with the same file. This reduces the redundancy of the file and decreases the number of operations made by EsMeCaTa.
* Modify how the log json files are created so if a run failed, a new log json file is created instead of erasing the previous ones.
* Move from `setup.py` and `setup.cfg` to `pyproject.toml`.
* Update readme and tutorial.
* Update license year.

## Remove

* Remove sub-commands `annotation_eggnog` and `workflow_eggnog` which are now the default sub-commands `annotation` and `workflow`.

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

