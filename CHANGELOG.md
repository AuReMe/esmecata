# Changelog

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

