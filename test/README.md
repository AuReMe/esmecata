# EsMeCaTa test

## Requirements

The test procedure was made with `pytest` and `pytest-mock`, so it requires these python packages to be installed:

```
pip install pytest pytest-mock
```

Furthermore, it requires `esmecata` and its following dependencies:

- [biopython](https://pypi.org/project/biopython/).
- [pandas](https://pypi.org/project/pandas/).
- [requests](https://pypi.org/project/requests/).
- [ete3](https://pypi.org/project/ete3/).
- [SPARQLwrapper](https://pypi.org/project/SPARQLWrapper/).
- [mmseqs2](https://github.com/soedinglab/MMseqs2).
- [datapane](https://github.com/datapane/datapane).
- [plotly](https://github.com/plotly/plotly.py).
- [kaleido](https://github.com/plotly/Kaleido)
- [ontosunburst](https://github.com/AuReMe/Ontology_sunburst).

## Files

This folder contains several files:

- `buchnera_workflow.tsv`: the test data file (that can be given as input to esmecata) used in most tests. It contains a taxonomic affiliation associated with the species `Buchnera aphidicola` (a symbiont of the pea aphid with a reduced genome, very useful to limit the time taken by the tests).
- `Example.tsv`: another example of input file for esmecata.
- `test_*.py`: scripts testing esmecata code.
- `annotation_input`: input folder for test on esmecata annotation part.
- `clustering_input`: input folder for test on esmecata clustering part.
- `annotation_expected`: eggnog-mapper predictions on test data for mocking the call to eggnog-mapper.
- `buchnera_database.zip`: precomputed database of esmecata for `Buchnera aphidicola` (used in `esmecata precomputed`).
- `uniprot_sprot.txt`: Swissprot annotation file for proteins of `Buchnera aphidicola` (used in `annotation_uniprot` with parameter `-annotation-files`).
- `uniprot_trembl.txt`: TrEMBL annotation file for proteins of `Buchnera aphidicola` (used in `annotation_uniprot` with parameter `-annotation-files`).


## Usage

As EsMeCaTa queries UniProt, some of the tests may fail due to UniProt updates.

To avoid potentially testing these failed tests, it is possible to select only the test that relies on already calculated results or uses a mock. This will not test all of esmecata, but will cover an important part.

To assert the tests associated with offline data and mocked tests, you can use the following command:

```
pytest -vv test*.py -k '_offline or _mocked'
```

The `-k` parameter of `pytest` will only select functions with `_offline` or `_mocked` in their names, i.e. functions that do not query UniProt.

And you can do the same to test only the test functions that query UniProt:

```
pytest -vv test*.py -k _online
```

Finally if you want to test all functions, just use:

```
pytest -vv test*.py
```