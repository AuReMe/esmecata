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