# EsMeCaTa test

## Requirements

The test procedure has been made with `pytest` and `pytest-mock`, so it requries these python packages that can be isntalled with:

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

Due to the fact that EsMeCaTa queries UniProt, several of the tests can fail due to update of UniProt.

To avoid testing potentially these failed tests, it is possible to select only the test relying on already computed results or using mock. This will not test all of esmecata but cover some major part.

To assert the tests associated with offline data and mocked tests, you can use the following command:

```
pytest -vv test*.py -k '_offline or _mocked'
```

Using the `-k` parameter of `pytest` will select only functions containing `_offline` in their names which are functions that do not query UniProt.

And you can do the same to test only the test functions that query UniProt:

```
pytest -vv test*.py -k _online
```

Finally if you want to test all functions, just use:

```
pytest -vv test*.py
```