# gene-harmony

[![image](https://img.shields.io/pypi/v/gene-harmony.svg)](https://pypi.python.org/pypi/gene-harmony)
[![image](https://img.shields.io/pypi/l/gene-harmony.svg)](https://pypi.python.org/pypi/gene-harmony)
[![image](https://img.shields.io/pypi/pyversions/gene-harmony.svg)](https://pypi.python.org/pypi/gene-harmony)
[![Actions status](https://github.com/cancervariants/gene-harmony/actions/workflows/checks.yaml/badge.svg)](https://github.com/cancervariants/gene-harmony/actions/checks.yaml)

<!-- description -->
This resource is useful for harmonizing ambiguous gene symbols powered by thousands of annotations characterizing the relationships between genes and their symbols.

Given a gene symbol, this package:
- finds all genes with this symbol as an alias or primary gene symbol
- provides relationships between the gene symbols and their genes, when available
- provides more detailed information when given the relationship type with the gene symbol

Resources used for annotations include:
[Ensembl](https://www.ensembl.org/index.html)
[HGNC](https://www.genenames.org/)
[NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/)
[Saccharomyces Genome Database](http://sgd-archive.yeastgenome.org/)
[FLJ Human cDNA Database](https://flj.lifesciencedb.jp/top/sys_info/02_about_database/accession_no/download_v032.html)
[OMIM](https://omim.org/downloads/)
[Uniprot](https://www.uniprot.org/)
[Mouse Genome Informatics](https://www.informatics.jax.org/)
<!-- /description -->

---

**[Documentation](https://gene-harmony.readthedocs.io/stable/)** · [Installation](https://gene-harmony.readthedocs.io/stable/install.html) · [Usage](https://gene-harmony.readthedocs.io/stable/usage.html) · [API reference](https://gene-harmony.readthedocs.io/stable/reference/index.html)

---

## Installation

Install from [PyPI](https://pypi.org/project/gene-harmony/):

```shell
python3 -m pip install gene-harmony
```

---

## Feedback and contributing

We welcome bug reports, feature requests, and code contributions from users and interested collaborators. The [documentation](https://gene-harmony.readthedocs.io/latest/contributing.html) contains guidance for submitting feedback and contributing new code.
