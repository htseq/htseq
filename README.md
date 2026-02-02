![CI](https://github.com/htseq/htseq/actions/workflows/ci.yml/badge.svg?branch=main)
[![Documentation Status](https://readthedocs.org/projects/htseq/badge/?version=latest)](https://htseq.readthedocs.io)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/htseq)](https://pypi.python.org/pypi/htseq)

# HTSeq
**DEVS**: https://github.com/htseq/htseq

**DOCS**: https://htseq.readthedocs.io

**CITATION** (please cite this new paper!): Putri et al. Analysing high-throughput sequencing data in Python with HTSeq 2.0. Bioinformatics, btac166, [https://doi.org/10.1093/bioinformatics/btac166](https://doi.org/10.1093/bioinformatics/btac166) (2022).

A Python library to facilitate programmatic analysis of data
from high-throughput sequencing (HTS) experiments. A popular component of ``HTSeq``
is ``htseq-count``, a script to quantify gene expression in bulk and single-cell RNA-Seq
and similar experiments.

## Requirements

To use ``HTSeq`` you need:

-  ``Python >= 3.10`` (**note**: ``Python 2`` support has been dropped)
-  ``numpy``
-  ``pysam``

To manipulate BigWig files, you also need:

- ``pyBigWig``

To run the ``htseq-qa`` script, you also need:

-  ``matplotlib``

To run ``htseq-count`` and ``htseq-count-barcodes`` with custom output formats for the counts table, you need:

- ``mtx`` file: ``scipy``
- ``h5ad`` file: ``anndata``
- ``loom`` file: ``loompy``

Both **Linux** and **OSX** are supported and binaries are provided on Pypi. We
would like to support **Windows** but currently lack the expertise to do so. If
you would like to take on the Windows release and maintenance, please open an
issue and we'll try to help.

A source package which should not require ``Cython`` nor ``SWIG`` is also
provided on Pypi.

To **develop** `HTSeq` you will **also** need:

-  ``Cython >=0.29.5``
-  ``SWIG >=3.0.8``

## Installation

### PIP

To install directly from PyPI:

```bash
pip install HTSeq
```

To install a specific version:

```bash
pip install 'HTSeq==2.0.0'
```

If this fails, please install all dependencies first:

```bash
pip install matplotlib
pip install Cython
pip install pysam
pip install HTSeq
```

### setup.py (distutils/setuptools)

Install the dependencies with your favourite tool (``pip``, ``conda``,
etc.).

To install ``HTSeq`` itself, run:

```bash
python setup.py build install
```

## Testing
To test locally, run

```bash
./test.sh
```

To test `htseq-count` alone, run it with the `-o` option.

A virtual environment is created in the `.venv` folder and `HTSeq` is installed inside it, including all modules and scripts.

## Authors
- 2021-: Givanna Putri ([ghar1821](https://github.com/ghar1821))
- 2016-: Fabio Zanini ([iosonofabio](https://github.com/iosonofabio))@ https://fabilab.org
- 2010-2015: Simon Anders ([simon-anders](https://github.com/simon-anders)), Wolfgang Huber
