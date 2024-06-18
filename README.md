# LG-TDGL

![testtesttest](docs/logo.png)

Time-dependent Ginzburg-Landau in Python

![PyPI](https://img.shields.io/pypi/v/tdgl)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/loganbvh/py-tdgl/lint-and-test.yml?branch=main)
[![Documentation Status](https://readthedocs.org/projects/py-tdgl/badge/?version=latest)](https://py-tdgl.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/loganbvh/py-tdgl/branch/main/graph/badge.svg?token=VXdxJKP6Ag)](https://codecov.io/gh/loganbvh/py-tdgl)
![GitHub](https://img.shields.io/github/license/loganbvh/py-tdgl)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/535746543.svg)](https://zenodo.org/badge/latestdoi/535746543)



## Motivation

`LG-TDGL` provides the platform to solve the light-matter interaction between structured light and supercinducting order parameters. The solving order parameters is based on the framework of "generalized time-depdendent Ginzburg-Landau (TDGL) equation" in `pyTDGL`*. <br>
The work "Structured light induced vorticity in superconductors I: Linearly polarized light"** is based on the calculation of `LG-TDGL`. <br>
In this notebook, we provide an example of linearly polarized light induced vortex and demonstrate the results for the light imprinted superflow***.

*Link of `pyTDGL`: https://py-tdgl.readthedocs.io/en/latest/ <br>
  This code has to be executed under the environment of `pyTDGL`.
**DOI:xxxxxxxxxxxxxxxxxxxxx. <br>
***This figure is corresponds to Fig. 5 in reference "xxxxxxxxxxxxxxx".


## Learn `pyTDGL`

The documentation for `pyTDGL` can be found at [py-tdgl.readthedocs.io](https://py-tdgl.readthedocs.io/en/latest/).

## Try `pyTDGL`

Click the badge below to try `pyTDGL` interactively online via [Google Colab](https://colab.research.google.com/):

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/loganbvh/py-tdgl/blob/main/docs/notebooks/quickstart.ipynb)

## Install `pyTDGL`

`pyTDGL` requires `python` `3.8`, `3.9`, `3.10`, or `3.11`. We recommend installing `pyTDGL` in a [`conda` environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), e.g.

```bash
conda create --name tdgl python="3.10"
conda activate tdgl
```

### Install via `pip`

From  [PyPI](https://pypi.org/project/tdgl/), the Python Package index:
    
```bash
pip install tdgl
```

From this [GitHub repository](https://github.com/loganbvh/py-tdgl/):

```bash
pip install git+https://github.com/loganbvh/py-tdgl.git
```

Editable installation:

```bash
git clone https://github.com/loganbvh/py-tdgl.git
cd py-tdgl
pip install -e ".[dev,docs]"
```
## About `pyTDGL`

### Authors

- Primary author and maintainer: [@TienTienYeh](https://github.com/loganbvh/).

### Citing `pyTDGL`

`LG-TDGL` is described in the following paper:

>* Structured light induced vorticity in superconductors I: Linearly polarized light (unpublished)

If you use `pyTDGL` in your research, please cite the paper linked above.

    % BibTeX citation
    @article{
    unpublished
    }


### Acknowledgments

This work is based on the `pyTDGL` developed by Logan Bishop-Van Horn.
The documentation for `pyTDGL` can be found at [py-tdgl.readthedocs.io](https://py-tdgl.readthedocs.io/en/latest/).

