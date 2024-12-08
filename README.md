# LG-TDGL

![testtesttest](docs/logo.png)

Time-dependent Ginzburg-Landau in Python

#![PyPI](https://img.shields.io/pypi/v/tdgl)
#![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/loganbvh/py-tdgl/lint-and-test.yml?branch=main)
#[![Documentation Status](https://readthedocs.org/projects/py-tdgl/badge/?version=latest)](https://py-tdgl.readthedocs.io/en/latest/?badge=latest)
#[![codecov](https://codecov.io/gh/loganbvh/py-tdgl/branch/main/graph/badge.svg?token=VXdxJKP6Ag)](https://codecov.io/gh/loganbvh/py-tdgl)
#![GitHub](https://img.shields.io/github/license/loganbvh/py-tdgl)
#[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
#[![DOI](https://zenodo.org/badge/535746543.svg)](https://zenodo.org/badge/latestdoi/535746543)


## Motivation

`LG-TDGL` provides a platform to explore the interaction between structured light and superconducting order parameters, enabling the realization of the *quantum printing* effect. The computation of order parameters is based on the framework of the generalized time-depdendent Ginzburg-Landau (TDGL) equation implemented in `pyTDGL` (see reference [1]). <br>
This script, `LG-TDGL`, builds upon the series of works titled ''*Structured light induced vorticity in superconductors*'' (see references [2],[3]). <br>
In this notebook, we showcase an example of a dynamics of superconducting vortices induced by Laguerre-Gaussian beam, and demonstrate the results for the light-imprinted superflow in reference.<br>
<br>

--
References:<br>
[1] This code requires the `pyTDGL` environment. See the documentation: https://py-tdgl.readthedocs.io/en/latest/, and publication DOI: https://doi.org/10.1016/j.cpc.2023.108799. <br>
[2] ''*Structured light and induced vorticity in superconductors I: Linearly polarized light.*'' DOI: https://arxiv.org/abs/2407.15834. <br>
[3] ''*Structured light and induced vorticity in superconductors II: Quantum Print with Laguerre-Gaussian beam.*'' DOI: https://arxiv.org/abs/2412.00935. <br>

## Try `LG-TDGL`

Click the badge below to try `LG-TDGL` interactively online via [Google Colab](https://colab.research.google.com/):

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TienTienYeh/lg-tdgl/blob/main/docs/quickstart.ipynb)


### Install via `pip`

From this [GitHub repository](https://github.com/TienTienYeh/lg-tdgl/):

```bash
pip install git+https://github.com/TienTienYeh/lg-tdgl.git
```

Editable installation:

```bash
git clone https://github.com/TienTienYeh/lg-tdgl.git
cd lg-tdgl
pip install -e ".[dev,docs]"
```
## About `pyTDGL`

### Authors

- Authors of works of *Structured light and induced vorticity in superconductors*: Tien-Tien Yeh, Hennadii Yerzhakov, Logan Bishop-Van Horn, Srinivas Raghu, Alexander Balatsky, 
- Primary author and maintainer of GitHub: [@TienTienYeh](https://github.com/TienTienYeh/).

### Citing `LG-TDGL`

`LG-TDGL` is described in the following papers:

>* ''*Structured light and induced vorticity in superconductors I: Linearly polarized light.*'' DOI: https://arxiv.org/abs/2407.15834. 
>* ''*Structured light and induced vorticity in superconductors II: Quantum Print with Laguerre-Gaussian beam.*'' DOI: https://arxiv.org/abs/2412.00935. 

If you use `LG-TDGL` in your research, please cite the paper linked above.

    % BibTeX citation

    @article{yeh2024structured,
    title={Structured light and induced vorticity in superconductors I: Linearly polarized light},
    author={Yeh, Tien-Tien and Yerzhakov, Hennadii and Horn, Logan Bishop-Van and Raghu, Srinivas and Balatsky, Alexander},
    journal={arXiv preprint arXiv:2407.15834},
    year={2024}
    }
    
    @article{yeh2024structured,
    title={Structured light and induced vorticity in superconductors II: Quantum Print with Laguerre-Gaussian beam},
    author={Yeh, Tien-Tien and Yerzhakov, Hennadii and Horn, Logan Bishop-Van and Raghu, Srinivas and Balatsky, Alexander},
    journal={arXiv preprint arXiv:2412.00935},
    year={2024}
    }


### Acknowledgments

This work is based on the `pyTDGL` developed by Logan Bishop-Van Horn.
The documentation for `pyTDGL` can be found at [py-tdgl.readthedocs.io](https://py-tdgl.readthedocs.io/en/latest/).

(Install `pyTDGL`)

`pyTDGL` requires `python` `3.8`, `3.9`, `3.10`, or `3.11`. We recommend installing `pyTDGL` in a [`conda` environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), e.g.

```bash
conda create --name tdgl python="3.10"
conda activate tdgl
```


