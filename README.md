PyMEMENTO
==============================
[//]: # (Badges)
[![pytest](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/run_pytest.yml/badge.svg)](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/run_pytest.yml)
[![codecov](https://codecov.io/github/simonlichtinger/PyMEMENTO/branch/master/graph/badge.svg?token=N663I9V356)](https://codecov.io/github/simonlichtinger/PyMEMENTO)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![CodeQL](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/codeql.yml/badge.svg)](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/codeql.yml)
[![Docs](https://img.shields.io/badge/pymemento.readthedocs.io-blueviolet)](https://pymemento.readthedocs.io)

PyMEMENTO is a simple python implementation of the MEMENTO method for generating paths between known protein conformations as inputs for umbrella sampling (Morphing Endstates by Modelling Ensembles with iNdependent TOpologies). The method is described in detail in the accompanying [publication](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00140). If MEMENTO or this repository have been useful in your work, please cite this paper.

### Installation

PyMEMENTO may be installed as the latest release from PyPI ( ``` pip install PyMEMENTO ``` ) or in the development version from this github repository. Gromacs, modeller and plumed are required but need to be installed separately. Detailed installation instructions can be found in the [documentation](https://pymemento.readthedocs.io/en/latest/installation.html).

### Usage

An exaplanation of the [scientific background](https://pymemento.readthedocs.io/en/latest/background.html) and [tutorials](https://pymemento.readthedocs.io/en/latest/examples.html) for example use cases including all required input files are in the documentation.

### Copyright

Copyright (c) 2022, Simon Lichtinger


#### Acknowledgements

Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
