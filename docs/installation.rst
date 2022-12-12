Installation
============

PyMEMENTO can be installed either from the PyPI repository or from the github source code.
It is recommended that you do this in a fresh conda environment:

.. code-block:: bash

    conda create -n PyMEMENTO python=3.9
    conda activate PyMEMENTO


PyPI installation
-----------------

PyMEMENTO is pip-installable, therefore the latest release can be installed with:

.. code-block:: bash

    pip install PyMEMENTO


Installation from source
------------------------


Alternatively, you may install the latest version using the source code on github.

.. code-block:: bash

    git clone https://github.com/simonlichtinger/PyMEMENTO.git
    cd PyMEMENTO
    pip install .


Additional packages
-------------------

PyMEMENTO depends on other packages that require some manual setup.

* **Modeller** 

    For template-based modelling, PyMEMENTO relies on `modeller <https://salilab.org/modeller/>`_ developed by the Sali lab.
    The package is available on conda (``conda install -c salilab modeller``), but requires a license to run and cannot be bundled into separate projects.
    Please obtain a license `here <https://salilab.org/modeller/registration.html>`_ and follow the installation instructions. PyMEMENTO
    checks for a correct modeller installation and will raise a warning if none is found.

* **Gromacs**

    PyMEMENTO runs solvation of boxes and energy minimizations via the gromacs MD engine. You can install it using the `instructions <https://manual.gromacs.org/current/install-guide/index.html>`_.
    Alternatively, you can use a pre-compiled binary available on conda (``conda install -c bioconda gromacs=2021.3``), which may be slower for running MD but will do the job for PyMEMENTO.
    PyMEMENTO checks for a correct gromacs installation and will raise a warning if none is found.

* **Plumed**

    If you wish to use PyMEMENTO to set up umbrella sampling for your equilibrated boxes, you will need to install plumed. This can be done using their `instructions <https://www.plumed.org/doc-v2.7/user-doc/html/_installation.html>`_
    or the conda package (``conda install -c conda-forge plumed``). A version of gromacs patched with plumed is not required to run this package, but is essential for running the subsequent umbrella sampling.
    PyMEMENTO checks for a correct plumed installation and will raise a warning if none is found.

