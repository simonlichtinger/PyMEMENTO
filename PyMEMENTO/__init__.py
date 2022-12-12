"""PyMEMENTO"""

import warnings

# Check whether Modeller is installed, if not, throw error and abort
try:
    import modeller
except:
    warnings.warn(
        "Dependency problem: Modeller not found and needs to be installed separately.\
     Try 'conda install -c salilab modeller' or follow the Modeller installation instructions here:\
        https://salilab.org/modeller/download_installation.html "
    )

# Check whether gromacs is found
import gromacs

gromacs.config.setup()
import os

devnull = open(os.devnull, "w")

try:
    gromacs.release()
except:
    warnings.warn(
        "Dependency problem: Gromacs not found and needs to be installed separately.\
        Follow the Gromacs installation instructions here:\
        https://manual.gromacs.org/current/install-guide/index.html or \
             install via 'conda install -c bioconda gromacs'."
    )


# Check whether plumed works
import subprocess

try:
    assert 0 == subprocess.call("plumed", shell=True, stdout=devnull)
except:
    warnings.warn(
        "Plumed not found, will not be able to set up umbrella sampling. Install via \
        the instructions on https://www.plumed.org/doc-v2.8/user-doc/html/index.html \
             or 'conda install -c conda-forge plumed'"
    )

# Add imports here
from .pymemento import MEMENTO

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions


# Ignore annoying MDAnalysis warnings specifically, which should be harmless in this context

import warnings

warnings.filterwarnings("ignore", message="UserWarning: missing dimension")
warnings.filterwarnings("ignore", message="UserWarning: Element information is missing")
warnings.filterwarnings("ignore", message="UserWarning: 1 A^3 CRYST1 record")
warnings.filterwarnings(
    "ignore", message="UserWarning: Unit cell dimensions not found."
)
warnings.filterwarnings("ignore", message="UserWarning: Reader has no dt information")
warnings.filterwarnings("ignore", message="UserWarning: Found missing chainIDs.")
warnings.filterwarnings("ignore", message="UserWarning: Found no information for attr:")
