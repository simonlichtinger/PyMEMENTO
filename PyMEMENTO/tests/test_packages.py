"""
Tests to check package requirements for MEMENTO package.
"""

import sys

import pytest

import PyMEMENTO


def test_PyMEMENTO_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "PyMEMENTO" in sys.modules


def test_modeller_found():
    """Test that modeller is correctly installed."""
    try:
        import modeller
    except:
        raise ModuleNotFoundError(
            "Dependency problem: Modeller not found and needs to be installed separately.\
        Try 'conda install Modeller' or follow the Modeller installation instructions here:\
        https://salilab.org/modeller/download_installation.html "
        )


# TODO: deprecate this once I'm using gromacs wrapper.
def test_gmx_working_shell():
    """Test that gromacs is working via shell."""
    import subprocess

    assert subprocess.call("gmx", shell=True) == 0


def test_gmx_wrapper():
    """Test that gromacs is working via the wrapper."""
    import gromacs

    release = gromacs.release
