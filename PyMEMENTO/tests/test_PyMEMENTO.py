"""
Unit and regression test for the PyMEMENTO package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import PyMEMENTO


def test_PyMEMENTO_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "PyMEMENTO" in sys.modules
