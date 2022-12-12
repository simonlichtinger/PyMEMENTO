"""
Unit and regression test for the PyMEMENTO package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

from PyMEMENTO import MEMENTO
import os
import shutil
from pathlib import Path
import MDAnalysis as mda
import py
from os.path import join

DATA_PATH = join(os.path.dirname(os.path.abspath(__file__)), "test_data/")

# TESTING_MDRUN_FLAGS = {"ntmpi": 1, "ntomp": 6, "gpu_id": "0"}
# TESTING_MDRUN_FLAGS = {"ntmpi": 1, "ntomp": 6, "nb": "cpu"}
# TESTING_MDRUN_FLAGS = {"nb": "cpu"}
TESTING_MDRUN_FLAGS = {}


#####   FULL WORKFLOW TESTS   #####

# Switch to these paths and change with tmpdir.as_cwd(): for with tmp2dir.as_cwd(): in order to
# run the tests locally instead of a tmp folder, for diagnostic purposes.
#    tmpdir = py.path.local("debugging_run/")
#    tmpdir.mkdir()


def test_workflow_nolipids(tmpdir):
    """Simplistic end-to-end test, which uses decaalanine to test
    whether the workflow can finish running -- as checked by the
    existance of the expected terminal file."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "decaalanine/helix.gro"),
            join(DATA_PATH, "decaalanine/extended.gro"),
            list(range(1, 11)),
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models()
        model.prepare_boxes()
        model.solvate_boxes(ion_concentration=0.15)

        model.minimize_boxes(TESTING_MDRUN_FLAGS)

        assert os.path.exists("testrun/boxes/sim1/em.gro")

        # copy in example data from running this, for just 100ps per box

        for sim in range(3):
            shutil.copyfile(
                join(DATA_PATH, f"decaalanine/equil_output/sim{sim}/prod_sc.gro"),
                f"testrun/boxes/sim{sim}/prod_sc.gro",
            )
            shutil.copyfile(
                join(DATA_PATH, f"decaalanine/equil_output/sim{sim}/prod_sc.xtc"),
                f"testrun/boxes/sim{sim}/prod_sc.xtc",
            )

        # Set up umbrella sampling folders
        model.prepare_umbrella_sampling(
            join(DATA_PATH, "decaalanine/plumed_distance/plumed_monitor.dat"),
            join(DATA_PATH, "decaalanine/plumed_distance/plumed.dat"),
            smoothen_ladder=0.5,
        )

        assert os.path.exists("testrun/umbrella/metadata.dat")

        # Also test umbrella sampling setup with an rmsd-based plumed file
        model.prepare_umbrella_sampling(
            join(DATA_PATH, "decaalanine/plumed_rmsd/plumed_monitor.dat"),
            join(DATA_PATH, "decaalanine/plumed_rmsd/plumed.dat"),
            extra_plumed_files=[
                join(DATA_PATH, "decaalanine/plumed_rmsd/ref_plumed.pdb")
            ],
            smoothen_ladder=0.5,
            folder_name="umbrella_rmsd",
        )

        assert os.path.exists("testrun/umbrella_rmsd/metadata.dat")


def test_multiprocessing_coverage():
    """pytest-cov 4.0.0 has a bug where it doesn't do multiprocessing coverage properly.
    It works for 3.0.0, which is now a requirement of this package. This test should hit the line in f(x)
    in pymemento.py. If it doesn't, then pytest-cov is still broken."""
    from PyMEMENTO.pymemento import test_multiprocessing_coverage

    test_multiprocessing_coverage()


def test_simplistic_ligand(tmpdir):
    """Use decalanine with an added ligand to test whether those codepaths work."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "decaalanine/helix_with_ligand.gro"),
            join(DATA_PATH, "decaalanine/extended_with_ligand.gro"),
            list(range(1, 11)),
            ligand="resname LIG",
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models()
        model.prepare_boxes(
            template_folder=join(DATA_PATH, "template_decaalanine_with_ligand")
        )
        model.solvate_boxes(ion_concentration=0.15)

        model.minimize_boxes(TESTING_MDRUN_FLAGS)

        assert os.path.exists("testrun/boxes/sim1/em.gro")


def test_simplistic_ligand_peptide(tmpdir):
    """Test that a peptide ligand, which gmx doesn't differentiate
    from the protein in index files, work correctly."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "decaalanine/helix_with_ligand_peptide.gro"),
            join(DATA_PATH, "decaalanine/extended_with_ligand_peptide.gro"),
            list(range(1, 11)),
            ligand="resid 12:13",
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models()
        model.prepare_boxes(
            template_folder=join(DATA_PATH, "template_decaalanine_with_ligand_peptide")
        )
        model.solvate_boxes(ion_concentration=0.15)

        model.minimize_boxes(TESTING_MDRUN_FLAGS)

        assert os.path.exists("testrun/boxes/sim1/em.gro")

        # check that the posre file has the expected length, ie 10 CA atoms
        with open("testrun/boxes/sim1/posre.itp") as f:
            assert len(f.readlines()) == 14


def test_setup_CHARMM(tmpdir):
    """Test that the CHARMM forcefield is handled correctly, using an
    deca-alanine test system."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "decaalanine/helix.gro"),
            join(DATA_PATH, "decaalanine/extended.gro"),
            list(range(1, 11)),
            forcefield="CHARMM36",
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models(cap_type="CHARMM")
        model.prepare_boxes()
        model.solvate_boxes(ion_concentration=0.15)

        model.minimize_boxes(TESTING_MDRUN_FLAGS)

        assert os.path.exists("testrun/boxes/sim1/em.gro")


def test_setup_protonation_states(tmpdir):
    """Test whether the non-standard protonation states are correctly assigned, by using a
    variant of the decaalanine test system with two copies of all protonatable residues."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "protonation_states_decaalanine/helix.gro"),
            join(DATA_PATH, "protonation_states_decaalanine/extended.gro"),
            list(range(1, 11)),
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models(
            his=True,
            his_protonation_states=(1, 2),
            glu=True,
            glu_protonation_states=(0, 1),
            asp=True,
            asp_protonation_states=(0, 1),
        )
        model.prepare_boxes()
        model.solvate_boxes(ion_concentration=0.15)

        # load a solvated - ie fuly processed - frame and check the correct protonation states.
        u = mda.Universe("testrun/boxes/sim1/solvated.gro")
        glu1 = u.select_atoms("resid 1")
        glu2 = u.select_atoms("resid 2")
        his1 = u.select_atoms("resid 3")
        his2 = u.select_atoms("resid 4")
        asp1 = u.select_atoms("resid 5")
        asp2 = u.select_atoms("resid 6")
        lengths = [len(glu1), len(glu2), len(his1), len(his2), len(asp1), len(asp2)]
        assert lengths == [15, 16, 17, 18, 12, 13]


def test_setup_protonation_states_CHARMM(tmpdir):
    """Test whether the non-standard protonation states are correctly assigned, by using a
    variant of the decaalanine test system with two copies of all protonatable residues,
    within the CHARMM forcefield."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "protonation_states_decaalanine/helix.gro"),
            join(DATA_PATH, "protonation_states_decaalanine/extended.gro"),
            list(range(1, 11)),
            forcefield="CHARMM36",
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models(
            his=True,
            his_protonation_states=(1, 2),
            glu=True,
            glu_protonation_states=(0, 1),
            asp=True,
            asp_protonation_states=(0, 1),
            cap_type="CHARMM",
        )
        model.prepare_boxes()
        model.solvate_boxes(ion_concentration=0.15)

        # load a solvated - ie fuly processed - frame and check the correct protonation states.
        u = mda.Universe("testrun/boxes/sim1/solvated.gro")
        glu1 = u.select_atoms("resid 1")
        glu2 = u.select_atoms("resid 2")
        his1 = u.select_atoms("resid 3")
        his2 = u.select_atoms("resid 4")
        asp1 = u.select_atoms("resid 5")
        asp2 = u.select_atoms("resid 6")
        lengths = [len(glu1), len(glu2), len(his1), len(his2), len(asp1), len(asp2)]
        assert lengths == [15, 16, 17, 18, 12, 13]


def test_multichain(tmpdir):
    """Test whether a simplistic multichain model is handled correctly."""
    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "multichain/multi.gro"),
            join(DATA_PATH, "multichain/multi.gro"),
            [1, 2, 1, 2],
            multiple_chains=["A", "A", "B", "B"],
        )
        # Perform morphing and modelling
        model.morph(3)
        model.make_models(2, include_residues=list(range(1, 11)))
        model.find_best_path()
        model.process_models(caps=False)
        model.prepare_boxes(template_folder=join(DATA_PATH, "template_multichain"))
        model.solvate_boxes(ion_concentration=0.15)

        model.minimize_boxes(TESTING_MDRUN_FLAGS)

        assert os.path.exists("testrun/boxes/sim1/em.gro")


@pytest.mark.slow
def test_lipid_setup(tmpdir):
    """Test whether a realisitc lipid membrane can be handled, using the LEUT example system.
    Perform only one EM step, and only check that the embedding EM ran correctly, to save time. This is
    a slow test and will be skipped unless --runslow is passed to pytest."""

    with tmpdir.as_cwd():
        model = MEMENTO(
            "testrun/",
            join(DATA_PATH, "LEUT/occ.gro"),
            join(DATA_PATH, "LEUT/if.gro"),
            list(range(11, 508)),
            forcefield="CHARMM36",
            lipid="resname POPE",
        )
        # Perform morphing and modelling
        model.morph(2)
        model.make_models(2, include_residues=list(range(1, 498)))
        model.find_best_path(poolsize=1, mc_replicates=1)
        model.process_models(
            caps=True,
            cap_type="CHARMM",
            his=True,
            his_protonation_states=[0, 0, 0, 0],
        )
        model.prepare_boxes(
            template_folder=join(DATA_PATH, "template_LEUT_if"),
            mdrun_flags=TESTING_MDRUN_FLAGS,
            embedding_steps=1,
        )

        assert os.path.exists("testrun/boxes/sim1/embed_1.gro")
