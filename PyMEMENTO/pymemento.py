""" This contains the core interactive functionality of the MEMENTO class. """

from multiprocessing import pool
from multiprocessing.sharedctypes import Value
import os
import py
import shutil
import multiprocessing
from matplotlib.pyplot import box
import numpy as np
from glob import glob
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis import transformations as trans
from .pdb_util import cap_termini, fix_residue_numbers, sed
from .modeller_util import create_ali_file, run_modeller
from .mc_path_sampling import MCPathSampler
from .gmx_util import (
    embed_in_lipids,
    pdb2gmx,
    solvate_ions,
    get_cubic_boxsize,
    edit_last_line,
    create_box,
    minimize,
    generate_posre,
    embed_in_lipids,
    generate_membrane_index,
)
from .plumed_util import get_monitor_value_from_xtc
from os.path import join
from pathlib import Path


class MEMENTO:
    """The MEMENTO class is the core handle for interacting with the modelling system. See the examples section\
        on how to sequentially use its functions."""

    def __init__(
        self,
        working_dir: str,
        start_file: str,
        target_file: str,
        residue_numbers: list,
        forcefield="AMBER14",
        forcefield_paths=None,
        lipid=None,
        ligand=None,
        ligand_type="rigid",
        PLUMED_PATH="plumed",
        multiple_chains=None,
        last_step_performed="",
    ):
        """Constructor for the MEMENTO class, which takes some initial information that every modelling run must provide.

        :param working_dir: Path to the directory in which all intermediate files should be written.
        :type start_file: str
        :param start_file: Path to the file containing the morphing starting structure.
        :type start_file: str
        :param end_file: Path to the file containing the morphing target structure.
        :type end_file: str
        :param residue_numbers:  A list of residue numbers which the protein residues should have after processing\
            (use this to take care of starting points other than 1, for example, and multichain overlapping residue numbers).
        :type residue_numbers: list
        :param forcefield: Which forcefield to use. 'AMBER14' will default to integrated amber14.sb, 'CHARMM36' to integrated charmm36-jul2021, 'AMBER_SLIPIDS' to integrated amber14.sb + slipids. \
        'Other' will read the forcefield_paths keyword argument. Defaults to 'AMBER14'.
        :type forcefield: str, optional
        :param forcefield_paths: Provice paths to a number of forcefield folders needed for running simulations, if not 'AMBER14' or 'CHARMM36', defaults to None.
        :type forcefield_paths: list, optional
        :param lipid: Specify the lipid selection group to be employed, None means no lipid, defaults to None
        :type lipid: str, optional
        :param ligand: Ligand selectio string for ligand in the binding pocket, which is to be moved during morphing, defaults to None
        :type ligand: str, optional
        :param ligand_type: Determines whether the ligand is treated as 'rigid' and only the COM is moved, or whether it should be morphed as 'single' atoms, defaults to 'rigid'
        :type ligand_type: str, optional
        :param PLUMED_PATH: Path to the plumed executable, defaults to 'plumed'
        :type PLUMED_PATH: str, optional
        :param multiple_chains: If the protein contains multiple chains, pass a list of lists of chain identifiers, \
            eg ["A", "A", "A", "B", "B", "B", "C", "C", "C"] for a trimer which has residue numbers 1-3 in chain A, 4-6 in chain B and 7-9 in chain C, defaults to None
        :type multiple_chains: list<str>, optional
        :param last_step_performed: If a previous run has already done some of the steps in the appropriate folder strucutre, pass the last step that was done\
            with this keyword argument to override order consistency checking. Can be 'modelling', 'pathfinding', 'processing', 'boxpreparation' or 'solvation'. Defaults to ''
        :type last_step_performed: str, optional
        """
        self.PLUMED_PATH = PLUMED_PATH

        # convert starting and end universes paths into absolute paths, as we'll be changing directories

        self.universe_start = mda.Universe(Path(start_file).absolute())
        self.universe_target = mda.Universe(Path(target_file).absolute())
        self.working_dir = working_dir
        self.ligand = ligand
        self.ligand_type = ligand_type
        self.lipid = lipid

        # Save the residue numbers which will be needed later.
        self.residue_numbers = residue_numbers

        # Assign forcefield folder
        if forcefield == "AMBER14":
            self.forcefield_paths = [
                join(os.path.dirname(os.path.abspath(__file__)), "data/amber14sb.ff")
            ]
        elif forcefield == "AMBER_SLIPIDS":
            self.forcefield_paths = [
                join(os.path.dirname(os.path.abspath(__file__)), "data/amber14sb.ff"),
                join(os.path.dirname(os.path.abspath(__file__)), "data/forcefield.ff"),
            ]
        elif forcefield == "CHARMM36":
            self.forcefield_paths = [
                join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "data/charmm36-jul2021.ff",
                )
            ]
        elif forcefield == "Other":
            self.forcefield_paths = forcefield_paths
        else:
            raise ValueError(f"Forcefield {forcefield} unknown.")

        self.forcefield = forcefield

        # Set up the correct forcefield path

        # We need to remember the box-size if we have lipids, as it's hard to reconstruct otherwise
        if lipid:
            with open(target_file) as f:
                self.boxsize_line = f.readlines()[-1]
        else:
            self.boxsize_line = None

        # if we have multiple chains, record this, and fix chains in the universe

        self.multiple_chains = multiple_chains

        if multiple_chains != None:
            self.universe_start.add_TopologyAttr("chainID")

            for r, res in enumerate(
                self.universe_start.select_atoms("protein").residues
            ):
                for atom in res.atoms:
                    atom.chainID = multiple_chains[r]

            self.universe_target.add_TopologyAttr("chainID")
            for r, res in enumerate(
                self.universe_target.select_atoms("protein").residues
            ):
                for atom in res.atoms:
                    atom.chainID = multiple_chains[r]

        # Create working directory if necessary
        os.makedirs(working_dir, exist_ok=True)

        # A few variables for checking which steps have been performed
        self.morph_done = False
        self.modelling_done = False
        self.pathfinding_done = False
        self.modelprocessing_done = False
        self.boxpreparation_done = False
        self.solvation_done = False

        # Continue from a previous run if appropriate

        if last_step_performed == "modelling":
            self.modelling_done = True
        if last_step_performed == "pathfinding":
            self.pathfinding_done = True
        if last_step_performed == "processing":
            self.modelprocessing_done = True
        if last_step_performed == "boxpreparation":
            self.boxpreparation_done = True
        if last_step_performed == "solvation":
            self.solvation_done = True

        if last_step_performed != "":
            # fetch number of morphing intermediates from folder structure
            self.number_of_intermediates = len(
                glob(join(self.working_dir, "morphing/frame*.pdb"))
            )
            self.number_of_models = len(
                glob(join(self.working_dir, "modeller/morph0/protein.B9999*.pdb"))
            )

    def morph(
        self, number_of_intermediates: int, fitting_selection="protein and name CA"
    ):
        """Calculate linear morphs between start and target structures which the class
        already holds (from constructor), with a given number of intermediates.

        :param number_of_intermediates: How many frames to create.
        :type number_of_intermediates: int
        :param fitting_selection: :class:`MDAnalysis` selection string for fitting, defaults to "protein and name CA"
        :type fitting_selection: str, optional
        """

        if number_of_intermediates < 2:
            raise ValueError("Need at least 2 frames to do morphing.")
        self.number_of_intermediates = number_of_intermediates

        # All the morphing happens in a separate folder
        with py.path.local(self.working_dir).as_cwd():
            local_path = "morphing/"
            os.makedirs(local_path, exist_ok=True)

            # Align the starting structure onto the target structure, and save files
            align.alignto(
                self.universe_start, self.universe_target, select=fitting_selection
            )
            self.universe_start.atoms.write(join(local_path, "start.pdb"))

            # Lipids come from target structure after aligning, if there's any
            if self.lipid:
                self.universe_target.select_atoms("not ( " + self.lipid + " )").write(
                    join(local_path, "target.pdb")
                )
                self.universe_target.select_atoms(self.lipid).write(
                    join(local_path, "lipids.gro")
                )
            else:
                self.universe_target.atoms.write(join(local_path, "target.pdb"))

            # Interpolate coordinates, (1-l)* original + l*(target-original)
            for n, l in enumerate(np.linspace(0, 1, number_of_intermediates)):
                intermediate_universe = self.universe_start.copy()

                # Exclude the ligand from the direct morphing
                if self.ligand == None:
                    intermediate_atoms = intermediate_universe.atoms
                else:
                    intermediate_atoms = intermediate_universe.select_atoms(
                        "not ( " + self.ligand + " )"
                    )

                # Morph all protein atoms directly
                for atom_id in range(len(intermediate_atoms)):
                    intermediate_atoms[atom_id].position *= 1 - l
                    intermediate_atoms[atom_id].position += (
                        l * self.universe_target.atoms[atom_id].position
                    )
                # Write protein atoms only
                intermediate_atoms.write(join(local_path, f"frame{n}"))

                # Process the ligand if we have one
                if self.ligand != None:
                    # Morph the ligand
                    ligand_morph = intermediate_universe.select_atoms(self.ligand)
                    ligand_target = self.universe_target.select_atoms(self.ligand)
                    for atom_id in range(len(ligand_morph)):
                        ligand_morph[atom_id].position *= 1 - l
                        ligand_morph[atom_id].position += (
                            l * ligand_target.atoms[atom_id].position
                        )

                    # Align a whole ligand, to keep structure, because we can't simply rebuild it
                    if self.ligand_type == "rigid":
                        ligand_template = self.universe_start.copy().select_atoms(
                            self.ligand
                        )
                        align.alignto(ligand_template, ligand_morph, select=self.ligand)
                        ligand_template.write(join(local_path, f"ligand{n}"))
                    elif self.ligand_type == "single":
                        ligand_morph.write(join(local_path, f"ligand{n}"))
                    else:
                        raise ValueError(
                            f"Ligand type {self.ligand_type} not supported."
                        )

        # remember what we did
        self.morph_done = True

    def make_models(self, number_of_models: int, include_residues=None, poolsize=12):
        """Use the modeller package to generate fixed models based on the morphs already
        present in the folder structure. Caps are removed at this stage.

        :param number_of_models: How many models to generate.
        :type number_of_models: int
        :param include_residues: List of all the residue_numbers to include, None means take all residues counting from 1, defaults to None
        :type include_residues: list, optional
        :param poolsize: How many multiprocessing tasks to run, defaults to 12
        :type poolsize: int, optional
        """

        if not self.morph_done:
            raise RuntimeError("Need to run morphing before making models.")

        if number_of_models < 2:
            raise ValueError(
                "Need to make at least 2 models per intermediate to have an ensemble."
            )

        self.number_of_models = number_of_models
        # Modelling happens in a separate folder
        with py.path.local(self.working_dir).as_cwd():
            local_path = "modeller/"

            os.makedirs(local_path, exist_ok=True)

            # Find appropriate residues to use
            if include_residues == None:
                include_residues = list(range(1, len(self.universe_start.residues) + 1))

            # Prepare folder structure
            frame_paths = []
            for n in range(self.number_of_intermediates):
                frame_path = join(local_path, f"morph{n}/")
                frame_paths.append(frame_path)
                os.makedirs(frame_path, exist_ok=True)

            # Prepare one ali file, copy for each intermediate

            if self.multiple_chains != None:
                use_all_chains = True
            else:
                use_all_chains = False

            create_ali_file(
                f"morphing/frame0.pdb",
                join(local_path, "morph0/"),
                include_residues,
                use_all_chains=use_all_chains,
            )
            for i in range(1, self.number_of_intermediates):
                shutil.copyfile(
                    join(local_path, f"morph0/morph->protein.ali"),
                    join(local_path, f"morph{i}/morph->protein.ali"),
                )
                sed(
                    join(local_path, f"morph{i}/morph->protein.ali"),
                    "frame0.pdb",
                    f"frame{i}.pdb",
                )

            # Run multiprocessed modeller for each intermediate
            # The normal with statement doesn't work with pytest-cov
            pool = multiprocessing.Pool(poolsize)
            pool.starmap(
                run_modeller, [(frame, number_of_models) for frame in frame_paths]
            )
            pool.close()
            pool.join()

        # remember what we did
        self.modelling_done = True

    def find_best_path(
        self, mc_starting_temp=50, mc_steps=10000, mc_replicates=12, poolsize=1
    ):
        """Call up an :class:`.MCPathSampler` instance to construct optimum RMSD paths
        through the space of the generated models.

        :param mc_starting_temp: 'Temperature' at which to start simulated annealing, defaults to 50
        :type mc_starting_temp: int, optional
        :param mc_steps: Number of MC steps to perform, defaults to 10000
        :type mc_steps: int, optional
        :param mc_replicates: Number of MC replicates to perform, defaults to 12
        :type mc_replicates: int, optional
        :param poolsize: Number of processes to spawn to speed up calculations, defaults to 1. **Note:** \
            With large universes this can cause memory issues because of an apparent bug in how MDAnalysis handles many open
            files in multiple processes, hence the default is 1.
        :type poolsize: int, optional
        """

        if not self.modelling_done:
            raise RuntimeError("Need to run modelling before finding best paths.")

        with py.path.local(self.working_dir).as_cwd():
            sampler = MCPathSampler(
                "modeller/",
                self.number_of_intermediates,
                self.number_of_models,
                mc_starting_temp,
                mc_steps,
                mc_replicates,
            )
            sampler.sample_path(poolsize=poolsize)

        # remember that we did path searching
        self.pathfinding_done = True

    def process_models(
        self,
        caps=True,
        cap_type="AMBER",
        proline_n_term=False,
        his=False,
        his_protonation_states=[],
        glu=False,
        glu_protonation_states=[],
        asp=False,
        asp_protonation_states=[],
    ):
        """Perform several processing steps on the identified minimal RMSD path:
        fixing residue numbers, aligning with starting structure, capping termini (optional) and adding H atoms
        via the use of 'gmx pdb2gmx'.

        :param caps: Whether or not to add ACE and NME caps, defaults to True. Currently not supported with multichain proteins.
        :type caps: bool, optional
        :param cap_type: Caps can be of the 'AMBER' or 'CHARMM' forcefield types, defaults to 'AMBER'
        :type caps_type: str, optional
        :param proline_n_term: Whether or not the N-term should be capped with \
        the dihedral angle preferences of a proline, defaults to False
        :type proline_n_term: bool, optional
        :param his: Whether or not his protonation states should be chosen manually, defaults to False
        :type his: bool, optional
        :param his_protonation_states: Those his protonation states which should be chosen (as integers prompted by pdb2gmx) , defaults to []
        :type his_protonation_states: list, optional
        :param glu: Whether or not glu protonation states should be chosen manually, defaults to False
        :type glu: bool, optional
        :param glu_protonation_states: Those glu protonation states which should be chosen (as integers prompted by pdb2gmx), defaults to []
        :type glu_protonation_states: list, optional
        :param asp: Whether or not asp protonation states should be chosen manually, defaults to False
        :type asp: bool, optional
        :param asp_protonation_states: Those asp protonation states which should be chosen (as integers prompted by pdb2gmx), defaults to []
        :type asp_protonation_states: list, optional
        """

        # Currently multichain is not supported with caps

        if caps and self.multiple_chains:
            raise RuntimeError(
                "Currently, combining caps and multichain proteins isn't supported."
            )

        if not self.pathfinding_done:
            raise RuntimeError(
                "Need to search for the best path before processing models."
            )

        # Do processing in separate folder

        with py.path.local(self.working_dir).as_cwd():
            local_path = "processing/"
            os.makedirs(local_path, exist_ok=True)

            # copy over the best path from before, align and fix residue numbers
            for n in range(self.number_of_intermediates):
                frame_path = join(local_path, f"frame{n}.pdb")
                shutil.copyfile(
                    f"modeller/morph{n}/best.pdb",
                    frame_path,
                )
                # aligning, just in case morph hadn't been called
                unaligned = mda.Universe(frame_path)
                select_string = "protein and name CA"

                # excude the ligand from the target selection, in case it fits the pattern
                select_target = select_string
                if self.ligand != None:
                    select_target = select_string + f" and not { self.ligand }"

                # There is a problem with MDAnalysis: it doesn't recognise the mass of "HID"-CA atoms. So
                # since these should all be CA atoms here,  I'm manually setting masses, so that the alignment
                # will work (since that checks for consistent masses).
                for atom in unaligned.select_atoms(select_string):
                    atom.mass = 12.011
                for atom in self.universe_target.select_atoms(select_target):
                    atom.mass = 12.011

                align.alignto(
                    unaligned,
                    self.universe_target,
                    select={"mobile": select_string, "reference": select_target},
                )
                unaligned.atoms.write(frame_path)

                fix_residue_numbers(
                    frame_path,
                    frame_path,
                    self.residue_numbers,
                )

            file_root = "frame"

            # Cap the protein termini if needed
            if caps:
                # Choose correct reference, do account for different dihedral preferences of a proline cap
                ref_folder = join(os.path.dirname(os.path.abspath(__file__)), "data")
                ref_file = (
                    join(ref_folder, "termini_ref_proline.pdb")
                    if proline_n_term
                    else join(ref_folder, "termini_ref.pdb")
                )
                for n in range(self.number_of_intermediates):
                    cap_termini(
                        join(local_path, file_root + f"{n}.pdb"),
                        join(local_path, file_root + f"capped{n}.pdb"),
                        ref_file,
                        self.residue_numbers[0],
                        self.residue_numbers[-1],
                    )

                    # Rename the cap residues to fit the protein naming convention
                    sed(
                        join(local_path, file_root + f"capped{n}.pdb"),
                        "ACE X   1",
                        "ACE A" + str(self.residue_numbers[0] - 1).rjust(4, " "),
                    )
                    sed(
                        join(local_path, file_root + f"capped{n}.pdb"),
                        "NME X   4",
                        "NME A" + str(self.residue_numbers[-1] + 1).rjust(4, " "),
                    )
                file_root = "framecapped"

            # Copy over the forcefield to our root directory for gmx to use
            for ff in self.forcefield_paths:
                shutil.copytree(ff, ff.split("/")[-1], dirs_exist_ok=True)

            # Now add hydrogen atoms and produce itp files
            for n in range(self.number_of_intermediates):
                pdb2gmx(
                    join(local_path, file_root + str(n) + ".pdb"),
                    local_path,
                    his=his,
                    his_protonation_states=his_protonation_states,
                    glu=glu,
                    glu_protonation_states=glu_protonation_states,
                    asp=asp,
                    asp_protonation_states=asp_protonation_states,
                    cap_type=cap_type,
                    multi_chain=(self.multiple_chains != None),
                )
                shutil.move(
                    join(local_path, file_root + str(n) + ".gro"),
                    join(local_path, f"processed{n}.gro"),
                )

                # see whether we'll have one protein itp file or several
                if self.multiple_chains:
                    unique_chain_ids = set(self.multiple_chains)
                    for chain_id in unique_chain_ids:
                        shutil.move(
                            join(
                                local_path,
                                file_root + str(n) + f"_Protein_chain_{chain_id}.itp",
                            ),
                            join(local_path, f"processed{n}_chain{chain_id}.itp"),
                        )
                else:
                    shutil.move(
                        join(local_path, file_root + str(n) + ".itp"),
                        join(local_path, f"processed{n}.itp"),
                    )

            # Remove the forcefield again from our root directory
            for ff in self.forcefield_paths:
                shutil.rmtree(ff.split("/")[-1])

        # remember what we did
        self.modelprocessing_done = True

    def prepare_boxes(
        self,
        template_folder: str = None,
        mdrun_flags={},
        embedding_starting_scale=1.15,
        embedding_end_scale=0.95,
        embedding_steps=5,
    ):
        """Aggregate all required files for starting to run and preprocess them in gromacs. Combine with ligand and embed in lipids 
        (by a procedure akin to inflate_gro) if appopriate.

        :param template_folder: Specify a template folder containing mdp inputs for minim, nvt, npt, prod_sc, and prod;\
        a template topology as well as a submission script. Defaults to None (which sources a folder form the data \
        suitable for a basic protein in water)
        :type template_folder: str, optional
        :param mdrun_flags: Extra flags to pass to gmx mdrun, related to gpu and cores etc (eg. {'ntomp': 6}), defaults to {}
        :type mdrun_flags: dict, optional
        :param starting_scale: Parameter for lipid embedding, factor by which the system should initially be stretched, defaults to 1.15
        :type starting_scale: float, optional
        :param starting_scale: Parameter for lipid embedding, factor by which the system should end up smaller if there's no clashes, defaults to 0.95
        :type starting_scale: float, optional
        :param steps: Parameter for lipid embedding, number of steps of energy minimisation, defaults to 5
        :type steps: int, optional
        """

        if not self.modelprocessing_done:
            raise RuntimeError("Need to process models before preparing boxes.")

        # The template folder should include all required files, except for the forcefields
        # which were specified earlier
        if template_folder == None:
            if self.forcefield == "AMBER14":
                template_folder = join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "data/template_nolipid_amber",
                )
            elif self.forcefield == "CHARMM36":
                template_folder = join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "data/template_nolipid_charmm",
                )

        if template_folder == None:
            raise ValueError(
                f"No template folder is included for forcefield {self.forcefield} and none is supplied\
                via the template_folder keyword argument."
            )

        # convert template folder into an absolute path if it isn't already
        template_folder = Path(template_folder).absolute()

        # Set up folder structure
        with py.path.local(self.working_dir).as_cwd():
            local_path = "boxes/"
            os.makedirs(local_path, exist_ok=True)

            for n in range(self.number_of_intermediates):
                shutil.copytree(
                    template_folder, join(local_path, f"sim{n}"), dirs_exist_ok=True
                )
                for ff in self.forcefield_paths:
                    shutil.copytree(
                        ff,
                        join(local_path, f"sim{n}/", ff.split("/")[-1]),
                        dirs_exist_ok=True,
                    )
                shutil.copyfile(
                    f"processing/processed{n}.gro",
                    join(local_path, f"sim{n}/processed.gro"),
                )
                if self.multiple_chains:
                    unique_chain_ids = set(self.multiple_chains)
                    for chain_id in unique_chain_ids:
                        shutil.copyfile(
                            join(f"processing/processed{n}_chain{chain_id}.itp"),
                            join(
                                local_path, f"sim{n}/topol_Protein_chain_{chain_id}.itp"
                            ),
                        )
                else:
                    shutil.copyfile(
                        f"processing/processed{n}.itp",
                        join(local_path, f"sim{n}/protein.itp"),
                    )

            # Add the ligand back in if that is required, asssuming the template topology has it already
            if self.ligand != None:
                for n in range(self.number_of_intermediates):
                    apo_universe = mda.Universe(
                        join(local_path, f"sim{n}/processed.gro")
                    )
                    ligand_universe = mda.Universe(f"morphing/ligand{n}.pdb")
                    holo_universe = mda.Merge(apo_universe.atoms, ligand_universe.atoms)
                    holo_universe.atoms.write(join(local_path, f"sim{n}/processed.gro"))

                # also determine which residues the ligand spans, if any
                self.ligand_residues = set([a.resid for a in ligand_universe.atoms])
            else:
                self.ligand_residues = []

            # if lipids are present, perform an inflate-gro like formalism to insert into available hole
            if self.lipid:
                for n in range(self.number_of_intermediates):
                    embed_in_lipids(
                        join(local_path, f"sim{n}/processed.gro"),
                        "morphing/lipids.gro",
                        join(local_path, f"sim{n}/"),
                        self.lipid,
                        self.boxsize_line,
                        mdrun_flags=mdrun_flags,
                        starting_scale=embedding_starting_scale,
                        end_scale=embedding_end_scale,
                        steps=embedding_steps,
                    )

        # remember what we did
        self.boxpreparation_done = True

    def solvate_boxes(self, ion_concentration=0.15):
        """This solvates the boxes with a consistent number of solvent molecules across the boxes.

        :param ion_concentration: Concentration of salt (NaCl) to use, defaults to 0.15
        :type ion_concentration: float, optional
        """

        if not self.boxpreparation_done:
            raise RuntimeError("Need to prepare boxes before solvating.")

        with py.path.local(self.working_dir).as_cwd():
            local_path = "boxes/"

            # if we don't have lipids, need to make boxes!
            if self.lipid == None:
                # Figure out the largest box size, as this will be used for all windows
                sizes, last_lines = [], []
                for n in range(self.number_of_intermediates):
                    size, last_line = get_cubic_boxsize(
                        join(local_path, f"sim{n}/processed.gro"),
                        join(local_path, f"sim{n}/"),
                        1,
                    )
                    sizes.append(size)
                    last_lines.append(last_line)
                largest_size = max(sizes)
                largest_size_line = last_lines[sizes.index(largest_size)]
                # Edit the box size
                for n in range(self.number_of_intermediates):
                    create_box(join(local_path, f"sim{n}/processed.gro"))
                    edit_last_line(
                        join(local_path, f"sim{n}/processed.gro"), largest_size_line
                    )

            # Solvate the boxes, intially without ions to get the minimal solvent count
            solvent_counts = []
            for n in range(self.number_of_intermediates):
                solvate_ions(
                    join(local_path, f"sim{n}/processed.gro"),
                    join(local_path, f"sim{n}/"),
                    0,
                )
                with open(join(local_path, f"sim{n}/topol.top"), "r") as f:
                    topology = f.readlines()
                    solvent_counts.append(int(topology[-1].split(" ")[-1]))

            # now solvate again, with consistent solvent count
            min_solvent_count = min(solvent_counts)
            for n in range(self.number_of_intermediates):
                solvate_ions(
                    join(local_path, f"sim{n}/processed.gro"),
                    join(local_path, f"sim{n}/"),
                    ion_concentration,
                    maxsol=min_solvent_count,
                )

        # remember what we did
        self.solvation_done = True

    def minimize_boxes(self, mdrun_flags={}, grompp_flags={}):
        """Perform an energy minimizations on the prepared boxes in the folder structure.

        :param mdrun_flags: Extra flags to pass to gmx mdrun, related to gpu and cores etc (eg. {'ntomp': 6}), defaults to {}
        :type mdrun_flags: dict, optional
        :param mdrun_flags: Extra flags to pass to gmx grompp, (eg. {'maxwarn': 1}), defaults to {}
        :type mdrun_flags: dict, optional
        """

        if not self.solvation_done:
            raise RuntimeError("Need to solvate boxes before minimizing.")

        with py.path.local(self.working_dir).as_cwd():
            # in case we are continuing from a previous step, we need to determine ligand residues again
            if self.ligand != None:
                ligand_universe = mda.Universe(f"morphing/ligand0.pdb")
                self.ligand_residues = set([a.resid for a in ligand_universe.atoms])
            else:
                self.ligand_residues = []

            local_path = "boxes/"

            for n in range(self.number_of_intermediates):
                minimize(
                    join(local_path, f"sim{n}/solvated.gro"),
                    join(local_path, f"sim{n}/"),
                    mdrun_flags=mdrun_flags,
                    grompp_flags=grompp_flags,
                )
                # Also generate backbone posre.itp files as we'll need them for equilibration
                generate_posre(
                    join(local_path, f"sim{n}/solvated.gro"),
                    join(local_path, f"sim{n}/"),
                    "C-alpha",
                    multiple_chains=self.multiple_chains,
                    exclude_res=self.ligand_residues,
                    residue_numbers=self.residue_numbers,
                )
                # For lipids, we also need to generate index files, as we'll be using two temperature coupling groups
                if self.lipid:
                    generate_membrane_index(
                        join(local_path, f"sim{n}/em.gro"),
                        join(local_path, f"sim{n}/"),
                    )

    def prepare_umbrella_sampling(
        self,
        plumed_monitor_file: str,
        plumed_umbrella_file: str,
        run_scripts: list = [],
        extra_plumed_files: list = [],
        smoothen_ladder: float = 0,
        folder_name: str = "umbrella/",
    ):
        """Prepare umbrella sampling simulations from equilibrated boxes.
        Appropriate input files for plumed must be provided. Also provides 'metadata.dat' which
        is an input file for the grossfield wham implementation.

        :param plumed_monitor_file: Path to a plumed input file that prints the CV of interest at stride 1 to file COLVAR_MONITOR.
        :type plumed_monitor_file: str
        :param plumed_umbrella_file: Path to a plumed input file that will be used for umbrella sampling. Must contain \
        the string $REPLICAS$, which is replaced with the determined CV values. Should have KAPPA set, which is read by this function. 
        :type plumed_umbrella_file: str
        :param run_scripts: List of file paths that should be copied into the umbrella sampling folder, eg. run scripts etc, defaults to []
        :type run_scripts: list, optional
        :param extra_plumed_files: List of file paths that are required for running plumed driver for monitoring CVs. \
        Also copied into umbrella path, defaults to []
        :type extra_plumed_files: list, optional
        :param smoothen_ladder: Smoothen the CV ladder for umbrella sampling by linear interpolation between the \
        observed CV values from the box equilibrations and an equidistant ladder. Defaults to 0, which means no smoothening.
        :type smoothen_ladder: float, optional
        :param folder_name: Name of the umbrella sampling folder to be created, defaults to "umbrella/"
        :type folder_name: str, optional
        """

        # make absolute paths out of plumed files

        plumed_monitor_file = os.path.abspath(plumed_monitor_file)
        plumed_umbrella_file = os.path.abspath(plumed_umbrella_file)

        # make absolute paths out of run_scripts and extra_plumed_files

        run_scripts = [os.path.abspath(script) for script in run_scripts]
        extra_plumed_files = [os.path.abspath(script) for script in extra_plumed_files]

        with py.path.local(self.working_dir).as_cwd():
            local_path = folder_name
            os.makedirs(local_path, exist_ok=True)

            # copy required temporary plumed files into the pwd.
            for required_file in extra_plumed_files:
                shutil.copy(required_file, "./")

            CV_values = []
            for n in range(self.number_of_intermediates):
                box_folder = f"boxes/sim{n}/"
                umbrella_folder = join(local_path, f"umbrella{n}/")
                os.makedirs(umbrella_folder, exist_ok=True)
                # Extract cenre of the CV for this umbrella window
                CV_values.append(
                    get_monitor_value_from_xtc(
                        join(box_folder, "prod_sc.xtc"),
                        plumed_monitor_file,
                        PLUMED_PATH=self.PLUMED_PATH,
                    )
                )

                # Copy all required simulation files over
                for ff in self.forcefield_paths:
                    shutil.copytree(
                        ff, join(umbrella_folder, ff.split("/")[-1]), dirs_exist_ok=True
                    )
                for subfile in (
                    glob(box_folder + "*.top")
                    + glob(box_folder + "*.itp")
                    + glob(box_folder + "*.mdp")
                ):
                    shutil.copy(subfile, umbrella_folder)
                shutil.copyfile(
                    join(box_folder, "prod_sc.gro"), join(umbrella_folder, "start.gro")
                )

                # For lipids, we also need to copy over the index file
                if self.lipid:
                    shutil.copyfile(
                        join(box_folder, "index.ndx"),
                        join(umbrella_folder, "index.ndx"),
                    )

            # If required, smoothen the CV-distribution.
            # This is fudging it a bit, but required if the chosen CV
            # doesn't correlate too well with the transition mechanism
            # found by modeller -> if this happens, there's danger of
            # hysteresis when smoothening. Consider a different CV!

            if smoothen_ladder != 0:
                uniform_CV_ladder = np.linspace(
                    CV_values[0], CV_values[-1], self.number_of_intermediates
                )

                print("Original CV ladder is: ", CV_values)
                print(f"Now smoothening with force constant {smoothen_ladder} ...")
                for n in range(len(CV_values)):
                    CV_values[n] = (
                        smoothen_ladder * uniform_CV_ladder[n]
                        + (1 - smoothen_ladder) * CV_values[n]
                    )
                print("New CV ladder is: ", CV_values)

            # Now make the plumed files and copy over the submit file
            shutil.copy(plumed_umbrella_file, local_path)
            plumed_copy = join(local_path, plumed_umbrella_file.split("/")[-1])
            sed(plumed_copy, "$REPLICAS$", ",".join([str(val) for val in CV_values]))
            with open(plumed_copy) as f:
                plumed_lines = f.readlines()
            plumed_lines = ["RESTART\n"] + plumed_lines
            with open(join(local_path, "plumed_restart.dat"), "w") as f:
                f.writelines(plumed_lines)
            for run_script in run_scripts + extra_plumed_files:
                shutil.copy(run_script, local_path)

            # Also make a metadata file for the grossfield wham later

            for line in plumed_lines:
                if (
                    "KAPPA" in line and "replicas" in line
                ):  # There may be other restraints active!
                    kappa_kjmol = float(line.split("KAPPA=")[-1].split(" ")[0])
            kappa_kcalmol = kappa_kjmol / 4.184

            metadata = []
            for n in range(self.number_of_intermediates):
                metadata.append(
                    f"convergence/umbrella{n}xyz {CV_values[n]} {kappa_kcalmol}\n"
                )
            with open(join(local_path, "metadata.dat"), "w") as f:
                f.writelines(metadata)


# This is taken directly from the pytest-cov documentation, to show whether multiprocessing is correctly recognised by
# pytest-cov.


def f(x):
    return x * x


def test_multiprocessing_coverage():
    from multiprocessing import Pool

    p = Pool(5)
    p.map(f, [1, 2, 3])
    p.close()
    p.join()
