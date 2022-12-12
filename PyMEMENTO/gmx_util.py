""" Some useful functions for running gromacs commands to be used by PyMEMENTO. Uses gromacswrapper."""

import shutil
import subprocess
import os
from time import sleep
from glob import glob
from matplotlib.pyplot import sca
import numpy as np
import MDAnalysis as mda
from .pdb_util import sed
import gromacs
import gromacs.cbook
from os.path import join
import MDAnalysis as mda


def crop_top_to_itp(file: str):
    """Crop a .top file to an .itp file, removing the original file.

    :param file: .top tile to be processed.
    :type file: str
    """
    new_file = file.split(".top")[-2] + ".itp"

    with open(file, "r") as f:
        data = f.readlines()
        itp = data[
            data.index("[ moleculetype ]\n") : data.index(
                "; Include Position restraint file\n"
            )
            - 1
        ]
    with open(new_file, "w") as f:
        f.writelines(itp)

    os.remove(file)


def pdb2gmx(
    file_to_process: str,
    folder_path: str,
    ignh=True,
    his=False,
    his_protonation_states=[],
    glu=False,
    glu_protonation_states=[],
    asp=False,
    asp_protonation_states=[],
    cap_type="AMBER",
    multi_chain=False,
):
    """Run gmx pdb2gmx on a given pdb file.

    :param file_to_process: Path to the pdb file.
    :type file_to_process: str
    :param folder_path: Path to the parent folder of the pdb file.
    :type folder_path: str
    :param ignh: Whether or not to pass the "-ignh" flag, ie rebuild hydrogen atoms, defaults to True
    :type ignh: bool, optional
    :param his: Whether or not to manually set his protonation states for consistency, defaults to False
    :type his: bool, optional
    :param his_protonation_states: If his states are set manually, this is a list of the desired ones, defaults to []
    :type his_protonation_states: list, optional
    :param glu: Whether or not to manually set glu protonation states, defaults to False
    :type glu: bool, optional
    :param glu_protonation_states: If glu states are set manually, this is a list of the desired ones, defaults to []
    :type glu_protonation_states: list, optional
    :param asp: Whether or not to manually set asp protonation states, defaults to False
    :type asp: bool, optional
    :param asp_protonation_states: If asp states are set manually, this is a list of the desired ones, defaults to []
    :type asp_protonation_states: list, optional
    :param cap_type: Can be 'AMBER' or 'CHARMM', defaults to 'AMBER'
    :type cap_type: str, optional
    :param multi_chain: Whether or not the protein has multiple chains, determines whether PyMEMENTO tried to \
        edit the generated top file into an itp, or keep the individual chain itp files. Defaults to False.
    :type multi_chain: bool, optional
    """
    # Set correct paths for everything
    root = file_to_process.split(".pdb")[-2]
    coordinate_output_file = root + ".gro"
    topol_output_file = root + ".top"

    pdb2gmx_input = [
        "1",
        "1",
    ]  # This will default to the forcefield present in the folder and TIP3P

    # Gather command-line options for pdb2gmx. Asp before glu before his before term
    pdb2gmx_input += [
        str(a) for a in asp_protonation_states
    ]  # Add asp protonation states if needed
    pdb2gmx_input += [
        str(g) for g in glu_protonation_states
    ]  # Add glu protonation states if needed
    pdb2gmx_input += [
        str(h) for h in his_protonation_states
    ]  # Add his protonation states if needed
    if cap_type == "CHARMM":
        pdb2gmx_input += ["6", "5"]  # choose 'None' caps if charmm

    # Construct options
    options = []
    if cap_type == "CHARMM":
        options.append("ter")
    if ignh:
        options.append("ignh")
    if his:
        options.append("his")
    if glu:
        options.append("glu")
    if asp:
        options.append("asp")

    # Call pdb2gmx
    gromacs.pdb2gmx(
        *options,
        f=file_to_process,
        o=coordinate_output_file,
        p=topol_output_file,
        input=pdb2gmx_input,
    )

    # Process output
    if not multi_chain:
        sed(topol_output_file, "Protein_chain_A", "Protein")
        crop_top_to_itp(topol_output_file)

        # remove temporary files
        os.remove("posre.itp")


def get_cubic_boxsize(file_to_process: str, folder_path: str, spacing: float):
    """Estimate the cubic boxsize including a given spacing from
    the box edge.

    :param file_to_process: Path to the file to be analysed.
    :type file_to_process: str
    :param folder_path: Path to the folder in which the processing takes place.
    :type folder_path: str
    :param spacing: Spacing from the box edge to be achieved in nm.
    :type spacing: float
    :return: Returns a tuble (size, last_line) of the box-size as number and the last line of a respective gro file.
    :rtype: tuple
    """
    output_file = join(folder_path, "temp_boxsize.gro")

    gromacs.editconf(
        f=file_to_process,
        o=output_file,
        c=True,
        d=spacing,
        bt="cubic",
    )

    with open(output_file, "r") as f:
        last_line = f.readlines()[-1][:-2]
    size = float(last_line.split(" ")[-1])
    os.remove(output_file)

    return (size, last_line)


def create_box(file_to_process: str, padding: float = 1.0):
    """Overwrite a coordinate file with one where a cubic box of defined padding around is created.

    :param file_to_process: Path to the coordinate file to process
    :type file_to_process: str
    :param file_to_process (optional): Padding for box generation, defaults to 1.0 nm
    :type file_to_process: float
    """

    gromacs.editconf(
        f=file_to_process,
        o="temp.gro",
        c=True,
        d=padding,
        bt="cubic",
    )
    shutil.move("temp.gro", file_to_process)


def edit_last_line(file_to_process: str, new_last_line: str):
    """Replace the last line of a file with a different last line."""
    with open(file_to_process, "r") as f:
        data = f.readlines()
    data[-1] = new_last_line
    with open(file_to_process, "w") as f:
        f.writelines(data)


def parse_indexkey_to_list(file_to_process: str):
    """This extracts an ordered list of the
    gmx index auto-generated labels for a given coordinate file.

    :param file_to_process: Path to the file to process.
    :type file_to_process: str
    :return: List of index handles which gmx would use for the structure.
    :rtype: list
    """
    # Make a temporary index file
    gromacs.make_ndx(
        f=file_to_process,
        o=join(os.path.dirname(file_to_process), "temp.ndx"),
        input=["", "q"],
    )

    raw_output = gromacs.cbook.get_ndx_groups(
        join(os.path.dirname(file_to_process), "temp.ndx")
    )

    output_list = []
    for item in raw_output:
        output_list.append(item["name"])

    os.remove(join(os.path.dirname(file_to_process), "temp.ndx"))

    return output_list


def solvate_ions(
    file_to_process: str,
    folder_path: str,
    ion_concentration: float,
    maxsol=1000000,
):
    """Solvate a box which has box size but no solvent present.

    :param file_to_process: Path to the file to solvate.
    :type file_to_process: str
    :param folder_path: Path to the folder in which to solvate.
    :type folder_path: str
    :param ion_concentration: Desired ion concentration (NaCl).
    :type ion_concentration: float
    :param maxsol: Number of solvent molecules to maximally insert, defaults to 1000000
    :type maxsol: int, optional
    """
    output_file = join(folder_path, "solvated.gro")
    shutil.copyfile(join(folder_path, "template.top"), join(folder_path, "topol.top"))
    topology_file = join(folder_path, "topol.top")

    gromacs.solvate(
        cp=file_to_process,
        cs="spc216.gro",
        o=output_file,
        p=topology_file,
        maxsol=maxsol,
    )

    # create ions if ion concetration >0 is specified
    if ion_concentration > 0:
        output_file_ions = join(folder_path, "solvated_ions.gro")

        gromacs.grompp(
            f=join(folder_path, "minim.mdp"),
            c=output_file,
            p=topology_file,
            o=join(folder_path, "ions.tpr"),
            maxwarn=2,
        )
        solvent_index = parse_indexkey_to_list(output_file).index("SOL")
        gromacs.genion(
            s=join(folder_path, "ions.tpr"),
            o=output_file_ions,
            p=topology_file,
            neutral=True,
            conc=ion_concentration,
            pname="NA",
            nname="CL",
            input=[str(solvent_index)],
        )
        shutil.move(output_file_ions, output_file)


def minimize(
    file_to_process: str,
    folder_path: str,
    name="em",
    grompp_flags={},
    mdrun_flags={},
):
    """Energy-minimize a box. A minim.mdp file is required in the folder.

    :param file_to_process: Path to the file to minimize.
    :type file_to_process: str
    :param folder_path: Path to the folder in which to minimize.
    :type folder_path: str
    :param name: Name for the minimization, defaults to "em"
    :type name: str, optional
    :param grompp_flags: Extra flags to pass to gmx grompp, eg {'maxwarn': 1}, defaults to {}
    :type grompp_flags: dict, optional
    :param mdrun_flags: Extra flags to pass to gmx mdrun, eg related to gpu and cores (eg. {'ntomp': 6}), defaults to {}
    :type mdrun_flags: dict, optional
    :return: ID of the overlapping atom, if failed. None if this worked.
    :rtype: int
    :raises ValueError: something went wrong with the minimizaion, but not atom overlap!
    """
    topology_file = join(folder_path, "topol.top")

    gromacs.grompp(
        f=join(folder_path, "minim.mdp"),
        c=file_to_process,
        p=topology_file,
        o=join(folder_path, name + ".tpr"),
        **grompp_flags,
    )

    gromacs.mdrun(deffnm=join(folder_path, name), v=True, **mdrun_flags)

    # Test for errors
    with open(join(folder_path, name + ".log")) as f:
        logfile = f.readlines()
    for line in logfile:
        if "inf on atom " in line:
            return int(line.split("inf on atom ")[-1])  # grab atom id
        if "ERROR" in line:
            raise ValueError(
                "something went wrong with the minimizaion, but not atom overlap!"
            )
    return None  # apparently, everything worked


def generate_posre(
    file_to_process: str,
    folder_path: str,
    group: str,
    exclude_res: list = [],
    multiple_chains: list = None,
    residue_numbers: list = None,
):
    """Generate a position restraint posre.itp file.

    :param file_to_process: Path to file to generate restraints from.
    :type file_to_process: str
    :param folder_path: Path to folder in which to place output.
    :type folder_path: str
    :param group: Name of group from index file to use.
    :type group: str
    :param exclude_res: List of residues to exclude from making restraints
    :type exclude_res: list<int>, optional
    :param multiple_chains: In case of a multichain protein, list of chain ID letters for each residue of the complete protein. Defaults to None.
    :type multiple_chains: list<str>, optional
    :param residue_numbers: In case of a multichain protein, list of residue numbers for the complete protein. Defaults to None.
    :type residue_numbers: list<int>, optional
    """

    if multiple_chains and residue_numbers == None:
        raise ValueError("Need to provided residue numbers for multichain proteins.")

    current_index_list = parse_indexkey_to_list(file_to_process)
    group_index = current_index_list.index(group)

    # for single-chain proteins, just take all protein atoms into account for posre generation.
    if multiple_chains == None:
        if len(exclude_res) > 0:
            exclude_string = " | ".join([f"r {res}" for res in exclude_res])

            gromacs.make_ndx(
                f=file_to_process,
                o=join(folder_path, "temp.ndx"),
                input=[
                    exclude_string,
                    f"!{len(current_index_list)} & {group_index}",
                    "q",
                ],
            )
            gromacs.genrestr(
                f=file_to_process,
                o=join(folder_path, "posre.itp"),
                n=join(folder_path, "temp.ndx"),
                input=[str(len(current_index_list) + 1)],
            )
            os.remove(join(folder_path, "temp.ndx"))
        else:
            gromacs.genrestr(
                f=file_to_process,
                o=join(folder_path, "posre.itp"),
                input=[str(group_index)],
            )
    else:
        # for multichain proteins, need to split into individual chains (while reassigning atom indices)
        for chainID in set(multiple_chains):
            # Make a selection string for the residues of a particular chain ID
            include_chainres = []
            for n, res_num in enumerate(residue_numbers):
                if multiple_chains[n] == chainID:
                    include_chainres.append(res_num)
            include_res_string = " | ".join([f"r {res}" for res in include_chainres])

            # Output a temporary coordinate file that only includes one chain

            if len(exclude_res) > 0:
                string_mda_sele = "protein and not (" + " or ".join(
                    [f"resnum {res}" for res in exclude_res] + ")"
                )
            else:
                string_mda_sele = "protein"

            temp_chain_file = join(folder_path, "temp.gro")

            universe_to_process = mda.Universe(file_to_process)
            atomgroups_to_merge = []
            for n, res in enumerate(
                universe_to_process.select_atoms(string_mda_sele).residues
            ):
                if multiple_chains[n] == chainID:
                    atomgroups_to_merge.append(res.atoms)

            universe_onechain = mda.Merge(*atomgroups_to_merge)
            universe_onechain.atoms.write(temp_chain_file)

            # Now restrain the desired atoms, updating the group index, in case that changed

            current_index_list = parse_indexkey_to_list(temp_chain_file)
            group_index = current_index_list.index(group)
            gromacs.genrestr(
                f=temp_chain_file,
                o=join(folder_path, f"posre_Protein_chain_{chainID}.itp"),
                input=[str(group_index)],
            )

            os.remove(temp_chain_file)


def generate_membrane_index(file_to_process: str, folder_path: str):
    """Generate an index file which separates the system in to Water_and_ions and Membrane.

    :param file_to_process: Path to file to generate index from.
    :type file_to_process: str
    :param folder_path: Path to folder in which to place output.
    :type folder_path: str
    """
    index_list = parse_indexkey_to_list(file_to_process)
    water_ions_index = index_list.index("Water_and_ions")

    gromacs.make_ndx(
        f=file_to_process,
        o=join(folder_path, "index.ndx"),
        input=[f"!{water_ions_index}", f"name {water_ions_index+1} Membrane", "q"],
    )


def stretch_relative(
    to_strech: mda.Universe, anchor: mda.Universe, scale_factor: float
):
    """Act on a Universe or AtomGroup to stretch it relative to the center of geometry of another Universe,
    in the x-y plane.

    :param to_strech: Universe or AtomGroup to stretch.
    :type to_strech: MDAnalysis.Universe or MDAnalysis.AtomGroup
    :param anchor: Universe or AtomGroup to use as the center of geometry.
    :type anchor: MDAnalysis.Universe or MDAnalysis.AtomGroup"""
    center = anchor.atoms.center(None)
    for atom in to_strech.atoms:
        atom.position = np.array(
            [
                center[0] + scale_factor * (atom.position[0] - center[0]),
                center[1] + scale_factor * (atom.position[1] - center[1]),
                atom.position[2],
            ]
        )


def stretch_boxsize(boxsize_line: str, scale_factor: float):
    """Scale box-size string in gro formatting by a scale_factor.

    :param boxsize_line: Box-size line in gro file format.
    :type boxsize_line: str
    :param scale_factor: Scale factor to apply to box-size.
    :type scale_factor: float"""
    scale = [scale_factor, scale_factor, 1]  # only in x-y plane!
    dimenions = [
        round(float(x) * scale[n], 5) if n < 3 else round(float(x), 5)
        for n, x in enumerate(boxsize_line.split())
    ]
    return "".join([str(x).rjust(10) for x in dimenions]) + "\n"


def change_boxsize(file_to_process: str, boxsize_line: str):
    """Change the dimensions of a file to values specified by a gro-formatted string.

    :param file_to_process: Path to the coordinate file to process.
    :type file_to_process: str
    :param boxsize_line: box size in format of last line of gro file
    :type boxsize_line: str"""
    dimensions = [float(x) * 10 for x in boxsize_line.split()]  # nm -> A
    u = mda.Universe(file_to_process)

    # Deal with hexagonal boxes, which will have additional box vectors in columns 4 to 6
    if len(dimensions) == 3:
        dimensions_out = dimensions + [90, 90, 90]
    elif (
        dimensions[3] == 0
        and dimensions[4] == 0
        and dimensions[5] > 59
        and dimensions[5] < 61
    ):
        dimensions_out = dimensions[:3] + [90, 90, 60]
    else:
        raise ValueError(
            f"The box vectors {boxsize_line} are not of a supported cubic or hexagonal type."
        )

    u.dimensions = np.array(dimensions_out)
    u.atoms.write(file_to_process)


def embed_in_lipids(
    protein_path: str,
    lipid_path: str,
    folder_path: str,
    lipid_selection: str,
    boxsize_line: str,
    starting_scale=1.15,
    end_scale=0.95,
    steps=5,
    mdrun_flags={},
):
    """Embed a protein in lipids (eg from a related conformational state), which already have the rough protein shape cut out, but there may be overlaps.

    :param protein_path: Path to the protein coordinate file.
    :type protein_path: str
    :param lipid_path: Path to the lipid coordinate file.
    :type lipid_path: str
    :param folder_path: Path to the folder in which to work.
    :type folder_path: str
    :param lipid_selection: MDAnalysis selecttion string for the lipids
    :type lipid_selection: str
    :param boxsize_line: Last line of the original lipid gro file.
    :type boxsize_line: str
    :param starting_scale: Factor by which the system should initially be stretched, defaults to 1.15
    :type starting_scale: float, optional
    :param starting_scale: Factor by which the system should end up smaller if there's no clashes, defaults to 0.95
    :type starting_scale: float, optional
    :param steps: Number of steps of energy minimisation, defaults to 5
    :type steps: int, optional
    :param mdrun_flags: Extra flags to pass to gmx mdrun, related to gpu and cores etc (eg. {'ntomp': 6}), defaults to {}
    :type mdrun_flags: dict, optional
    """
    boxsize = boxsize_line

    # make an intitial embed
    protein_universe = mda.Universe(protein_path)
    lipid_universe = mda.Universe(lipid_path)
    stretch_relative(lipid_universe, protein_universe, starting_scale)
    boxsize = stretch_boxsize(boxsize, starting_scale)
    combined_universe = mda.Merge(protein_universe.atoms, lipid_universe.atoms)
    combined_universe.atoms.write(join(folder_path, "embed_0.gro"))
    change_boxsize(join(folder_path, "embed_0.gro"), boxsize)

    stepwise_compression = np.power(end_scale / starting_scale, 1 / steps)
    for n in range(1, steps + 1):
        # Compress in multiple steps ...
        combined_universe = mda.Universe(join(folder_path, f"embed_{n-1}.gro"))
        protein = combined_universe.select_atoms("not ( " + lipid_selection + " )")
        lipid = combined_universe.select_atoms(lipid_selection)
        stretch_relative(lipid, protein, stepwise_compression)
        boxsize = stretch_boxsize(boxsize, stepwise_compression)
        combined_universe.atoms.write(join(folder_path, f"embed_raw{n}.gro"))
        change_boxsize(join(folder_path, f"embed_raw{n}.gro"), boxsize)

        # ... while energy minimizing
        return_code = -1
        while return_code != None:
            shutil.copyfile(
                join(folder_path, "template.top"), join(folder_path, "topol.top")
            )
            return_code = minimize(
                join(folder_path, f"embed_raw{n}.gro"),
                folder_path,
                name=f"embed_{n}",
                mdrun_flags=mdrun_flags,
                grompp_flags={"maxwarn": 1},
            )
            if return_code != None:
                # There was a problem of overlapping atoms!
                # Need to repeat with small rocking about of atoms,
                # this tends to happen at periodic boundary faces.
                sel = combined_universe.select_atoms("bynum " + str(return_code))
                sel[0].position = sel[0].position + np.array([1, 0, 0])
                combined_universe.atoms.write(folder_path + f"embed_raw{n}.gro")
                change_boxsize(join(folder_path, f"embed_raw{n}.gro"), boxsize)

    # Cleanup and move final file
    shutil.copyfile(join(folder_path, f"embed_{steps}.gro"), protein_path)
