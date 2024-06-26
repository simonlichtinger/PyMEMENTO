""" This contains useful functions which are to do with handling pdb files. """

import MDAnalysis as mda
from MDAnalysis.analysis import align


def cap_termini(
    inpath: str, outpath: str, reference: str, first_res: int, last_res: int, cap_chain:str= None
):
    """Patch a protein with termini, as specified in a reference structure.

    :param inpath: Path to the pdb file to be processed.
    :type inpath: str
    :param outpath: Path to the pdb file to be written to output.
    :type outpath: str
    :param reference: Path to the reference pdb file. This contains a 4-amino acid sequence: cap-nterm-cterm-cap.
    :type reference: str
    :param first_res: Number of the first residue in the inpath file (onto which to attach the cap).
    :type first_res: int
    :param last_res: Number of the last residue in the inpath file (onto which to attach the cap).
    :type last_res: int
    :param cap_chain: Chain ID to cap for multichain proteins. If None, the entire protein without regard to chains is called. Defaults to None.
    :type cap_chain: str
    """
    # Load the pdb files for protein and reference termini
    term_ace = mda.Universe(reference)
    term_nme = mda.Universe(reference)
    prot = mda.Universe(inpath)

    # Chain select string if needed, otherwise empty string
    chain_select = f" and chainid {cap_chain}" if cap_chain else ""

    print(f"Processing chain {cap_chain}")

    # Align termini to protein
    align.alignto(
        term_ace,
        prot,
        select=(
            "resid 2 and (name N or name CA or name C)",
            f"resid {first_res} and (name N or name CA or name C)" + chain_select,
        ),
    )

    # get the first atom id of the first residue in the relevant chain
    first_atom_id = min(prot.select_atoms(f"resid {first_res} and (name N or name CA or name C)" + chain_select)[x].id for x in range(3)) - 1

    align.alignto(
        term_nme,
        prot,
        select=(
            "resid 3 and (name N or name CA or name C)",
            f"resid {last_res} and (name N or name CA or name C)" + chain_select,
        ),
    )

    # get the last atom id of the last residue in the relevant chain
    last_atom_id = prot.select_atoms(f"resid {last_res} and (name OXT)" + chain_select)[0].id + 1

    # Construct new molecule

    # this is ugly but MDAnalysis doesns't let me merge empty selections or None
    if len(prot.select_atoms(f"id 0-{first_atom_id}")) == 0 and len(prot.select_atoms(f"id {last_atom_id}-9999999")) > 0: 
        patched_prot = mda.Merge(
            term_ace.select_atoms("resid 1"),
            prot.select_atoms("not name OXT "+chain_select),
            term_nme.select_atoms("resid 4"),
            prot.select_atoms(f"id {last_atom_id}-9999999")
        )
    elif len(prot.select_atoms(f"id {last_atom_id}-9999999")) == 0 and len(prot.select_atoms(f"id 0-{first_atom_id}")) > 0:
        patched_prot = mda.Merge(
            prot.select_atoms(f"id 0-{first_atom_id}"),
            term_ace.select_atoms("resid 1"),
            prot.select_atoms("not name OXT "+chain_select),
            term_nme.select_atoms("resid 4")
        )
    elif len(prot.select_atoms(f"id {last_atom_id}-9999999")) == 0 and len(prot.select_atoms(f"id 0-{first_atom_id}")) == 0:
        patched_prot = mda.Merge(
            term_ace.select_atoms("resid 1"),
            prot.select_atoms("not name OXT "+chain_select),
            term_nme.select_atoms("resid 4"),
        )
    else:
        patched_prot = mda.Merge(
            prot.select_atoms(f"id 0-{first_atom_id}"),
            term_ace.select_atoms("resid 1"),
            prot.select_atoms("not name OXT "+chain_select),
            term_nme.select_atoms("resid 4"),
            prot.select_atoms(f"id {last_atom_id}-9999999")
        )

    # write output
    patched_prot.atoms.write(outpath)


def sed(path, replace, by):
    """Find and replace within a file.

    :param path: Path to the file to be modified.
    :type path: str
    :param replace: String to be replaced.
    :type replace: str
    :param by: String to replace with.
    :type by: str"""
    with open(path, "r") as f:
        lines = f.readlines()
    lines_out = []
    for l in lines:
        lines_out.append(l.replace(replace, str(by)))
    with open(path, "w") as f:
        f.writelines(lines_out)


def fix_residue_numbers(inpath: str, outpath: str, corrected_residue_numbers: list):
    """Modify the residue numbers in a pdb file to reach a target sequence.
    This function will also clean up any caps that might have been propagated through
    MODELLER, because they will not have proper geometry.

    :param inpath: Path to the origin pdb file
    :type inpath: str
    :param outpath: Where to write the modifed pdb file to.
    :type outpath: str
    :param corrected_residue_numbers: List with the desired residue numbers.
    :type corrected_residue_numbers: list<int>
    """
    with open(inpath) as f:
        data = f.readlines()

    dataout = []
    dropped_residue_positions = []

    previous_resnum = -100
    counter = -1
    for line in data:
        # drop caps here
        if "ACE" in line or "NME" in line:
            continue
        if line[:4] == "ATOM":
            res_num = int(line[22:26])
            if res_num != previous_resnum:
                previous_resnum = res_num
                counter += 1

            dataout.append(
                line[:22]
                + str(corrected_residue_numbers[counter]).rjust(4, " ")
                + line[26:]
            )
        else:
            dataout.append(line)

    with open(outpath, "w") as f:
        f.writelines(dataout)

    # In case this was freshly capped, we can now adjust the names of the caps
    sed(
        outpath,
        "ACE X   1",
        "ACE A" + str(corrected_residue_numbers[0] - 1).rjust(4, " "),
    )
    sed(
        outpath,
        "NME X   4",
        "NME X" + str(corrected_residue_numbers[-1] + 1).rjust(4, " "),
    )
