""" This contains methods to interact with the modeller python package. """

from shutil import move
import os

# Modeller imports needs to lie in functions, because otherwise the docs can't build without a modeller license.
# from modeller import Environ, Alignment, Model
# from modeller.automodel import AutoModel


from glob import glob
from os.path import join


def create_ali_file(
    pdb_file: str,
    out_path: str,
    include_residues: list,
    ali_file_name: str = "morph->protein.ali",
    use_all_chains=False,
):
    """Create the *.ali file which is necessary for modelling. In this case,
    the only differnce between knowns and sequence is removing the caps.

    :param pdb_file: Path to file with coordinates.
    :type pdb_file: str
    :param out_path: Path into which to write the ali file.
    :type out_path: str
    :param include_residues: List of all residues which should be included in the model.
    :type include_residues: list<int>
    :param ali_file_name: Name of the ali file to be written, defaults to "morph->protein.ali"
    :type ali_file_name: str, optional
    :param use_all_chains: Whether to use all chains of a mutiple-domain protein. If False, only chain 'X' (as written by MDAnalysis\
        for a protein which didn't have any chain identifiers) is used. If True, use all residues that were present in the original coordinate file. Defaults to False.
    :type ali_file_name: bool, optional
    """

    from modeller import Environ, Alignment, Model
    from modeller.automodel import AutoModel

    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb_file)
    aln.append_model(mdl, align_codes="morph", atom_files=pdb_file)

    # This is needed to avoid including the caps in the knowns, because
    # modeller doesn't recognise them! Will add them back
    start_count = min(include_residues)
    end_count = max(include_residues)

    if use_all_chains:
        mdl = Model(env, file=pdb_file)
    else:
        mdl = Model(
            env, file=pdb_file, model_segment=(f"{start_count}:X", f"{end_count}:X")
        )

    aln.append_model(mdl, align_codes="protein", atom_files=pdb_file)
    aln.align2d()
    aln.write(file=join(out_path, ali_file_name), alignment_format="PIR")


def run_modeller(
    path: str, number_of_models: int, ali_file_name: str = "morph->protein.ali"
):
    """Run modeller on the ali file contained in a specified directory.

    :param path: Directory in which to run modeller.
    :type path: str
    :param number_of_models: How many models to generate.
    :type number_of_models: int
    :param ali_file_name: Name of the ali file to be used, defaults to "morph->protein.ali"
    :type ali_file_name: str, optional
    """

    from modeller import Environ, Alignment, Model
    from modeller.automodel import AutoModel

    # Remember current working directory, workaround to get modeller files in correct place
    working_dir = os.getcwd()
    os.chdir(path)

    # Edit the ali file to fix paths given that now we're in the modeller folder
    out_lines = []
    with open(ali_file_name, "r") as f:
        data = f.readlines()
        for line in data:
            if "morphing" in line:
                out_lines.append(
                    line.split(":")[0] + ":../../morphing" + line.split("morphing")[-1]
                )
            else:
                out_lines.append(line)
    with open(ali_file_name, "w") as f:
        f.writelines(out_lines)

    # Do the actual modelling based on an ali file
    env = Environ()
    a = AutoModel(env, alnfile=ali_file_name, knowns="morph", sequence="protein")
    a.starting_model = 1
    a.ending_model = number_of_models
    a.make()

    # Return to previous working directory
    os.chdir(working_dir)
