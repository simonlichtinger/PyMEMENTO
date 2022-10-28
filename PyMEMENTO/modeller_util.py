""" This contains methods to interact with the modeller python package. """

from shutil import move
import os
from modeller import Environ, Alignment, Model
from modeller.automodel import AutoModel
from glob import glob
from os.path import join


def create_ali_file(pdb_file: str, out_path: str, include_residues: list, caps=True):
    """This function creates the *.ali file which is necessary for modelling. In this case,
    the only differnce between knowns and sequence is removing the caps though.

    :param pdb_file: Path to file with morphed coordinates.
    :type pdb_file: str
    :param out_path: Path into which to write the ali file.
    :type out_path: str
    :param include_residues: List of all residues which should be included in the model.
    :type include_residues: int
    """

    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb_file)
    aln.append_model(mdl, align_codes="morph", atom_files=pdb_file)

    # This is needed to avoid including the caps in the knowns, because
    # modeller doesn't recognise them! Will add them back
    start_count = min(include_residues)
    end_count = max(include_residues)

    mdl = Model(
        env, file=pdb_file, model_segment=(f"{start_count}:X", f"{end_count}:X")
    )
    aln.append_model(mdl, align_codes="protein", atom_files=pdb_file)
    aln.align2d()
    aln.write(file=join(out_path, "morph->protein.ali"), alignment_format="PIR")


def run_modeller(path: str, number_of_models: int):
    """This runs modeller on the ali file contained in a specified directory.

    :param path: Directory in which to run modeller.
    :type path: str
    :param number_of_models: How many models to generate.
    :type number_of_models: int
    """

    # Remember current working directory, workaround to get modeller files in correct place
    working_dir = os.getcwd()
    os.chdir(path)

    # Edit the ali file to fix paths given that now we're in the modeller folder
    out_lines = []
    with open("morph->protein.ali", "r") as f:
        data = f.readlines()
        for line in data:
            if "morphing" in line:
                out_lines.append(
                    line.split(":")[0] + ":../../morphing" + line.split("morphing")[-1]
                )
            else:
                out_lines.append(line)
    with open("morph->protein.ali", "w") as f:
        f.writelines(out_lines)

    # Do the actual modelling based on an ali file
    env = Environ()
    a = AutoModel(env, alnfile="morph->protein.ali", knowns="morph", sequence="protein")
    a.starting_model = 1
    a.ending_model = number_of_models
    a.make()

    # Return to previous working directory
    os.chdir(working_dir)