""" This file contains functions which interface with plumed. """

import subprocess
import os
import pandas as pd
import numpy as np


def get_monitor_value_from_xtc(
    PLUMED_PATH: str, file_to_process: str, plumed_monitor_file: str
):
    """This uses plumed to extract the value of a CV from the last frame of an xtc file.

    :param PLUMED_PATH: Path to the plumed executable.
    :type PLUMED_PATH: str
    :param file_to_process: Path to the xtc file to be analysed.
    :type file_to_process: str
    :param plumed_monitor_file: Path to the plumed input file, should print to COLVAR_MONITOR
    :type plumed_monitor_file: str
    :return: Value of the CV for the xtc file, last frame.
    :rtype: float
    """
    subprocess.call(
        f"{PLUMED_PATH} driver --plumed {plumed_monitor_file} --mf_xtc {file_to_process}",
        shell=True,
    )
    data = pd.read_table(
        "COLVAR_MONITOR", delim_whitespace=True, comment="#", header=None
    )
    value = np.mean(data[1].to_numpy())
    os.remove("COLVAR_MONITOR")
    return float(value)
