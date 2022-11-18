""" This file contains functions which interface with plumed. """

import subprocess
import os
import pandas as pd
import numpy as np


def get_monitor_value_from_xtc(
    file_to_process: str,
    plumed_monitor_file: str,
    PLUMED_PATH: str = "plumed",
):
    """Use plumed to extract the average value of a CV from an xtc file.


    :param file_to_process: Path to the xtc file to be analysed.
    :type file_to_process: str
    :param plumed_monitor_file: Path to the plumed input file, should print to COLVAR_MONITOR
    :type plumed_monitor_file: str
    :param PLUMED_PATH: Path to the plumed executable.
    :type PLUMED_PATH: str, optional

    :return: Average value of the CV in the xtc file.
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
