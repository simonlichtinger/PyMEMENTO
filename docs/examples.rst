Examples and tutorial
=====================

This section gives several examples of how to use the PyMEMENTO package to perform MEMENTO. A detailed analysis of the results from these
runs will be presented in our paper. We assume in this tutorial that the user is familiar with running MD simulations using Gromacs. Experience
with the PLUMED package is required for setting up umbrella sampling simulations.


Simple use case --- deca-alanine
--------------------------------

The simplest system we have used to validate PyMEMENTO is deca-alanine, a capped peptide composed 
of 10 alanine residues that can exist as either a helix or an unfolded, extended state.

For the most basic use case, we need two end-state coordinate files: ``helix.gro`` and ``extended.gro``
which can be found under the `test data directory <https://github.com/simonlichtinger/PyMEMENTO/tree/master/PyMEMENTO/tests/test_data/decaalanine>`_.
Copy these into a location that is convenient for you, and then run the following python script.

.. code-block:: python

    from PyMEMENTO import MEMENTO

    # Initiliase an instance of the MEMENTO class
    model = MEMENTO(
            "example_decaalanine/",
            "helix.gro", # Change to wherever you have saved the file
            "extended.gro", # Change to wherever you have saved the file
            list(range(1, 11))
        )

    # Perform coordinate interpolation in 24 windows
    model.morph(24)

    # Make 5 models per intermediate to fix morph, using all residues except for caps (give a list of 
    # residue numbers like in the class initialisation if you do not want the default, which would be 
    # a continuous range from the minimum to the maximum residue given in the constructor 
    # --- for deca-alanine, that would have been fine).
    model.make_models(5, include_residues=list(range(1, 11)))

    # Find the smoothest path through the space of models
    model.find_best_path()

    # Process the models for MD, including caps by default
    model.process_models()

    # Make MD boxes for the models
    model.prepare_boxes()

    # Solvate and add ions to the MD boxes
    model.solvate_boxes(ion_concentration=0.15)

    # Perform an energy minimisation on the MD boxes
    model.minimize_boxes()


.. note::  
    It is a known issue that the pathfinding can crash if too many models have to be loaded by MDAnalysis, which can be limited by the system as MDAnalysis sometimes doesn't close files after loading if the universe is preserved.
    If you're making a lot of models per intermediate and encounter an OSError, cannot open file, then try running ``` ulimit -n 10000 ``` (or whatever number of models you need) to increase the number of files that can be open at the same time.

The script first initialises a :class:`PyMEMENTO.MEMENTO` instance, which holds some core information: a working path which will later
contain the created directory tree, an initial coordinate file, a final coordinate file and a list with the residue numbers that should
be assigned to the protein on completion. This is useful if your protein has dis-continuous residue numbering or residue numbers do not start at 1. In the case of deca-alanine, ``list(range(1, 11))`` simply tells PyMEMENTO that we want a list of residue numbers as [1, 2, 3, 4, 5, 6, 7, 8, 9, 10].  Note that every step PyMEMENTO performs is saved into the working directory, so MEMENTO can be restarted after
only some of the following commands have run, provided that the ``last_step_performed`` keyword argument is passed.

We then perform the steps of the MEMENTO procedure, which are presented as individual methods of the class (here, for simplicity
we choose to make 5 models per window, in our paper we use 50). The order in which these methods
are to be called is fixed (there will be an error if one attempts otherwise), however each step can be customised via the keyword arguments
as described in the documentation of the individual methods. Note that PyMEMENTO only ships template topology files for simple systems and the
amber14.sb and CHARMM36 (July 2021) forcefields, so you will need to provide your own template folders should you wish to extend your work beyond
this scope.

The main output of the procedure lies in the ``boxes`` subdirectory of the working path. These are now ready for running equilibration MD, for 
which appropriate gromacs mdp files (``nvt.mdp``, ``nvt.mdp`` and ``prod_sc.mdp`` --- to further equilibrate the side chains) will already be in the folders if you are using the built-in template --- as in this example.
These equilibrations are by default 200ps NVT ensemble MD at a timestep of 1fs (v-rescale thermostat), then 1ns NPT at a time-step of 2fs (berendsen barostat and v-rescale thermostat)
followed by further 10ns NPT (Parrinello-Rahman barostat, v-rescale thermostat). All of these are run with position restraints on the C-alpha atoms of the protein.

.. note:: 
    Running this will be most efficient on a distributed cluster, so you will need to provide your own submission script for this purpose.

Setting up umbrella sampling
----------------------------

The previous example covers the core functionality of PyMEMENTO, and you may wish to process the output and set up free energy simulations as is appropriate for your test system.
There is however further code available if you intend to set up 1-dimensional umbrella sampling with gromacs and plumed.

After the equilibrations are complete and have been transferred to the original directory structure, the following python code will set up umbrella sampling.

.. code-block:: python

    from PyMEMENTO import MEMENTO

    # Keep the initialisation line from the previous example
    # as PyMEMENTO needs to know where to expect files
    model = MEMENTO(
            "example_decaalanine/",
            "helix.gro",
            "extended.gro",
            list(range(1, 11))
        )
    
    # Set up umbrella sampling
    model.prepare_umbrella_sampling(
        "plumed_monitor.dat",
        "plumed.dat",
        smoothen_ladder=0.5
    )

Running this script requires two further input files, ``plumed_monitor.dat`` and ``plumed.dat``, examples of which can be found `here <https://github.com/simonlichtinger/PyMEMENTO/tree/master/PyMEMENTO/tests/test_data/decaalanine/plumed_distance>`_.

.. code-block::
    :caption: plumed_monitor.dat

    d: DISTANCE ATOMS=9,99

    PRINT ARG=d STRIDE=1  FILE=COLVAR_MONITOR

PyMEMENTO uses plumed to extract the value of an arbitrary collective variable as the average across the prod_sc equilibration trajectories. The ``plumed_monitor.dat``
input file defines the collective variable of interest. You can modify this file to fit your purposes, as long as it prints a single column of values into a file called
``COLVAR_MONITOR`` at every frame of the trajectory (stride 1). In this case, it computes the end-to-end distance of deca-alanine.

The average values of each collective variable will be inserted into a provided ``plumed.dat`` file for running MD. You can use any plumed file for this, PyMEMENTO only replaces the
``$REPLICAS$`` string.

.. code-block::
    :caption: plumed.dat

    d: DISTANCE ATOMS=9,99

    restraint: RESTRAINT ARG=d AT=@replicas:$REPLICAS$ KAPPA=1000

    PRINT ARG=* STRIDE=5000  FILE=../COLVAR_MULTI

This should be passed via the ``-plumed plumed.dat`` flag to gromacs to run umbrella sampling MD.

.. note::
    The method :meth:`PyMEMENTO.MEMENTO:setup_umbrella_sampling` takes the keyword argument ``smoothen_ladder`` which is set to 0.5 for deca-alanine.
    This variable determines to what extent the collective variable spacing between umbrella sampling windows should be uniform (smoothen_ladder=1) or follow
    the values found in the simulations trajectories (smoothen_ladder=0). Because deca-alanine is such as small system with a drastic conformational change,
    the modeller step doesn't produce a uniform spacing of collective variable values along the intermediates if the end-to-end distance collective variable is
    used, which can lead to poor histogram overlap in some cases.
    By trial and error, we found that a smoothening of 0.5 works well for deca-alanine. For all of our larger systems this was not necessary, however, so unless you have good reason 
    you should keep the default value of zero here.



Handling ligands
----------------

A note on limitations of ligand morphing
........................................

MEMENTO fixes unphysical morphing intermediates using modeller to reconstruct sensible intermediate states that perform well in umbrella sampling.
This works only for proteins at this stage: fixing arbitrary ligand topologies is not supported and will not be without significant future effort. 
Nonetheless, PyMEMENTO can morph ligands if they are present in the starting and end coordinate files.

The way in which this is achieved is handled by the ``ligand_type`` keyword argument in the MEMENTO constructor :meth:`PyMEMENTO.MEMENTO.__init__`. The *'rigid'* value (default)
will translate the centre-of-mass of the ligand as a linear morph, but retain all internal ligand strucutre. This option is useful for complex organic ligands. In all all of our validation cases,
the equilibration runs were sufficient to generate a reasonable ligand conformation for each intermediate. **It is essential that you check this manually to avoid unphysical intermediates.**

The *'single'* option instructs PyMEMENTO to linearly interpolate ligand coordinates. There is no fixing of interactions with the protein or of ligand topology, so this will only work for organic molecules
if the movement of the ligand is very small and an energy minimisation is sufficient to restore structure. The option is thus mainly useful for cases like multiple ion binding sites in a protein, where
the independent displacement of ligand molecules is key.

A real-world example --- adenylate kinase (ADK)
...............................................

As described in our publication, adenylate kinase has a large ligand-dependent conformational change between open and closed states. The conformational equilibrium
as a function of the ligand has been studied in previous publications by various authors (which we review in our paper), and it is known that the open
state is the global minimum of apo-ADK, while holo-ADK favours the closed state.

The ligand is a rather complex inhibitor and its position is known for the closed state (pdb 1AKE), while the open state (pdb 4AKE) was solved without the inhibitor.
We can run MEMENTO of holo ADK only by assuming a similar ligand position in the open state and running MD for equilibration. Suitably prepared input files for
MEMENTO can be found `here <https://github.com/simonlichtinger/PyMEMENTO/tree/master/PyMEMENTO/tests/test_data/ADK>`_.

.. code-block:: python

    from PyMEMENTO import MEMENTO

    # Initiliase an instance of the MEMENTO class, change all paths to where you've saved the files
    model = MEMENTO(
            "example_ADK_holo/",
            "holo_closed.gro",
            "holo_open.gro",
            list(range(1, 215)),
            ligand="resname AP5"
        )

    # Perform morphing and modelling with 24 windows and 5 models per window
    model.morph(24)
    model.make_models(5)
    model.find_best_path()
    model.process_models(caps=False, his=True, his_protonation_states=[0,1,0])
    model.prepare_boxes(template_folder="template_ADK_holo")

    # Solvate and energy-minimise the boxes
    model.solvate_boxes(ion_concentration=0.15)
    model.minimize_boxes()

Running this script will prepare gromacs simulation inputs for restrained equilibrations. The subsequent setup of umbrella sampling 
works the same as for deca-alanine above, though ``plumed.dat`` and ``plumed_monitor.dat`` need to be modified to fit the system.

.. note::
    ADK contains three histidine residues. By default, these would be assigned by gmx pdb2gmx, however this is almost never a good idea
    because different histidine prototation states between umbrella sampling windows will crash replica exchange (because the atom ordering
    and connectivity are not identical). It is therefore recommended to first run pdb2gmx on one of the end states and use the resulting 
    histidine protonation states for all intermediates, as was done in this example. In order to avoid problems, with pre-assigned histidine
    protonation states in your input coordinate files (which come with residue names like HISE, depending on the forcefield), it is often safest to 
    rename all histidines to 'HIS' before passing to PyMEMENTO and reassigning them later.


Handling lipids
---------------

How MEMENTO handles lipids
...........................

One of the strengths of MEMENTO (for the purpose of which we have indeed developed the method) is the ability to handle membrane proteins, where
techniques based for example on normal-mode analysis will fail and hysteresis often becomes more pronounced and problematic for approaches such
as metadynamics.

Much like for ligands molecules, reconstructing lipid conformations and interactions is beyond the capabilities of the modeller package. If however, 
the change in protein--lipid interactions is small enough to be sampled with a feasible amount of equilibration time (we find that 100ns per window
usually are sufficient for our systems), MEMENTO can treat the membrane as essentially static.

Usually, the shape presented by the protein to the membrane changes in the process of a conformational change. PyMEMENTO therefore expects a
lipid membrane as input surrounding whatever conformation is deemed 'wider' (this needs to be the target state by definition). It first stretches the
membrane in the x-y plane to allow for new intermediate side-chain conformations to fit, followed by a sequential compression and energy minimisations to
pack the membrane snugly to the protein. This is essentially a python reimplementation of the inflate-gro method, and is
done in the :meth:`PyMEMENTO.MEMENTO:prepare_boxes` method. The parameters of the embedding procedure like initial expansion coefficient, number
of reduction rounds and size reduction coefficient can be set manually, but the defaults usually work well.


A membrane protein example - LeuT
................................

LeuT is a bacterial leucine transport protein that is thought to follow an alternating-access conformational cycle for which several
conformations are known experimentally. For our paper, we have chosen the occluded-outwards-facing (OCC/OF, pdb 3F3E) <-> inwards-facing (IF, pdb 3TT3) transition
as a prospective application. Coordinate files of the OCC/OF and POPE-embedded IF states can be found `here <https://github.com/simonlichtinger/PyMEMENTO/tree/master/PyMEMENTO/tests/test_data/LEUT>`_
and a template folder for MD equlibration here `here <https://github.com/simonlichtinger/PyMEMENTO/tree/master/PyMEMENTO/tests/test_data/template_LEUT_if>`_.

.. note::  
    When running equilibration MD of the CHARMM-GUI-embedded boxes, we noticed that the IF state is not stable in unbiased
    MD. Transmembrane helix 1 closes up over a few hundred nanoseconds, which we address in detail in our publication. For the purposes of illustrating
    the use of lipids with PyMEMENTO in this
    tutorial, just bear in mind that the IF-like states will not be stable without position restraints.

The following python script can be used to set up MEMENTO boxes for equilibration.

.. code-block:: python

    from PyMEMENTO import MEMENTO
    
    # Set up MEMENTO class, replace paths as appropriate for your file tree
    model = MEMENTO(
        "example_LEUT/",
        "occ.gro",
        "if.gro",
        list(range(11, 508)),
        forcefield="CHARMM36",
        lipid="resname POPE",
    )

    # Perform morphing and modelling
    model.morph(24)
    model.make_models(5, include_residues=list(range(1, 498)))
    
    # using poolsize=1 as multiprocessing can trigger an MDAnalysis memory bug
    # here if the system is as large as LEUT
    model.find_best_path(poolsize=1)
    model.process_models(
        caps=True,
        cap_type="CHARMM",
        his=True,
        his_protonation_states=[0, 0, 0, 0],
    )
    model.prepare_boxes("template_LEUT_if")

    model.solvate_boxes(ion_concentration=0.15)
    model.minimize_boxes()
