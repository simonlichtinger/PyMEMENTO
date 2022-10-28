#!/bin/bash

# Generic options:

#SBATCH --account=bdhbs02  # Run job under project <project>
#SBATCH --time=48:0:0         # Run for a max of 48 hour

# Node resources:
# (choose between 1-4 gpus per node)

#SBATCH --partition=gpu    # Choose either "gpu" or "infer" node type
#SBATCH --nodes=1          # Resources from a single node
#SBATCH --gres=gpu:1       # One GPU per node (plus 25% of node CPU and RAM per GPU)

module load hecbiosim
module add gromacs/2020.4-plumed-2.6.2

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
bede-mpirun --bede-par 1ppg mdrun_mpi -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -r em.gro -p topol.top -o npt.tpr -maxwarn 1
bede-mpirun --bede-par 1ppg mdrun_mpi -deffnm npt

gmx grompp -f prod_sc.mdp -c npt.gro -r em.gro -p topol.top -o prod_sc.tpr
bede-mpirun --bede-par 1ppg mdrun_mpi -deffnm prod_sc

