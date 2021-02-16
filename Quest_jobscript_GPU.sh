#!/bin/bash
#SBATCH -A p30009               # Allocation
#SBATCH -p gengpu               # Queue: short, normal, long
#SBATCH -t 11:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --ntasks-per-node=2     # Number of Cores (Processors)
#SBATCH --gres gpu:2            # Extra line for gpu jobs specifying number of gpus
#SBATCH --mem=10G               # Memory request per node
#SBATCH --mail-user=<email_id>  # Designate email address for job communications
#SBATCH --mail-type=END         # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --output=<file_path>    # Path for output must already exist - e.g. /path/job.out
#SBATCH --error=<file_path>     # Path for errors must already exist - e.g. /path/job.err
#SBATCH --job-name="<name>"     # Name of job

# unload any modules that carried over from your command line session
module purge

# add a project directory to your PATH (if needed)
# export PATH=$PATH:/projects/p30009/<dirname>

# load modules you need to use
module load gromacs/5.0.4-cuda

## A command you actually want to execute: Example showed for running gmx

## You may want to create the tpr for mdrun in the cluster core itself for compatibility reasons
gmx_mpi grompp -f sim_eq.mdp -p PA.top -c sim_min.gro -o sim_eq.tpr -maxwarn 2

## Run mdrun on single core. Throw the output to out.log
# mpirun -np 1 gmx_mpi mdrun -deffnm sim_eq -v -cpt 30 &> out.log

## Run mdrun on multiple cores. Throw the output to out.log
mpirun -np 2 gmx_mpi mdrun -deffnm sim_eq -v -cpt 30 &> out.log

# Using a quick benchmark, 
# GROMACS 2020.3 runs fastest using GPU with the following command on Quest
# mpirun -np 4 gmx_mpi mdrun -deffnm md -v -cpt 5 -rdd 2 -nb gpu -pme gpu -npme 1 # using 4 GPUs
# gmx mdrun -deffnm md -v -cpt 5 -nb gpu -pme gpu -update gpu # using 1 GPU


