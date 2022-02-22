#!/bin/bash
#SBATCH -A p30009               # Allocation
#SBATCH -p normal               # Queue: short, normal, long
#SBATCH -t 11:00:00             # Walltime/duration of the job
#SBATCH -nodes 1                # Number of Nodes
#SBATCH --ntasks-per-node=2     # Number of Cores (Processors)
#SBATCH --mail-user=<email_id>  # Designate email address for job communications
#SBATCH --mail-type=END         # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --output=<file_path>    # Path for output must already exist - e.g. /path/job.out
#SBATCH --error=<file_path>     # Path for errors must already exist - e.g. /path/job.err
#SBATCH --job-name="<name>"     # Name of job


# add a project directory to your PATH (if needed)
# export PATH=$PATH:/projects/p30009/<dirname>

# load modules you need to use
module purge
module use /software/spack_production/spack/share/spack/modules/linux-rhel7-x86_64/
module load gcc/10.3.0-gcc
module load intel-oneapi-compilers/2021.3.0-gcc
module load intel-mkl/2020.4.304-intel
module load cuda/11.4.0-intel
source /software/spack_production/spack/opt/spack/linux-rhel7-x86_64/intel-2021.3.0/gromacs-2021.5/bin/GMXRC

export OMP_NUM_THREADS=13

## A command you actually want to execute: Example showed for running gmx

## You may want to create the tpr for mdrun in the cluster core itself for compatibility reasons
# gmx grompp -f sim_eq.mdp -p PA.top -c sim_min.gro -o sim_eq.tpr -maxwarn 2

## Run mdrun on single core. Throw the output to out.log
# gmx mdrun -deffnm sim_eq -v -cpt 30 &> out.log

## Run mdrun on multiple cores. Throw the output to out.log
gmx -ntmpi 2 mdrun -deffnm sim_eq -v -cpt 30 &> out.log
