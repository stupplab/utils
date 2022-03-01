#!/bin/bash
#SBATCH -A p30009               # Allocation
#SBATCH -p gengpu               # Queue: short, normal, long
#SBATCH -t 11:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --ntasks-per-node=4     # Number of Cores (CPU Processors). 2 GPU nodes have 52 available CPU cores
#SBATCH --gres=gpu:2            # Extra line for gpu jobs specifying number of gpus
#SBATCH --mem=10G               # Memory request per node
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


# run
gmx grompp -f em.mdp -c *_water.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v -ntmpi 4 -ntomp 13
gmx grompp -f md_anneal.mdp -c em.gro -p topol.top -o md.tpr
gmx mdrun -ntmpi 4 -ntomp 13 -deffnm md -v -cpt 5 -nb gpu -pme gpu -npme 1 &> log

# continue initial run
#gmx mdrun -ntmpi 4 -ntomp 13 -deffnm md -v -cpt 5 -cpi md.cpt -append -nb gpu -pme gpu -npme 1 &> log

