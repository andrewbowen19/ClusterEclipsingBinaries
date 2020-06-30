#!/bin/bash
#SBATCH -A p30137               # Allocation
#SBATCH -p short                # Queue
#SBATCH -t 04:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=18G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --ntasks-per-node=6     # Number of Cores (Processors)
#SBATCH --mail-user=andrewbowen2020@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=END     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="EBLSST Histogram run"       # Name of job

# unload any modules that carried over from your command line session
module purge

# add a project directory to your PATH (if needed)
export PATH=$PATH:/projects/p30137/projects/abowen/CEB/testing/eblsst/

# load modules you need to use
module load python/anaconda3.6

# A command you actually want to execute:

# Another command you actually want to execute, if needed:
python runPlotter.py
