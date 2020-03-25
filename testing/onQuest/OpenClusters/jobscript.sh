#!/bin/bash
#SBATCH -A p30137               # Allocation
#SBATCH -p short                # Queue
#SBATCH -t 04:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH --mem=18G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --ntasks-per-node=6     # Number of Cores (Processors)
#SBATCH --mail-user=andrewbowen2020@u.northqestern.edu # Designate email address for job communications
#SBATCH --mail-type=END     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --output=/projects/p30137/abowen/CEB/testing/OpenClusters/    # Path for output must already exist
#SBATCH --error=/projects/p30137/abowen/CEB/testing/OpenClusters/     # Path for errors must already exist
#SBATCH --job-name="analyse test"       # Name of job

# unload any modules that carried over from your command line session
module purge

# add a project directory to your PATH (if needed)
export PATH=$PATH:/projects/p30137/abowen/CEB/testing/OpenClusters/

# load modules you need to use
module load python/anaconda3.6

python analyseEBLSSTcluster-OC.py
