#!/bin/sh
#
#SBATCH --account=free      # The account name for the job.
#SBATCH --job-name=CPU    # The job name.
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=00:30              # The time the job will take to run.
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.

./main
 
# End of script