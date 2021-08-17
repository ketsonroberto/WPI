#!/bin/sh
#
#SBATCH --account=free      # The account name for the job.
#SBATCH --job-name=CPU    # The job name.
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=2:00:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=12gb        # The memory the job will use per cpu core.

./main2dof
 
# End of script
