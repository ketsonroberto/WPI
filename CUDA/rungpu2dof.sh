#!/bin/sh
#
#SBATCH --account=free      # The account name for the job.
#SBATCH --job-name=GPU2DOF    # The job name.
#SBATCH --gres=gpu:1             # Request 1 gpu (Up to 4 on K80s, or up to 2 on P100s are valid).
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=10:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=12gb        # The memory the job will use per cpu core.
 
module load cuda80/toolkit
./array2D
 
# End of script

