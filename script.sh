#!/bin/bash -l

#SBATCH --partition=workq
#SBATCH --account=q97
#SBATCH --job-name=uwgeo_demo
#SBATCH --nodes=20
#SBATCH --time=20:00:00

echo "PRINTING ENVIRONMENT"
env

# load Singularity
module load singularity

srun -N 20 -n 480 singularity exec \
        /group/q97/ak5446p/images2/v2.9.x_latest.sif \
	bash -c "python3 file3.py"
