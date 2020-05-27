#!/bin/bash -l

#SBATCH --partition=workq
#SBATCH --account=q97
#SBATCH --job-name=uwgeo_demo
#SBATCH --nodes=2
#SBATCH --time=00:20:00

echo "PRINTING ENVIRONMENT"
env

# load Singularity
module load singularity

srun -N 2 -n 48 singularity exec \
/group/q97/ak5446p/images3/uwgeodynamics_v2.9.5.sif \
bash -c "python3 geo_moresi_2014_subduction.py"
