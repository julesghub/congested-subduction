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

SINGULARITYENV_LD_LIBRARY_PATH=/opt/cray/pe/mpt/7.7.0/gni/mpich-gnu-abi/4.9/lib:/opt/cray/xpmem/2.2.15-6.0.7.1_5.16__g7549d06.ari/lib64:/opt/cray/ugni/6.0.14.0-6.0.7.1_3.18__gea11d3d.ari/lib64:/opt/cray/udreg/2.3.2-6.0.7.1_5.18__g5196236.ari/lib64:/opt/cray/pe/pmi/5.0.13/lib64:/opt/cray/alps/6.6.43-6.0.7.1_5.54__ga796da32.ari/lib64:/opt/cray/wlm_detect/1.3.3-6.0.7.1_5.8__g7109084.ari/lib64:/usr/lib64 srun -N 2 -n 48 singularity exec \
/group/q97/ak5446p/images3/uwgeodynamics_v2.9.5.sif \
bash -c "python3 geo_moresi_2014_subduction.py"
