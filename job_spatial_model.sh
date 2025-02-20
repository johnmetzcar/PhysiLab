#!/bin/bash

#SBATCH --mail-user=jpmetzca@iu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=spatial000000
#SBATCH -p general
#SBATCH -o spatial_model_%j.txt
#SBATCH -e spatial_model_%j.err
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH -A r00241


module load python/3.9.8
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -A r00241 --cpu-bind=sockets python RunSimulations_testing.py leukemia_spatial_model_files
# srun -A r00241 python scripts/loadDataParallel.py leukemia_spatial_output dataframes aggregated_live_cells_spatial.csv

#SBATCH --mem=240G
