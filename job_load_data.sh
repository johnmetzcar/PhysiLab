#!/bin/bash

#SBATCH --mail-user=jpmetzca@iu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=spatial000000
#SBATCH -p general
#SBATCH -o load_data_%j.txt
#SBATCH -e load_data_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH -A r00241


module load python/3.9.8
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -A r00241 --cpu-bind=sockets python scripts/create_population_df.py

#SBATCH --mem=240G
