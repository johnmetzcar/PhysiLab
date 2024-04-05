#!/bin/bash

#SBATCH --mail-user=jpmetzca@iu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=spatial000000
#SBATCH -p general
#SBATCH -o load_data_%j.txt
#SBATCH -e load_data_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -A r00241
#SBATCH --mem=240G

module load python/3.9.8
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -A r00241 --cpu-bind=sockets python scripts/loadDataParallel.py


