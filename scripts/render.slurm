#!/bin/bash
#SBATCH --job-name=stars_render
#SBATCH --output=render_%j.out
#SBATCH -N 2
#SBATCH -n 16
#SBATCH -c 4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3500MB
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=pablo.rey.ansemil@udc.es

module load paraview/5.11.0
# paths
INIT_DIR=
FINAL_DIR=/path/to/chunks_final
OUT_DIR=/path/to/output_frames
NFRAMES=200

srun pvbatch render.py $INIT_DIR $FINAL_DIR $OUT_DIR $NFRAMES
