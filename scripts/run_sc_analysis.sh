#!/bin/bash
#SBATCH --job-name=sc_harmony
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load python
python run_single_cell_analysis.py
