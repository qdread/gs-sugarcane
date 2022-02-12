#!/bin/bash
#SBATCH --job-name=gs_md
#SBATCH --ntasks=64
#SBATCH --mem=32gb
#SBATCH --partition=medium
#SBATCH --time=7-00:00:00

cd /home/quentin.read/GitHub/gs-sugarcane
module load r/4.1.2
Rscript2 GS_sugarcane_2022_markerdensity.R
