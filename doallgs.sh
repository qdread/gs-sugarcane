#!/bin/bash
#SBATCH --job-name=gs
#SBATCH --ntasks=72
#SBATCH --mem=80gb
#SBATCH --partition=medium
#SBATCH --time=7-00:00:00

cd /home/quentin.read/GitHub/gs-sugarcane
module load r/4.1.2
Rscript2 GS_sugarcane_2022_runGS.R
