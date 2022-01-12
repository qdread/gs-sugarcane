#!/bin/bash
#SBATCH --job-name=do_all_gs
#SBATCH --ntasks=16
#SBATCH --mem=16gb
#SBATCH --time=3-00:00:00

cd /home/quentin.read/GitHub/gs-sugarcane
module load r/4.1.2
Rscript2 GS_sugarcane_2022.R
