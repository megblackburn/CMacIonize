#!/bin/bash
#SBATCH --job-name=tdi_tc_f185
#SBATCH --partition=large-long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=90
#SBATCH --output=%x_%j.log
#SBATCH --time=24:00:00
# Load required modules

conda activate mgbenv # insert your env name here 

echo date
../../../CMacIonize --params bursty_sfr_tallbox.param --task-based-rhd --threads 90
echo date
