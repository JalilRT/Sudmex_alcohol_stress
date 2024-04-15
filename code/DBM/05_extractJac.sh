#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab
### Usage example: sbatch code/05_extractJac.sh ###

cd /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo/DBM/
Rscript ../code/05_extractJac.R
