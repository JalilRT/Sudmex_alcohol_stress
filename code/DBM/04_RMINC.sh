#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab
### Usage example: sbatch code/04_RMINC.sh ###

cd /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo/DBM/
Rscript ../code/04_RMINC.R

for i in $(ls -d tmaps/*mnc); do su=$(basename $i .mnc); mnc2nii -nii $i tmaps/${su}.nii.gz; done

chmod 777 -R *
