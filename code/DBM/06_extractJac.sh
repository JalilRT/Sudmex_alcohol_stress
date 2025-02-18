#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab
### Usage example: sbatch code/05_extractJac.sh smooth###
smooth=$1

path="/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo/DBM/"
cd $path

Rscript ../code/06_extractJac.R ${smooth} --save
