#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab
### Usage example: sbatch code/07_plotTmaps.sh smooth###
smooth=$1

path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo
cd $path

Rscript code/08_VolumesROIs.R ${smooth} --save
