#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab

# Reading variables
path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri

cd ${path}/Data

read -p "which session do you want to run? (i.e. T5): " session
read -p "which subject do you want to run? (i.e. sub-42): " subject

for i in $(ls -d ${subject}/ses-${session}); do 
 sub=$(dirname $i)
 ses=$(basename $i)

 mkdir -p ../data/Session-${session}/${sub}/ses-${session}
 cp -r ${path}/Data/${sub}/ses-${session}/* ${path}/data/Session-${session}/${sub}/ses-${session}/
 rm ${path}/data/Session-${session}/${sub}/ses-${session}/anat/*.mnc

done
