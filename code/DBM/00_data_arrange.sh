#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

bids=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/BIDS
path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo
input=${path}/data


for ses in T1 T2 T3 T5; do
 for i in $(cat code/subjid); do 
  mkdir -p ${input}/$i/ses-${ses}/anat/

  cp -r ${bids}/${i}/ses-${ses}/anat/* ${input}/$i/ses-${ses}/anat/
  
  rm -d ${input}/$(basename $i)/ses-${ses}/anat
  rm -d ${input}/$(basename $i)/ses-${ses}/
  rm -d ${input}/$(basename $i)/

  chmod -R 777 ${input}/

 done
done
