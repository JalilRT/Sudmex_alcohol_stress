#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab

## Usage example: sbatch code/01_rat-average.sh T1 sub-*

ses=$1
sub=$2
path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo
input=${path}/data
preproc=${path}/preproc

for i in ${input}/${sub}; do 
 mkdir -p ${preproc}/$(basename $i)/ses-${ses}/anat/

 bash ${path}/code/rat-invivo-average.sh $(basename $i) ses-${ses} run-02_T2w

 chmod -R 777 ${preproc}/

done
