#!/bin/bash

cd /scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo/preproc
path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo

rm *txt
for i in $(ls -d sub-66/ses-??); do
 sub=$(dirname $i)
 ses=$(basename $i)
 echo $PWD/${sub}/${ses}/anat/${sub}_${ses}_pp.nii
done >> out.txt
 
for i in $(wc -l out.txt); do
 head -n 1 out.txt >> aut.txt
#line=$(cat out.txt) 
#li=$(cat aut.txt)
 grep -vxFf aut.txt out.txt >> meh.txt
 echo ${line},
done

for i in $(ls -d sub-*/); do echo $PWD/$(basename $i)/ses-T1/anat/$(basename $i)_ses-T1_pp.nii,$PWD/$(basename $i)/ses-T2/anat/$(basename $i)_ses-T2_pp.nii,$PWD/$(basename $i)/ses-T3/anat/$(basename $i)_ses-T3_pp.nii,$PWD/$(basename $i)/ses-T5/anat/$(basename $i)_ses-T5_pp.nii >> ../code/input_jacobians.csv; done
