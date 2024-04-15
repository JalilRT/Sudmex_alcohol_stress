#!/bin/bash

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas
ROIs=${path}/ROIs_atlas
ROIs_fmri=${path}/ROIs_fmri
mkdir -p ${ROIs}
mkdir -p ${ROIs_fmri}
module load cobralab

# How to use
# For example we have insular values: 11, 21, 121, 451 and 12, 22, 122, 452. Tag to 11 or 12
# run:
# bash rois_mask_joining.sh

echo "Enter your tag"
read tag
echo "Enter the numbers of each ROI to merge"
read numbers

echo "#### Running left hemisphere ####"
echo "tag: "${tag}" "
echo "labels: "${numbers}" "

 fslmaths ${path}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii \
  -mul 0 ${ROIs_fmri}/ROI_${tag}.nii.gz

 for i in ${numbers}; do

  fslmaths ${ROIs_fmri}/ROI_${tag}.nii.gz -add ${ROIs}/ROI_${i}.nii.gz ${ROIs_fmri}/ROI_${tag}.nii.gz

 done

 fslmaths ${ROIs_fmri}/ROI_${tag}.nii.gz -bin ${ROIs_fmri}/ROI_${tag}_bin.nii.gz

 fslmaths ${ROIs_fmri}/ROI_${tag}_bin.nii.gz -mul ${tag} ${ROIs_fmri}/ROI_${tag}_mul

 rm ${ROIs_fmri}/ROI_${tag}.nii.gz ${ROIs_fmri}/ROI_${tag}_bin.nii.gz

 mv ${ROIs_fmri}/ROI_${tag}_mul.nii.gz ${ROIs_fmri}/ROI_${tag}.nii.gz


tag2=$(
for t in ${tag}; do 
 sum=$((${t} + 1))
 echo $sum
done
)

numbers2=$(
for s in ${numbers}; do
 suma=$((${s} + 1))
 echo $suma
done
)

echo ""
echo "#### Running right hemisphere ####"
echo "tag: "${tag2}" "
echo "labels: "${numbers2}" "

 fslmaths ${path}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii \
  -mul 0 ${ROIs_fmri}/ROI_${tag2}.nii.gz

 for i in ${numbers2}; do

  fslmaths ${ROIs_fmri}/ROI_${tag2}.nii.gz -add ${ROIs}/ROI_${i}.nii.gz ${ROIs_fmri}/ROI_${tag2}.nii.gz

 done

 fslmaths ${ROIs_fmri}/ROI_${tag2}.nii.gz -bin ${ROIs_fmri}/ROI_${tag2}_bin.nii.gz

 fslmaths ${ROIs_fmri}/ROI_${tag2}_bin.nii.gz -mul ${tag2} ${ROIs_fmri}/ROI_${tag2}_mul

 rm ${ROIs_fmri}/ROI_${tag2}.nii.gz ${ROIs_fmri}/ROI_${tag2}_bin.nii.gz

 mv ${ROIs_fmri}/ROI_${tag2}_mul.nii.gz ${ROIs_fmri}/ROI_${tag2}.nii.gz


