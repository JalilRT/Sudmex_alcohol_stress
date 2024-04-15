#!/bin/bash

#set -x

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas
ROIs=${path}/ROIs_atlas
ROIs_fmri=${path}/ROIs_fmri
mkdir -p ${ROIs_fmri}
module load cobralab

cp ${ROIs}/rois_mask.txt ${ROIs_fmri}/rois_mask.txt

for i in $(cat ${ROIs_fmri}/rois_mask.txt); do cp ${ROIs}/ROI_${i}.nii.gz ${ROIs_fmri}/ ; done

for i in $(ls ${ROIs_fmri}/ROI_*gz); do fslmaths ${i} -ero ${i}; done

cp ${ROIs}/ROI_VTA_* ${ROIs_fmri}/
cp ${ROIs}/ROI_6[1,2].* ${ROIs_fmri}/
cp ${ROIs}/ROI_17[1,2].* ${ROIs_fmri}/
cp ${ROIs}/ROI_18[1,2].* ${ROIs_fmri}/

fslmaths ${path}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii -mul 0 ${ROIs_fmri}/Mask_erode

for i in $(ls ${ROIs_fmri}/ROI_*gz); do roi=$(basename ${i} .nii.gz); fslmaths ${ROIs_fmri}/Mask_erode.nii.gz -add ${ROIs_fmri}/${roi}.nii.gz ${ROIs_fmri}/Mask_erode.nii.gz ; done

nii2mnc ${ROIs_fmri}/Mask_erode.nii.gz ${ROIs_fmri}/Mask_erode.mnc

fslmaths ${ROIs_fmri}/Mask_erode.nii.gz -bin ${ROIs_fmri}/Mask_erode_bin.nii.gz

nii2mnc ${ROIs_fmri}/Mask_erode_bin.nii.gz ${ROIs_fmri}/Mask_erode_bin.mnc

chmod -R 777 ${ROIs_fmri}
