#!/bin/bash

set -x

module load cobralab

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis
ROIs_folder=${path}/fmri/Atlas/ROIs
#Template_fmri=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Functional_Imaging/SIGMA_EPI_Brain_Template_Masked.nii
Template_fmri=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template.nii

fslmaths ${Template_fmri} -mul 0 ${ROIs_folder}/ROI_fmri_mask
for i in $(cat ${ROIs_folder}/rois_seed_combined.txt); do fslmaths ${ROIs_folder}/ROI_fmri_mask.nii.gz \
 -add ${ROIs_folder}/${i} ${ROIs_folder}/ROI_fmri_mask.nii.gz ; done

cp ${ROIs_folder}/ROI_fmri_mask.nii.gz ${ROIs_folder}/../ROI_fmri_mask.nii.gz
gunzip ${ROIs_folder}/../ROI_fmri_mask.nii.gz

chmod -R 777 ${ROIs_folder}
