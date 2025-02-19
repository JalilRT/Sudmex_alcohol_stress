#!/bin/bash

### The following scripts is to perform group level analysis using multivariate approach ###
###			Using longitudinal sessions of melodic dual regression		 ###
## How to use bash code/06_3dLMEr_melodic.sh ##

module load cobralab
# Setting paths
lv_path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/fmri
path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/fmri
atlas_model=${path}/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Functional_Imaging/SIGMA_EPI_Brain_Mask.nii 
cd $lv_path

#sbatch code/LME_SB/t2z_transform.sh

#Rscript code/LME_SB/arranging_datasetDR.R

if [ ! -f ${path}/Dual_reg/SIGMA_mask_registered.nii.gz ]; then

 echo "sending mask to subjects space"

 antsRegistration_affine_SyN.sh \
  --float --fast \
  --clobber \
  ${atlas_model} \
  ${path}/rabies_smri/preproc_out/Session-T5/bold_datasink/commonspace_labels/_scan_info_subject_id42.sessionT5_split_name_sub-42_ses-T5_run-02_T2w/_run_1/sub-42_ses-T5_task-rest_run-01_bold_RAS_EPI_anat_labels.nii.gz \
  ${path}/Dual_reg/SIGMA_mask_registered

 antsApplyTransforms -d 3 -i ${atlas_model} -t ${path}/Dual_reg/SIGMA_mask_registered0GenericAffine.mat \
  -t ${path}/Dual_reg/SIGMA_mask_registered1InverseWarp.nii.gz \
  -r ${path}/rabies_smri/preproc_out/Session-T5/bold_datasink/commonspace_labels/_scan_info_subject_id42.sessionT5_split_name_sub-42_ses-T5_run-02_T2w/_run_1/sub-42_ses-T5_task-rest_run-01_bold_RAS_EPI_anat_labels.nii.gz \
  -o ${path}/Dual_reg/SIGMA_mask_registered.nii.gz -n GenericLabel --verbose

else

 echo "mask already in subjects space"

fi

bash code/LME_DR/model4H.sh
