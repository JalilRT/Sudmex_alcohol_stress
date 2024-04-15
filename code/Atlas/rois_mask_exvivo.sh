
#!/bin/bash

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas
ROIs=${path}/ROIs_atlas_exvivo
mkdir -p ${ROIs}
module load cobralab

for i in {1..1172}; do fslmaths ${path}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/ExVivo_Atlas/SIGMA_ExVivo_Anatomical_Brain_Atlas.nii -thr $i -uthr $i ${ROIs}/ROI_${i}.nii.gz; done

fslmaths ${path}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/ExVivo_Atlas/SIGMA_ExVivo_Anatomical_Brain_Atlas.nii -mul 0 ${ROIs}/ROI_mask

for i in $(cat ${ROIs}/rois_mask.txt); do fslmaths ${ROIs}/ROI_mask.nii.gz -add ${ROIs}/ROI_${i}.nii.gz ${ROIs}/ROI_mask.nii.gz ; done

nii2mnc ${ROIs}/ROI_mask.nii.gz ${ROIs}/ROI_mask.mnc

fslmaths ${ROIs}/ROI_mask.nii.gz -bin ${ROIs}/ROI_mask_bin.nii.gz

nii2mnc ${ROIs}/ROI_mask_bin.nii.gz ${ROIs}/ROI_mask_bin.mnc

chmod -R 777 ${ROIs}
