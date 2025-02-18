

for i in {1..1162}; do fslmaths /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii -thr $i -uthr $i ROI_${i}.nii.gz; done

fslmaths /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii -mul 0 DBM_results/ROI_mask

for i in $(cat DBM_results/rois_mask.txt); do  fslmaths DBM_results/ROI_mask.nii.gz  -add ROIs_atlas/ROI_${i}.nii.gz DBM_results/ROI_mask.nii.gz ; done
