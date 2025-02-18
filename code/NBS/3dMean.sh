#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri
container=${path}/container/afni.sif
mkdir -p ${path}/Seed_based/3dMean
#structures=(left_granule_cell_level_of_the_cerebellum_left_2_Pmod4L left_granule_cell_level_of_the_cerebellum_left_Pmod4L left_molecular_layer_of_the_cerebellum_left_Pmod4L left_pre_limbic_system_left_Pmod4H left_presubiculum_left_Pmod6 right_cornu_ammonis_2_Pmod4H right_granule_cell_level_of_the_cerebellum_Pmod6 right_hypothalamic_region_Pmod4L right_olfactory_bulb_Pmod4H)

for roi in left_granule_cell_level_of_the_cerebellum_left_2_Pmod4L left_granule_cell_level_of_the_cerebellum_left_Pmod4L left_molecular_layer_of_the_cerebellum_left_Pmod4L left_pre_limbic_system_left_Pmod4H left_presubiculum_left_Pmod6 right_cornu_ammonis_2_Pmod4H right_granule_cell_level_of_the_cerebellum_Pmod6 right_hypothalamic_region_Pmod4L right_olfactory_bulb_Pmod4H; do
 # obtain mean by group
 for gp in Alc Ctrl; do
  #obtain mean by session
  for ses in T1 T2 T3 T5; do

   singularity exec -B /scratch:/scratch ${container} 3dMean \
    -prefix ${path}/Seed_based/3dMean/${ses}_${gp}_${roi}.nii.gz \
    ${path}/rabies/analysisDR_out/Session-${ses}/${gp}/analysis_datasink/seed_correlation_maps/_split_name_sub-*_ses-${ses}_task-rest_run-01_bold/_seed_name_${roi}_peaks_sphere_bin/sub-*_ses-${ses}_task-rest_run-01_bold_RAS_combined_cleaned_${roi}_peaks_sphere_bin_corr_map.nii.gz

  done
 done
done
