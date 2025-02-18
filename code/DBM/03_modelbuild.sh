#!/bin/bash

module load cobralab 

export QBATCH_CHUNKSIZE=8
export QBATCH_CORES=80
export QBATCH_NODES=1

path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo
atlas_model=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template_Masked.nii
atlas_mask=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Mask.nii
jacobians=${path}/jacobians
mkdir -p ${jacobians}

cd ${jacobians}

${path}/code/optimized_antsMultivariateTemplateConstruction/twolevel_modelbuild.sh ${path}/code/input_jacobians.csv \
 --walltime-nonlinear 05:00:00 \
 --masks ${path}/code/input_jacobians_mask.csv \
 --iterations 2 \
 --final-target \
 ${atlas_model} \
 --final-target-mask \
 ${atlas_mask}
 
done
chmod -R 777 ${jacobians}
