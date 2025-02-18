#!/bin/bash

module load cobralab && module load minc-toolkit ANTs

path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis/smri/DBM_invivo
jacobians=${path}/jacobians

cd ${jacobians}

${path}/code/optimized_antsMultivariateTemplateConstruction/twolevel_dbm.sh --walltime 05:00:00 \
 --jacobian-smooth 4vox,1mm ${path}/code/input_jacobians.csv \
 --mask ${jacobians}/output/secondlevel/final/average/mask_shapeupdate.nii.gz

chmod -R 777 ${jacobians}
