#!/bin/bash

module load cobralab && module load minc-toolkit ANTs

export PYTHONUNBUFFERED=TRUE
export PATH=/home/m/mchakrav/egarza/twolevel_ants_dbm:$PATH

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo
jacobians=${path}/jacobians
mkdir -p ${jacobians}

cd ${jacobians}

twolevel_dbm.py --jacobian-sigmas 0.12 0.17 \
 --modelbuild-command ${QUARANTINE_PATH}/twolevel-dbm/niagara2-antsMultivariateTemplateConstruction2.sh \
 --cluster-type=slurm \
 --rigid-model-target /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template_Masked.nii \
 2level /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo/code/input_jacobians.csv 2>&1 | tee -a dbm_logfile.log


chmod -R 777 ${jacobians}
