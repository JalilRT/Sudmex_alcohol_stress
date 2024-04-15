#!/bin/bash
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

singularity exec -B /scratch:/scratch /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/container/afni.sif  3dLMEr -prefix /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/Seed_based/models/M4L_fc_dataset_right_olfactory_bulb.nii.gz -jobs 1 	-mask /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/Seed_based/mask_registered.nii.gz         -model 'poly(Age,2)*Intake*Sex+Batch+(1|Subj)'         -qVars 'Age'         -gltCode interaction 'Intake : 1*Low -1*Ctrl Sex : 1*male -1*female'         -dataTable @/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/code/LME_SB/fc_dataset_right_olfactory_bulb.txt

