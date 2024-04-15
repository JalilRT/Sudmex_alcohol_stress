#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

# load necessary modules
module load cobralab

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri
mkdir -p ${path}/rabies/melodic
mask=${path}/rabies/preproc_out/Session-T1/bold_datasink/commonspace_mask/_scan_info_subject_id135.sessionT1_split_name_sub-135_ses-T1_run-02_T2w/_run_1/sub-135_ses-T1_task-rest_run-01_bold_RAS_EPI_brain_mask.nii.gz

melodic -i code/inputs_melodic.txt -m ${mask} \
 --outdir=${path}/rabies/melodic_prueba --dim=30 --tr=1 --report --seed=1
