#!/bin/bash
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

#export QBATCH_NODES=1
#export QBATCH_SYSTEM="slurm"
module load cobralab

#################################
## 2023-09-14
## Created by JalilRT

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri
mkdir -p ${path}/code/LME_SB/models
mkdir -p ${path}/Seed_based/models
container=${path}/container/afni.sif
cd ${path}
fc=$(ls ${path}/code/LME_SB/fc_dataset_*)

#################################
for i in ${fc}; do
rois=$(basename ${i} .txt)

cat << EOF > ${path}/code/LME_SB/models/M4L_${rois}.sh
#!/bin/bash
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

singularity exec -B /scratch:/scratch ${container} \
 3dLMEr -prefix /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/Seed_based/models/M4L_${rois}.nii.gz -jobs 1 \
	-mask /scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/Seed_based/mask_registered.nii.gz \
        -model 'poly(Age,2)*Intake*Sex+Batch+(1|Subj)' \
        -qVars 'Age' \
        -gltCode interaction 'Intake : 1*Low -1*Ctrl Sex : 1*male -1*female' \
        -dataTable @/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri/code/LME_SB/${rois}.txt

EOF

sbatch ${path}/code/LME_SB/models/M4L_${rois}.sh

done
