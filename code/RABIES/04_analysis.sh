#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --account=rrg-mchakrav-ab

# load necessary modules
module load cobralab

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/fmri
atlas_folder=${path}/Atlas/
TR=1.0

#atlas=${atlas_folder}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii
#atlas=${atlas_folder}/Mask_erode.nii
atlas=${atlas_folder}/ROI_fmri_mask.nii

container=${path}/container/rabies_v047.sif

### define input/output folders ###
ses=$1
group=$2

data_input=${path}/data/Session-${ses}
preproc_output=${path}/rabies/preproc_out/Session-${ses}
confound_out=${path}/rabies/confound_out/Session-${ses}
analysis_out=${path}/rabies/analysis_out/Session-${ses}

mkdir -p ${analysis_out}/${group}

### Creating list to use for scan_list ###

for i in $(cat code/list_${group}.txt); do 
 cuca=$(echo $(ls data/Session-${ses}/${i}/ses-${ses}/func/${i}_ses-${ses}_task-rest_run-01_bold.nii.gz))
 coco=$(basename ${cuca})
 echo $coco >> ${group}_${ses}.txt
 sed -i '/^$/d' ${group}_${ses}.txt
done

# If you need to remove a subject from analysis write the subject (i.e. sub-34) #
if [ $ses == "T1" ]
then
  echo "removing subject"
  sed -i '/sub-34/d' ${group}_${ses}.txt
else
  echo "without any exclusion"
fi

if [ $ses == "T5" ]
then
  echo "removing subject"
  sed -i '/sub-43/d' ${group}_${ses}.txt
else
  echo "without any exclusion"
fi

### RABIES call ###
singularity run -B /scratch:/scratch -B ${data_input}:/data_input:ro -B ${atlas_folder}:/Atlas \
-B ${preproc_output}:/preproc_output -B ${confound_out}:/confound_out -B ${analysis_out}/${group}:/analysis \
${container} -p Linear \
analysis /confound_out /analysis \
--data_diagnosis \
--group_ica apply=true,dim=30,random_seed=1 \
--scan_list ${path}/${group}_${ses}.txt

chmod -R 777 ${analysis_out}

rm ${group}_${ses}.txt
