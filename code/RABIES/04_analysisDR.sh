#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

# load necessary modules
module load cobralab

path=/scratch/m/mchakrav/jrasgado/sudmex_b45/func
atlas_folder=${path}/Atlas/

atlas=${atlas_folder}/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Brain_Atlases/SIGMA_Anatomical_Atlas/InVivo_Atlas/SIGMA_InVivo_Anatomical_Brain_Atlas.nii

container=${path}/container/rabies_v047.sif

### define input/output folders ###
ses=$1
group=$2

data_input=${path}/data/Session-${ses}
preproc_output=${path}/rabies/preproc_out/Session-${ses}
confound_out=${path}/rabies/confound_out/Session-${ses}
analysis_out=${path}/rabies/analysis_out/Session-${ses}
analysisDR_out=${path}/rabies/analysisDR_out

mkdir -p ${analysisDR_out}/${group}

### Creating list to use for scan_list ###

for i in $(cat code/list_${group}.txt); do 
 cuca=$(echo $(ls data/Session-${ses}/${i}/ses-${ses}/func/${i}_ses-${ses}_task-rest_run-01_bold.nii.gz))
 coco=$(basename $cuca)
 echo $coco >> ${group}_${ses}.txt
 sed -i '/^$/d' ${group}_${ses}.txt
done

#RABIES call
singularity run -B /scratch:/scratch -B ${data_input}:/data_input:ro -B ${atlas_folder}:/Atlas \
-B ${preproc_output}:/preproc_output -B ${confound_out}:/confound_out -B ${analysisDR_out}/${group}:/analysis \
${container} -p Linear \
analysis /confound_out /analysis \
--DR_ICA --prior_map ${analysis_out}/${group}/analysis_main_wf/analysis_wf/group_ICA/group_melodic.ica/melodic_IC.nii.gz \
--scan_list ${path}/${group}_${ses}.txt
