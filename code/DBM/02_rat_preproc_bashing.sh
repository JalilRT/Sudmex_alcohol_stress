#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab/2019b
module load  minc-toolkit-extras
#module load eigen/3.3.7
#module load fftw/3.3.8
#module load mrtrix/3.0.0

export PSILANTRO_PATH=/scratch/m/mchakrav/dangeles/PSILANTRO

path=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/analysis/smri/DBM_invivo
input=${path}/data
preproc=${path}/preproc

echo "### Running preprocessing script v.7 to in-vivo MRI ###"
echo ""
read -p "Do you want to run an only one volume or do the average? (i.e one or average): " vol
echo ""
read -p "What session do you want to run? (i.e. T1, T2, T3 or T5): " ses
echo ""
read -p "Do you want to run all subjects or a particular one? (i.e all or one): " subs

if [[ $vol == "average" ]]; then
 if [[ $subs == "all" ]]; then

   bash ${path}/code/02_rat_preproc.sh ${ses} sub-*

 elif [[ $subs == "one" ]]; then

   read -p "which subject do you want to run? (i.e 42): " sub
   bash ${path}/code/02_rat_preproc.sh	${ses} sub-${sub}

 else
  echo "wrong definition"
 fi

elif [[ $vol == "one" ]]; then

 read -p "which subject do you want to run? (i.e 42): " sub

 fslsplit ${input}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_run-02_T2w.nii.gz \
  ${preproc}/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_run-02_T2w_
 
 for i in $(ls -d ${preproc}/sub-${sub}); do 

  mkdir -p ${preproc}/$(basename $i)/ses-${ses}/anat/

  nii2mnc ${preproc}/$(basename $i)/ses-${ses}/anat/sub-${sub}_ses-${ses}_run-02_T2w_0000.nii.gz \
   ${preproc}/$(basename $i)/ses-${ses}/anat/sub-${sub}_ses-${ses}_run-02_T2w_0000.mnc

  bash ${path}/code/rat-preprocessing-v7.sh \
   ${preproc}/$(basename $i)/ses-${ses}/anat/sub-${sub}_ses-${ses}_run-02_T2w_0000.mnc \
   ${preproc}/$(basename $i)/ses-${ses}/anat/$(basename $i)_ses-${ses}_pp.mnc

  mnc2nii -nii ${preproc}/$(basename $i)/ses-${ses}/anat/$(basename $i)_ses-${ses}_pp.mnc \
   ${preproc}/$(basename $i)/ses-${ses}/anat/$(basename $i)_ses-${ses}_pp.nii

  mincpik -scale 20 -t ${preproc}/$(basename $i)/ses-${ses}/anat/$(basename $i)_ses-${ses}_pp.mnc \
   ${preproc}/$(basename $i)/ses-${ses}/anat/$(basename $i)_ses-${ses}_pp.jpg

  chmod -R 777 ${preproc}/

 done

fi
