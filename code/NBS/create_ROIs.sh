#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=rrg-mchakrav-ab

module load cobralab

size=$1

path=/scratch/m/mchakrav/jrasgado/sudmex_alc_stress/analysis
ROIs_folder=${path}/fmri/Atlas/ROIs
mkdir -p $ROIs_folder

ROIs_path=${path}/fmri/Atlas/ROIs/ROI_${size}
#Template_fmri=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Functional_Imaging/SIGMA_EPI_Brain_Template_Masked.nii
Template_fmri=/scratch/m/mchakrav/jrasgado/sudmex_alcohol_rat/Atlas/SIGMA_Wistar_Rat_Brain_TemplatesAndAtlases_Version1.2.1/SIGMA_Rat_Anatomical_Imaging/SIGMA_Rat_Anatomical_InVivo_Template/SIGMA_InVivo_Brain_Template.nii
mkdir -p ${ROIs_path}

CELESTE="\e[;36m"
RESET='\e[0m'
YELLOW='\e[33m'

echo -e "${YELLOW} ## Creating ROIs based on smri in-vivo results ## ${RESET}"
	
vector_x=$(cut -d ',' -f 4 ${path}/smri/DBM_invivo/DBM/smooth_1mm/Trayectories/ROIs.csv)
vector_y=$(cut -d ',' -f 5 ${path}/smri/DBM_invivo/DBM/smooth_1mm/Trayectories/ROIs.csv)
vector_z=$(cut -d ',' -f 6 ${path}/smri/DBM_invivo/DBM/smooth_1mm/Trayectories/ROIs.csv)
vector_names=$(cut -d ',' -f 11 ${path}/smri/DBM_invivo/DBM/smooth_1mm/Trayectories/ROIs.csv)

nrow=$(cat ${path}/smri/DBM_invivo/DBM/smooth_1mm/Trayectories/ROIs.csv | wc -l)

for i in $(seq 2 ${nrow} ); do

names_o=$(echo $vector_names | cut -d ' ' -f $i)
echo -e "${CELESTE} $names_o ${RESET}"

echo -e "## Converting mm2RAS ##"
x_o=$(echo $vector_x | cut -d ' ' -f $i)
y_o=$(echo $vector_y| cut -d ' ' -f $i)
z_o=$(echo $vector_z | cut -d ' ' -f $i)

x=`echo "${x_o} * (6.66) + 62.93" | bc`
y=`echo "${y_o} * (6.66) + 118.73" | bc`
z=`echo "${z_o} * (6.66) + 45.10" | bc`

echo -e "## Rounding ##"
round_x=$(printf "%.${precision}f" $x)
round_y=$(printf "%.${precision}f" $y)
round_z=$(printf "%.${precision}f" $z)

echo -e "## Creating the ROIs ##"

echo fslmaths $Template_fmri -mul 0 -add 1 -roi $x_o 1 $y_o 1 $z_o 1 0 1 ${ROIs_path}/${names_o} -odt float
echo fslmaths $Template_fmri -mul 0 -add 1 -roi $round_x 1 $round_y 1 $round_z 1 0 1 ${ROIs_path}/${names_o} -odt float

fslmaths $Template_fmri -mul 0 -add 1 -roi $round_x 1 $round_y 1 $round_z 1 0 1 ${ROIs_path}/${names_o} -odt float

echo -e "## Inflate to an sphere ##"
fslmaths ${ROIs_path}/${names_o} -kernel sphere ${size} -fmean ${ROIs_path}/${names_o}_sphere -odt float

fslmaths ${ROIs_path}/${names_o}_sphere -bin ${ROIs_path}/${names_o}_sphere_bin 

fslmaths ${ROIs_path}/${names_o}_sphere_bin -mul ${i} ${ROIs_path}/${names_o}_sphere_bin_tagged

done

## Select only the ROIs for NBR ##
fslmaths ${Template_fmri} -mul 0 ${ROIs_folder}/ROI_fmri_mask_${size}
for i in $(cat ${ROIs_folder}/rois_seed.txt); do fslmaths ${ROIs_folder}/ROI_fmri_mask_${size}.nii.gz \
 -add ${ROIs_path}/${i}_sphere_bin_tagged.nii.gz ${ROIs_folder}/ROI_fmri_mask_${size}.nii.gz ; done

chmod -R 777 ${ROIs_path}

for file in Atlas/ROIs/ROI_0.5/*; do
    mv "$file" "$(echo $file | sed 's/\^M//g')"
done