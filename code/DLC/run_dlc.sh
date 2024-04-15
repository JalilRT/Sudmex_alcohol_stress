#!bin/bash
## SGE batch file - dlc
#$ -S /bin/bash
#$ -N deeplabcut
#$ -V
#$ -l mem_free=15G
#$ -pe openmp 8
#$ -j y
#$ -wd /mnt/MD1200B/egarza/jrasgado/Alcohol_model/Behavior/logs

module load anaconda3/2021.05
source activate deeplabcut

path='/mnt/MD1200B/egarza/jrasgado/Alcohol_model/Behavior'

python ${path}/dlc_options.py
