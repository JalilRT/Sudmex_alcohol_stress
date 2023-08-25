#!/bin/bash
# @Author: Jalil Rasgado
# Date: 2022-12-24
# Description: This script is intended to crete a conda environment with speficic dependencies to reproduce the code of sudmex-alcohol. 
# 	Before running it make sure to have an installation of Anaconda or Miniconda in your system.


### a comand to export the envieronment. DO Uncomment this line. 


 echo -e "\e[0;35m"

 echo 
 echo " ++ Creating the environment" 
 echo

# If conda exist start environment creation

 if [ -z $(which conda) ] ; then 

 echo -e "\e[0;31m" 
 echo " ++ Error: No conda installation was found. Program exits!! "

	exit 1
 else

	conda create  --force -y -n sudmex_alcohol_stress  -c conda-forge python=3.8 jupyter jupyterlab notebook \
		r-base=4.2 r-tidyverse r-rlang r-pacman r-biocmanager r-rstatix  r-devtools r-viridis r-irkernel
	echo
	echo " ++ Base environment created "
	echo " ++ Installing additional dependencies"
	echo

	conda install -y -n sudmex_alcohol_stress -c  bioconda bioconductor-summarizedexperiment 


 fi


echo -e "\e[0m"

exit
