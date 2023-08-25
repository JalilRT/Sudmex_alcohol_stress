# Enviroment installation and configuration
---

## Instruction


1. Launch a conda or bash terminal

2. Clone this repository and navigate into it:

  ```{shell}
  # cloning
  git clone https://github.com/psilantrolab/Sudmex-alcohol-rat
 
  # get into the directory
  cd Sudmex-alcohol-rat
  ```
  
3. Create the Conda environment:

  ```{shell}
 conda create --force -y -v --name sudmex_alcohol --file utils/sudmex_alcohol_spec-file.txt
  ```
  
4. Activate the Conda environment:
  
  ```{shell}
  conda activate sudmex_alcohol
  conda install -c r r-irkernel
  ```
  
5. Launch the Jupyter notebook: 

  ```{shell}
  # example
  jupyter-lab code/Alcohol_consumption.ipynb
  ```
