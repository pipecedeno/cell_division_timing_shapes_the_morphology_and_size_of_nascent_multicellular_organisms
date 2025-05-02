# Emergence of coordinated cell division during the evolution of multicellularity

Things to add:
- Conda environment and versions (also upload environment yml file, and add command to create environment with the file)
- Explanation of how to run the codes (adding them to path) (explanation on how to run the matlab simulations)
- document with which simulations code were used for which figures
- List of commands used to run each of the simulation results
- R codes:
    - Clean the codes to only have important parts
    - Change the input files (save single file and just have commented code on how to load the codes when they are in different files)
- simulation codes (change name of files to something easier to understand)

-Environment:
    -Python 3.10.12
    -matplotlib=3.8.2
    -Networkx=3.2.1
    -pyvis=0.3.2
    -pandas=2.2.0
    -numpy=1.26.4
    -parallel=20240122

Commands to create the conda environment: 
```
# Create a new conda environment with Python 3.10.12
conda create -n network_sim python=3.10.12

# Activate the new environment
conda activate network_sim

# Install the required packages 
conda install matplotlib=3.8.2 networkx=3.2.1 pandas=2.2.0 numpy=1.26.4 pyvis=0.3.2
conda install -c conda-forge parallel=20240122
```

Alternatively, the conda environment can be crated by executing ```conda env create -f environment.yml``` with the environment.yml file being either ```network_sim_macos.yml``` or ```network_sim_ubuntu.yml```.

Matlab version used: R2023b

## How to execute the code?

All the code was executed in macOS Sequoia 15.4.1 and Linux (Ubuntu 22.04.2).

In order to execute the code used 

