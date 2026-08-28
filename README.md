# Cell division timing shapes the morphology and size of nascent multicellular organisms


## Conda environment creation:

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

In order to use the code from this repository you need to give execution permission to all the python and bash scripts, this can be done with the following two commands:

```
chmod +x *.py
chomd +x *.sh
```

Next, the folder Code/simulation_code needs to be added to the PATH variable of the terminal, that can be done by adding the complete path of the folder into the .zshrc or .bashrc depending if you are in macOS or Linux, respectively. In the file the line needs to be added if not already there: ```export PATH="$PATH:<whole path>/Code/simulation_code/```

### Simulations with matlab

For the settling selection simulations, given that this ones use matlab executed from the terminal, the path also needs to be added within the matlab app so that the code can be executed from anywhere in the terminal. First, make sure that matlab can be opened in the terminal, this can be done by typing ```matlab``` in the terminal, if this command doesn't work then the executable file needs to be added to the PATH variable. Then, within matlab the following commands need to be executed to add the folder Code/simulation_code in the matlab path:

```
addpath('/path/to/Code/simulation_code');
savepath;
```

Now the bash scripts should be capable of executing the matlab commands.

## Github Structure

The github has the following structure:

```
emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/
├── Code
│   ├── R_code
│   ├── timelapse_analysis
│   └── simulation_code
│       ├── table_1
│       └── visualization_codes
├── Data
│   ├── fig_3_network_growth_no_fragmentation
│   ├── fig_4_network_growth_with_fragmentation
│   ├── fig_5_delay_variation_in_cluster_properties
│   ├── fig_8_fast_first_division
│   ├── Images
│   ├── network_files
│   ├── settling_simulation_data
│   ├── supp_fig_3_size_difference_physics_sim
│   └── supp_fig_4_growth_with_fragmentation_all
├── Paper_figures
└── Supplementary_videos
```

In the Code folder the simulation_code includes all the python, bash, and matlab scripts required to generate all the data of the results. The folder R_code includes all the R scripts that were used to analyze the data from the simulations and generate all the plots of the article. And the timelapse_analysis folder contain the python scripts needed to generate the dataset of doubling time distribution and cell synchrony. Each of these folders has their own README.md file with instructions of which codes were used for each figure and the commands used in the simulations.  




