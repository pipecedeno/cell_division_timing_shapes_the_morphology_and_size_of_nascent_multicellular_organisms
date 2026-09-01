
## Codes used to generate results


### Code for figure 3: Network growth without fragmentation

```
sim_clust_no_fragmentation_15nov2024.py -d petite_doubling_time_dist_2023may30.csv -n 300 -o petite_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d petite_only_second_doubling_8feb2024.csv -n 300 -o petite-second-only_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d grande_doubling_time_dist_2024dec2.csv -n 300 -o grande_300n_200m -m 200
```


### Code for figure 4: Network growth with fragmentation keeping track of all cluster fractures

```
par_sim_frag_all_clusters_exponential_11mar2024.sh -d petite_doubling_time_dist_2023may30.csv -n 100 -o petite_15e_100g_all-clusters -g 9 -e 15 -t 10

par_sim_frag_all_clusters_exponential_11mar2024.sh -d petite_only_second_doubling_8feb2024.csv -n 100 -o petite-second-only_15e_100g_all-clusters -g 9 -e 15 -t 10

par_sim_frag_all_clusters_exponential_11mar2024.sh -d grande_doubling_time_dist_2024dec2.csv -n 100 -o grande_15e_100g_all-clusters -g 9 -e 15 -t 10
```
Note: this bash script calls the following scripts:  
- sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py: This code executes each simulation of network growth
- concat_files_output_18feb2024.py: simple code to concatenate all the files of each simulation into only one file

### Code for figure 5: Effects of delay and variation in cluster properties
To execute this code the file params_4diffs_4stds_23sep2025.csv is needed from causes_of_asynchrony_1oct2024.R, as this file contains the parameters of the log-normal distributions that are going to be used for the simulations.

```
exec_sim_frag_clust_from_file_params.sh -n 30 -g 50 -e 15 -c random -i ~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_5_delay_variation_in_cluster_properties/params_4diffs_4stds_23sep2025.csv
```
The bash script executes the following python script:  
- sim_frag_clust_diff_mean_dists_24apr2024.py  

### Figure 6: Settling selection simulations

All the following codes are needed for figure 6, and supplementary figures 7 and 9.  

```
# Testing different carrying capacities
par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10k_50sim -e 15 -c 50:50 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 10000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_100k_50sim -e 15 -c 500:500 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 100000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_1mill_50sim -e 15 -c 5000:5000 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 1000000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10mill_50sim -e 15 -c 50000:50000 -r 20 -t 5 -m 50 -g 0.1 -p 0.1 -k 10000000 -h 12 -i 50 -y random
#Note: this code consume up to 20GB per thread, so adjust -t accordingly to the amount of RAM in your systems

# Different carrying capacities using bootstrap to increase population size
par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10mill-boot_50sim -e 15 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 10000000 -h 12 -i 50 -u 30 -y random

par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_100mill-boot_50sim -e 15 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 100000000 -h 12 -i 50 -u 30 -y random

par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_1bill-boot_50sim -e 15 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 1000000000 -h 12 -i 50 -u 30 -y random

# Simulations with equal growth rate for grande and petite strains
par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_plus_15min_28june2024.csv -n petite:second-15min -o petite_second-15min_1mill_50sim -e 15 -c 5000:5000 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 1000000 -h 12 -i 50 -y random
```

Notes: 
- to perform the top-settling simulations the flag ```-y top``` needs to be added in the commands.
- par_set_sim_v0_17feb2024.sh executes the following scripts:
    - settling_selection_sim_v0_15feb2024.sh
    - initial_cluster_growth_15mar2024.py
    - growth_phase_clusters_v0_9feb2024.py
    - settling_selection_phase_v0_15feb2024.m
    - write_settling_selection_output_v0_16feb2024.py
- par_boot_set_sim_v0_16may2024.sh executes the following scripts:
    - bootstrap_settling_sim_13may2024.sh
    - initial_cluster_growth_15mar2024.py
    - growth_phase_clusters_v0_9feb2024.py
    - settling_selection_bootstrap_v0_10may2024


### Figure 8: Fast first division simulations


```
exec_sim_fast_first_div_frag_from_file_params.sh -n 30 -g 50 -u 5 -x 15 -c random -i ../params_distributions_50percent_15var_7july2025.csv

exec_sim_fast_first_div_no_frag_from_file_params_parallel.sh -n 300 -m 1300 -i ../params_distributions_50percent_15var_7july2025.csv -t 8
```
The following python scripts are executed inside of the bash scripts:  
- sim_frag_clust_diff_mean_dists_24apr2024.py  
- sim_clust_no_fragmentation_lognormal_27may2025.py  

### Supplementary Table 3

1.- Find the threshold that will obtain the desired volume for the grande strain
The matlab simulations were done using the scripts:
```
snowflake_asp_ratio_volume_aprox.m
snowflake_cell_diam_volume_aprox.m
snowflake_synchrony_volume_aprox.m
```
With the data created from this scripts, the results were loaded to an R script (volume_predictions_8jan2023.R, all the analysis are in this script), and a linear regression was fitted to the data to obtain what threshold value would be the best for the following step.

Commands used to create the network files for the synchrony predictions
```
sim_edges_network_null_1nov2023.py -d grande_doubling_time_dist_2024dec2.csv -m 3000 -n 50 -o grande_3000m_50n

sim_edges_network_null_1nov2023.py -d petite_doubling_time_dist_2023may30.csv -m 3000 -n 50 -o petite_3000m_50n
```

2.- Do simulation with constant threshold and vary the desired variable
This matlab simulations were done using the following scripts:
```
snowflake_different_aspect_ratios.m
snowflake_different_cell_diameters.m
snowflake_synchrony_volume_prediction.m
```
With this data a linear regression was fitted to the data, and a the difference from the predictions was obtained by using the actual parameter for each of the strains.

The following code was executed to obtain the result of all variabled tested combined:
```
snowflake_synchrony_volume_prediction_combined_effect.m
```

### Supplementary Figure 3: Finding edge degree threshold for simulations

```
diff_edge_degrees_grande_28may2025.sh
```
Note: this code executes the python script ```sim_frag_clust_edge_degree_after_all_cells_22feb2024.py``` same as supplementary figure 6.


### Supplementary figure 4: Validation of the network model with biophysical simulations

First some networks need to be created for the matlab code to be executed:
```
sim_edges_network_null_1nov2023.py -d grande_doubling_time_dist_2024dec2.csv -m 1000 -n 500 -o test_grande_1000m_500n

sim_edges_network_null_1nov2023.py -d petite_doubling_time_dist_2023may30.csv -m 1000 -n 500 -o test_petite_1000m_500n
```

Then the matlab code ``` size_diff_sync_vs_async_2024sep26_v3.m ``` is executed to calculate the size at fragmentation and also the overlaps between the cells.

To generate the statistics of the overlap data the following python code was used:
```
process_overlap_pos_physics_sim_1nov2024.py -i . -f _overlap_pos_30sim_500n_1.2aspr_10attempts_70overlap.csv -o stats_overlap_pos_10jun2025.csv
```

Finally, the code matlab code ``` size_diff_overlap_acumulation_2024dec2.m ``` is executed to calculate how the overlap changes as a cell gets added.

Notes: both matlab codes need have their input parameters in the first lines of the code. For the overlap acumulation code only 100 networks were created instead of 500.



### Supplementary Figure 5: Normalization method for network diameter

The code used to generate the results is the following:  
```
test_network_diam_normalization.ipynb
```
The value of the linear regression obtained here is what was used in the for the normalization, even though an R code was used to generate the plots.  


### Supplementary figure 6: Network growth with fragmentation tracking only one cluster after fragmentation

```
sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_doubling_time_dist_2023may30.csv -n 500 -g 100 -e 15 -o petite_15e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_only_second_doubling_8feb2024.csv -n 500 -g 100 -e 15 -o petite-second-only_15e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d grande_doubling_time_dist_2024dec2.csv -n 500 -g 100 -e 15 -o grande_15e_100g_random -c random
```


### Supplementary Figure 10 & 11

```
exec_death_prob_no_frag_parallel_delays.sh
```
Bash script runs the following python script:  
sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py  


### Codes to create visualizations

- cell_overlap_viz_26sep2024.m: was used to create the 3D visualizations shown in supplementary figure 4.  
- network_plot_no_fragmentation_27may2025.ipynb: code used to create the networks in figure 3.  
- syn_first_div_filaments_3oct2024.ipynb: code used to create the networks in figure 8.  
- visualization_delay_constant.ipynb: code used to create the networks in supplementary figure 10.  


