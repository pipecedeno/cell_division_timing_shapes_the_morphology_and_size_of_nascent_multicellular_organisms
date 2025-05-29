
## Codes used to generate results

### Code for figure 3: Network growth without fragmentation

```
sim_clust_no_fragmentation_15nov2024.py -d petite_doubling_time_dist_2023may30.csv -n 300 -o petite_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d petite_only_second_doubling_8feb2024.csv -n 300 -o petite-second-only_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d grande_doubling_time_dist_2024dec2.csv -n 300 -o grande_300n_200m -m 200
```

### updated code for figure 3:
```
sim_clust_no_fragmentation_27may2025.py -d petite_doubling_time_dist_2023may30.csv -n 300 -o petite_300n_1200m -m 1200

sim_clust_no_fragmentation_27may2025.py -d grande_doubling_time_dist_2024dec2.csv -n 300 -o grande_300n_1200m -m 1200

sim_clust_no_fragmentation_27may2025.py -d petite_only_second_doubling_8feb2024.csv -n 300 -o petite-second-only_300n_1200m -m 1200 
```

### Code for figure 4: Network growth with fragmentation keeping track of all cluster fractures

```
par_sim_frag_all_clusters_exponential_11mar2024.sh -d petite_doubling_time_dist_2023may30.csv -n 100 -o petite_14e_100g_all-clusters -g 9 -e 14 -t 10

par_sim_frag_all_clusters_exponential_11mar2024.sh -d petite_only_second_doubling_8feb2024.csv -n 100 -o petite-second-only_14e_100g_all-clusters -g 9 -e 14 -t 10

par_sim_frag_all_clusters_exponential_11mar2024.sh -d grande_doubling_time_dist_2024dec2.csv -n 100 -o grande_14e_100g_all-clusters -g 9 -e 14 -t 10
```
Note: this bash script calls the following scripts:  
- sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py: This code executes each simulation of network growth
- concat_files_output_18feb2024.py: simple code to concatenate all the files of each simulation into only one file


### Code for figure 5: Effects of delay and variation in cluster properties
To execute this code the file params_4diffs_4stds_17oct2024.csv is needed from causes_of_asynchrony_1oct2024.R, as this file contains the parameters of the log-normal distributions that are going to be used for the simulations.

```
exec_sim_frag_clust_from_file_params.sh -n 30 -g 50 -e 14 -c random -i ~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_5_delay_variation_in_cluster_properties/params_4diffs_4stds_17oct2024.csv
```

Note: this same data is used for Supplementary Figure 6



### Figure 6: Settling selection simulations

All the following codes are needed for figure 6, and supplementary figures 7 and 9.  

```
# Testing different carrying capacities
par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10k_50sim -e 14 -c 50:50 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 10000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_100k_50sim -e 14 -c 500:500 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 100000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_1mill_50sim -e 14 -c 5000:5000 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 1000000 -h 12 -i 50 -y random

par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10mill_50sim -e 14 -c 50000:50000 -r 20 -t 5 -m 50 -g 0.1 -p 0.1 -k 10000000 -h 12 -i 50 -y random
#Note: this code consume up to 20GB per thread, so adjust -t accordingly to the amount of RAM in your systems

#Different carrying capacities using bootstrap to increase population size
par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_10mill-boot_50sim -e 14 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 10000000 -h 12 -i 50 -u 30

par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_100mill-boot_50sim -e 14 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 100000000 -h 12 -i 50 -u 30

par_boot_set_sim_v0_16may2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_8feb2024.csv -n petite:second -o petite_second_1bill-boot_50sim -e 14 -c 5000:5000 -m 50 -t 10 -g 0.1 -p 0.1 -k 1000000 -b 1000000000 -h 12 -i 50 -u 30
```

Notes: 
- to perform the top-settling simulations the flag ```-y top``` needs to be added in the commands.
- par_set_sim_v0_17feb2024.sh executes the following scripts:
    - settling_selection_sim_v0_14feb2024.sh
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
sim_frag_clust_diff_mean_dists_24apr2024.py -n 100 -o synchronized_strain -g 100 -e 14 -i 4.4979,0.0821 -j 4.4979,0.0821 -c random

sim_frag_clust_diff_mean_dists_24apr2024.py -n 100 -o fast_first_div -g 100 -e 14 -i 4.0938,0.0809 -j 4.4979,0.0821 -c random

sim_frag_clust_diff_mean_dists_24apr2024.py -n 100 -o slow_first_div -g 100 -e 14 -i 4.7872,0.0843 -j 4.4979,0.0821 -c random
```

Note: this is the same code used for figure 5, just with different parameters.


### Table 1

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

### Supplementary figure 3: Validation of the network model with biophysical simulations

First some networks need to be created for the matlab code to be executed:
```
sim_edges_network_null_1nov2023.py -d grande_doubling_time_dist_2024dec2.csv -m 1000 -n 500 -o test_grande_1000m_500n

sim_edges_network_null_1nov2023.py -d petite_doubling_time_dist_2023may30.csv -m 1000 -n 500 -o test_petite_1000m_500n
```

Then the matlab code ``` size_diff_sync_vs_async_2024sep26_v3.m ``` is executed to calculate the size at fragmentation and also the overlaps between the cells.

Finally, the code matlab code ``` size_diff_overlap_acumulation_2024dec2.m ``` is executed to calculate how the overlap changes as a cell gets added.

Notes: both matlab codes need have their input parameters in the first lines of the code. For the overlap acumulation code only 100 networks were created instead of 500.

### Supplementary Figure 4:

The code used to generate the results is the following:

test_network_diam_normalization.ipynb

The value of the linear regression obtained here is what was used in the for the normalization, even though an R code was used to generate the plots.

### Supplementary figure 5: Network growth with fragmentation tracking only one cluster after fragmentation

```
sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_doubling_time_dist_2023may30.csv -n 500 -g 100 -e 14 -o petite_14e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_only_second_doubling_8feb2024.csv -n 500 -g 100 -e 14 -o petite-second-only_14e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d grande_doubling_time_dist_2024dec2.csv -n 500 -g 100 -e 14 -o grande_14e_100g_random -c random
```


### Supplementary figure 8: Settling selection simulations with same growth rates

```
par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_plus_15min_28june2024.csv -n petite:second-15min -o petite_second-15min_1mill_50sim -e 14 -c 5000:5000 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 1000000 -h 12 -i 50 -y random
```


### Codes to create visualizations

- cell_overlap_viz_26sep2024.m: was used to create the 3D visualizations shown in supplementary figure 3.
- network_plot_no_fragmentation.ipynb: code used to create the networks in figure 3.
- syn_first_div_filaments_3oct2024.ipynb: code used to create the networks in figure 8.




