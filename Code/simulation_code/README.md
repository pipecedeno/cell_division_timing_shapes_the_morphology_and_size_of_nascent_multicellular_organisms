
## Codes used to generate results

### Code for figure 3: Network growth without fragmentation

```
sim_clust_no_fragmentation_15nov2024.py -d petite_doubling_time_dist_2023may30.csv -n 300 -o petite_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d petite_only_second_doubling_8feb2024.csv -n 300 -o petite-second-only_300n_200m -m 200

sim_clust_no_fragmentation_15nov2024.py -d grande_doubling_time_dist_2024dec2.csv -n 300 -o grande_300n_200m -m 200
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
exec_sim_frag_clust_from_file_params.sh
```

Note: this same data is used for Supplementary Figure 5

# Note: update this code so that it can be executed with flags for the other parameters and for the input file


### Figure 6: Settling selection simulations

All the following codes are needed for figure 6, and supplementary figures 6 and 8.  

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


### Supplementary figure 5: Validation of the network model with biophysical simulations

First some networks need to be created for the matlab code to be  
```

```

### Supplementary figure 4: Network growth with fragmentation tracking only one cluster after fragmentation

```
sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_doubling_time_dist_2023may30.csv -n 500 -g 100 -e 14 -o petite_14e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d petite_only_second_doubling_8feb2024.csv -n 500 -g 100 -e 14 -o petite-second-only_14e_100g_random -c random

sim_frag_clust_edge_degree_after_all_cells_22feb2024.py -d grande_doubling_time_dist_2024dec2.csv -n 500 -g 100 -e 14 -o grande_14e_100g_random -c random
```


### Supplementary figure 7: Settling selection simulations with same growth rates

```
par_set_sim_v0_17feb2024.sh -f petite_doubling_time_dist_2023may30.csv -s petite_only_second_doubling_plus_15min_28june2024.csv -n petite:second-15min -o petite_second-15min_1mill_50sim -e 14 -c 5000:5000 -r 20 -t 10 -m 50 -g 0.1 -p 0.1 -k 1000000 -h 12 -i 50 -y random
```



