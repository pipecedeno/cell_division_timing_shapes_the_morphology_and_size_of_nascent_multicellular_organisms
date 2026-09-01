#!/bin/bash

# Configuration variables for edge degree range
min_edge_degree=5
max_edge_degree=20

# Fixed parameters
data_file="grande_doubling_time_dist_2024dec2.csv"
num_simulations=300
num_generations=100
propagule_selection="random"

echo "Starting simulation with edge degree range: $min_edge_degree to $max_edge_degree"
echo "================================================"

# Iterate through edge degree values
for edge_degree in $(seq $min_edge_degree $max_edge_degree); do
    # Create output filename based on edge degree
    output_name="grande_${edge_degree}e_${num_generations}g_${propagule_selection}"
    
    echo "Running simulation with edge degree: $edge_degree"
    echo "Output: $output_name"
    
    # Execute the Python script
    sim_frag_clust_edge_degree_after_all_cells_22feb2024.py \
        -d "$data_file" \
        -n $num_simulations \
        -g $num_generations \
        -e $edge_degree \
        -o "$output_name" \
        -c "$propagule_selection"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully completed edge degree $edge_degree"
    else
        echo "✗ Error occurred with edge degree $edge_degree"
    fi
    
    echo "------------------------------------------------"
done

echo "All simulations completed!"
