#!/bin/bash

# Check if the CSV file exists
if [ ! -f "params_4diffs_4stds_17oct2024.csv" ]; then
    echo "Error: CSV file not found!"
    exit 1
fi

# Skip the header line and process each row of the CSV file
tail -n +2 "params_4diffs_4stds_17oct2024.csv" | while IFS=',' read -r delay_level var_level col3 col4 log_mean1 log_std1 col7 col8 log_mean2 log_std2 col11
do
    # Construct the command
    cmd="sim_frag_clust_diff_mean_dists_24apr2024.py -n 30 -o test_${var_level}_var_${delay_level}_diff -g 50 -e 14 -i ${log_mean1},${log_std1} -j ${log_mean2},${log_std2} -c random"
    
    # Execute the command
    echo "Executing: $cmd"
    eval $cmd
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Command executed successfully."
    else
        echo "Error: Command failed."
    fi
    
    echo "----------------------------------------"
done

echo "All commands have been executed."
