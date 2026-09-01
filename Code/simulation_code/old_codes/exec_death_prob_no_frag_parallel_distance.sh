 #!/bin/bash
# Code to execute simulations with different death probabilities for cluster growth without fragmentation.
# This script runs simulations across a range of death probabilities in parallel.
#
# Inputs:
# -n: Number of simulations to be executed for each death probability
# -m: Maximum cluster size (number of cells the clusters are allowed to have)
# -i: First division parameters (mean,variance) for log-normal distribution
# -j: Second division parameters (mean,variance) for log-normal distribution  
# -p: Minimum death probability (0.0 to 1.0)
# -P: Maximum death probability (0.0 to 1.0)
# -s: Step size for death probability increments
# -t: Number of threads/parallel processes to use

while getopts "n:m:i:j:p:P:s:t:" option
do
case "${option}"
in
n) sim_number=${OPTARG};;
m) max_clust_size=${OPTARG};;
i) first_div_params=${OPTARG};;
j) second_div_params=${OPTARG};;
p) min_death_prob=${OPTARG};;
P) max_death_prob=${OPTARG};;
s) death_prob_step=${OPTARG};;
t) num_threads=${OPTARG};;
*) echo "Invalid option or missing option argument for -$OPTARG" >&2; exit 1;;
esac
done

# Check if required arguments are provided
if [ -z "$sim_number" ] || [ -z "$max_clust_size" ] || [ -z "$first_div_params" ] || [ -z "$second_div_params" ] || [ -z "$min_death_prob" ] || [ -z "$max_death_prob" ] || [ -z "$death_prob_step" ] || [ -z "$num_threads" ]; then
    echo "Error: Missing required arguments!"
    echo "Usage: $0 -n <number_of_simulations> -m <max_cluster_size> -i <first_div_params> -j <second_div_params> -p <min_death_prob> -P <max_death_prob> -s <death_prob_step> -t <number_of_threads>"
    echo ""
    echo "Example: $0 -n 100 -m 500 -i '2.5,0.5' -j '2.0,0.3' -p 0.0 -P 0.5 -s 0.1 -t 4"
    exit 1
fi


echo "Starting simulations with parameters:"
echo "Number of simulations per death probability: ${sim_number}"
echo "Maximum cluster size: ${max_clust_size}"
echo "First division parameters: ${first_div_params}"
echo "Second division parameters: ${second_div_params}"
echo "Death probability range: ${min_death_prob} to ${max_death_prob} (step: ${death_prob_step})"
echo "Number of parallel threads: ${num_threads}"
echo "========================================"

# Function to execute a single simulation
run_simulation() {
    local death_prob=$1
    local sim_number=$2
    local max_clust_size=$3
    local first_div_params=$4
    local second_div_params=$5
    
    # Format death probability for directory name using printf for consistent formatting
    # This ensures consistent formatting like 0.00, 0.05, 0.10, etc.
    local death_prob_formatted=$(printf "%.2f" "$death_prob" | sed 's/\./_/g')
    
    # Create output directory name
    output_dir="test_${death_prob_formatted}"
    
    # Construct the command
    cmd="sim_clust_no_fragmentation_lognormal_death_prob_distance_5aug2025.py -n ${sim_number} -o ${output_dir} -m ${max_clust_size} -i ${first_div_params} -j ${second_div_params} -p ${death_prob}"
    
    # Execute the command
    echo "Executing: $cmd"
    eval $cmd
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Command executed successfully for death_prob=${death_prob}"
    else
        echo "Error: Command failed for death_prob=${death_prob}"
        return 1
    fi
    echo "----------------------------------------"
}

# Export the function and variables so they're available to parallel processes
export -f run_simulation
export sim_number
export max_clust_size
export first_div_params
export second_div_params

# Generate death probability values and create a temporary file
temp_file=$(mktemp)

# Use bc for floating point arithmetic to generate the range
current_prob=$min_death_prob
while (( $(echo "$current_prob <= $max_death_prob" | bc -l) )); do
    echo "$current_prob" >> "$temp_file"
    current_prob=$(echo "$current_prob + $death_prob_step" | bc -l)
done

# Count how many death probabilities we'll test
num_death_probs=$(wc -l < "$temp_file")
echo "Testing ${num_death_probs} different death probabilities..."
echo ""

# Use GNU parallel for execution
echo "Using GNU parallel for execution..."
parallel -j ${num_threads} run_simulation {} ${sim_number} ${max_clust_size} ${first_div_params} ${second_div_params} :::: ${temp_file}

# Check if parallel command was successful
if [ $? -eq 0 ]; then
    echo ""
    echo "All simulations completed successfully!"
    echo "Output directories created: test_0_0, test_0_1, test_0_2, etc."
else
    echo ""
    echo "Some simulations may have failed. Check the output above for details."
fi

# Clean up temporary file
rm ${temp_file}

echo "========================================"
echo "Script execution completed."